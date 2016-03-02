#pragma once
#ifndef _NLMEANS_H_
#define _NLMEANS_H_

//mem leak checking
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>


#include "core.h"
#include "util.h"

#include <omp.h>
#include <algorithm>
#include <vector>
#include <assert.h>
#include <iostream>
#include <ctime>

template<typename T> class DenoiserInput {
	// Generic class to encapsulate denoiser input - useful since different denoisers have
	// different set of necessary bitmaps
public:
	DenoiserInput(std::string filename = "./renderoutput") :m_basefilename(filename) {}

	/// Serialize the filter to a binary data stream
	inline bool isInitialized() { return m_isinitialized; }
	inline std::string getBaseFilename() const { return m_basefilename; }
	inline std::vector<TImageBlock<T>*> getImageBlocks() { return m_imageblocks; }
	inline const TImageBlock<T>* getImageBlock(int i = 0) const { 
		if (i < m_imageblocks.size())
			return m_imageblocks[i];
		else
			return NULL;
	}
	void addImageBlock(TImageBlock<T>* imageblock) { m_imageblocks.push_back(imageblock); }

protected:
	~DenoiserInput() {}
protected:
	bool m_isinitialized;
	std::string m_basefilename; //base path to dump bitmaps
	std::vector<TImageBlock<T>*> m_imageblocks; // input imageblocks
};

template<typename T> class DenoiserOutput {
	// Generic class to encapsulate denoiser output including statistics for next iteration
	// in case of iterative denoising
public:
	long m_denoiseduration; // time in milliseconds to denoise the input
	DenoiserOutput(TImageBlock<T> *denoisedoutput = NULL) : m_denoisedimage(denoisedoutput) {}
	template<typename O> DenoiserOutput(DenoiserInput<O> *denoiserinput) {
		m_denoisedimage = new TImageBlock<T>(denoiserinput->getImageBlocks()[0]);
		m_denoisedimage->clear();
	}
	inline TImageBlock<T>* getDenoisedImage() { return m_denoisedimage; }
	inline long getDenoiseDuration() { return m_denoiseduration; }
protected:
	virtual ~DenoiserOutput() {}
protected:
	TImageBlock<T>* m_denoisedimage;
};

template<typename I, typename O> class NLMeansDenoiser {
public:
	// calculate parameters automatically
	NLMeansDenoiser(Float sigma, bool dumpm = true, bool dumpparams = true) : m_sigma(sigma), m_dumpmaps(dumpm) {
		if (sigma > 0.0f && sigma <= 25.0f) {
			m_f = 1;
			m_r = 10;
			m_k = 0.55f;
		}
		else if (sigma > 25.0f && sigma <= 55.0f) {
			m_f = 2;
			m_r = 17;
			m_k = 0.4f;
		}
		else if (sigma <= 100.0f) {
			m_f = 3;
			m_r = 17;
			m_k = 0.35f;
		}
		else {
			LOG(EWarn, "Variance cannot be above 100.f - reverting to default parameters!");
			m_f = 1;
			m_r = 10;
			m_k = 0.55f;
			m_sigma = 1.0f;
		}

		m_windowsize = 2 * m_r + 1;
		m_patchsize = 2 * m_f + 1;
		m_patch2 = static_cast<Float>(m_patchsize * m_patchsize);
		m_sigma2 = m_sigma * m_sigma * m_patch2; // variance (sigma squared)
		m_h2 = m_k * m_k * m_sigma2; // filter parameter squared and normalized with patch size
		m_isinitialized = false;

		if (dumpparams)
			dump();
	}
	NLMeansDenoiser(int r = 7, int f = 3, Float k = 0.45f, Float sigma = 1.0f, bool dumpm = true, bool dumpparams = true) : m_r(r), m_f(f), m_k(k),
				 m_sigma(sigma), m_dumpmaps(dumpm) {

		m_windowsize = 2 * m_r + 1;
		m_patchsize = 2 * m_f + 1;
		m_patch2 = static_cast<Float>(m_patchsize * m_patchsize);
		m_sigma2 = m_sigma * m_sigma * m_patch2; // variance (sigma squared)
		m_h2 = m_k * m_k * m_sigma2; // filter parameter squared and normalized with patch size
		m_isinitialized = false;

		if (dumpparams)
			dump();
	}

	void inline dump() {
		LOG(ECustom, "Denoiser: Non-Local Means Denoiser");
		LOG(ECustom, "Window Radius(r) = %d", m_r);
		LOG(ECustom, "Window size = %d", m_windowsize);
		LOG(ECustom, "Patch Radius(f) = %d", m_f);
		LOG(ECustom, "Patch size = %d", m_patchsize);
		LOG(ECustom, "Filter parameter(k) = %f", m_k);
		LOG(ECustom, "Standard deviation (sigma) = %f", m_sigma);
	}

	DenoiserOutput<O>* denoise(DenoiserInput<I> *in, bool patchbased = true, bool progressbar = true) {
		LOG(EInfo, "Denoising started");
		if (patchbased)
			return patchBasedDenoise<I, O>(in, progressbar);
		else
			return pixelBasedDenoise(in, progressbar);
	}


	DenoiserOutput<O>* pixelBasedDenoise(DenoiserInput<I> *in, bool progressbar = true) {
		bool check = false;
		int inputcount = in->getImageBlocks().size();
		assert(inputcount > 0);
		ImageBlockF* denoised = convert<I,Float>(in->getImageBlocks()[0]);
		denoised->clear();

		if (!m_isinitialized)
			check = initialize(in);
		assert(check);

		std::string basefilename = in->getBaseFilename();
		BitmapF *output = denoised->getBitmap();
		Float *outdata = output->getData();
		int outchannelcount = output->getChannelCount();

		// progress reporting!
		Float progress = 0.f, inc = 100.f / static_cast<Float>(m_size(0)*m_size(1));
		int channelcountA = m_imageblockA->getChannelCount();
		std::clock_t start = std::clock();
			// for each pixel compute weighted average from pixels in window
		for (int y = 0; y < m_size(1); ++y) {// iterate through each pixel in imageblock ignoring boundaries
	#pragma omp parallel for
			{
				for (int x = 0; x < m_size(0); ++x) {
					// pixel being filtered
					const Point2i p(x, y);
					Float weightSumA = 0.f;
					Float *filteredValueA = new Float[channelcountA];
					for (auto i = 0; i < channelcountA; ++i)
						filteredValueA[i] = 0.f;

					int width = computeSquareWindowWidth(x, y, m_size, m_r);
					const Point2i windowmin(x - width, y - width), windowmax(x + width, y + width);

					const Float *baseimageblockAdata = m_imageblockA->getBitmap()->getData();
					int offsetA = (windowmin(0) + m_size(0) * windowmin(1)) * channelcountA;
					int xincrementA = m_size(0) * channelcountA;

					// trivial case in case of pixels on edges
					if (width == 0) {
						for (int k = 0; k < outchannelcount; ++k)
							*outdata++ = *(baseimageblockAdata + offsetA++);
						continue;
					}

					Float maxweight = 0.f; // max weight of neighbor pixlels in window to be set for denoisee pixel
					// loop through each pixel in window and compute filter weight on equivalent pixel in other imageblock
					// and apply it to pixel in this imageblock - cross filtering
					for (int j = windowmin(1); j <= windowmax(1); ++j) {
						const Float *imageblockAdata = baseimageblockAdata + offsetA;
						for (int i = windowmin(0); i <= windowmax(0); ++i) {
							const Point2i q(i, j);
							Float weightA = 0.f;
							if (p(0) != i && p(1) != j) {
								weightA = computeNeighborWeightContribution(m_imageblockA->getBitmap(), p, q);
								if (weightA > maxweight) {
									maxweight = weightA;
								}
							}
							weightSumA += weightA;
							for (int k = 0; k < channelcountA; ++k)
								filteredValueA[k] += *imageblockAdata++ * weightA;
						}
						offsetA += xincrementA;
					}
					// account for weight contribution manually set for denoisee pixel (as the max of other weights)
					const Float *imageblockAdata = baseimageblockAdata + (x + m_size(0) * y) * channelcountA;
					for (int k = 0; k < channelcountA; ++k)
						filteredValueA[k] += *imageblockAdata++ * maxweight;
					weightSumA += maxweight;

					for (int k = 0; k < outchannelcount; ++k)
						*outdata++ = filteredValueA[k] / weightSumA;
					if (progressbar) {
						progress += inc;
						std::cout << progress << " % denoised \r";
					}
				}
			}
		}

		TImageBlock<O> *finaldenoised = convert<Float,O>(denoised);
		if (finaldenoised == NULL)
			return NULL;
		DenoiserOutput<O> *finaloutput = new DenoiserOutput<O>(finaldenoised);
		finaloutput->m_denoiseduration = (std::clock() - start) / CLOCKS_PER_SEC;
		LOG(EInfo, "Denoising finished. Time taken = %d", finaloutput->m_denoiseduration);
		return finaloutput;
	}

	template<typename I, typename O> DenoiserOutput<O>* patchBasedDenoise(DenoiserInput<I> *in, bool progressbar = true) {
		bool check = false;
		int inputcount = in->getImageBlocks().size();
		assert(inputcount > 0);
				if (!m_isinitialized)
			check = initialize(in);
		assert(check);
		ImageBlockF *denoised = convert<I, Float>(in->getImageBlocks()[0]);
		LOG(EInfo, "Data range of input after conversion to floating point representation = %f", denoised->getBitmap()->getDataRange());
		denoised->clear();

		std::string basefilename = in->getBaseFilename();
		BitmapF *output = denoised->getBitmap();
		// use spp bitmap to store denoised values per pixel
		BitmapF *valuesPerPixel = denoised->getSppBitmap();
		Float *outdata = output->getData();
		Float *sppdata = valuesPerPixel->getData();
		int outchannelcount = output->getChannelCount();
		const Float *bitmapData = m_imageblockA->getBitmap()->getData();

		// progress reporting!
		Float progress = 0.f, inc = 100.f / static_cast<Float>(m_size(0)*m_size(1));
		int channelcountA = m_imageblockA->getChannelCount();
		std::clock_t start = std::clock();

		// for each pixel compute weighted average from pixels in window
		for (int y = 0; y < m_size(1); ++y) {// iterate through each pixel in imageblock ignoring boundaries
//#pragma omp parallel for
			{
				// patch to be denoised
				Patch<Float> *denoisedPatch = new Patch<Float>(m_patchsize, m_patchsize, channelcountA);
				Float *denoisedData = denoisedPatch->getData();

				for (int x = 0; x < m_size(0); ++x) {
					// pixel being filtered
					const Point2i p(x, y);
					Float weightSumA = 0.f, maxWeightA = 0.f;

					// clear denoised patch
					denoisedPatch->clear(0.f);

					// reduce the size of the comparison window if we are near the boundary
					int m_f0 = min(m_f, min(m_size(0) - 1 - x, min(m_size(1) - 1 - y, min(x, y))));

					// research zone depending on the boundary and the size of the window
					int imin = max(x - m_r, m_f0);
					int jmin = max(y - m_r, m_f0);

					int imax = min(x + m_r, m_size(0) - 1 - m_f0);
					int jmax = min(y + m_r, m_size(1) - 1 - m_f0);

					for (int j = jmin; j <= jmax; ++j)
						for (int i = imin; i <= imax; ++i)
							if (i != x && j != y) {
								const Point2i q(i, j);

								Float weightA = computeNeighborWeightContribution(m_imageblockA->getBitmap(), p, q, m_f0);

								if (weightA > maxWeightA)
									maxWeightA = weightA;
								weightSumA += weightA;

								for (int is = -m_f0; is <= m_f0; ++is) {
									int aiindex = ((m_f + is) * m_patchsize + m_f) * channelcountA;
									int ail = ((j + is)*m_size(0) + i) * channelcountA;

									for (int ir = -m_f0 * channelcountA; ir <= m_f0 * channelcountA; ++ir) {
										int iindex = aiindex + ir;
										int il = ail + ir;
										denoisedData[iindex] += weightA * bitmapData[il];
									}
								}
							}

					// current patch with maxWeight
					for (int is = -m_f0; is <= m_f0; ++is) {
						int aiindex = ((m_f + is) * m_patchsize + m_f) * channelcountA;
						int ail = ((y + is)*m_size(0) + x) * channelcountA;

						for (int ir = -m_f0 * channelcountA; ir <= m_f0 * channelcountA; ++ir) {
							int iindex = aiindex + ir;
							int il = ail + ir;
							denoisedData[iindex] += maxWeightA * bitmapData[il];
						}
					}
					weightSumA += maxWeightA;

					//TODO : FIX BUG CAUSING 1 PIXEL VERTICAL BAR ON THE RIGHT END!

					// normalize average value when fTotalweight is not near zero
					if (weightSumA > FLOAT_EPSILON) {
						Float invWeightSumA = 1.f / weightSumA;
						for (int is = -m_f0; is <= m_f0; ++is) {
							int aiindex = ((m_f + is) * m_patchsize + m_f) * channelcountA;
							int ail = ((y + is)*m_size(0) + x) * channelcountA;

							for (int ir = -m_f0 * channelcountA; ir <= m_f0 * channelcountA; ++ir) {
								int iindex = aiindex + ir;
								int il = ail + ir;
								outdata[il] += denoisedData[iindex] * invWeightSumA;
								sppdata[il] += 1;
							}
						}

					}

					if (progressbar) {
						progress += inc;
						std::cout << progress << " % denoised \r";
					}
				}
				delete denoisedPatch;
			}
		}

		LOG(EInfo, "Data range of denoised output before spp normalization = %f", output->getDataRange());
		// normalize output according to no. of denoised values used per pixel
		if (!normalizeBitmap(output, valuesPerPixel, output)) {
			LOG(EError, "Error : Normalization of denoised bitmap failed!");
		}
		LOG(EInfo, "Data range of denoised output after spp normalization = %f", output->getDataRange());
		TImageBlock<O> *finaldenoised = convert<Float,O>(denoised);
		if (finaldenoised == NULL)
			return NULL;
		LOG(EInfo, "Data range of denoised output after conversion to output representation = %f", finaldenoised->getBitmap()->getDataRange());
		DenoiserOutput<O> *finaloutput = new DenoiserOutput<O>(finaldenoised);
		finaloutput->m_denoiseduration = (std::clock() - start) / CLOCKS_PER_SEC;
		return finaloutput;
	}


protected:
	// check validity of input, allocate space for temporary bitmaps and adjust parameters
	template<typename I> bool initialize(DenoiserInput<I> *in) {
		int inputcount = in->getImageBlocks().size();
		if (inputcount > 0) {
			Float minValue = 0.f, maxValue = 0.f;
			m_imageblockA = convert<I, Float>(in->getImageBlocks()[0]);
			Float scale = static_cast<Float>(m_imageblockA->getBitmap()->getDataRange(minValue, maxValue));
			LOG(ECustom, "Original scale of data = %f; MIN_VALUE = %f, MAX_VALUE = %f", scale, minValue, maxValue); 
			scale /= 255.f;
			m_size = m_imageblockA->getSize();
			m_sigma *= scale;
			m_sigma2 = m_sigma * m_sigma * m_imageblockA->getChannelCount();
			m_h2 *= m_imageblockA->getChannelCount() * scale * scale; // normalizing the parameters
			return true;
		}
		m_imageblockA = NULL;
		return false;
	}

	Vector2i getSize() const { return m_size; }

	bool isInitialized() const { return m_isinitialized; }

	// for pixel based weight computation
	template<typename I> Float computeNeighborWeightContribution(const TBitmap<I> *bitmap/*, const TBitmap<I> *varbitmap*/, const Point2i &p, const Point2i &q) const {
		Float patchDistance = 0.f;
		const Vector2i size = bitmap->getSize();

		int pwidth = computeSquareWindowWidth(p(0), p(1), size, m_f);
		int qwidth = computeSquareWindowWidth(q(0), q(1), size, m_f);
		int width = std::min(pwidth, qwidth);

		patchDistance = computePatchDistance(bitmap, p, q, width);
		return exp(-std::max(0.f, patchDistance - 2.f * m_sigma2)/m_h2);
	}

	// for patch based weight computation
	template<typename I> Float computeNeighborWeightContribution(const TBitmap<I> *bitmap/*, const TBitmap<I> *varbitmap*/, const Point2i &p, const Point2i &q, int width, bool print = false) const {
		Float patchDistance = computePatchDistance(bitmap, p, q, width);
		return exp(-std::max(0.f, patchDistance - 2.f * m_sigma2) / m_h2);
	}

	template<typename I> inline Float computePatchDistance(const TBitmap<I> *bitmap,// const TBitmap<I> *varbitmap,
		const Point2i &p, const Point2i &q, int width = -1) const {
		const Float *baseimageblockAdata = bitmap->getData();
		Vector2i size = bitmap->getSize();
		int stride = bitmap->getChannelCount();
		Float patchDistance = 0.f;
		if (width == -1)
			width = m_f;

		else if (width == 0) {
			int offsetP = (p(0) + size(0) * p(1)) * stride;
			int offsetQ = (q(0) + size(0) * q(1)) * stride;
			const Float *pvalue = baseimageblockAdata + offsetP;
			const Float *qvalue = baseimageblockAdata + offsetQ;
			return perPixelsquaredDistance(*pvalue, *qvalue);
		}

		const Point2i patchminP(p(0) - width, p(1) - width), patchmaxP(p(0) + width, p(1) + width);
		const Point2i patchminQ(q(0) - width, q(1) - width), patchmaxQ(q(0) + width, q(1) + width);

		int offsetP = (patchminP(0) + size(0) * patchminP(1)) * stride;
		int offsetQ = (patchminQ(0) + size(0) * patchminQ(1)) * stride;
		int xincrement = m_size(0) * stride;

		for (int j = patchminP(1); j <= patchmaxP(1); ++j) {
			const Float *pvalue = baseimageblockAdata + offsetP;
			const Float *qvalue = baseimageblockAdata + offsetQ;
			for (int i = patchminP(0); i <= patchmaxP(0); ++i) {
				for (int k = 0; k < stride; ++k)
					patchDistance += perPixelsquaredDistance(*pvalue++, *qvalue++);
			}
			offsetP += xincrement;
			offsetQ += xincrement;
		}

		return patchDistance;
	}

	/* Computes original nl means per-pixel squared distance assuming uniform variance */
	template<typename T> inline T perPixelsquaredDistance(const T &p, const T &q) const {
		return (p - q) * (p - q);
	}

protected:
	int m_r, m_f;
	Float m_k, m_sigma, m_sigma2, m_h2, m_patch2;
	int m_windowsize, m_patchsize; // m_windowsize = 2r + 1, m_patchsize = 2f + 1
	const ImageBlockF *m_imageblockA;
	Vector2i m_size;
	bool m_isinitialized, m_dumpmaps;
};

#endif