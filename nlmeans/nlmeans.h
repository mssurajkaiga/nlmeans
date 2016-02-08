#ifndef _NLMEANS_H_
#define _NLMEANS_H_

#include "core.h"
#include "util.h"

#include <omp.h>
#include <algorithm>
#include <vector>
#include <assert.h>
#include <iostream>
#include <ctime>

template<typename T> class NLMeansDenoiser;

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
	friend class NLMeansDenoiser<T>;
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
	DenoiserOutput(TImageBlock<T> *denoisedoutput = NULL) : m_denoisedimage(denoisedoutput) {}
	DenoiserOutput(DenoiserInput<T> *denoiserinput) {
		m_denoisedimage = new TImageBlock<T>(denoiserinput->getImageBlocks()[0]);
		m_denoisedimage->clear();
	}
	inline TImageBlock<T>* getDenoisedImage() { return m_denoisedimage; }
	inline long getDenoiseDuration() { return m_denoiseduration; }
protected:
	friend class NLMeansDenoiser<T>;
	virtual ~DenoiserOutput() {}
protected:
	TImageBlock<T>* m_denoisedimage;
	long m_denoiseduration; // time in milliseconds to denoise the input
};

template<typename T> class NLMeansDenoiser {
public:
	NLMeansDenoiser(int r = 7, int f = 3, float k = 0.45, float sigma= 1.0, int vr = 1, bool dumpm=true) : m_r(r), m_f(f), m_k(k),
				 m_sigma(sigma), m_dumpmaps(dumpm) {

		m_windowsize = 2 * m_r + 1;
		m_patchsize = 2 * m_f + 1;
		m_sigma2 = m_sigma * m_sigma; // variance (sigma squared)
		m_h2 = m_k * m_k * m_sigma2; // filteri parameter squared
		m_isinitialized = false;

		dump();
	}

	void inline dump() {
		Log(ECustom, "Modified NL Means Denoiser Parameters");
		Log(ECustom, "Window Radius(r) = %d", m_r);
		Log(ECustom, "Window size = %d", m_windowsize);
		Log(ECustom, "Patch Radius(f) = %d", m_f);
		Log(ECustom, "Patch size = %d", m_patchsize);
		Log(ECustom, "Damping Factor(k) = %f", m_k);
		Log(ECustom, "Standard deviation (sigma) = %f", m_sigma);
		//Log(ECustom, "NOTE: DENOISER USES CROSS FILTERING ONLY WITH SAMPLE BITMAP AND USES SCALED MEAN SAMPLE VARIANCE BITMAP");
	}

	DenoiserOutput<T>* denoise(DenoiserInput<T> *in, bool progressbar = true) {
		bool check = false;
		int inputcount = in->getImageBlocks().size();
		assert(inputcount > 0);
		TImageBlock<T>* denoised = new TImageBlock<T>(in->getImageBlocks()[0]);
		denoised->clear();
		DenoiserOutput<T>* out = new DenoiserOutput<T>(denoised);

		if (!m_isinitialized)
			check = initialize(in);
		assert(check);

		std::string basefilename = in->getBaseFilename();
		TBitmap<T> *output = out->getDenoisedImage()->getBitmap();
		T *outdata = output->getData();
		int outchannelcount = output->getChannelCount();

		// progress reporting!
		Float progress = 0.f, inc = 100.f / static_cast<Float>(m_size(0)*m_size(1));
		int channelcountA = m_imageblockA->getChannelCount();
		std::clock_t start = std::clock();
			// for each pixel compute weighted average from pixels in window
		for (size_t y = 0; y < m_size(1); ++y) {// iterate through each pixel in imageblock ignoring boundaries
	#pragma omp parallel for
			{
				for (size_t x = 0; x < m_size(0); ++x) {
					// pixel being filtered
					const Point2i p(x, y);
					Float weightSumA = 0.f;
					Float *filteredValueA = new Float[channelcountA];
					for (auto i = 0; i < channelcountA; ++i)
						filteredValueA[i] = 0.f;

					//const Point2i tempmin(x - std::max((int)x - m_r, 0), y - (std::max((int)y - m_r, 0)));
					//const Point2i tempmax(std::min((int)x + m_r, m_size(0) - 1) - x, std::min((int)y + m_r, m_size(1) - 1) - y);
					//int width = std::min(std::min(tempmin(0), tempmin(1)), std::min(tempmax(0), tempmax(1)));
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
			out->m_denoiseduration = (std::clock() - start)/CLOCKS_PER_SEC;
		return out;
	}

	DenoiserOutput<T>* patchBasedDenoise((DenoiserInput<T> *in, bool progressbar = true) {
		bool check = false;
		int inputcount = in->getImageBlocks().size();
		assert(inputcount > 0);
		TImageBlock<T>* denoised = new TImageBlock<T>(in->getImageBlocks()[0]);
		denoised->clear();
		DenoiserOutput<T>* out = new DenoiserOutput<T>(denoised);

		if (!m_isinitialized)
			check = initialize(in);
		assert(check);

		std::string basefilename = in->getBaseFilename();
		TBitmap<T> *output = out->getDenoisedImage()->getBitmap();
		T *outdata = output->getData();
		int outchannelcount = output->getChannelCount();

		// progress reporting!
		Float progress = 0.f, inc = 100.f / static_cast<Float>(m_size(0)*m_size(1));
		int channelcountA = m_imageblockA->getChannelCount();
		std::clock_t start = std::clock();
		// for each pixel compute weighted average from pixels in window
		for (size_t y = 0; y < m_size(1); ++y) {// iterate through each pixel in imageblock ignoring boundaries
#pragma omp parallel for
			{
				// patch to be denoised
				Patch<T> *denoisedPatch = new Patch<T>(m_f, m_f, channelcountA);
				const T *bitmapData = m_imageblockA->getBitmap()->getData();
				T *denoisedData = denoisedPatch->getData();

				for (size_t x = 0; x < m_size(0); ++x) {
					// pixel being filtered
					const Point2i p(x, y);
					Float weightSumA = 0.f, maxWeightA = 0.f;
					Float *filteredValueA = new Float[channelcountA];
					for (int i = 0; i < channelcountA; ++i)
						filteredValueA[i] = 0.f;

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
							if (i != x || j != y) {
								const Point2i q(i, j);

								//float fDif = fiL2FloatDist(fpI, fpI, x, y, i, j, m_f0, iChannels, m_size(0), m_size(0));

								//// dif^2 - 2 * fSigma^2 * N      dif is not normalized
								//fDif = MAX(fDif - 2.0f * (float)icwl *  fSigma2, 0.0f);
								//fDif = fDif / fH2;

								//float fWeight = wxSLUT(fDif, fpLut);

								Float weightA = computeNeighborWeightContribution(m_imageblockA->getBitmap(), p, q);

								if (weightA > maxWeightA)
									maxWeight = weightA;
								weightSumA += weightA;

								for (int is = -m_f0; is <= m_f0; ++is) {
									int aiindex = ((m_f + is) * m_patchsize + m_f ) * channelcountA;
									int ail = ((j + is)*m_size(0) + i) * channelcountA;

									for (int ir = -m_f0 * channelcountA; ir <= m_f0 * channelcountA; ++ir) {
										int iindex = aiindex + ir;
										int il = ail + ir;
										denoisedData[iindex] += weightA * bitmapData[il];W
									}
								}
							}


					// current patch with fMaxWeight
					for (int is = -m_f0; is <= m_f0; ++is) {
						int aiindex = ((m_f + is) * m_patchsize + m_f) * channelcountA;
						int ail = ((y + is)*m_size(0) + x) * channelcountA;

						for (int ir = -m_f0 * channelcountA; ir <= m_f0 * channelcountA; ++ir) {
							int iindex = aiindex + ir;
							int il = ail + ir;
							denoisedData[iindex] += maxWeightA * bitmapData[il]; W
						}
					}
					weighSumA += maxWeightA;

					// normalize average value when fTotalweight is not near zero
					if (weightSumA > FLOAT_EPSILON) {

						for (int is = -m_f0; is <= m_f0; ++is) {
							int aiindex = ((m_f + is) * m_patchsize + m_f) * channelcountA;
							int ail = ((y + is)*m_size(0) + x) * channelcountA;

							for (int ir = -m_f0 * channelcountA; ir <= m_f0 * channelcountA; ++ir) {
								int iindex = aiindex + ir;
								int il = ail + ir;
								denoisedData[iindex] += maxWeightA * bitmapData[il]; W
							}
						}



						for (int is = -m_f0; is <= m_f0; is++) {
							int aiindex = (m_f + is) * m_patchsize + m_f;
							int ail = (y + is)*m_size(0) + x;

							for (int ir = -m_f0; ir <= m_f0; ir++) {
								int iindex = aiindex + ir;
								int il = ail + ir;

								fpCount[il]++;

								for (int ii = 0; ii < iChannels; ii++) {
									fpO[ii][il] += fpODenoised[ii][iindex] / fTotalWeight;

								}

							}
						}


					}


			}
		}
	}


private:
	// check validity of input, allocate space for temporary bitmaps and generate pre-denoising bitmaps
	bool initialize(DenoiserInput<T> *in) {
		int inputcount = in->getImageBlocks().size();
		if (inputcount > 0) {
			m_imageblockA = in->getImageBlocks()[0];
			m_size = m_imageblockA->getSize();
			m_sigma2 *= m_imageblockA->getChannelCount(); // normalizing the parameters
			m_h2 *= m_imageblockA->getChannelCount(); // normalizing the parameters
			return true;
		}
		m_imageblockA = NULL;
		return false;
	}

	Vector2i getSize() const { return m_size; }

	bool isInitialized() const { return m_isinitialized; }

	

	//loop through each pixel in patch and compute neighbor pixel contribution
	template<typename T> Float computeNeighborWeightContribution(const TBitmap<T> *bitmap/*, const TBitmap<T> *varbitmap*/, const Point2i &p, const Point2i &q) const {
		Float patchDistance = 0.f;
		const Vector2i size = bitmap->getSize();

		int pwidth = computeSquareWindowWidth(p(0), p(1), size, m_f);
		int qwidth = computeSquareWindowWidth(q(0), q(1), size, m_f);
		int width = std::min(pwidth, qwidth);

		patchDistance = computePatchDistance(bitmap, p, q, width);
		
		//return exp(-std::max(0.f, patchDistance));
		//return exp(-std::max(0.f, patchDistance - 2.f * m_sigma2)/m_h2);
		return exp(-patchDistance/m_h2);
	}

	template<typename T> inline Float computePatchDistance(const TBitmap<T> *bitmap,// const TBitmap<T> *varbitmap,
		const Point2i &p, const Point2i &q, int width = -1) const {
		const Float *baseimageblockAdata = bitmap->getData();
		Vector2i size = bitmap->getSize();
		int stride = bitmap->getChannelCount();
		Float patchDistance = 0.f;
		if (width == -1)
			width = m_f;
		int patchwidth = 2 * width + 1;

		if (width == 0) {
			int offsetP = (p(0) + size(0) * p(1)) * stride;
			int offsetQ = (q(0) + size(0) * q(1)) * stride;
			const Float *pvalue = baseimageblockAdata + offsetP;
			const Float *qvalue = baseimageblockAdata + offsetQ;
			return perPixelsquaredDistance(*pvalue, *qvalue) / static_cast<Float>(stride);
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

		patchDistance /= static_cast<Float>(stride * patchwidth * patchwidth);
		return patchDistance;
	}

	/* Computes original nl means per-pixel squared distance assuming uniform variance */
	template<typename T> inline T perPixelsquaredDistance(const T &p, const T &q) const {
		return (p - q) * (p - q);
	}

private:
	int m_r, m_f;
	Float m_k, m_sigma, m_sigma2, m_h2;
	int m_windowsize, m_patchsize; // m_windowsize = 2r + 1, m_patchsize = 2f + 1
	const TImageBlock<T> *m_imageblockA, *m_imageblockB;
	Vector2i m_size;
	bool m_isinitialized, m_dumpmaps;
};

#endif