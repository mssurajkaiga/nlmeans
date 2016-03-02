#pragma once
#ifndef _MLT_NLMEANS_H_
#define _MLT_NLMEANS_H_

//mem leak checking
#define _CRTDBG_MAP_ALLOC

#include "modifiednlmeans.h"

/* A modifed version of Rouselle.et.al's modified NL-means denoiser (which works well only for MC renders)
to denoise multiple instance renders (works for MC and MCMC)
*/

template<typename I, typename O> class MLTNLMeansDenoiser : public ModifiedNLMeansDenoiser < I, O > {
public:
	MLTNLMeansDenoiser(int r = 7, int f = 3, Float k = 0.45, Float sigma = 1.0, int numinstances = 1, bool dumpm = true)
		: ModifiedNLMeansDenoiser<I, O>(r, f, k, sigma, dumpm, false), m_numinstances(numinstances) {
		m_isinitialized = false;
		dump();
	}

	void inline dump() {
		LOG(ECustom, "DENOISER: Parallel MLT NL Means Denoiser");
		LOG(ECustom, "Window Radius(r) = %d", m_r);
		LOG(ECustom, "Window size = %d", m_windowsize);
		LOG(ECustom, "Patch Radius(f) = %d", m_f);
		LOG(ECustom, "Patch size = %d", m_patchsize);
		LOG(ECustom, "Damping Factor (k) = %f", m_k);
		LOG(ECustom, "Variance cancellation (alpha) = %f", m_sigma);
		LOG(ECustom, "Number of render instances(samples for variance estimation) = %d", m_numinstances);
	}

	DenoiserOutput<O>* denoise(DenoiserInput<I> *in, bool patchbased = true, bool progressbar = true) {
		LOG(EInfo, "Denoising started");
		if (patchbased)
			return patchBasedDenoise<I, O>(in, progressbar);
		else
			return pixelBasedDenoise(in, progressbar);
	}

	// NOT IMPLEMENTED - disabled for now
	DenoiserOutput<O>* pixelBasedDenoise(DenoiserInput<I> *in, bool progressbar = true) {
		return NULL;
	}

	template<typename I, typename O> DenoiserOutput<O>* patchBasedDenoise(DenoiserInput<I> *in, bool progressbar = true) {
		bool check = false;
		ImageBlockF *denoised = convert<I, Float>(in->getImageBlocks()[0]);
		denoised->clear();

		if (!m_isinitialized)
			check = initialize(in);
		assert(check);

		std::string basefilename = in->getBaseFilename();
		dumpIntermediateBitmaps(basefilename + "_intermediate");
		BitmapF *output = denoised->getBitmap();

		// use spp bitmap to store denoised values per pixel
		BitmapF *valuesPerPixel = denoised->getSppBitmap();
		Float *outdata = output->getData();
		Float *sppdata = valuesPerPixel->getData();
		int outchannelcount = output->getChannelCount();
		const Float *bitmapData = m_inputbitmap->getData();

		// progress reporting!
		Float progress = 0.f, inc = 100.f / static_cast<Float>(m_size(0)*m_size(1));
		int channelcount = m_inputbitmap->getChannelCount();
		std::clock_t start = std::clock();

		// for each pixel compute weighted average from pixels in window
		for (int y = 0; y < m_size(1); ++y) {// iterate through each pixel in imageblock ignoring boundaries
			//#pragma omp parallel for
			{
				// patches to be denoised
				Patch<Float> *denoisedPatch = new Patch<Float>(m_patchsize, m_patchsize, channelcount);
				Float *denoisedData = denoisedPatch->getData();

				for (int x = 0; x < m_size(0); ++x) {
					// pixel being filtered
					const Point2i p(x, y);
					Float weightSum = 0.f, maxWeight = 0.f;
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

								Float weight = computeNeighborWeightContribution(m_inputbitmap, m_samplevariance, p, q, m_f0);

								if (weight > maxWeight)
									maxWeight = weight;
								weightSum += weight;
								
								for (int is = -m_f0; is <= m_f0; ++is) {
									int aiindex = ((m_f + is) * m_patchsize + m_f) * channelcount;
									int ail = ((j + is)*m_size(0) + i) * channelcount;

									for (int ir = -m_f0 * channelcount; ir <= m_f0 * channelcount; ++ir) {
										int iindex = aiindex + ir;
										int il = ail + ir;
										denoisedData[iindex] += weight * bitmapData[il];
									}
								}
							}

					// current patch with maxWeight
					for (int is = -m_f0; is <= m_f0; ++is) {
						int aiindex = ((m_f + is) * m_patchsize + m_f) * channelcount;
						int ail = ((y + is)*m_size(0) + x) * channelcount;

						for (int ir = -m_f0 * channelcount; ir <= m_f0 * channelcount; ++ir) {
							int iindex = aiindex + ir;
							int il = ail + ir;
							denoisedData[iindex] += maxWeight * bitmapData[il];
						}
					}
					weightSum += maxWeight;

					Float invWeightSum = 0.f;
					// normalize average value when fTotalweight is not near zero
					if (weightSum > FLOAT_EPSILON) {
						invWeightSum = 1.f / weightSum;
							for (int is = -m_f0; is <= m_f0; ++is) {
								int aiindex = ((m_f + is) * m_patchsize + m_f) * channelcount;
								int ail = ((y + is)*m_size(0) + x) * channelcount;

								for (int ir = -m_f0 * channelcount; ir <= m_f0 * channelcount; ++ir) {
									int iindex = aiindex + ir;
									int il = ail + ir;
									outdata[il] += denoisedData[iindex] * invWeightSum;
									sppdata[il] += 1; // a sample from each of the two patches
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


		// normalize output according to no. of denoised values used per pixel
		if (!normalizeBitmap<Float, Float>(output, valuesPerPixel, output)) {
			LOG(EError, "Normalization of denoised bitmap failed!\n");
		}

		TImageBlock<O> *finaldenoised = convert<Float, O>(denoised);
		if (finaldenoised == NULL)
			return NULL;
		DenoiserOutput<O> *finaloutput = new DenoiserOutput<O>(finaldenoised);
		finaloutput->m_denoiseduration = (std::clock() - start) / CLOCKS_PER_SEC;
		LOG(EInfo, "Denoising finished. Time taken = %d", finaloutput->m_denoiseduration);
		return finaloutput;
	}


protected:
	// check validity of input, allocate space for temporary bitmaps and generate pre-denoising bitmaps
	template<typename I> bool initialize(DenoiserInput<I> *in) {
		std::vector<TImageBlock<I>*> imageblocks = in->getImageBlocks();
		std::vector<BitmapF*> convertedbitmaps;
		int inputcount = imageblocks.size();

		m_size = imageblocks[0]->getSize();
		int channelcount = imageblocks[0]->getChannelCount();
		m_inputbitmap = new BitmapF(imageblocks[0]->getBitmap());
		m_inputbitmap->clear();
		Float scale = 1.f / static_cast<Float>(inputcount); // normalization factor while averaging input bitmaps
		//check if all input bitmaps have same sizes and average the bitmaps to generate final bitmap to be denoised
		for (auto i = imageblocks.begin(); i != imageblocks.end(); ++i) {
			Assert((*i)->getSize() == m_size);
			Assert((*i)->getChannelCount() == channelcount);
			BitmapF *converted = convert<I, Float>((*i)->getBitmap());
			convertedbitmaps.push_back(converted);
			m_inputbitmap->accumulate(converted, scale);
		}

		Float minValue = 0.f, maxValue = 0.f;
		scale = m_inputbitmap->getDataRange(minValue, maxValue) / 255.f;
		
		LOG(EInfo, "Normalization scale of data = %f", scale);
		m_sigma *= scale;
		m_sigma2 = m_sigma * m_sigma * channelcount;
		m_h2 *= channelcount * scale * scale; // normalizing the parameters

		m_samplevariance = new BitmapF(m_inputbitmap);
		m_samplevariance->clear();
		Assert(updateSampleVariance(convertedbitmaps, m_samplevariance, m_inputbitmap), "updateSampleVariance failed!");

		return true;
	}

	Vector2i getSize() const { return m_size; }

	bool isInitialized() const { return m_isinitialized; }

	template<typename I> bool updateSampleVariance(std::vector<TBitmap<I>*> inputbitmaps, TBitmap<I> *samplevariance, TBitmap<I> *meanbitmap = NULL) {
		int spp = inputbitmaps.size(); // each bitmap essentially represents a sample per pixel for variance estimation
		std::vector<const I*> inputdatas;
		if (meanbitmap == NULL) {
			LOG(EInfo, "Generating mean bitmap for sample variance computation");
			meanbitmap = new TBitmap<I>(inputbitmaps[0]);
			meanbitmap->clear();
			Float scale = 1.f / static_cast<Float>(spp);
			for (auto i = inputbitmaps.begin(); i != inputbitmaps.end(); ++i) {
				meanbitmap->accumulate(*i, scale);
				inputdatas.push_back((*i)->getData());
			}
		}
		else {
			for (auto i = inputbitmaps.begin(); i != inputbitmaps.end(); ++i) {
				inputdatas.push_back((*i)->getData());
			}
		}
			

		const int &channels = inputbitmaps[0]->getChannelCount();
		const Vector2i &size = inputbitmaps[0]->getSize();
		if (samplevariance->getSize() != size || samplevariance->getChannelCount() != channels)
			return false;

		const I *destmean = meanbitmap->getData();
		I *destsv = samplevariance->getData();

		if (spp > 1) {
			for (int i = 0; i < size(0) * size(1); ++i) {
				for (int k = 0; k < channels; ++k) {
					I mean = *destmean;
					I sum = static_cast<I>(0);
					for (int l = 0; l < spp; ++l) {
						const I* samplevalue = inputdatas[l];
						sum += (*samplevalue - mean) * (*samplevalue - mean);
						++inputdatas[l];
					}
					sum /= static_cast<Float>(spp - 1);
					*destsv = sum;
					++destsv;
					++destmean;
				}
			}
		}
		else
			for (int i = 0; i < size(0) * size(1) * channels; ++i, ++destsv)
				*destsv = 0.f; // no sufficient samples in this pixel to estimate variance

		return true;
	}

	// for patch based weight computation
	template<typename I> Float computeNeighborWeightContribution(const TBitmap<I> *bitmap, const TBitmap<I> *varbitmap,
		const Point2i &p, const Point2i &q, int width) const {
		Float patchDistance = computePatchDistance(bitmap, varbitmap, p, q, width);
		return exp(-std::max(0.f, patchDistance));
		//return exp(-std::max(0.f, patchDistance - 2.f * m_sigma2) / m_h2);
		//return exp(-patchDistance/m_h2);
	}

	template<typename I> inline Float computePatchDistance(const TBitmap<I> *bitmap, const TBitmap<I> *varbitmap,
		const Point2i &p, const Point2i &q, int width = -1) const {
		const Float *baseimageblockAdata = bitmap->getData();
		const Float *varbaseimageblockAdata = varbitmap->getData();
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
			const Float *varpvalue = varbaseimageblockAdata + offsetP;
			const Float *varqvalue = varbaseimageblockAdata + offsetQ;
			return perPixelsquaredDistance(*pvalue, *qvalue, *varpvalue, *varqvalue);
		}

		const Point2i patchminP(p(0) - width, p(1) - width), patchmaxP(p(0) + width, p(1) + width);
		const Point2i patchminQ(q(0) - width, q(1) - width), patchmaxQ(q(0) + width, q(1) + width);

		int offsetP = (patchminP(0) + size(0) * patchminP(1)) * stride;
		int offsetQ = (patchminQ(0) + size(0) * patchminQ(1)) * stride;
		int xincrement = m_size(0) * stride;

		for (int j = patchminP(1); j <= patchmaxP(1); ++j) {
			const Float *pvalue = baseimageblockAdata + offsetP;
			const Float *qvalue = baseimageblockAdata + offsetQ;
			const Float *varpvalue = varbaseimageblockAdata + offsetP;
			const Float *varqvalue = varbaseimageblockAdata + offsetQ;
			for (int i = patchminP(0); i <= patchmaxP(0); ++i) {
				for (int k = 0; k < stride; ++k)
					patchDistance += perPixelsquaredDistance(*pvalue++, *qvalue++, *varpvalue++, *varqvalue++);
			}
			offsetP += xincrement;
			offsetQ += xincrement;
		}

		return patchDistance;
	}

	/* Computes modified nl means per-pixel squared distance */
	template<typename T> inline T perPixelsquaredDistance(const T &p, const T &q,
		const T &varp, const T &varq) const {
		//return (p - q) * (p - q);
		return ((p - q) * (p - q) - m_sigma * (varp + min(varp, varq))) /
			(static_cast<T>(FLT_EPSILON)+m_k * m_k * (varp + varq));
	}

	void dumpIntermediateBitmaps(std::string basefilename = "") {
		dumpMap(m_samplevariance, basefilename + "_sample_variance", EHDR);
		dumpMap(m_inputbitmap, basefilename + "_mean", EHDR);
	}

protected:
	int m_numinstances; //no. of parallel instances of MLT used for denoising
	TBitmap<I> *m_inputbitmap, *m_samplevariance;
};

#endif