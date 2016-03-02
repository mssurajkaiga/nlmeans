#pragma once
#ifndef _MODIFIED_NLMEANS_H_
#define _MODIFIED_NLMEANS_H_

//mem leak checking
#define _CRTDBG_MAP_ALLOC

#include "nlmeans.h"

/* Modified NL Means Denoiser as per Rouselle et. al's paper
some parameters with same name are different from original nl means implementation in nlmeans.h
k = damping factor (similar in pratice to original k)
sigma = strength of variance cancellation and not variance value
*/

template<typename I, typename O> class ModifiedNLMeansDenoiser: public NLMeansDenoiser<I, O> {
public:
	ModifiedNLMeansDenoiser(int r = 7, int f = 3, Float k = 0.45, Float sigma = 1.0, bool dumpm = true, bool dumpparams = true)
		: NLMeansDenoiser<I, O>(r, f, k, sigma, dumpm, false) {

		m_vR = M_R;
		m_vF = M_F;
		m_vSigma = M_ALPHA * 10.f;
		m_vK = M_K;

		//rescale parameters since input values are according to the paper
		if (m_sigma < 10.f) {
			m_sigma *= 100.f;
			LOG(EInfo, "Rescaled variance cancellation (alpha) by 100.0");
		}
		/*else if (m_sigma >= 1.f && m_sigma < 10.f){
			m_sigma *= 10.f;
			LOG(EInfo, "Rescaled variance cancellation (alpha) by 10.0");
		}*/
		if (m_sigma >= 500.f ){
			m_sigma = 500.f;
			LOG(EInfo, "Clamped variance cancellation (alpha) to 500.0");
		}
		m_sigma2 = m_sigma * m_sigma *m_patch2; // variance (sigma squared)
		m_h2 = m_k * m_k * m_sigma2; // filter parameter squared and normalized with patch size
		m_isinitialized = false;

		if (dumpparams)
			dump();
	}

	void inline dump() {
		LOG(ECustom, "DENOISER: Modified NL Means Denoiser");
		LOG(ECustom, "Window Radius(r) = %d", m_r);
		LOG(ECustom, "Window size = %d", m_windowsize);
		LOG(ECustom, "Patch Radius(f) = %d", m_f);
		LOG(ECustom, "Patch size = %d", m_patchsize);
		LOG(ECustom, "Damping Factor (k) = %f", m_k);
		LOG(ECustom, "Variance cancellation (alpha) = %f", m_sigma);
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
		int inputcount = in->getImageBlocks().size();
		assert(inputcount == 2);
		in->getImageBlocks()[0]->getBitmap()->logToFile("original_input.txt");
		ImageBlockF *denoised = convert<I, Float>(in->getImageBlocks()[0]);
		denoised->clear();

		if (!m_isinitialized)
			check = initialize(in);
		assert(check);
		const BitmapF *converted = m_imageblockA->getBitmap();
		converted->logToFile("converted_input.txt");

		std::ofstream logFile;
		logFile.open("log.txt");

		std::string basefilename = in->getBaseFilename();
		dumpIntermediateBitmaps(basefilename+"_intermediate");
		BitmapF *output = denoised->getBitmap();
		// use spp bitmap to store denoised values per pixel
		BitmapF *valuesPerPixel = denoised->getSppBitmap();
		Float *outdata = output->getData();
		Float *sppdata = valuesPerPixel->getData();
		int outchannelcount = output->getChannelCount();
		const Float *bitmapAData = m_imageblockA->getBitmap()->getData();
		const Float *bitmapBData = m_imageblockB->getBitmap()->getData();

		// progress reporting!
		Float progress = 0.f, inc = 100.f / static_cast<Float>(m_size(0)*m_size(1));
		int channelcountA = m_imageblockA->getChannelCount();
		std::clock_t start = std::clock();

		// for each pixel compute weighted average from pixels in window
		for (int y = 0; y < m_size(1); ++y) {// iterate through each pixel in imageblock ignoring boundaries
			//#pragma omp parallel for
			{
				// patches to be denoised
				Patch<Float> *denoisedPatchA = new Patch<Float>(m_patchsize, m_patchsize, channelcountA);
				Float *denoisedDataA = denoisedPatchA->getData();
				Patch<Float> *denoisedPatchB = new Patch<Float>(m_patchsize, m_patchsize, channelcountA);
				Float *denoisedDataB = denoisedPatchB->getData();

				for (int x = 0; x < m_size(0); ++x) {
					// pixel being filtered
					const Point2i p(x, y);
					Float weightSumA = 0.f, maxWeightA = 0.f, weightSumB = 0.f, maxWeightB = 0.f;
					// clear denoised patch
					denoisedPatchA->clear(0.f);
					denoisedPatchB->clear(0.f);

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

								Float weightA = computeNeighborWeightContribution(m_imageblockB->getBitmap(), m_filteredbuffervariance, p, q, m_f0);
								Float weightB = computeNeighborWeightContribution(m_imageblockA->getBitmap(), m_filteredbuffervariance, p, q, m_f0);

								if (weightA > maxWeightA)
									maxWeightA = weightA;
								weightSumA += weightA;
								if (weightB > maxWeightB)
									maxWeightB = weightB;
								weightSumB += weightB;

								for (int is = -m_f0; is <= m_f0; ++is) {
									int aiindex = ((m_f + is) * m_patchsize + m_f) * channelcountA;
									int ail = ((j + is)*m_size(0) + i) * channelcountA;

									for (int ir = -m_f0 * channelcountA; ir <= m_f0 * channelcountA; ++ir) {
										int iindex = aiindex + ir;
										int il = ail + ir;
										denoisedDataA[iindex] += weightA * bitmapAData[il];
										denoisedDataB[iindex] += weightB * bitmapBData[il];
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
							denoisedDataA[iindex] += maxWeightA * bitmapAData[il];
							denoisedDataB[iindex] += maxWeightB * bitmapBData[il];
						}
					}
					weightSumA += maxWeightA;
					weightSumB += maxWeightB;

					Float invWeightSumA = 0.f, invWeightSumB = 0.f;
					// normalize average value when fTotalweight is not near zero
					if (weightSumA > FLOAT_EPSILON) {
						invWeightSumA = 1.f / weightSumA;
						if (weightSumB > FLOAT_EPSILON) {
							invWeightSumB = 1.f / weightSumB;
							for (int is = -m_f0; is <= m_f0; ++is) {
								int aiindex = ((m_f + is) * m_patchsize + m_f) * channelcountA;
								int ail = ((y + is)*m_size(0) + x) * channelcountA;

								for (int ir = -m_f0 * channelcountA; ir <= m_f0 * channelcountA; ++ir) {
									int iindex = aiindex + ir;
									int il = ail + ir;
									outdata[il] += denoisedDataA[iindex] * invWeightSumA
										+ denoisedDataB[iindex] * invWeightSumB;
									sppdata[il] += 2; // a sample from each of the two patches
								}
							}
						}
						else {
							for (int is = -m_f0; is <= m_f0; ++is) {
								int aiindex = ((m_f + is) * m_patchsize + m_f) * channelcountA;
								int ail = ((y + is)*m_size(0) + x) * channelcountA;

								for (int ir = -m_f0 * channelcountA; ir <= m_f0 * channelcountA; ++ir) {
									int iindex = aiindex + ir;
									int il = ail + ir;
									outdata[il] += denoisedDataA[iindex] * invWeightSumA;
									sppdata[il] += 1;
								}

							}
						}
					}
					else if (weightSumB > FLOAT_EPSILON) {
						invWeightSumB = 1.f / weightSumB;
						for (int is = -m_f0; is <= m_f0; ++is) {
							int aiindex = ((m_f + is) * m_patchsize + m_f) * channelcountA;
							int ail = ((y + is)*m_size(0) + x) * channelcountA;

							for (int ir = -m_f0 * channelcountA; ir <= m_f0 * channelcountA; ++ir) {
								int iindex = aiindex + ir;
								int il = ail + ir;
								outdata[il] += denoisedDataB[iindex] * invWeightSumB;
								sppdata[il] += 1;
							}
						}
					}

					if (progressbar) {
						progress += inc;
						std::cout << progress << " % denoised \r";
					}
				}
				delete denoisedPatchA;
				delete denoisedPatchB;
			}
		}


		// normalize output according to no. of denoised values used per pixel
		if (!normalizeBitmap<Float, Float>(output, valuesPerPixel, output)) {
			LOG(EError, "Normalization of denoised bitmap failed!\n");
		}

		logFile.close();
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
		int inputcount = in->getImageBlocks().size();
		if (inputcount != 2) {
			LOG(EError, "Need dual buffer samples!");
			return false;
		}

		Float minValue = 0.f, maxValue = 0.f;
		m_imageblockA = convert<I, Float>(in->getImageBlocks()[0]);
		Float scaleA = static_cast<Float>(m_imageblockA->getBitmap()->getDataRange(minValue, maxValue)) / 255.f;
		m_size = m_imageblockA->getSize();
		int channelcount = m_imageblockA->getChannelCount();
		m_imageblockB = convert<I, Float>(in->getImageBlocks()[1]);
		Float scaleB = static_cast<Float>(m_imageblockB->getBitmap()->getDataRange(minValue, maxValue)) / 255.f;
		Assert(m_size == m_imageblockB->getSize(), "size mismatch between input imageblocks");
		Assert(channelcount == m_imageblockB->getChannelCount(), "channel count mismatch between input imageblocks");
		Float scale = (scaleA + scaleB) / 2.f; // needs to be redone later!
		std::cout << "Normalization scale of data = " << scale << "\n";
		m_sigma *= scale;
		m_sigma2 = m_sigma * m_sigma * channelcount;
		m_h2 *= channelcount * scale * scale; // normalizing the parameters

		const BitmapF *var = m_imageblockA->getVarBitmap();
		m_samplevarianceA = new BitmapF(var->getSize(), var->getChannelCount());
		m_buffervariance = new BitmapF(var->getSize(), var->getChannelCount());
		m_meansamplevariance = new BitmapF(var->getSize(), var->getChannelCount());
		m_samplevariancevariance = new BitmapF(var->getSize(), var->getChannelCount());
		m_filteredbuffervariance = new BitmapF(var->getSize(), var->getChannelCount());
		Assert(updateSampleVariance(m_imageblockA, m_samplevarianceA), "updateSampleVariance failed for image block A!");

		m_samplevarianceB = new BitmapF(var->getSize(), var->getChannelCount());
		Assert(updateSampleVariance(m_imageblockB, m_samplevarianceB), "updateSampleVariance failed for image block A!");
		//Assert(updateBufferVariance(m_samplevarianceA, m_samplevarianceB, m_buffervariance), "updateBufferVariance failed!");
		Assert(updateBufferVariance(m_imageblockA->getBitmap(), m_imageblockB->getBitmap(), m_buffervariance), "updateBufferVariance failed!");
		Assert(updateSampleVarianceVariance(), "updateSampleVarianceVariance failed!");
		Assert(getWeightedMean(m_samplevarianceA, m_samplevarianceB, m_imageblockA->getSppBitmap(),
			m_imageblockB->getSppBitmap(), m_meansamplevariance), "getWeightedMean failed!");
		//delete m_filteredbuffervariance;
		//m_filteredbuffervariance = m_meansamplevariance;
		Assert(crossfilter(m_buffervariance, m_meansamplevariance, m_samplevariancevariance, m_filteredbuffervariance));
		return true;
	}

	Vector2i getSize() const { return m_size; }

	bool isInitialized() const { return m_isinitialized; }

	// denoise buffer variance estimate by cross filtering with weights from empirical sample variance estimate
	template<typename I> bool crossfilter(TBitmap<I> *buffervariance, TBitmap<I> *samplevariance, TBitmap<I> *samplevariancevariance, TBitmap<I> *filteredoutput) {
		const Vector2i size = buffervariance->getSize();
		int channelcount = buffervariance->getChannelCount();
		Assert(size == samplevariance->getSize() && channelcount == samplevariance->getChannelCount());
		Assert(size == samplevariancevariance->getSize() && channelcount == samplevariancevariance->getChannelCount());
		Assert(size == filteredoutput->getSize() && channelcount == filteredoutput->getChannelCount());

		//Create temporary bitmap for storing values per pixel
		TBitmap<I> *sppbitmap = new TBitmap<I>(buffervariance);
		sppbitmap->clear();

		// Temporarily swapping variables for the purpose of cross-filtering.. I swear, I'll return them to original values!
		int temp_r = m_r;
		int temp_f = m_f;
		Float temp_sigma = m_sigma;
		Float temp_k = m_k;
		Float temp_sigma2 = m_sigma2;
		Float temp_h2 = m_h2;
		int temp_patchsize = m_patchsize;
		int temp_patch2 = m_patch2;

		m_r = m_vR;
		m_f = m_vF;
		m_sigma = m_vSigma;
		m_k = m_vK;
		m_patchsize = 2 * m_f + 1;
		m_patch2 = m_patchsize * m_patchsize;
		m_sigma2 = m_sigma * m_sigma *m_patch2; // variance (sigma squared)
		m_h2 = m_k * m_k * m_sigma2; // filter parameter squared and normalized with patch size
		Float scale = buffervariance->getDataRange();
		m_sigma *= scale;
		m_sigma2 = m_sigma * m_sigma * channelcount;
		m_h2 *= channelcount * scale * scale; // normalizing the parameters

		const Float *bitmapData = buffervariance->getData();
		Float *outdata = filteredoutput->getData();
		Float *sppdata = sppbitmap->getData();

		// for each pixel compute weighted average from pixels in window
		for (int y = 0; y < size(1); ++y) {// iterate through each pixel in imageblock ignoring boundaries
			// create patch to be filtered
			Patch<Float> *denoisedPatch = new Patch<Float>(m_patchsize, m_patchsize, channelcount);
			Float *denoisedData = denoisedPatch->getData();
			for (int x = 0; x < size(0); ++x) {
				// pixel being filtered
				const Point2i p(x, y);
				Float weightSum = 0.f, maxWeight = 0.f;

				denoisedPatch->clear(0.f);
				//compute window size
				int m_f0 = min(m_f, min(size(0) - 1 - x, min(size(1) - 1 - y, min(x, y))));

				// research zone depending on the boundary and the size of the window
				int imin = max(x - m_r, m_f0);
				int jmin = max(y - m_r, m_f0);

				int imax = min(x + m_r, size(0) - 1 - m_f0);
				int jmax = min(y + m_r, size(1) - 1 - m_f0);

				// loop through each pixel in window and compute filter weight
				for (int j = jmin; j <= jmax; ++j)
					for (int i = imin; i <= imax; ++i)
						if (i != x && j != y) {
							const Point2i q(i, j);

							Float weight = computeNeighborWeightContribution(samplevariance, samplevariancevariance, p, q, m_f0);

							if (weight > maxWeight)
								maxWeight = weight;
							weightSum += weight;

							for (int is = -m_f0; is <= m_f0; ++is) {
								int aiindex = ((m_f + is) * m_patchsize + m_f) * channelcount;
								int ail = ((j + is)*size(0) + i) * channelcount;

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
					int ail = ((y + is)*size(0) + x) * channelcount;

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
						int ail = ((y + is)*size(0) + x) * channelcount;

						for (int ir = -m_f0 * channelcount; ir <= m_f0 * channelcount; ++ir) {
							int iindex = aiindex + ir;
							int il = ail + ir;
							outdata[il] += denoisedData[iindex] * invWeightSum;
							sppdata[il] += 1;
						}
					}
				}
			}
			delete denoisedPatch;
		}

		Assert(normalizeBitmap(filteredoutput, sppbitmap, filteredoutput), "Normalizing cross-filtered bitmap failed!");

		// As promised, I deliver the swap ;)
		m_r = temp_r;
		m_f = temp_f;
		m_sigma = temp_sigma;
		m_k = temp_k;
		m_sigma2 = temp_sigma2;
		m_h2 = temp_h2;
		m_patchsize = temp_patchsize;
		m_patch2 = temp_patch2;

		return true;
	}

	template<typename I> bool updateSampleVariance(const TImageBlock<I> *imageblock, TBitmap<I> *samplevariance) {
		const TBitmap<I> *varbitmap = imageblock->getVarBitmap(),
			*varsbitmap = imageblock->getVarsBitmap();
		const TBitmap<I> *sppbitmap = imageblock->getSppBitmap();

		if (varbitmap == NULL || varsbitmap == NULL)
			return false;

		const int &channels = varbitmap->getChannelCount();
		const Vector2i &size = varbitmap->getSize();
		if (samplevariance->getSize() != size)
			return false;

		const I *destvar = varbitmap->getData(),
			*destvars = varsbitmap->getData();

		if (sppbitmap == NULL) {
			for (int i = 0; i < size(0) * size(1) * channels; ++i, ++destsv) {
				*destsv = 0.f; // no sufficient samples in this pixel to estimate variance
			}
			return true;
		}

		const I *destspp = sppbitmap->getData();
		I *destsv = samplevariance->getData();

		for (int i = 0; i < size(0) * size(1); ++i) {
			if (*destspp > 1) {
				for (int k = 0; k < channels; ++k, ++destvar, ++destvars, ++destspp, ++destsv) {
					I mean = *destvar / static_cast<I>(*destspp);
					*destsv = (*destvars - *destvar * mean) / static_cast<I>(*destspp - 1);
					if (*destsv < 0.f) {
						//LOG(EError, "Sample variance cannot be negative!! - check variance map computation in renderer at position %d", i);
						*destsv = std::abs(*destsv);
					}
				}
			}
			else {
				for (int k = 0; k < channels; ++k, ++destvar, ++destvars, ++destspp, ++destsv) {
					*destsv = 0.f; // no sufficient samples in this pixel to estimate variance
				}
			}
		}

		return true;
	}

	template<typename I> bool getWeightedMean(const TBitmap<I> *input1, const TBitmap<I> *input2,
		const TBitmap<I> *weight1, const TBitmap<I> *weight2, TBitmap<I> *output) {

		const int &channels = input1->getChannelCount();
		const Vector2i &size = input1->getSize();
		if (input2->getSize() != size || weight1->getSize() != size || weight2->getSize() != size || output->getSize() != size)
			return false;
		if (input2->getChannelCount() != channels || output->getChannelCount() != channels)
			return false;

		const I *dest1 = input1->getData(),
			*dest2 = input2->getData();
		I *destO = output->getData();

		const I *destW1 = weight1->getData(),
			*destW2 = weight2->getData();

		for (int i = 0; i < size(0) * size(1); ++i) {
			for (int k = 0; k < channels; ++k, ++destW1, ++destW2, ++dest1, ++dest2, ++destO) {
				if (*destW1 == 0 && *destW2 == 0) {
					*destO == static_cast<I>(0);
					continue;
				}
				*destO = (*dest1 * *destW1 + *dest2 * *destW2) / static_cast<I>(*destW1 + *destW2);
			}
		}

		return true;
	}

	template<typename I> bool clampAndScale(TBitmap<I> *samplevariance, TBitmap<I> *buffervariance, I clampmin = 0.f, I clampmax = 10.f) {
		I *destsv = samplevariance->getData();
		const I *destbv = buffervariance->getData();
		I avgvarianceratio[3];
		for (int i = 0; i < 3; ++i)
			avgvarianceration[i] = 0.f;

		const int &channels = buffervariance->getChannelCount();
		const Vector2i &size = samplevariance->getSize();
		int pixelscounted = 0;

		for (int i = 0; i < size.x * size.y; ++i) {
			for (int k = 0; k < channels; ++k) {
				if (*destsv > FLT_EPSILON) {
					avgvarianceratio[k] += clamp((*destbv / *destsv), clampmin, clampmax);
					++pixelscounted;
				}
				++destsv;
				++destbv;
			}
		}
		avgvarianceratio /= static_cast<I>(pixelscounted);

		destsv = samplevariance->getData();
		//scale the sample variance bitmap with the average variance ratio
		for (int i = 0; i < size.x * size.y; ++i) {
			for (int k = 0; k < channels; ++k) {
				*destsv++ *= avgvarianceratio[k];
			}
		}

		LOG(ECustom, "Average Variance Ratio = %f %f %f", avgvarianceratio[0], avgvarianceratio[1],
			avgvarianceratio[2]);
		return true;
	}

	template<typename I> bool updateBufferVariance(const TBitmap<I> *inputA, const TBitmap<I> *inputB, TBitmap<I> *output) {
		return computeBitmapDifference<2,I>(inputA, inputB, output, 0.5);
	}

	bool updateSampleVarianceVariance() {
		return computeBitmapDifference<1,Float>(m_samplevarianceA, m_samplevarianceB, m_samplevariancevariance, 1.f, true);
	}

	// for pixel based weight computation - DISABLED FOR NOW
	template<typename I> Float computeNeighborWeightContribution(const TBitmap<I> *bitmap,
		const TBitmap<I> *varbitmap, const Point2i &p, const Point2i &q) const {
		//Float patchDistance = 0.f;
		//const Vector2i size = bitmap->getSize();

		//int pwidth = computeSquareWindowWidth(p(0), p(1), size, m_f);
		//int qwidth = computeSquareWindowWidth(q(0), q(1), size, m_f);
		//int width = std::min(pwidth, qwidth);

		//patchDistance = computePatchDistance(bitmap, p, q, width);

		////return exp(-std::max(0.f, patchDistance));
		//return exp(-std::max(0.f, patchDistance - 2.f * m_sigma2) / m_h2);
		////return exp(-patchDistance/m_h2);
		return 0.f;
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
			(static_cast<T>(FLT_EPSILON) + m_k * m_k * (varp + varq));
	}

	void dumpIntermediateBitmaps(std::string basefilename = "") {
		dumpMap(m_samplevarianceA, basefilename + "_sample_variance_A", EHDR);
		dumpMap(m_samplevarianceB, basefilename + "_sample_variance_B", EHDR);
		dumpMap(m_buffervariance, basefilename + "_buffer_variance", EHDR);
		dumpMap(m_meansamplevariance, basefilename + "_mean_sample_variance", EHDR);
		dumpMap(m_samplevariancevariance, basefilename + "_sample_variance_variance", EHDR);
		dumpMap(m_filteredbuffervariance, basefilename + "_filtered_buffer_variance", EHDR);
	}

protected:
	int m_vR, m_vF;
	Float m_vSigma, m_vK;
	const ImageBlockF *m_imageblockA, *m_imageblockB;
	BitmapF *m_samplevarianceA, *m_samplevarianceB,
		*m_buffervariance, // estimate of pixel variance using squared difference of buffers
		*m_meansamplevariance, // spp weighted average of sample variances of buffers A and B
		*m_samplevariancevariance, // variance of sample variance - needed for smoothing buffer variance with cross filtering
		*m_filteredbuffervariance;
};

#endif