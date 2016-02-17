// nlmeans.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "nlmeans.h"
#include <string>
#include <iostream>

template<typename I, typename O> int denoiseMain(std::string inputfile, std::string outputfile, EBitmapType type, int r = 3, int f = 1, Float k = 0.7, Float sigma = 0.1, bool dump = true) {
	NLMeansDenoiser<I, O> *denoiser = new NLMeansDenoiser<I, O>(r, f, k, sigma, dump);
	TBitmap<I> *input = new TBitmap<I>();
	if (!input->loadBitmap(inputfile, type)){
		std::cout << "Input bitmap " << inputfile << " failed to load!\n";
		std::cin.get();
		return 0;
	}
	I mi, ma;
	I value = input->getDataRange(mi, ma);
	DenoiserInput<I> *dInput = new DenoiserInput<I>(outputfile);
	dInput->addImageBlock(new TImageBlock<I>(Point2i(0., 0.), input->getSize(), 1, false, input));

	std::cout << "Denoising started\n";
	DenoiserOutput<O> *dOutput = denoiser->denoise(dInput);
	std::cout << "Denoising finished in " << dOutput->getDenoiseDuration() << " seconds \n";
	dumpMap(dOutput->getDenoisedImage()->getBitmap(), outputfile, type);
	input->unloadBitmap();
	std::cin.get();
	_CrtDumpMemoryLeaks(); // prints mem leaks
	return 0;
}

int _tmain(int argc, char* argv[])
{
	int r = 4, f = 2;
	Float k = 0.45f, sigma = 40.f;
	bool dump = true;
	//std::string inputfile="test\\input_1_sigma_40.png" , outputfile="test\\output_1_0@4_5_10_denoised_3_4";
	std::string inputfile = "test\\cbox_pssmlt_denoise_1st_stage_A", outputfile = "test\\cbox_pssmlt_denoise_1st_stage_A_denoised_4";
	//std::string inputfile = "test\\noisy_alley_small.png", outputfile = "test\\denoised_alley_small2";
	//std::string inputfile = "test\\matlab\\input_gaussian_ones_0_0-3.hdr", outputfile = "test\\matlab\\output_gaussian_ones_0_0-3";

	return denoiseMain<Float, Float>(inputfile, outputfile, EHDR, r, f, k, sigma, dump);
}