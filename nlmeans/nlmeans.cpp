// nlmeans.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "nlmeans.h"
#include "modifiednlmeans.h"
#include <string>
#include <iostream>

template<typename I, typename O> int denoiseModifiedNLMeans(std::string inputfile, std::string outputfile, EBitmapType type, int r = 3, int f = 1, Float k = 0.7, Float sigma = 0.1, bool dump = true) {
	ModifiedNLMeansDenoiser<I, O> *denoiser = new ModifiedNLMeansDenoiser<I, O>(r, f, k, sigma, dump);
	TBitmap<I> *inputA = new TBitmap<I>;
	BitmapI *inputsppA = new BitmapI;
	TBitmap<I> *inputvarA = new TBitmap<I>;
	TBitmap<I> *inputvarsA = new TBitmap<I>;
	TBitmap<I> *inputB = new TBitmap<I>;
	BitmapI *inputsppB = new BitmapI;
	TBitmap<I> *inputvarB = new TBitmap<I>;
	TBitmap<I> *inputvarsB = new TBitmap<I>;
	Assert(inputA->loadBitmap(inputfile + "_A", type), "Input bitmap A failed to load!");
	Assert(inputsppA->loadBitmap(inputfile + "_spp_A", type), "Input spp bitmap A failed to load!");
	Assert(inputvarA->loadBitmap(inputfile + "_variance_A", type), "Input variance bitmap A failed to load!");
	Assert(inputvarsA->loadBitmap(inputfile + "_variance_square_A", type), "Input variance square bitmap A failed to load!");
	Assert(inputB->loadBitmap(inputfile + "_B", type), "Input bitmap B failed to load!");
	Assert(inputsppB->loadBitmap(inputfile + "_spp_B", type), "Input spp bitmap B failed to load!");
	Assert(inputvarsB->loadBitmap(inputfile + "_variance_square_B", type), "Input variance square bitmap B failed to load!");
	Assert(inputvarB->loadBitmap(inputfile + "_variance_B", type), "Input variance bitmap B failed to load!");

	DenoiserInput<I> *dInput = new DenoiserInput<I>(outputfile);
	dInput->addImageBlock(new TImageBlock<I>(Point2i(0., 0.), inputA->getSize(), 1, false, inputA, inputsppA, inputvarA, inputvarsA));
	dInput->addImageBlock(new TImageBlock<I>(Point2i(0., 0.), inputB->getSize(), 1, false, inputB, inputsppB, inputvarB, inputvarsB));

	std::cout << "Denoising started\n";
	DenoiserOutput<O> *dOutput = denoiser->denoise(dInput);
	std::cout << "Denoising finished in " << dOutput->getDenoiseDuration() << " seconds \n";
	dumpMap(dOutput->getDenoisedImage()->getBitmap(), outputfile, type);
	inputA->unloadBitmap();
	inputsppA->unloadBitmap();
	inputvarA->unloadBitmap();
	inputvarsA->unloadBitmap();
	inputB->unloadBitmap();
	inputsppB->unloadBitmap();
	inputvarB->unloadBitmap();
	inputvarsB->unloadBitmap();
	std::cin.get();
	_CrtDumpMemoryLeaks(); // prints mem leaks
	return 0;
}

template<typename I, typename O> int denoiseNLMeans(std::string inputfile, std::string outputfile, EBitmapType type, int r = 3, int f = 1, Float k = 0.7, Float sigma = 0.1, bool dump = true) {
	NLMeansDenoiser<I, O> *denoiser = new NLMeansDenoiser<I, O>(r, f, k, sigma, dump);
	TBitmap<I> *input = new TBitmap<I>();
	if (!input->loadBitmap(inputfile, type)){
		std::cout << "Input bitmap " << inputfile << " failed to load!\n";
		std::cin.get();
		return 0;
	}
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
	//std::string inputfile = "test\\cbox_pssmlt_denoise_1st_stage_A", outputfile = "test\\cbox_pssmlt_denoise_1st_stage_A_denoised_4";
	//std::string inputfile = "test\\noisy_alley_small.png", outputfile = "test\\denoised_alley_small2";
	//std::string inputfile = "test\\matlab\\input_gaussian_ones_0_0-3.hdr", outputfile = "test\\matlab\\output_gaussian_ones_0_0-3";

	//return denoiseNLMeans<Float, Float>(inputfile, outputfile, EHDR, r, f, k, sigma, dump);
	std::string inputfile = "test\\cbox_output\\cbox_pssmlt_denoise_1st_stage", outputfile = "test\\cbox_output\\cbox_pssmlt_denoise_1st_stage_denoised";
	return denoiseModifiedNLMeans<Float, Float>(inputfile, outputfile, EHDR, r, f, k, sigma, dump);
}