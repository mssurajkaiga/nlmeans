// nlmeans.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "nlmeans.h"
#include <string>
#include <iostream>


int _tmain(int argc, char* argv[])
{
	int r = 5, f = 2;
	Float k = 0.55, sigma = 0.4;
	bool dump = true;
	//std::string inputfile="test\\input_1_sigma_40.hdr" , outputfile="test\\output_1_4_17_10_denoised_3", extension=".hdr";
	std::string inputfile = "test\\noisy_alley_small.hdr", outputfile = "test\\denoised_alley_small", extension = ".hdr";

	NLMeansDenoiser<Float> *denoiser = new NLMeansDenoiser<Float>(r, f, k, sigma, dump);
	BitmapF *input = new BitmapF();
	if (!input->loadBitmap(inputfile)){
		std::cout << "Input bitmap " << inputfile << " failed to load!\n";
		std::cin.get();
		return 0;
	}
	DenoiserInput<Float> *dInput = new DenoiserInput<Float>(outputfile);
	dInput->addImageBlock(new TImageBlock<Float>(Point2i(0., 0.), input->getSize(), 1, false, input));

	std::cout << "Denoising started\n";
	DenoiserOutput<Float> *dOutput = denoiser->denoise(dInput);
	std::cout << "Denoising finished in " << dOutput->getDenoiseDuration() << " seconds \n";
	dumpMap<Float>(dOutput->getDenoisedImage()->getBitmap(), outputfile, EHDR);
	input->unloadBitmap();
	std::cin.get();
	_CrtDumpMemoryLeaks(); // prints mem leaks
	return 0;
}