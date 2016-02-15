// nlmeans.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "nlmeans.h"
#include <string>
#include <iostream>


int _tmain(int argc, char* argv[])
{
	int r = 3, f = 1;
	Float k = 0.7, sigma = 0.1;
	bool dump = true;
	//std::string inputfile="test\\input_1_sigma_40.hdr" , outputfile="test\\output_1_0@4_5_10_denoised_3_3", extension=".hdr";
	std::string inputfile = "test\\cbox_pssmlt_denoise_1st_stage_A.hdr", outputfile = "test\\cbox_pssmlt_denoise_1st_stage_A_denoised_2", extension = ".png";
	//std::string inputfile = "test\\noisy_alley_small.hdr", outputfile = "test\\denoised_alley_small", extension = ".hdr";
	//std::string inputfile = "test\\matlab\\input_gaussian_ones_0_0-3.hdr", outputfile = "test\\matlab\\output_gaussian_ones_0_0-3", extension = ".hdr";

	NLMeansDenoiser<Float, Uchar> *denoiser = new NLMeansDenoiser<Float, Uchar>(r, f, k, sigma, dump);
	BitmapF *input = new BitmapF();
	if (!input->loadBitmap(inputfile)){
		std::cout << "Input bitmap " << inputfile << " failed to load!\n";
		std::cin.get();
		return 0;
	}
	DenoiserInput<Float> *dInput = new DenoiserInput<Float>(outputfile);
	dInput->addImageBlock(new TImageBlock<Float>(Point2i(0., 0.), input->getSize(), 1, false, input));

	std::cout << "Denoising started\n";
	DenoiserOutput<Uchar> *dOutput = denoiser->denoise(dInput);
	std::cout << "Denoising finished in " << dOutput->getDenoiseDuration() << " seconds \n";
	dumpMap<Uchar>(dOutput->getDenoisedImage()->getBitmap(), outputfile, EPNG);
	input->unloadBitmap();
	std::cin.get();
	_CrtDumpMemoryLeaks(); // prints mem leaks
	return 0;
}