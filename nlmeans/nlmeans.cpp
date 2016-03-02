// nlmeans.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "nlmeans.h"
#include "modifiednlmeans.h"
#include "mltnlmeans.h"
#include <string>
#include <iostream>

template<typename I, typename O> int denoiseMLTNLMeans(std::string inputfile, EBitmapType inputtype, std::string outputfile, EBitmapType outputtype, int r = 3, int f = 1, Float k = 0.7, Float sigma = 0.1, int numinstances = 1, bool dump = true) {
	MLTNLMeansDenoiser<I, O> *denoiser = new MLTNLMeansDenoiser<I, O>(r, f, k, sigma, numinstances, dump);

	DenoiserInput<I> *dInput = new DenoiserInput<I>(outputfile);
	// load bitmaps into DenoiserInput
	for (int i = 1; i <= numinstances; ++i) {
		TBitmap<I> *inputBitmap = new TBitmap<I> ;
		std::string stri = std::to_string(i);
		Assert(inputBitmap->loadBitmap(inputfile + "_" + stri, inputtype), "Input bitmap no."+ stri +" failed to load!");
		dInput->addImageBlock(new TImageBlock<I>(Point2i(0., 0.), inputBitmap->getSize(), 1, false, inputBitmap));
	}

	DenoiserOutput<O> *dOutput = denoiser->denoise(dInput);
	LOG(EInfo, "Denoising finished. Time taken = %d seconds", dOutput->getDenoiseDuration());
	dumpMap(dOutput->getDenoisedImage()->getBitmap(), outputfile, outputtype);

	for (int i = 0; i < numinstances; ++i) {
		dInput->getImageBlocks()[i]->getBitmap()->unloadBitmap();
	}
	std::cin.get();
	_CrtDumpMemoryLeaks(); // prints mem leaks
	return 0;
}


template<typename I, typename O> int denoiseModifiedNLMeans(std::string inputfile, EBitmapType inputtype, std::string outputfile, EBitmapType outputtype, int r = 3, int f = 1, Float k = 0.7, Float sigma = 0.1, bool dump = true) {
	ModifiedNLMeansDenoiser<I, O> *denoiser = new ModifiedNLMeansDenoiser<I, O>(r, f, k, sigma, dump);
	TBitmap<I> *inputA = new TBitmap<I>;
	//BitmapI *inputsppA = new BitmapI;
	TBitmap<I> *inputsppA = new TBitmap<I>;
	TBitmap<I> *inputvarA = new TBitmap<I>;
	TBitmap<I> *inputvarsA = new TBitmap<I>;
	TBitmap<I> *inputB = new TBitmap<I>;
	//BitmapI *inputsppB = new BitmapI;
	TBitmap<I> *inputsppB = new TBitmap<I>;
	TBitmap<I> *inputvarB = new TBitmap<I>;
	TBitmap<I> *inputvarsB = new TBitmap<I>;
	Assert(inputA->loadBitmap(inputfile + "_A", inputtype), "Input bitmap A failed to load!");
	Assert(inputsppA->loadBitmap(inputfile + "_spp_A", inputtype), "Input spp bitmap A failed to load!");
	Assert(inputvarA->loadBitmap(inputfile + "_variance_A", inputtype), "Input variance bitmap A failed to load!");
	Assert(inputvarsA->loadBitmap(inputfile + "_variance_square_A", inputtype), "Input variance square bitmap A failed to load!");
	Assert(inputB->loadBitmap(inputfile + "_B", inputtype), "Input bitmap B failed to load!");
	Assert(inputsppB->loadBitmap(inputfile + "_spp_B", inputtype), "Input spp bitmap B failed to load!");
	Assert(inputvarsB->loadBitmap(inputfile + "_variance_square_B", inputtype), "Input variance square bitmap B failed to load!");
	Assert(inputvarB->loadBitmap(inputfile + "_variance_B", inputtype), "Input variance bitmap B failed to load!");

	DenoiserInput<I> *dInput = new DenoiserInput<I>(outputfile);
	dInput->addImageBlock(new TImageBlock<I>(Point2i(0., 0.), inputA->getSize(), 1, false, inputA, inputsppA, inputvarA, inputvarsA));
	dInput->addImageBlock(new TImageBlock<I>(Point2i(0., 0.), inputB->getSize(), 1, false, inputB, inputsppB, inputvarB, inputvarsB));
	DenoiserOutput<O> *dOutput = denoiser->denoise(dInput);
	LOG(EInfo, "Denoising finished. Time taken = %d seconds", dOutput->getDenoiseDuration());
	dumpMap(dOutput->getDenoisedImage()->getBitmap(), outputfile, outputtype);
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

template<typename I, typename O> int denoiseNLMeans(std::string inputfile, EBitmapType inputtype, std::string outputfile, EBitmapType outputtype, int r = 3, int f = 1, Float k = 0.7, Float sigma = 0.1, bool dump = true) {
	NLMeansDenoiser<I, O> *denoiser = new NLMeansDenoiser<I, O>(sigma, dump);
	TBitmap<I> *input = new TBitmap<I>();
	if (!input->loadBitmap(inputfile, inputtype)){
		std::cout << "Input bitmap " << inputfile << " failed to load!\n";
		std::cin.get();
		return 0;
	}
	DenoiserInput<I> *dInput = new DenoiserInput<I>(outputfile);
	dInput->addImageBlock(new TImageBlock<I>(Point2i(0., 0.), input->getSize(), 1, false, input));
	DenoiserOutput<O> *dOutput = denoiser->denoise(dInput);
	LOG(EInfo, "Denoising finished. Time taken = %d seconds", dOutput->getDenoiseDuration());
	dumpMap(dOutput->getDenoisedImage()->getBitmap(), outputfile, outputtype);
	input->unloadBitmap();
	std::cin.get();
	_CrtDumpMemoryLeaks(); // prints mem leaks
	return 0;
}


int _tmain(int argc, char* argv[])
{
	int r = 10, f = 4;
	Float k = 0.45f, sigma = 20.f, alpha = 1.f;
	bool dump = true, ltf = true, ltc = true;
	std::string inputfile, outputfile;
	// parse arguments
	for (int count = 1; count < argc; ++count) {
		std::string arg(argv[count]);
		if ( arg == "-i" || arg == "--input") {
			inputfile = std::string(argv[++count]);
		}
		if (arg == "-o" || arg == "--output") {
			outputfile = std::string(argv[++count]);
		}
		if (arg == "-l" || arg == "--log") {
			ltf = true;
		}
		if (arg == "-v" || arg == "--verbose") {
			ltc = true;
		}
	}

	//inputfile = "test\\input_1_sigma_40"; outputfile = "test\\output_1_0@4_5_10_denoised_3_5";
	//inputfile = "test\\cbox_pssmlt_denoise_1st_stage_A"; outputfile = "test\\cbox_pssmlt_denoise_1st_stage_A_denoised_5";
	//inputfile = "test\\noisy_alley_small"; outputfile = "test\\denoised_alley_small2";
	//inputfile = "test\\matlab\\input_gaussian_40"; outputfile = "test\\matlab\\output_gaussian_40";
	//inputfile = "results\\nlmeans\\dice_scene\\noisy"; outputfile = "results\\nlmeans\\dice_scene\\denoised";
	//inputfile = "results\\nlmeans\\white_gaussian_scene\\noisy"; outputfile = "results\\nlmeans\\white_gaussian_scene\\denoised";
	//inputfile = "results\\modified_nlmeans\\veach_door_diffuse\\veach_door_diffuse_1st_stage"; outputfile = "results\\modified_nlmeans\\veach_door_diffuse\\veach_door_diffuse_1st_stage_denoised";
	//Logger::initialize(ltf, ltc, inputfile + ".txt");

	//return denoiseNLMeans<Uchar, Uchar>(inputfile, EHDR, outputfile, EHDR, r, f, k, sigma, dump);
	//inputfile = "test\\cbox_output\\cbox_pssmlt_denoise_1st_stage"; outputfile = "test\\cbox_output\\cbox_pssmlt_denoise_1st_stage_denoised11";
	//inputfile = "results\\modified_nlmeans\\veach_door_diffuse\\veach_door_diffuse_1st_stage"; outputfile = "results\\modified_nlmeans\\veach_door_diffuse\\veach_door_diffuse_1st_stage_denoised";

	//Logger::initialize(ltf, ltc, inputfile + ".log");
	//return denoiseModifiedNLMeans<Float, Float>(inputfile, EHDR, outputfile, EHDR, r, f, k, alpha, dump);

	inputfile = "results\\mlt_nlmeans\\veach_door_diffuse\\parallel_renders\\veach_door_diffuse"; outputfile = "results\\mlt_nlmeans\\veach_door_diffuse\\veach_door_diffuse_denoised";
	Logger::initialize(ltf, ltc, inputfile + ".txt");
	return denoiseMLTNLMeans<Float, Float>(inputfile, EHDR, outputfile, EHDR, r, f, k, alpha, 8, dump);
}