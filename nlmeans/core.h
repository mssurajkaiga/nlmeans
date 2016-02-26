#pragma once
#ifndef _CORE_H_
#define _CORE_H_

#define _CRT_SECURE_NO_WARNINGS
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include <fstream>
#include <vector>
#include <string>
#include "Eigen/Dense"
using namespace Eigen;

#define STBI_FAILURE_USERMSG

#define FLOAT_EPSILON 0.00000001f
#define INT_EPSILON 1e-10
#define IS_WINDOWS defined(_MSC_VER) && _MSC_VER >= 1400
#define IS_DEBUG false
#ifdef _DEBUG //Visual studio compiler
	#define IS_DEBUG _DEBUG
#else
	#define IS_DEBUG false // insert environment specific debug flag
#endif
typedef ArrayXf GeneralSpectrum;
typedef Array3f RGBSpectrum;
typedef RGBSpectrum Spectrum; /* modify as per need to choose appropriate representations */
typedef float Float;
typedef unsigned char Uchar;

typedef Vector2i Point2i;
typedef Vector2f Point2f;



/// Available Log message types
enum ELogLevel {
	ETrace,   // Trace message, for extremely verbose debugging
	EDebug, // Debug message, usually turned off
	EInfo,  // More relevant debug / information message
	ECustom, // Custom Log level for custom appenders
	EWarn,  // Warning message
	EError  // Error message, causes an exception to be thrown
};

std::vector<std::string> ELogLevelString = { "ETrace", "EDebug", "EInfo", "ECustom", "EWarn", "EError" };

enum EBitmapType {
	EPNG, EBMP, ETGA, EHDR
};

#endif