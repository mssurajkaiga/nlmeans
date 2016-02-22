#pragma once
#ifndef _CORE_H_
#define _CORE_H_

#define _CRT_SECURE_NO_WARNINGS
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include <fstream>
#include "Eigen/Dense"
using namespace Eigen;

#define STBI_FAILURE_USERMSG

#define FLOAT_EPSILON 0.00000001f
#define INT_EPSILON 1e-10
#define IS_WINDOWS defined(_MSC_VER) && _MSC_VER >= 1400
typedef ArrayXf GeneralSpectrum;
typedef Array3f RGBSpectrum;
typedef RGBSpectrum Spectrum; /* modify as per need to choose appropriate representations */
typedef float Float;
typedef unsigned char Uchar;

typedef Vector2i Point2i;
typedef Vector2f Point2f;

/// Available Log message types
enum ELogLevel {
	ETrace = 0,   // Trace message, for extremely verbose debugging
	EDebug = 100, // Debug message, usually turned off
	EInfo = 200,  // More relevant debug / information message
	ECustom = 250, // Custom Log level for custom appenders
	EWarn = 300,  // Warning message
	EError = 400  // Error message, causes an exception to be thrown
};

enum EBitmapType {
	EPNG, EBMP, ETGA, EHDR
};

#endif