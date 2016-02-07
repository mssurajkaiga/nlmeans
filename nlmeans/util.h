#ifndef _UTIL_H_
#define _UTIL_H_

#include "core.h"
#include <iostream>
#include <varargs.h>

// Default parameters for cross filtering the buffer estimate - values as per paper.
#define M_R 1
#define M_F 3
#define M_ALPHA 4
#define M_K 0.45

// Utility functions for generic min and max computations
template<typename T> inline T min(const T &a, const T &b) {
	return std::min(a, b);
}

template<> Spectrum inline min(const Spectrum &a, const Spectrum &b) {
	Spectrum retval;
	for (int i = 0; i < retval.size(); ++i)
		retval[i] = std::min(a[i], b[i]);
	return retval;
}

template<typename T> inline T max(const T &a, const T &b) {
	return std::max(a, b);
}

template<> Spectrum inline max(const Spectrum &a, const Spectrum &b) {
	Spectrum retval;
	for (int i = 0; i < retval.size(); ++i)
		retval[i] = std::max(a[i], b[i]);
	return retval;
}

template<typename T> inline T clamp(const T a, const T tmin, const T tmax) {
	return min(max(a, tmin), tmax);
}

inline int computeSquareWindowWidth(int xpos, int ypos, Vector2i bounds, int radius) {
	int xwidth = 0, ywidth = 0;
	if (xpos < radius)
		xwidth = xpos;
	else if (xpos > bounds(0) - radius)
		xwidth = bounds(0) - xpos;
	else
		xwidth = radius;

	if (ypos < radius)
		ywidth = ypos;
	else if (ypos > bounds(0) - radius)
		ywidth = bounds(0) - ypos;
	else
		ywidth = radius;

	return std::max(0, std::min(xwidth, ywidth));
}

// utility function to compute difference between two bitmaps
// absolute - compute absolute difference
template<int C = 1, typename T> inline bool computeBitmapDifference(TBitmap<T> *inputA, TBitmap<T> *inputB, TBitmap<T> *output, const Float scale = 1.f, bool absolute = false) {
	const int &channels = inputA->getChannelCount();
	const Vector2i &size = inputA->getSize();
	if (inputB->getSize() != size || output->getSize() != size ||
		inputB->getChannelCount() != channels || output->getChannelCount() != channels)
		return false;

	T *destA = inputA->getData(),
		*destB = inputB->getData(),
		*destO = output->getData();

	if (C == 1)
		for (size_t i = 0; i < size(0) * size(1); ++i) {
			for (int k = 0; k < channels; ++k) {
				*destO++ = absolute ? abs(*destA - *destB) * scale : (*destA - *destB) * scale;
				++destA;
				++destB;
			}
		}
	else
		for (size_t i = 0; i < size(0) * size(1); ++i) {
			for (int k = 0; k < channels; ++k) {
				*destO++ = absolute ? pow(abs(*destA - *destB), C) * scale : pow(*destA - *destB, C) * scale;
				++destA;
				++destB;
			}
		}

	return true;
}

void Log(ELogLevel level, const char *filename, const char *fmt, ...) {
	FILE * pFile;
	pFile = fopen(filename, "w");

	va_list args;
	va_start(args, fmt);
	vfprintf(pFile, fmt, args);
	va_end(args);
	fclose(pFile);
}

void Log(ELogLevel level, const char *fmt, ...) {
	va_list args;
	va_start(args, fmt);
	Log(level, "log.txt", fmt, args);
	va_end(args);
}


#endif