#pragma once
#ifndef _UTIL_H_
#define _UTIL_H_

#include "datastructures.h"
#include <iostream>
#include <varargs.h>

// Default parameters for cross filtering the buffer estimate - values as per paper.
#define M_R 1
#define M_F 3
#define M_ALPHA 4
#define M_K 0.45

#define LOG(X, ...) Logger::getInstance()->log(X, ##__VA_ARGS__) 


class Logger;
static Logger *logger = NULL;

class Logger {
public:
	static void initialize(bool ltf = true, bool ltc = true, std::string filename = "") {
		if (ltf && filename == "")
			filename = "renderlog.txt";
		// clear existing log file with same filename
		FILE * pFile;
		pFile = fopen(filename.c_str(), "w");
		fclose(pFile);
		logger = new Logger(ltf, ltc, filename);
	}

	static Logger* getInstance() {
		if (logger == NULL)
			logger = new Logger();
		return logger;
	}

	void log(ELogLevel level, const char *fmt, ...) {
		va_list args;
		va_start(args, fmt);
		if (logToFile)
			LogToFile(level, m_filename, fmt, args);
		if (logToConsole)
			LogToConsole(level, fmt, args);
		va_end(args);
	}

	bool getLtf() { return logToFile; }
	bool getLtc() { return logToConsole; }
	std::string getFilename() { return m_filename;  }

protected:
	Logger(bool ltf = true, bool ltc = true, std::string filename = "") : logToFile(ltf), logToConsole(ltc), m_filename(filename) {
		if (ltf && filename == "")
			m_filename = "renderlog.txt";
	}

	Logger(Logger *logger) {
		logToFile = logger->getLtf();
		logToConsole = logger->getLtc();
		m_filename = logger->getFilename();
	}

	void LogToConsole(ELogLevel level, const char *fmt, va_list args) {
#if IS_WINDOWS
		printf("%s: ",ELogLevelString[level]);
		vprintf_s(fmt, args);
#else
		printf("%s: ", ELogLevelString[level].c_str());
		vprintf(fmt, args);
#endif
		printf("\n");
	}

	void LogToFile(ELogLevel level, std::string filename, const char *fmt, va_list args) {
		FILE * pFile;
		pFile = fopen(filename.c_str(), "a+");
		fprintf(pFile, "%s: ", ELogLevelString[level].c_str());
		vfprintf(pFile, fmt, args);
		fprintf(pFile, "\n");
		fclose(pFile);
	}

protected:
	std::string m_filename;
	bool logToFile, logToConsole;
};


void Assert(bool expression, std::string customError= "") {
	if (!expression) {
		LOG(EError, std::string("Error: " + customError).c_str());
		//assert(false);
		std::cin.get();
		exit(0);
	}
	//assert(true);
}

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
template<int C = 1, typename T> inline bool computeBitmapDifference(const TBitmap<T> *inputA, const TBitmap<T> *inputB, TBitmap<T> *output, const Float scale = 1.f, bool absolute = false) {
	const int &channels = inputA->getChannelCount();
	const Vector2i &size = inputA->getSize();
	if (inputB->getSize() != size || output->getSize() != size ||
		inputB->getChannelCount() != channels || output->getChannelCount() != channels)
		return false;

	const T *destA = inputA->getData(),
		*destB = inputB->getData();
	T *destO = output->getData();

	if (C == 1)
		for (int i = 0; i < size(0) * size(1); ++i) {
			for (int k = 0; k < channels; ++k) {
				*destO++ = absolute ? abs(*destA - *destB) * scale : (*destA - *destB) * scale;
				++destA;
				++destB;
			}
		}
	else
		for (int i = 0; i < size(0) * size(1); ++i) {
			for (int k = 0; k < channels; ++k) {
				*destO++ = absolute ? pow(abs(*destA - *destB), C) * scale : pow(*destA - *destB, C) * scale;
				++destA;
				++destB;
			}
		}

	return true;
}


//utility function to scale Bitmap using another Bitmap
template<typename T, typename I> bool scaleBitmap(const TBitmap<T> *input, const TBitmap<I> *scale, TBitmap<T> *output) {

	const int &channels = input->getChannelCount();
	const int &sChannels = scale->getChannelCount();
	const Vector2i &size = input->getSize();
	if (scale->getSize() != size || output->getSize() != size)
		return false;
	if (output->getChannelCount() != channels)
		return false;

	const T *dest1 = input->getData(),
	const I	*dest2 = scale->getData();
	T *destO = output->getData();

	if (sChannels == channels) {
		for (int i = 0; i < size(0) * size(1); ++i)
			for (int k = 0; k < channels; ++k)
				*destO++ = (*dest1++ * static_cast<T>(*dest2++));
	}
	else if (sChannels == 1) {
		for (int i = 0; i < size(0) * size(1); ++i, ++dest2)
			for (int k = 0; k < channels; ++k)
				*destO++ = (*dest1++ * static_cast<T>(*dest2));
	}
	else
		return false;

	return true;
}

//utility function to normalize Bitmap using another Bitmap
template<typename T, typename I> bool normalizeBitmap(const TBitmap<T> *input, const TBitmap<I> *normal, TBitmap<T> *output) {

	const int &channels = input->getChannelCount();
	const int &sChannels = normal->getChannelCount();
	const Vector2i &size = input->getSize();
	if (normal->getSize() != size || output->getSize() != size)
		return false;
	if (output->getChannelCount() != channels)
		return false;

	const T *dest1 = input->getData();
	const I	*dest2 = normal->getData();
	T *destO = output->getData();

	if (sChannels == channels) {
		for (int i = 0; i < size(0) * size(1); ++i)
			for (int k = 0; k < channels; ++k) {
				T den = static_cast<T>(*dest2++);
				if (den > 0.f)
					*destO++ = *dest1++ / den;
				else
					*destO++ = *dest1++;
			}
	}
	else if (sChannels == 1) {
		for (int i = 0; i < size(0) * size(1); ++i, ++dest2) {
			T den = static_cast<T>(*dest2);
			for (int k = 0; k < channels; ++k) {
				if (den > 0.f)
					*destO++ = *dest1++ / den;
				else
					*destO++ = *dest1++;
			}
		}
	}
	else
		return false;

	return true;
}


// UTILITY functions for bitmap data inter-conversion
// allocate memory and convert data from one type to another
template<typename I, typename O> O* convert(const I *input, int size) {
	O *output = new O[size];
	for (int i = 0; i < size; ++i)
		output[i] = static_cast<O>(input[i]);
	return output;
}

// convert data from one type to another assuming memory is already allocated
template<typename I, typename O> void convert(const I *input, int size, O *output) {
	for (int i = 0; i < size; ++i)
		output[i] = static_cast<O>(input[i]);
}

// UTILITY functions for bitmap inter-conversion
template<typename I, typename O> TBitmap<O>* convert(TBitmap<I> *input) {
	if (std::is_same<I, O>::value)
		return new TBitmap<O>(input);
	Vector2i size = input->getSize();
	int channels = input->getChannelCount();
	TBitmap<O> *output = new TBitmap<O>(size, channels, NULL);
	convert(input->getData(), size(0)*size(1)*channels, output->getData());
	return output;
}

// special case for handling float to unsigned char conversion
template<> BitmapC* convert(BitmapF *input) {
	BitmapC *output = new BitmapC(input->getSize(), input->getChannelCount(), NULL);
	Float minValue = 0.f, maxValue = 0.f;
	Float invScale = 1.f / input->getDataRange(minValue, maxValue);
	Uchar *outputdata = output->getData();
	Float *inputdata = input->getData();
	Vector2i size = input->getSize();
	int channels = input->getChannelCount();
	for (int i = 0; i < size(0)*size(1)*channels; ++i) {
		*outputdata++ = static_cast<Uchar>((*inputdata++ - minValue) * invScale * 255);
	}
	return output;
}

// UTILITY Functions for imageblock data inter-conversion
template<typename I, typename O> TImageBlock<O>* convert(TImageBlock<I> *input) {
	//BitmapI *outputsppbitmap = new BitmapI(input->getSppBitmap());
	TBitmap<O> *outputsppbitmap = convert<I, O>(input->getSppBitmap());
	TBitmap<O> *outputbitmap = convert<I, O>(input->getBitmap());
	TBitmap<O> *outputvarbitmap = convert<I, O>(input->getVarBitmap());
	TBitmap<O> *outputvarsbitmap = convert<I, O>(input->getVarsBitmap());
	if (outputbitmap == NULL || outputvarbitmap == NULL || outputvarsbitmap == NULL) {
		std::cout << "Conversion between imageblocks of these data types failed!\n";
		return NULL;
	}
	TImageBlock<O> *output = new TImageBlock<O>(input->getOffset(), input->getSize(), input->getBorderSize(),
		input->getWarn(), outputbitmap, outputsppbitmap, outputvarbitmap, outputvarsbitmap);
	return output;
}


template <typename T> int dumpMap(const TBitmap<T> *bitmap, std::string filename, EBitmapType format = EPNG, std::string bitmapname = "") {
	if (bitmapname == ""){
		bitmapname = filename;
	}
	LOG(EInfo, "Writing bitmap %s", bitmapname.c_str());
	int x = bitmap->getSize()(0);
	int y = bitmap->getSize()(1);
	int comp = bitmap->getChannelCount();
	const T *data = bitmap->getData();
	Uchar *chardata = NULL; float *floatdata = NULL;
	int ret;

	switch (format) {
	case EPNG:
		chardata = convert<T, Uchar>(data, x * y * comp);
		ret = stbi_write_png(std::string(filename + ".png").c_str(), x, y, comp, chardata, 0); // pass 0 to let library calculate row_stride
		break;

	case EBMP:
		chardata = convert<T, Uchar>(data, x * y * comp);
		ret = stbi_write_bmp(std::string(filename + ".bmp").c_str(), x, y, comp, chardata);
		break;

	case ETGA:
		chardata = convert<T, Uchar>(data, x * y * comp);
		ret = stbi_write_tga(std::string(filename + ".tga").c_str(), x, y, comp, chardata);
		break;

	case EHDR:
		floatdata = convert<T, float>(data, x * y * comp);
		ret = stbi_write_hdr(std::string(filename + ".hdr").c_str(), x, y, comp, floatdata);
		break;

	default:
		ret = -1;
		LOG(EError, "Invalid format!");
		break;
	}
	if (ret != -1) {
		LOG(EInfo, "Bitmap written to %s", filename.c_str());
	}
	return ret;
}

#endif