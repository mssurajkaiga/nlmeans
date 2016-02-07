#ifndef _CORE_H_
#define _CORE_H_

#define _CRT_SECURE_NO_WARNINGS
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include "Eigen/Dense"
using namespace Eigen;

#define STBI_FAILURE_USERMSG

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


//generic bitmap data class
template<typename T> class TBitmap {
public:
	TBitmap(int w, int h, int comp, const T *d = NULL) :w(w), h(h), channels(comp) {
		size = Vector2i(w, h);
		if (d != NULL)
			data = (void*)d;
		else
			data = new T[w * h * channels];
	}
	TBitmap(TBitmap<T> *input = NULL) {
		if (input == NULL) {
			w = 0;
			h = 0;
			size = Vector2i(w, h);
			channels = 0;
			data = NULL;
			return;
		}
		w = input->getSize()(0);
		h = input->getSize()(1);
		size = input->getSize();
		channels = input->getChannelCount();
		T* inputdata = input->getData();
		data = static_cast<void *>(new T[w*h*channels]);
		T* thisdata = static_cast<T*>(data);
		for (size_t i = 0; i < w*h*channels; ++i) {
			*thisdata++ = *inputdata++;
		}
	}
	Vector2i getSize() { return size; }
	const Vector2i getSize() const { return size; }
	int getChannelCount() { return channels; }
	const int getChannelCount() const { return channels; }
	T* getData() { return static_cast<T*>(data); }
	const T* getData() const { return static_cast<T*>(data); }
	void clear() {
		T *thisdata = static_cast<T*>(data);
		for (size_t i = 0; i < w*h*channels; ++i) {
			*thisdata++ = 0.f;
		}
	}

	bool loadBitmap(std::string filename) {
		bool is_hdr = stbi_is_hdr(filename.c_str());
		int x, y, n;
		if (is_hdr)
			data = static_cast<void*>(stbi_loadf(filename.c_str(), &x, &y, &n, 0));
		else
			data = static_cast<void*>(stbi_load(filename.c_str(), &x, &y, &n, 0));
		if (data != NULL) {
			w = x;
			h = y;
			size = Vector2i(w, h);
			channels = n;
			return true;
		}
		return false;
	}

	void unloadBitmap() {
		stbi_image_free(data);
	}

	~TBitmap() {
		unloadBitmap();
	}
protected:
	int w, h;
	int channels;
	Vector2i size;
	void *data;
};

typedef TBitmap<Float> BitmapF;
typedef TBitmap<int> BitmapI;
typedef TBitmap<Uchar> BitmapC;


template<typename T> class TImageBlock {
public:
	TImageBlock(Point2i offset, Vector2i size, int border, bool warn = false, TBitmap<T> *bitmap = NULL/*,
		BitmapI *sppbitmap = NULL, TBitmap<T> *varbitmap = NULL, TBitmap<T> *varsbitmap = NULL*/) : m_offset(offset),
		m_size(size), m_borderSize(border), m_warn(warn), m_bitmap(bitmap)//, m_sppbitmap(sppbitmap),
		//m_varbitmap(varbitmap), m_varsbitmap(varsbitmap) 
	{}

	TImageBlock(TImageBlock<T> *input) {
		m_offset = input->getOffset();
		m_size = input->getSize();
		m_borderSize = input->getBorderSize();
		m_warn = input->getWarn();
		m_bitmap = new TBitmap<T>(input->getBitmap());
		/*m_sppbitmap = new TBitmap<int>(input->getSppBitmap());
		m_varbitmap = new TBitmap<T>(input->getVarBitmap());
		m_varsbitmap = new TBitmap<T>(input->getVarsBitmap());*/
	}

	/// Set the current block offset
	inline void setOffset(const Point2i &offset) { m_offset = offset; }

	/// Return the current block offset
	inline const Point2i &getOffset() const { return m_offset; }

	/// Set the current block size
	inline void setSize(const Vector2i &size) { m_size = size; }

	/// Return the current block size
	inline const Vector2i &getSize() const { return m_size; }

	/// Return the bitmap's width in pixels
	inline int getWidth() const { return m_size(0); }

	/// Return the bitmap's height in pixels
	inline int getHeight() const { return m_size(1); }

	/// Warn when writing bad sample values?
	inline bool getWarn() const { return m_warn; }

	/// Warn when writing bad sample values?
	inline void setWarn(bool warn) { m_warn = warn; }

	/// Return the border region used by the reconstruction filter
	inline int getBorderSize() const { return m_borderSize; }

	/// Return the number of channels stored by the image block
	inline int getChannelCount() const { return m_bitmap->getChannelCount(); }

	/// Return a pointer to the underlying bitmap representation
	inline TBitmap<T>* getBitmap() { return m_bitmap; }
	inline const TBitmap<T>* getBitmap() const { return m_bitmap; }
/*
	/// Return a pointer to the underlying bitmap representation
	inline TBitmap<int>* getSppBitmap() { return m_sppbitmap; }
	inline const TBitmap<int>* getSppBitmap() const{ return m_sppbitmap; }

	/// Return a pointer to the underlying variance bitmap representation
	inline TBitmap<T>* getVarBitmap() { return m_varbitmap; }
	inline const TBitmap<T>* getVarBitmap() const { return m_varbitmap; }

	/// Return a pointer to the underlying variance square bitmap representation
	inline TBitmap<T>* getVarsBitmap() { return m_varsbitmap; }
	inline const TBitmap<T>* getVarsBitmap() const { return m_varsbitmap; }*/

	/// Clear everything to zero
	inline void clear() {
		m_bitmap->clear();
		//m_sppbitmap->clear();
		//m_varbitmap->clear();
		//m_varsbitmap->clear();
	}

protected:
	Point2i m_offset;
	Vector2i m_size;
	int m_borderSize;
	Float *m_weightsX, *m_weightsY;
	bool m_warn;
	/* Stores samples-per-pixel, per-pixel variance, squared variance
	*/
	TBitmap<T> *m_bitmap;// , *m_varbitmap, *m_varsbitmap;
	//TBitmap<int> *m_sppbitmap;
};

typedef TImageBlock<float> ImageBlockF;
typedef TImageBlock<float> ImageBlockC;

template <typename T> void dumpMap(const TBitmap<T> *bitmap, std::string filename, EBitmapType format, std::string bitmapname="") {
	if (bitmapname == "") {
		bitmapname = filename;
	}
	//Log(EInfo, "Writing bitmap %s to \"%s\" ..", bitmapname.c_str(), filename.c_str());
	int x = bitmap->getSize()(0);
	int y = bitmap->getSize()(1);
	int comp = bitmap->getChannelCount();
	const T *data = bitmap->getData();

	switch (format) {
		case EPNG:
			stbi_write_png(std::string(filename+".png").c_str(), x, y, comp, data, 16);
			break;

		case EBMP:
			stbi_write_bmp(std::string(filename + ".bmp").c_str(), x, y, comp, data);
			break;

		case ETGA:
			stbi_write_tga(std::string(filename + ".tga").c_str(), x, y, comp, data);
			break;

		case EHDR:
			stbi_write_hdr(std::string(filename + ".hdr").c_str(), x, y, comp, data);
			break;
		default:
			//Log(EError, "Invalid format!");
			break;
	}
}


#endif