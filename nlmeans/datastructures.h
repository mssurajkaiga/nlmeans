#pragma once
#ifndef _DATASTRUCTURES_H_
#define _DATASTRUCTURES_H_

#include "core.h"

// represents a pixel
template<typename T> class Pixel {
public:
	Pixel(T *data, int channels = 3) : m_channels(channels){
		m_data = new T[channels];
		for (int i = 0; i < channels; ++i)
			m_data[i] = data[i];
	}

	Pixel(T r, T g, T b, T a = -1) {
		m_channels = (a == -1) ? 3 : 4;
		m_data = new T[m_channels];
		m_data[0] = r;
		m_data[1] = g;
		m_data[2] = b;
		if (m_channels == 4)
			m_data[3] = a;
	}

	void clear(const T value = 0) {
		for (int i = 0; i < m_channels; ++i)
			m_data[i] = static_cast<T>(value);
	}

	~Pixel() {
		delete[] m_data;
	}

protected:
	int m_channels; // no. of color channels
	T* m_data;
};

// represents a patch to be denoised
template<typename T> class Patch {
public:
	Patch(int width, int height, int channels, T *data = NULL) : m_width(width), m_height(height), m_channels(channels), m_data(data) {
		m_size = Vector2i(width, height);
		m_datasize = width * height * channels;
		if (data == NULL) {
			m_data = new T[m_datasize];
			clear();
		}
	}

	Patch(const Patch &input) {
		m_width = input.getSize()(0);
		m_height = input.getSize()(1);
		m_channels = input.getChannelCount();
		m_datasize = m_width * m_height * m_channels;
		m_data = new T[m_datasize];
		T *data = input.getData();
		for (int i = 0; i < m_datasize; ++i) {
			*m_data++ = *data++;
		}
	}

	inline void clear(const T value = 0.f) {
		T *data = m_data;
		for (int i = 0; i < m_datasize; ++i) {
			*data++ = static_cast<T>(value);
		}
	}

	inline Pixel<T>* getPixel(int x, int y) {
		T *pos = m_data[(y * m_width + x) * m_channels];
		return new Pixel<T>(pos, m_channels);
	}

	inline T* getData() {
		return static_cast<T*>(m_data);
	}

	inline const T* getData() const {
		return static_cast<T*>(m_data);
	}

	Vector2i getSize() { return m_size; }
	const Vector2i getSize() const { return m_size; }
	Vector2i getChannelCount() { return m_channels; }
	const Vector2i getChannelCount() const { return m_channels; }

	~Patch() {
		if (m_data != NULL)
			delete[] m_data;
	}

private:
	int m_width, m_height, m_channels, m_datasize;
	Vector2i m_size;
	T *m_data;
};


//generic bitmap data class
template<typename T> class TBitmap {
public:
	TBitmap(int w, int h, int comp, const T *d = NULL) : w(w), h(h), channels(comp) {
		size = Vector2i(w, h);
		if (d != NULL)
			data = convert<T, T>(d, w*h*comp);
		else {
			data = new T[w * h * channels];
			clear();
		}
	}

	TBitmap(Vector2i s, int comp, const T *d = NULL) : size(s), channels(comp) {
		w = size(0);
		h = size(1);
		if (data == NULL) {
			data = new T[w * h * channels];
			clear();
		}
		else {
			data = convert<T, T>(d, w*h*comp);
		}
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
		data = new T[w*h*channels];
		T* thisdata = data;
		for (int i = 0; i < w*h*channels; ++i) {
			*thisdata++ = *inputdata++;
		}
	}

	template<typename O> TBitmap<T>(TBitmap<O> *input) {
		w = input->getSize()(0);
		h = input->getSize()(1);
		size = input->getSize();
		channels = input->getChannelCount();
		O* inputdata = input->getData();
		data = convert<O, T>(inputdata, w*h*channels);
	}

	TBitmap(const TBitmap<T> &input) {
		w = input.getSize()(0);
		h = input.getSize()(1);
		size = input.getSize();
		channels = input.getChannelCount();
		T* inputdata = input.getData();
		data = new T[w*h*channels];
		T* thisdata = data;
		for (int i = 0; i < w*h*channels; ++i) {
			*thisdata++ = *inputdata++;
		}
	}

	Vector2i getSize() { return size; }
	const Vector2i getSize() const { return size; }
	int getChannelCount() { return channels; }
	const int getChannelCount() const { return channels; }
	T* getData() { return static_cast<T*>(data); }
	const T* getData() const { return static_cast<T*>(data); }
	void clear(const T value = 0) {
		T *thisdata = static_cast<T*>(data);
		for (int i = 0; i < w*h*channels; ++i) {
			*thisdata++ = static_cast<T>(value);
		}
	}

	//accumulate data from another bitmap after scaling it into current bitmap
	void accumulate(TBitmap<T> *input, T scale = 1) {
		T *tdata = static_cast<T*>(data);
		const T *idata = input->getData();
		for (int i = 0; i < size(0) * size(1) * channels; ++i, ++tdata, ++idata) {
			*tdata += *idata * scale;
		}
	}

	//stores min value of data and max value of data and returns the range of data
	T getDataRange(T &tmin, T &tmax) const {
		const T *tdata = getData();
		T maxValue = static_cast<T>(0), minValue = static_cast<T>(0);
		for (int i = 0; i < size(0) * size(1) * channels; ++i, ++tdata) {
			if (maxValue < *tdata)
				maxValue = *tdata;
			if (minValue > *tdata)
				minValue = *tdata;
		}
		tmin = minValue;
		tmax = maxValue;
		return static_cast<T>(maxValue - minValue);
	}

	//get just the data range
	T getDataRange() const {
		const T *tdata = getData();
		T maxValue = static_cast<T>(0), minValue = static_cast<T>(0);
		for (int i = 0; i < size(0) * size(1) * channels; ++i, ++tdata) {
			if (maxValue < *tdata)
				maxValue = *tdata;
			if (minValue > *tdata)
				minValue = *tdata;
		}
		return static_cast<T>(maxValue - minValue);
	}

	bool loadBitmap(std::string filename, EBitmapType format = EHDR) {
		switch (format) {
		case EPNG:
			filename += ".png";
			break;

		case EBMP:
			filename += ".bmp";
			break;

		case ETGA:
			filename += ".tga";
			break;

		case EHDR:
		default:
			filename += ".hdr";
			break;
		}

		bool is_hdr = stbi_is_hdr(filename.c_str());
		int x, y, n;
		if (is_hdr){
			float *fdata = stbi_loadf(filename.c_str(), &x, &y, &n, 0);
			data = convert<float, T>(fdata, x*y*n);
		}
		else {
			Uchar *udata = stbi_load(filename.c_str(), &x, &y, &n, 0);
			data = convert<Uchar, T>(udata, x*y*n);
		}
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

	void logToFile(const std::string filename = "log.txt") const {
		std::ofstream logFile;
		logFile.open(filename);
		T *thisdata = static_cast<T*>(data);
		for (int j = 0; j < h; ++j)
			for (int i = 0; i < w; ++i) {
				logFile << "x = " << i << " y = " << j << " bitmap value = " << static_cast<Float>(thisdata[0]) << ", " << static_cast<Float>(thisdata[1]) << ", " << static_cast<Float>(thisdata[2]) << "\n";
				thisdata += 3;
			}
		logFile.close();
	}

	~TBitmap() {
		unloadBitmap();
	}
protected:
	int w, h;
	int channels;
	Vector2i size;
	T *data;
};

typedef TBitmap<Float> BitmapF;
typedef TBitmap<int> BitmapI;
typedef TBitmap<Uchar> BitmapC;


template<typename T> class TImageBlock {
public:
	TImageBlock(Point2i offset, Vector2i size, int border, bool warn = false, TBitmap<T> *bitmap = NULL,
		TBitmap<T> *sppbitmap = NULL, TBitmap<T> *varbitmap = NULL, TBitmap<T> *varsbitmap = NULL) : m_offset(offset),
		m_size(size), m_borderSize(border), m_warn(warn), m_bitmap(bitmap), m_sppbitmap(sppbitmap),
		m_varbitmap(varbitmap), m_varsbitmap(varsbitmap)
	{
		if (bitmap == NULL)
			m_bitmap = new TBitmap<T>(size(0), size(1), 3);
		int channelcount = m_bitmap->getChannelCount();
		int x = m_bitmap->getSize()(0);
		int y = m_bitmap->getSize()(1);
		if (sppbitmap == NULL)
			m_sppbitmap = new TBitmap<T>(x, y, channelcount);
		if (varbitmap == NULL)
			m_varbitmap = new TBitmap<T>(x, y, channelcount);
		if (varsbitmap == NULL)
			m_varsbitmap = new TBitmap<T>(x, y, channelcount);
	}

	TImageBlock(TImageBlock<T> *input) {
		m_offset = input->getOffset();
		m_size = input->getSize();
		m_borderSize = input->getBorderSize();
		m_warn = input->getWarn();
		m_bitmap = new TBitmap<T>(input->getBitmap());
		m_sppbitmap = new TBitmap<T>(input->getSppBitmap());
		m_varbitmap = new TBitmap<T>(input->getVarBitmap());
		m_varsbitmap = new TBitmap<T>(input->getVarsBitmap());
	}

	TImageBlock(const TImageBlock<T> &input) {
		m_offset = input.getOffset();
		m_size = input.getSize();
		m_borderSize = input.getBorderSize();
		m_warn = input.getWarn();
		m_bitmap = new TBitmap<T>(input.getBitmap());
		m_sppbitmap = new TBitmap<T>(input.getSppBitmap());
		m_varbitmap = new TBitmap<T>(input->getVarBitmap());
		m_varsbitmap = new TBitmap<T>(input->getVarsBitmap());
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

	/// Return a pointer to the underlying bitmap representation
	inline TBitmap<T>* getSppBitmap() { return m_sppbitmap; }
	inline const TBitmap<T>* getSppBitmap() const{ return m_sppbitmap; }

	/// Return a pointer to the underlying variance bitmap representation
	inline TBitmap<T>* getVarBitmap() { return m_varbitmap; }
	inline const TBitmap<T>* getVarBitmap() const { return m_varbitmap; }

	/// Return a pointer to the underlying variance square bitmap representation
	inline TBitmap<T>* getVarsBitmap() { return m_varsbitmap; }
	inline const TBitmap<T>* getVarsBitmap() const { return m_varsbitmap; }

	/// Clear everything to zero
	inline void clear() {
		m_bitmap->clear();
		m_sppbitmap->clear();
		m_varbitmap->clear();
		m_varsbitmap->clear();
	}
	~TImageBlock() {
		delete m_bitmap;
		delete m_varbitmap;
		delete m_varsbitmap;
		delete m_sppbitmap;
	}

protected:
	Point2i m_offset;
	Vector2i m_size;
	int m_borderSize;
	Float *m_weightsX, *m_weightsY;
	bool m_warn;
	/* Stores samples-per-pixel, per-pixel variance, squared variance
	*/
	TBitmap<T> *m_bitmap, *m_varbitmap, *m_varsbitmap;
	TBitmap<T> *m_sppbitmap;
};

typedef TImageBlock<Float> ImageBlockF;
typedef TImageBlock<Uchar> ImageBlockC;

#endif