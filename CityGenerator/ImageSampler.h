#pragma once
#include "FreeImage.h"
#include "Polygon.h"

class ImageSampler
{
	enum sampleType {density, elevation, valid};
public:
	FIBITMAP* bitmapDensity;
	FIBITMAP* bitmapElevation;
	FIBITMAP* bitmapLand;
	int width, height;

	ImageSampler(void);
	~ImageSampler(void);
	ImageSampler(char* dFileName, char* eFileName, char* lFileName);
	double getDensity(const Polygon p) const;
	double getElevation(const Polygon p) const;
	bool isValid(const Polygon p) const;
	double getDensity(double x, double y) const;
	double getElevation(double x, double y) const;
	bool isValid(double x, double y) const;
	bool isInside(Polygon p, int x, int y) const;
	void done();
	double sample(const Polygon p, sampleType type) const;
	double samplePoint(int x, int y, sampleType type) const;
	int getHeight() const;
	int getWidth() const;
};

