#include "StdAfx.h"
#include "ImageSampler.h"
#include "FreeImage\FreeImage.h"
#include "limits.h"
#include "Polygon.h"
#include <algorithm>
#include <iostream>
using namespace std;

ImageSampler::ImageSampler(void) {
}

ImageSampler::~ImageSampler(void) {
}

ImageSampler::ImageSampler(char* dFileName, char* eFileName, char* lFileName) {
	bitmapDensity = FreeImage_Load(FIF_PNG,dFileName,PNG_DEFAULT);
	bitmapElevation = FreeImage_Load(FIF_PNG,eFileName,PNG_DEFAULT);
	bitmapLand = FreeImage_Load(FIF_PNG,lFileName,PNG_DEFAULT);
	width = FreeImage_GetWidth(bitmapDensity);
	height = FreeImage_GetHeight(bitmapDensity);

	if(!bitmapDensity) {
		cout << "Failed to load " << dFileName << endl;
		exit(1);
	}
	if(!bitmapElevation) {
		cout << "Failed to load " << eFileName << endl;
		exit(1);
	}
	if(!bitmapLand) {
		cout << "Failed to load " << lFileName << endl;
		exit(1);
	}
	if(width != FreeImage_GetWidth(bitmapElevation) && width != FreeImage_GetWidth(bitmapLand) && 
		height != FreeImage_GetHeight(bitmapElevation) && height != FreeImage_GetHeight(bitmapLand)) {
		cout << "Maps are of different size" << endl;
		exit(1);
	}
}

//"http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html#The C Method"
bool ImageSampler::isInside(const Polygon p, int x, int y) const {
	int size = p.vertices.size();
	int i,j;
	bool c = false;
	for(i = 0, j = size - 1; i < size; j = i++) {
		if(((p.vertices[i][1] > y) != (p.vertices[j][1] > y)) &&
			(x < 
			(p.vertices[j][0] - p.vertices[i][0]) * (y - p.vertices[i][1]) / 
			(p.vertices[j][1] - p.vertices[i][1]) + p.vertices[i][0])) {
			c = !c;
		}
	}
	return c;
}

double ImageSampler::samplePoint(int x, int y, sampleType type) const {
	FIBITMAP* bitmap;
	if(type == density) {
		bitmap = bitmapDensity;
	} else if(type == elevation) {
		bitmap = bitmapElevation;
	} else if(type == valid) {
		bitmap = bitmapLand;
	} else {
		return 0;
	}

	RGBQUAD c;
	int r,g,b;
	if(x < width && y < height && x >= 0 && y >= 0) {
		FreeImage_GetPixelColor(bitmap,x,y,&c);
		r = c.rgbRed;
		g = c.rgbGreen;
		b =	c.rgbBlue;
		//cout << "in " << x << "," << y << ":" << r << " " << g << " " << b << endl;
	} else {
		r = 0;
		g = 0;
		b =	0;
	}

	if(type == density || type == elevation || type == valid) {
		return (r+g+b)/3;
	} else {
		return 0;
	}
}

double ImageSampler::sample(const Polygon p, sampleType type) const {
	FIBITMAP* bitmap;
	if(type == density) {
		bitmap = bitmapDensity;
	} else if(type == elevation) {
		bitmap = bitmapElevation;
	} else if(type == valid) {
		bitmap = bitmapLand;
	} else {
		return 0;
	}
	
	double maxX = -1, maxY = -1;
	double minX = INT_MAX, minY = INT_MAX;
	for(unsigned i = 0; i < p.vertices.size(); i++) {
		maxX = max(maxX, p.vertices[i][0]);
		maxY = max(maxY, p.vertices[i][1]);
		minX = min(minX, p.vertices[i][0]);
		minY = min(minY, p.vertices[i][1]);
	}


	double total = 0, ctr = 0;
	RGBQUAD c;
	int r,g,b;
	bool isValid = true;
	for(int x = (int)minX; x < (int)maxX; x++) {
		for(int y = (int)minY; y < (int)maxY; y++) {
			if(isInside(p,x,y)) {
				if(x < width && y < height && x >= 0 && y >= 0) {
					FreeImage_GetPixelColor(bitmap,x,y,&c);
					r = c.rgbRed;
					g = c.rgbGreen;
					b =	c.rgbBlue;
				} else {
					r = 0;
					g = 0;
					b =	0;
				}
				//cout << "in " << x << "," << y << ":" << r << " " << g << " " << b << endl;
				total = total + (r+g+b)/3; //assuming grayscale input, r,g,b values should be the same
				ctr++;
				if((r+g+b)/3 == 0) {
					isValid = false;
				}
			} else {
				//cout << "out " << x << "," << y << endl;
			}
		}
	}
	
	if(type == density) {
		printf("\t\tdensity sample\n");
		printf("\t\t total = %f\n", total);
		printf("\t\t ctr = %i\n", ctr);
		return total/ctr;
	} else if(type == elevation) {
		return total/ctr;
	} else if(type == valid) {
		return isValid;
	} else {
		return 0;
	}
}

double ImageSampler::getDensity(const Polygon p) const {
	return sample(p,density);
}

double ImageSampler::getElevation(const Polygon p) const {
	return sample(p,elevation);
}

bool ImageSampler::isValid(const Polygon p) const {
	return sample(p,valid) != 0;
}

double ImageSampler::getDensity(double x, double y) const {
	//printf("\t\tsampling at %f, %f\n", x, y);
	return samplePoint((int)(x*width),(int)(y*height),density);
}

double ImageSampler::getElevation(double x, double y) const {
	return samplePoint((int)(x*width),(int)(y*height),elevation);
}

bool ImageSampler::isValid(double x, double y) const {
	//printf("got sample %i\n", sample);
	return samplePoint((int)(x*width),(int)(y*height),valid) != 0;
}

int ImageSampler::getHeight() const {
	return height;
}

int ImageSampler::getWidth() const {
	return width;
}

void ImageSampler::done() {
	FreeImage_Unload(bitmapDensity);
	FreeImage_Unload(bitmapElevation);
	FreeImage_Unload(bitmapLand);
}




