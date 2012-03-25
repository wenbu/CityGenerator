#pragma once
#include "Lot.h"
#include "Polygon.h"
#include "ImageSampler.h"
#include <vector>
#include "Eigen/StdVector"
#include <iostream>
#include <fstream>

//todo: clean up implementation
// Mesh implementation + extrude operator?
using namespace Eigen;
using namespace std;

class BuildingGenerator {

public:
	vector<Lot, aligned_allocator<Lot>> lots; 
	
	vector<vector<Polygon>> buildings;

	BuildingGenerator(void);
	~BuildingGenerator(void);
	BuildingGenerator(vector<Lot> &l);
	void generateObjFile();
	void getBuildings(double maxHeight, ImageSampler &is);
	void getPolygons(vector<Polygon> &polygons) const;
	double jitter(double c, double r) const;
};

