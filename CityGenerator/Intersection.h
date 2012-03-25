#pragma once
#include "stdafx.h"
#include "Road.h"
#include <vector>

using namespace std;
using namespace Eigen;

class Road;

class Intersection {
	vector<int> roads;
	//intersectiontype?
public:
	Vector3d position;
	Intersection(void);
	~Intersection(void);
	Intersection(double x, double y);
	Intersection(Vector3d &position);
	void getPosition(double &x, double &y, double &z) const;
	void printPosition() const {
		printf("Intersection: [%f,%f]\n",position[0],position[1]);
		printf("           %f at 0x%p\n", position[0],&position[0]);
		printf("           %f at 0x%p\n", position[1],&position[1]);
		printf("           %f at 0x%p\n", position[2],&position[2]);
	}
	double operator[] (unsigned int dim) const;
	void getRoads(vector<int> &roads) const;
	int numRoads() const {return roads.size();}
	void addRoad(int roadIndex);
	bool removeRoad(int roadIndex);
	bool isEndpoint() const;
};

