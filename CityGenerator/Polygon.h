#pragma once
#include "stdafx.h"
#include <vector>
using namespace Eigen;
using namespace std;

class Polygon {
public:
	vector<Vector3d> vertices;
	vector<int> cutPoints;
	Polygon(void);
	~Polygon(void);
	//implement
	void insertCutPoints(vector<int> &indices, vector<Vector3d> &points);
};

