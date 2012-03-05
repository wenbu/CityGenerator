#pragma once
#include "stdafx.h"
#include <vector>
using namespace Eigen;
using namespace std;

class Polygon {
public:
	vector<Vector3d> vertices;
	Polygon(void);
	~Polygon(void);
};

