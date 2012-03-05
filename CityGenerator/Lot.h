#pragma once
#include "stdafx.h"
#include "Polygon.h"
using namespace Eigen;
using namespace std;

class Lot {
public:
	Polygon lot;
	Lot(void);
	~Lot(void);
	Lot(Polygon &p);
};

