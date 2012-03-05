#pragma once
#include "stdafx.h"
#include "Intersection.h"
#include <vector>

using namespace Eigen;
using namespace std;

class Intersection; //stupid c++

enum RoadType {
	RT_AXIOM,
	RT_GRID_MINOR,
	RT_GRID_MAJOR,
	RT_HIGHWAY
};

struct ruleAttr {
	RoadType type;

	//angle road is oriented in, relative to x+
	double angle;

	//angle of major axis of road, relative to x+
	double majorAxis;
};

class Road {
	int intersections[2];
	ruleAttr r;
	//double angle;
public:
	Road(void);
	~Road(void);
	//Road(double x1, double y1, double x2, double y2);
	//Road(Intersection* endpoint1, Intersection* endpoint2, double angle);
	//Road(int i1, int i2, double angle);
	Road(int i1, int i2, ruleAttr &r);
	Road(int i1, int i2, RoadType type, double angle, double majorAxis);
	int getIntersection(int index) const;
	double getAngle() const;
	double getAxis() const;
	RoadType getType() const;
	void setIntersection(int intersectionIndex, int newIndex);
};

class ProposedRoad {
public:
	int t;
	bool usesExistingIntersections[2];
	int existingIntersection[2];
	vector<Intersection> newIntersection;
	ruleAttr r;
	bool splitRoad; //specifies whether an existing road needs
					//to be split
	int splitRoadIndex; //specifies road to be split, if necessary

	ProposedRoad(void);
	ProposedRoad(int t, RoadType type, double angle, double majorAxis);
	void setExistingIntersection(int index, int existingIntersection);
	void setNewIntersection(int index, Intersection &newIntersection);
	void setSplitRoad(int index);
};

class ProposedRoadComparator {
public:
	bool operator() (ProposedRoad &r1, ProposedRoad &r2) {
		return r1.t > r2.t;
	}
};

