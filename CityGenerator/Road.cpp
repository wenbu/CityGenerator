#include "StdAfx.h"
#include "Road.h"

Road::Road(void) {
	intersections[0] = 0;
	intersections[1] = 0;
	r.angle = 0;
	r.majorAxis = 0;
	r.type = RT_GRID_MINOR;
}

Road::~Road(void) {}

Road::Road(int i1, int i2, ruleAttr &rule) {
	intersections[0] = i1;
	intersections[1] = i2;

	//set up ruleAttr
	r.type = rule.type;
	r.angle = rule.angle;
	r.majorAxis = rule.majorAxis;
}

Road::Road(int i1, int i2, RoadType t, double a, double mA) {
	intersections[0] = i1;
	intersections[1] = i2;
	r.type = t;
	r.angle = a;
	r.majorAxis = mA;
}

int Road::getIntersection(int n) const {
	return intersections[n];
}

double Road::getAngle() const {
	return r.angle;
}

double Road::getAxis() const {
	return r.majorAxis;
}

RoadType Road::getType() const {
	return r.type;
}

void Road::setIntersection(int i, int newindex) {
	intersections[i] = newindex;
}

ProposedRoad::ProposedRoad(void) {
	t = 0;
	r.type = RT_GRID_MINOR;
	r.angle = 0;
	r.majorAxis = 0;
	usesExistingIntersections[0] = 1;
	usesExistingIntersections[1] = 1;
	existingIntersection[0] = 0;
	existingIntersection[1] = 0;
	newIntersection = vector<Intersection>(2);
	splitRoad = 0;
	splitRoadIndex = 0;
}

ProposedRoad::ProposedRoad(int time, RoadType type, double a, double mA) {
	t = time;
	r.type = type;
	r.angle = a;
	r.majorAxis = mA;
	usesExistingIntersections[0] = 1;
	usesExistingIntersections[1] = 1;
	existingIntersection[0] = 0;
	existingIntersection[1] = 0;
	newIntersection = vector<Intersection>(2);
	splitRoad = 0;
	splitRoadIndex = 0;
}

void ProposedRoad::setExistingIntersection(int i, int nI) {
	usesExistingIntersections[i] = 1;
	existingIntersection[i] = nI;
}

void ProposedRoad::setNewIntersection(int i, Intersection &nI) {
	usesExistingIntersections[i] = 0;
	newIntersection[i] = nI;
}

void ProposedRoad::setSplitRoad(int i) {
	splitRoad = 1;
	splitRoadIndex = i;
}