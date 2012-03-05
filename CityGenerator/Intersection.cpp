#include "StdAfx.h"
#include "Intersection.h"

Intersection::Intersection(void) {
	position[0] = 0.0;
	position[1] = 0.0;
	position[2] = 0.0;
}

Intersection::~Intersection(void) {
	roads.clear();
}

Intersection::Intersection(Vector3d &pos) {
	for(unsigned i = 0; i < 3; i++) {
		position[i] = pos[i];
	}
}

Intersection::Intersection(double x, double y) {
	//position[0] = x;
	//position[1] = y;
	//position[2] = 0.0;
	position = Vector3d(x, y, 0);
	//printf("Intersection: init:\n");
	//printf("              [%f, %f, %f]\n",position[0],position[1],position[2]);
}

void Intersection::getRoads(vector<int> &r) const {
	r = roads;
}

void Intersection::addRoad(int r) {
	roads.push_back(r);
}

void Intersection::getPosition(double &X, double &Y, double &Z) const {
	//printf("Intersection: position: [%f, %f, %f]\n", position[0],position[1],position[2]);
	X = position[0];
	Y = position[1];
	Z = position[2];
}

double Intersection::operator[] (unsigned int d) const {
	return position[d];
}

bool Intersection::isEndpoint() const {
	return roads.size() <= 1;
}

//removes r from roads, if it's there
//returns true if r was successfully removed, false otherwise
//i.e. if r isn't in roads
bool Intersection::removeRoad(int r) {
	for(unsigned i = 0; i < roads.size(); i++) {
		if(roads[i] == r) {
			roads.erase(roads.begin()+i);
			return true;
		}
	}
	return false;
}