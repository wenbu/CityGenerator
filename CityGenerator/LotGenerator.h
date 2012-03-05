#pragma once
#include "Lot.h"
#include "Polygon.h"
#include "ImageSampler.h"
#include <math.h>
#include <map>
#include <queue>
#include <deque>
#include <iostream>
#include <vector>

class LotGenerator {
public:
	double maxWidth;
	double maxDeviance;
	double minArea;
	int skippedArea;
	int skippedVertex;
	ImageSampler is;
	vector<Polygon> polygons;
	vector<Lot> lots;

	LotGenerator(void);
	~LotGenerator(void);
	LotGenerator(vector<Polygon> b, double width, double deviance, double area, ImageSampler &img);
	void getLots(vector<Lot> &lots);
	void subdivide(); 
	void subdivideBlock(Polygon &lot, vector<Lot> &lots);

	static bool intersection3d(Vector3d& x1, Vector3d& x2, Vector3d& y1, Vector3d& y2, Vector3d& intersection);
	void splitPolygon(Polygon& p, Vector3d& slStart, Vector3d& slEnd, vector<Polygon>& polys);
	void split(int startIndex, vector<vector<Vector3d>>& polys, 
		vector<Vector3d>& points, 
		vector<Vector3d>& cross, vector<Vector3d>& done);
	double getArea(Polygon &p);
	double getDistance(Vector3d& a, Vector3d& b);	
	double getLongestEdge(Polygon &lot, Vector3d &start, Vector3d &end);
	double getEdgeLength(Vector3d &start, Vector3d &end);
	double getPerpendicularBisectorForLongestSide(Polygon &p, Vector3d& start, Vector3d& end);
	double cross2d(Vector3d &v0, Vector3d &v1);
	bool isOverlapLine(Vector3d& x1, Vector3d& x2, Vector3d& y1, Vector3d& y2);
	bool isOverlapLine(Polygon &lot, Vector3d& x, Vector3d& y);
	bool isOverlapLine(Polygon &lot, Polygon &p);
};

