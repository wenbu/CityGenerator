#pragma once
#include "Intersection.h"
#include <GL/glut.h>
#include <GL/glu.h>
#include <set>

class QuadPoint {
public:
	//Intersection* intersection;
	Vector3d position;
	int index;
	QuadPoint(Intersection* isecn, int i) {position = Vector3d(isecn->position[0],isecn->position[1],isecn->position[2]); index = i;}
	QuadPoint(Vector3d &v, int i) {position = Vector3d(v); index = i;}
	QuadPoint(double x, double y, double z, int i) {position = Vector3d(x,y,z); index = i;}
};

class QuadNode {
	static const int MAX_SIZE;

	QuadNode* children[4];
	//set<QuadPoint*> contents;
	vector<QuadPoint*> contents;
	bool isLeaf; //isLeaf == !(children are defined)
	Vector3d SplitPoint;
	Vector3d bounds[2];
	double mostRecentQuery[4];
public:
	string name;

	QuadNode(void);
	~QuadNode(void);
	QuadNode(const Vector3d &minBounds, const Vector3d &maxBounds);
	QuadNode(double x0, double y0, double x1, double y1);

	void insert(Intersection *isecn, int index);
	void insert(double x, double y, int i);
	void remove(double x, double y, int i);
	void query(const Vector3d &qMin, const Vector3d &qMax, vector<QuadPoint*> &result);
	void query(const Vector3d &q0, const Vector3d &q1,
			   const Vector3d &q2, const Vector3d &q3,
			   vector<QuadPoint*> &result);
	void draw() const;
	void reset();
private:
	void setName(string n);
	void insert(QuadPoint* qp);
	void insertIntoChild(Intersection *isecn, int index);
	void insertIntoChild(double x, double y, int i);
	void removeFromChild(double x, double y, int i);
	int getCorrectChild(double x, double y) const;
	bool contains(const Vector3d &v) const;
	bool contained(const Vector3d &v1, const Vector3d &v2, int i) const;
	int intersectQueryBox(const Vector3d &qMin, const Vector3d &qMax) const;
	int intersectQueryBox(const Vector3d &q0, const Vector3d &q1,
						  const Vector3d &q2, const Vector3d &q3) const;
	bool contains2d(const Vector2d &q0, const Vector2d &q1,
					const Vector2d &q2, const Vector2d &q3,
					const Vector2d &p) const;
	bool contains3d(const Vector3d &q0, const Vector3d &q1,
					const Vector3d &q2, const Vector3d &q3,
					const Vector3d &p) const;
	int whichSide(const Vector2d &q0, const Vector2d &q1,
				  const Vector2d &q2, const Vector2d &q3,
				  const Vector2d &v0, const Vector2d &v1) const;
	void dumpContents(vector<QuadPoint*> &result) const;
	void split();
	void drawBox(const Vector3d &min, const Vector3d &max) const;
};

