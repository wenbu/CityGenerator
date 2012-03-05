#include "StdAfx.h"
#include "QuadNode.h"

using namespace std;

const int QuadNode::MAX_SIZE = 20;

QuadNode::QuadNode(void) {
	for(unsigned i = 0; i < 4; i++)
		children[i] = 0;
	isLeaf = 1;
	bounds[0] = Vector3d(-1, -1, 0);
	bounds[1] = Vector3d(1, 1, 0);
	name = ".";
}

QuadNode::~QuadNode(void){
	for(unsigned i = 0; i < 4; i++)
		delete children[i];
	contents.clear();
}

void QuadNode::reset() {
	for(unsigned i = 0; i < 4; i++)
		delete children[i];
	contents.clear();
}

QuadNode::QuadNode(const Vector3d &min, const Vector3d &max) {
	for(unsigned i = 0; i < 4; i++)
		children[i] = 0;
	isLeaf = 1;
	bounds[0] = Vector3d(min);
	bounds[1] = Vector3d(max);
	name = ".";
}

QuadNode::QuadNode(double x0, double y0, double x1, double y1) {
	for(unsigned i = 0; i < 4; i++)
		children[i] = 0;
	isLeaf = 1;
	bounds[0] = Vector3d(x0, y0, 0);
	bounds[1] = Vector3d(x1, y1, 0);
	name = ".";
}

void QuadNode::insert(QuadPoint *qp) {
	if(!isLeaf) {
		printf("QuadNode: tried to insert QuadPoint into non-leaf\n");
		exit(1);
	}
	if(contents.size() >= MAX_SIZE) {
		printf("QuadNode: tried to insert QuadPoint into full node\n");
		exit(1);
	}
	contents.push_back(qp);
}

void QuadNode::insert(Intersection *isecn, int i) {
	//double x, y, z;
	//isecn->getPosition(x,y,z);
	//printf("QuadNode: inserting [%f,%f]\n",x,y);
	//printf("QuadNode: inserting with index %i\n", i);
	if(isLeaf) {
		if(contents.size() < MAX_SIZE) {
			//case 1: node is leaf and has room: insert directly
			QuadPoint *qp =  new QuadPoint(isecn, i);
			double x = qp->position[0];
			double y = qp->position[1];
			double z = qp->position[2];
			//printf("QuadNode: inserting:\n");
			//printf("                 qp: %p\n", qp);
			//printf("          [%f, %f]\n", x,y);
			contents.push_back(qp);
			//printf("QuadNode: inserted:\n");
			//printf("          [%f, %f]\n", contents[contents.size()-1]->position[0],contents[contents.size()-1]->position[1]);
		} else {
			//case 2: node is leaf and has no room: split & insert into child
			split();
			insertIntoChild(isecn, i);
		}
	} else {
		//case 3: node is not leaf: insert recursively into child
		insertIntoChild(isecn, i);
	}
}

void QuadNode::insert(double x, double y, int i) {
	if(isLeaf) {
		if(contents.size() < MAX_SIZE) {
			QuadPoint *qp = new QuadPoint(x, y, 0, i);
			contents.push_back(qp);
		} else {
			split();
			insertIntoChild(x, y, i);
		}
	} else
		insertIntoChild(x, y, i);
}

void QuadNode::remove(double x, double y, int i) {
	if(isLeaf) {
		//scan through contents, look for quadpoint w/ index i
		for(unsigned j = 0; j < contents.size(); j++) {
			if(contents[j]->index == i) {
				contents.erase(contents.begin() + i);
				return;
			}
		}
		printf("QuadNode: cannot delete nonexistent point\n");
	} else
		removeFromChild(x, y, i);
}

void QuadNode::insertIntoChild(Intersection *isecn, int i) {
	double x, y, z;
	isecn->getPosition(x, y, z);
	int j = getCorrectChild(x, y);
	if(j >= 0) {
		//printf("QuadNode: insertIntoChild: %i -> [%i]\n",i,j);
		children[j]->insert(isecn, i);
	} else {
		printf("QuadNode: tried to get child of leaf\n");
		exit(1);
	}
}

void QuadNode::insertIntoChild(double x, double y, int i) {
	int j = getCorrectChild(x, y);
	if(j >= 0)
		children[j]->insert(x, y, i);
	else {
		printf("QuadNode: tried to get child of leaf\n");
		exit(1);
	}
}

void QuadNode::removeFromChild(double x, double y, int i) {
	int j = getCorrectChild(x, y);
	if(j >= 0)
		children[j]->remove(x, y, i);
	else {
		printf("QuadNode: tried to get child of leaf\n");
		exit(1);
	}
}

void QuadNode::query(const Vector3d &qMin, const Vector3d &qMax, vector<QuadPoint*> &result) {
	mostRecentQuery[0] = qMin[0];
	mostRecentQuery[1] = qMin[1];
	mostRecentQuery[2] = qMax[0];
	mostRecentQuery[3] = qMax[1];
	int i = intersectQueryBox(qMin, qMax);
	//disable optimization for now
	i = 1;
	//printf("QuadNode::query: case %i\n", i);
	switch(i) {
	case 0:
		//this node is not in query box; skip
		printf("QuadNode: query %s: skipping\n", name);
		return;
	case 1:
		//this node is partially in query box
		if(isLeaf) {
			//check each datum for containment
			double qx0 = qMin[0], qy0 = qMin[1];
			double qx1 = qMax[0], qy1 = qMax[1];
			for(unsigned j = 0; j < contents.size(); j++) {
					double x = contents[j]->position[0];
					double y = contents[j]->position[1];
					//printf("QuadNode: query: considering %i\n",contents[j]->index);
					if(qx0 < x && x < qx1 && qy0 < y && y < qy1) {
						//point is in query box
						//printf("                   accepted\n");
						result.push_back(contents[j]);
					}
			}
			return;
		} else {
			//query each child
			for(unsigned j = 0; j < 4; j++)
				children[j]->query(qMin, qMax, result);
		}
		break;
	case 2:
		//this node is covered by query box
		//dump contents
		dumpContents(result);
		break;
	default:
		printf("QuadNode: got illegal intersectQueryBox result: %i\n", i);
		exit(1);
	}
}

void QuadNode::query(const Vector3d &q0, const Vector3d &q1,
					 const Vector3d &q2, const Vector3d &q3,
					 vector<QuadPoint*> &result) {
	//int i = intersectQueryBox(q0, q1, q2, q3);
	//printf("QuadNode::query: case %i\n", i);
	int i = 1;
	switch(i) {
	case 0:
		//this node is not in query box; skip
		//printf("QuadNode: query %s: skipping\n", name);
		return;
	case 1:
		//this node is partially in query box
		
		if(isLeaf) {
			//check each datum for containment
			for(unsigned j = 0; j < contents.size(); j++) {
					// borrowed from http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
					bool contained = contains3d(q0, q1, q2, q3, contents[j]->position);
					if(contained) {
						//point is in query box
						result.push_back(contents[j]);
					}
			}
			return;
		} else {
			//query each child
			for(unsigned j = 0; j < 4; j++)
				children[j]->query(q0, q1, q2, q3, result);
		}
		break;
		
	case 2:
		//this node is covered by query box
		//dump contents
		dumpContents(result);
		break;
	default:
		printf("QuadNode: got illegal intersectQueryBox result: %i\n", i);
		exit(1);
	}
}

void QuadNode::split() {
	if(!isLeaf) {
		printf("QuadNode:: tried to split non-leaf\n");
		exit(1);
	}
	isLeaf = 0;
	//get split point
	double Sx = 0, Sy = 0;
	//set<QuadPoint*>::iterator it;
	//unsigned i;
	//for(it = contents.begin(), i = 0; it != contents.end(); it++, i++) {
	for(unsigned i = 0; i < contents.size(); i++) {
		//double x, y, z;
		//QuadPoint* qp = (*it);
		//Intersection* j = qp->intersection;
		//j->getPosition(x, y, z);
		double x = contents[i]->position[0];
		double y = contents[i]->position[1];
		//calculate mean position of points iteratively
		Sx += (1.0/(i+1)) * (x-Sx);
		Sy += (1.0/(i+1)) * (y-Sy);
	}
	//printf("QuadNode: splitting at [%f, %f]\n",Sx,Sy);
	SplitPoint = Vector3d(Sx, Sy, 0);

	//initialize children
	double x0 = bounds[0][0], y0 = bounds[0][1];
	double x1 = bounds[1][0], y1 = bounds[1][1];
	children[0] = new QuadNode(Sx, Sy, x1, y1);
	children[0]->setName(name+"0");
	children[1] = new QuadNode(x0, Sy, Sx, y1);
	children[1]->setName(name+"1");
	children[2] = new QuadNode(x0, y0, Sx, Sy);
	children[2]->setName(name+"2");
	children[3] = new QuadNode(Sx, y0, x1, Sy);
	children[3]->setName(name+"3");

	//need to sort nodes into children
	for(unsigned i = 0; i < contents.size(); i++) {
		double x = contents[i]->position[0];
		double y = contents[i]->position[1];
		int j = getCorrectChild(x,y);
		//printf("QuadNode: split: %i -> [%i]\n",contents[i]->index,j);
		children[j]->insert(contents[i]);
	}
}

//returns the index of the child that (x,y) belongs in
int QuadNode::getCorrectChild(double x, double y) const {
	if(isLeaf)
		return -1;
	//return index of correct child
	double Sx = SplitPoint[0];
	double Sy = SplitPoint[1];

	if(y >= Sy) {
		if(x >= Sx) return 0;
		else return 1;
	} else {
		if(x >= Sx) return 3;
		else return 2;
	}
}

//checks whether v is within this node's bounds
inline bool QuadNode::contains(const Vector3d &v) const {
	bool result = v[0] > bounds[0][0] && v[1] > bounds[0][1] && v[0] < bounds[1][0] && v[1] < bounds[1][1];
	return result;
}

//checks whether bounds[i] is contained within the axis-aligned box defined
//by v1, v2
inline bool QuadNode::contained(const Vector3d &v1, const Vector3d &v2, int i) const {
	double x = bounds[i][0];
	double y = bounds[i][1];
	bool result = v1[0] < x && x < v2[0] && v1[1] < y && y < v2[1];
	return result;
}

//returns 0 for no intersection with query rectangle, 1 for partial intersection
//with query rectangle, 2 for complete encapsulation within query rectangle
int QuadNode::intersectQueryBox(const Vector3d &qMin, const Vector3d &qMax) const {
	int result = 0;
	if(contained(qMin, qMax, 0)) result += 1;
	if(contained(qMin, qMax, 1)) result += 1;
	if(result == 0 && contains(qMin)) result = 1;
	return result;
}

int QuadNode::intersectQueryBox(const Vector3d &q0, const Vector3d &q1, 
								const Vector3d &q2, const Vector3d &q3) const {
	Vector2d c0[4] = { q0.head<2>(), q1.head<2>(),
					   q2.head<2>(), q3.head<2>() };
	Vector2d c1[4] = { Vector2d(bounds[0][0], bounds[0][1]),
					   Vector2d(bounds[1][0], bounds[0][1]),
					   Vector2d(bounds[1][0], bounds[1][1]),
					   Vector2d(bounds[0][0], bounds[1][1]) };
	int numContained1;
	for(unsigned i = 0; i < 4; i++) {
		if(contains2d(c0[0], c0[1], c0[2], c0[3], c1[i]))
			numContained1++;
	}
	int numContained2;
	for(unsigned i = 0; i < 4; i++) {
		if(contains2d(c1[0], c1[1], c1[2], c1[3], c0[i]))
			numContained2++;
	}

	if(numContained1 == 4 || numContained2 == 4) {
		//total encapsulation
		return 2;
	}
	else if(numContained1 != 0 || numContained2 != 0) {
		//partial containment
		return 1;
	}
	else {
		//test query
		unsigned i0, i1;
		for(i0 = 0, i1 = 3; i0 < 4; i1 - i0, i0++) {
			Vector2d d = c0[i0] - c0[i1];
			Vector2d D = Vector2d(-d[1], d[0]);
			if(whichSide(c1[0], c1[1], c1[2], c1[3], D, c0[i0]) > 0)
				return 0;
		}

		//test bounding
		for(i0 = 0, i1 = 3; i0 < 4; i1 = i0, i0++) {
			Vector2d d = c1[i0] - c1[i1];
			Vector2d D = Vector2d(-d[1], d[0]);
			if(whichSide(c0[0], c0[1], c0[2], c0[3], D, c1[i0]) > 0)
				return 0;
		}
		return 1;
	}
}

bool QuadNode::contains2d(const Vector2d &q0, const Vector2d &q1,
						  const Vector2d &q2, const Vector2d &q3,
						  const Vector2d &p) const {
	bool c = 0;
	double x = p[0];
	double y = p[1];
	double qx[4] = { q0[0], q1[0], q2[0], q3[0] };
	double qy[4] = { q0[1], q1[1], q2[1], q3[1] };
	unsigned i, j;
	for (i = 0, j = 3; i < 4; j = i++) {
		if ( ((qy[i] > y) != (qy[j] > y)) &&
			(x < (qx[j]-qx[i]) * (y-qy[i]) / (qy[j]-qy[i]) + qx[i]) )
			c = !c;
	}
	return c;
}

bool QuadNode::contains3d(const Vector3d &q0, const Vector3d &q1,
						  const Vector3d &q2, const Vector3d &q3,
						  const Vector3d &p) const {
	return contains2d(q0.head<2>(), q1.head<2>(), q2.head<2>(),
					  q3.head<2>(), p.head<2>());
}

//returns 0 if the polygon defined by vn straddles the line defined
//by p + t*d
int QuadNode::whichSide(const Vector2d &v0, const Vector2d &v1,
						const Vector2d &v2, const Vector2d &v3,
						const Vector2d &p, const Vector2d &d) const {
	int positive = 0, negative = 0;
	Vector2d S[4] = { v0, v1, v2, v3 };
	Vector2d D = d.head<2>();
	Vector2d P = p.head<2>();
	for(unsigned i = 0; i < 3; i++) {
		double t = D.dot(S[i] - P);
		if(t > 0) positive++;
		else if(t < 0) negative++;
		if(positive && negative) return 0;
	}
	return (positive? 1 : -1);
}

void QuadNode::dumpContents(vector<QuadPoint*> &v) const {
	if(!isLeaf) {
		for(unsigned i = 0; i < 4; i++)
			children[i]->dumpContents(v);
	} else {
		for(unsigned i = 0; i < contents.size(); i++) {
				//printf("QuadNode: query: dumping %i\n",contents[i]->index);
				v.push_back(contents[i]);
		}
	}
}

void QuadNode::drawBox(const Vector3d &min, const Vector3d &max) const {
	printf("QuadNode::drawBox\n");
	glColor3f(0.7, 0.3, 0.3);
	glBegin(GL_QUADS);
		glVertex3f(min[0], min[1], 0);
		glVertex3f(max[0], min[1], 0);
		glVertex3f(max[0], max[1], 0);
		glVertex3f(max[0], min[1], 0);
	glEnd();
}

void QuadNode::setName(string s) {
	name = s;
}

void QuadNode::draw() const {
	//draw bounds
	glBegin(GL_LINE_STRIP);
		glColor3f(0.2,0.2,0.2);
		glVertex3f(bounds[0][0],bounds[0][1],0);
		glVertex3f(bounds[1][0],bounds[0][1],0);
		glVertex3f(bounds[1][0],bounds[1][1],0);
		glVertex3f(bounds[0][0],bounds[1][1],0);
	glEnd();
	if(name == ".") {
		//draw most recent query
		glBegin(GL_QUADS);
			glColor4f(1,0,0,0.1);
			glVertex3f(mostRecentQuery[0],mostRecentQuery[1],0);
			glVertex3f(mostRecentQuery[2],mostRecentQuery[1],0);
			glVertex3f(mostRecentQuery[2],mostRecentQuery[3],0);
			glVertex3f(mostRecentQuery[0],mostRecentQuery[3],0);
		glEnd();
	}
	if(!isLeaf) {
		for(unsigned i = 0; i < 4; i++) {
			children[i]->draw();
		}
	}
}