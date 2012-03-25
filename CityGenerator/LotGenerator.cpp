#include "StdAfx.h"
#include "LotGenerator.h"
using namespace std;

LotGenerator::LotGenerator(void) {
	vector<Polygon> p;
	polygons = p;
	maxWidth = 0;
	maxDeviance = 0;
	minArea = 0;
}

LotGenerator::~LotGenerator(void) {
}

LotGenerator::LotGenerator(vector<Polygon> b, double width, 
						   double deviance, double area,
						   ImageSampler &img) {
	polygons = b;
	maxWidth = width;
	if(deviance > 0.25 || deviance < 0.0) {
		deviance = 0.0;
	}
	maxDeviance = deviance;
	minArea = area;
	skippedArea = 0;
	skippedVertex = 0;
	is = img;
}

void LotGenerator::getLots(vector<Lot> &l) {
	//lots.clear();
	subdivide();
	printf("got %i lots\n", lots.size());
	for(unsigned i = 0; i < lots.size(); i++) {
		l.push_back(lots[i]);
	}
}

void LotGenerator::subdivide() {
	//printf("lg: polygons.size: %i\n", polygons.size());
	int numBlocks = polygons.size();
	for(unsigned i = 0; i < numBlocks; i++) {
		//printf("subdividing polygon %i\n", i);
		subdivideBlock(polygons[i],lots);
	}
}

void LotGenerator::draw() {
	for(unsigned i = 0; i < drawBuffer.size(); i++) {
		glColor3f(1, 0.8, 0.8);
		glBegin(GL_POLYGON);
		for(unsigned j = 0; j < drawBuffer[i].vertices.size(); j++) {
			glVertex3f(drawBuffer[i].vertices[j][0],
					   drawBuffer[i].vertices[j][1],
					   drawBuffer[i].vertices[j][2]);
		}
		glEnd();
		glColor3f(0, 0, 0);
		for(unsigned j = 0; j < drawBuffer[i].vertices.size(); j++) {
			glBegin(GL_LINE_STRIP);
			glVertex3f(drawBuffer[i].vertices[j][0],
					   drawBuffer[i].vertices[j][1],
					   drawBuffer[i].vertices[j][2]);
			if(i != drawBuffer[i].vertices.size()-1) {
				glVertex3f(drawBuffer[i].vertices[j+1][0],
						   drawBuffer[i].vertices[j+1][1],
						   drawBuffer[i].vertices[j+1][2]);
			} else {
				glVertex3f(drawBuffer[i].vertices[0][0],
						   drawBuffer[i].vertices[0][1],
						   drawBuffer[i].vertices[0][2]);
			}
			glEnd();
		}
	}
	glColor3f(1, 0, 0);
	for(unsigned i = 0; i < edgeBuffer.size(); i += 2) {
		glBegin(GL_LINE_STRIP);
		glVertex3f(edgeBuffer[i][0], edgeBuffer[i][1], 0);
		glVertex3f(edgeBuffer[i+1][0], edgeBuffer[i+1][1], 0);
		glEnd();
	}
}

void LotGenerator::subdivideBlock(Polygon &lot, vector<Lot> &lots) {
	deque<Polygon> q;
	q.push_back(lot);
	double density = is.getDensity(lot.vertices[0][0], lot.vertices[0][1])/255;
	density = (density/4)+0.75;
	double width = density*maxWidth;

	//iterate through list of polygons
	while(!q.empty()) {
		Polygon block = q.front();

		//checks: try to skip degenerate polygons
		if(block.vertices.size() <= 2) {
			++skippedVertex;
			q.pop_front();
			continue;
		}

		double area = getArea(block);

		if(area <= 1e-20) {
			++skippedArea;
			q.pop_front();
			continue;
		}

		Vector3d start,end;

		//get longest edge & perp bisector
		double longestEdge = getPerpendicularBisectorForLongestSide(block,start,end);

		if(longestEdge <= width && area >= minArea) {
			//done splitting this lot, put in lots (if meets criteria)
			if(isOverlapLine(lot, block))
				lots.push_back(block);
			else {
				printf("LOT:\n");
				for(unsigned i = 0; i < lot.vertices.size(); i++)
					printf("\t[%f, %f]\n", lot.vertices[i][0], lot.vertices[i][1]);
				printf("BLOCK:\n");
				for(unsigned i = 0; i < block.vertices.size(); i++)
					printf("\t[%f, %f]\n", block.vertices[i][0], lot.vertices[i][1]);
			}
		} else if(longestEdge > width && area >= minArea) {
			//lot needs further splitting
			vector<Polygon> polys;
			splitPolygon(block,start,end,polys);
			for(unsigned i = 0; i < polys.size(); i++) {
				if(isOverlapLine(lot,polys[i]))
					q.push_back(polys[i]);
			}
		}
		q.pop_front();
	}
}

double LotGenerator::getPerpendicularBisectorForLongestSide(Polygon &p, Vector3d& start, Vector3d& end) {
	Vector3d lStart, lEnd;
	double distance = getLongestEdge(p,lStart,lEnd);

	double deviance = maxDeviance * ((rand() % 10) / 10.0);
	Vector3d dirVector = (lEnd - lStart).normalized(); //normalized
	Vector3d midPt = lStart + (0.5 * distance * dirVector);

	Vector3d perpVector(-dirVector[1], dirVector[0], 0);
	start = midPt - perpVector;
	end = midPt + perpVector;

	return distance;
}

double LotGenerator::getLongestEdge(Polygon &lot, Vector3d &start, Vector3d &end) {
	if(lot.vertices.size() <= 2) return -1;
	double maxlength = -1;
	for(unsigned i = 0; i < lot.vertices.size(); i++) {
		double length = getEdgeLength(lot.vertices[i],lot.vertices[(i+1)%lot.vertices.size()]);
		if(length > maxlength) {
			start = lot.vertices[i];
			end = lot.vertices [(i+1)%lot.vertices.size()];
			maxlength = length;
		}
	}
	if(maxlength == -1) {
		printf("LotGenerator::getLongestEdge, edge not found\n");
	}
	return maxlength;
}

//http://alienryderflex.com/polygon_area/
double LotGenerator::getArea(Polygon &p) {
	double area = 0;
	int j = p.vertices.size() - 1;
	for(unsigned i = 0; i < p.vertices.size(); i++) {
		area += (p.vertices[i][0] + p.vertices[j][0]) * (p.vertices[j][1] - p.vertices[i][1]);
		j = i;
	}
	return -(area * 0.5);
}

inline double LotGenerator::getEdgeLength(Vector3d &start, Vector3d &end) {
	return (start-end).norm();
}

class CompareDistanceFromPt {
	Vector3d pt;
public:
	CompareDistanceFromPt(const double& a=0, const double& b=0, const double& c=0) {
		pt = Vector3d(a,b,c);
	}

	double getDistance(Vector3d& a, Vector3d& b) {
		return (b-a).norm();
	}

	bool operator() (Vector3d& lhs, Vector3d& rhs) {
		double dlhs = getDistance(pt,lhs);
		double drhs = getDistance(pt,rhs);
		return dlhs > drhs;
	}
};

double LotGenerator::getDistance(Vector3d& a, Vector3d& b) {
	return (b-a).norm();
}

//algo: http://stackoverflow.com/questions/1775457/generate-new-polygons-from-a-cut-polygon-2d
//pq advice: http://comsci.liu.edu/~jrodriguez/cs631sp08/c++priorityqueue.html
//slStart/End: split line defn
//assembles vertex/intersection lists (with intersections inserted into vertex list)
void LotGenerator::splitPolygon(Polygon& p, Vector3d& slStart, Vector3d& slEnd, vector<Polygon>& results) {
	edgeBuffer.push_back(slStart);
	edgeBuffer.push_back(slEnd);
	vector<Vector3d> vertices = p.vertices;
	vector<Vector3d> points;
	vector<Vector3d> crossPts;

	points.push_back(vertices[0]);
	for(unsigned i = 0; i < vertices.size(); i++) {
		//want to go around polygon starting at vert i
		int startI = i;
		int endI = (i+1) % vertices.size();
		Vector3d start = vertices[startI];
		Vector3d end = vertices[endI];
		Vector3d intersection;
		bool intersects = intersection3d(start, end, slStart, slEnd, intersection);
		if(intersects) {
			if(end[0] == intersection[0] && end[1] == intersection[1] && end[2] == intersection[2]) {
			} else {
				crossPts.push_back(intersection);
				//if intersection is not in points, put it in points
				if(find(points.begin(),points.end(),intersection) == points.end()) {
					points.push_back(intersection);
				}
			}
		}
		if((i + 1) % vertices.size() != 0) {
			points.push_back(end);
		}
	}

	//sorts crossPts with a priority queue
	typedef priority_queue<Vector3d, vector<Vector3d>, CompareDistanceFromPt> my_pq;
	my_pq pq(CompareDistanceFromPt(slStart[0], slStart[1], slStart[2]));

	for(unsigned i = 0; i < crossPts.size(); i++) {
		pq.push(crossPts[i]);
	}
	crossPts.clear();
	while(!pq.empty()) {
		crossPts.push_back(pq.top());
		pq.pop();
	}

	//crossPtsIds now sorts crossPts by closest distance to start pt of split line
	vector<vector<Vector3d>> splitPolys;
	vector<Vector3d> done;
	split(0,splitPolys,points,crossPts,done);

	for(unsigned i = 0; i < splitPolys.size(); i++) {
		Polygon p;
		p.vertices = splitPolys[i];
		results.push_back(p);
	}
}

//couldn't get right values from *crossBackItr, somehow pointing to wrong thing, hence var goal to store value
void LotGenerator::split(int startIndex, vector<vector<Vector3d>>& polys, vector<Vector3d>& points, 
	vector<Vector3d>& cross, vector<Vector3d>& done) {

	vector<Vector3d> cut;
	int current = startIndex;
	vector<Vector3d>::iterator crossBackItr = cross.end();
	Vector3d goal;
	vector<Vector3d> result;
	result.push_back(points[current]);
	current = (current+1) % points.size();

	//traverse vertex list
	while(current != startIndex) {
		Vector3d& currentPt = points[current];
		vector<Vector3d>::iterator itr = find(cross.begin(),cross.end(),currentPt);
		bool isCrossPt = itr != cross.end();
		if(crossBackItr == cross.end() && !isCrossPt) {
			//current vert isn't a split point
			result.push_back(currentPt);
		} else if(crossBackItr == cross.end() && isCrossPt) {
			//current vert is a split point
			result.push_back(currentPt);

			//what's the position of this split point in cross?
			int crossIndex = distance(cross.begin(),itr);
			if(crossIndex % 2 == 0) {
				if((itr+1) == cross.end()) {
					crossBackItr = cross.begin();
				} else {
					crossBackItr = itr + 1;
				}
			} else {
				if((itr-1) == cross.end()) {
					crossBackItr = cross.begin();
				} else {
					crossBackItr = itr - 1;
				}
			}
			goal = *crossBackItr;
			cut.push_back(currentPt);
			done.push_back(goal);
		} else if (crossBackItr != cross.end() && goal[0] == currentPt[0] 
			&& goal[1] == currentPt[1] && goal[2] == currentPt[2]) {
			result.push_back(currentPt);
			crossBackItr = cross.end();
		}
		current = (current + 1) % points.size();
	}
	polys.push_back(result);
	for(unsigned i = 0; i < cut.size(); i++) {
		if(find(done.begin(),done.end(),cut[i]) == done.end()) {
			int d = distance(points.begin(),find(points.begin(),points.end(),cut[i]));
			split(d,polys,points,cross,done);
		}
	}
}

bool LotGenerator::intersection3d(Vector3d& x1, Vector3d& x2, 
	Vector3d& y1, Vector3d& y2, Vector3d& intersection) {

	Vector2d c1 = x1.head<2>();
	Vector2d c2 = x2.head<2>();
	Vector2d c3 = y1.head<2>();
	Vector2d c4 = y2.head<2>();

	Vector2d s43 = c4-c3;
	Vector2d s13 = c1-c3;
	Vector2d s21 = c2-c1;
	double numer1 = (s43[0]*s13[1]) - (s43[1]*s13[0]);
	double numer2 = (s21[0]*s13[1]) - (s21[1]*s13[0]);
	double denom  = (s43[1]*s21[0]) - (s43[0]*s21[1]);
	double u1, u2;

	if(denom <= 1e-5 && denom >= -1e-5) {
		return false;
	}

	u1 = numer1/denom;
	u2 = numer2/denom;

	double x = c1[0]+u1*(c2[0]-c1[0]);
	double y = c1[1]+u1*(c2[1]-c1[1]);
	intersection = Vector3d(x, y, 0);
	return 0 <= u1 && u1 <= 1 && 0 <= u2 && u2 <= 1;
}

bool LotGenerator::isOverlapLine(Polygon &lot, Vector3d& x, Vector3d& y) {
	for(unsigned i = 0; i < lot.vertices.size(); i++) {
		Vector3d s = lot.vertices[i];
		Vector3d e = lot.vertices[(i+1)%lot.vertices.size()];
		if(isOverlapLine(x,y,s,e)) {
			return true;
		}
	}
	return false;
}

bool LotGenerator::isOverlapLine(Vector3d& x1, Vector3d& x2, Vector3d& y1, Vector3d& y2) {
	Vector3d r = (x2 - x1).normalized();
	Vector3d q = (y1 - x1).normalized();
	double cross = cross2d(r, q);
	return (cross < 1e-7) && (cross > -1e-7);
}

bool LotGenerator::isOverlapLine(Polygon &lot, Polygon &p) {
	for(unsigned i = 0; i < p.vertices.size(); i++) {
		Vector3d s = p.vertices[i];
		Vector3d e = p.vertices[(i+1)%p.vertices.size()];
		if(isOverlapLine(lot,s,e)) return true;
	}
	return false;
}

inline double LotGenerator::cross2d(Vector3d &v1, Vector3d &v2) {
	return (v1[0]*v2[1]) - (v1[1]*v2[0]);
}