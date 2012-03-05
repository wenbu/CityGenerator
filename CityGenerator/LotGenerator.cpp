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
	for(unsigned i = 0; i < polygons.size(); i++) {
		//printf("subdividing polygon %i\n", i);
		subdivideBlock(polygons[i],lots);
	}
}

void LotGenerator::subdivideBlock(Polygon &lot, vector<Lot> &lots) {
	deque<Polygon> q;
	q.push_back(lot);
	double density = is.getDensity(lot.vertices[0][0], lot.vertices[0][1])/255;
	density = (density/2)+0.5;
	double width = density*maxWidth;
	while(!q.empty()) {
		//printf("loop\n");
		Polygon block = q.front();
		//printf("\tblock: %i verts\n", block.vertices.size());
		//if(block.vertices.size() == 7)
			//printf("!!!\n");
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
		double longestEdge = getPerpendicularBisectorForLongestSide(block,start,end);

		if(longestEdge <= width && area >= minArea) {
			//printf("if 1\n");
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
			//printf("else 1\n");
			vector<Polygon> polys;
			splitPolygon(block,start,end,polys);
			for(unsigned i = 0; i < polys.size(); i++) {
				/**
				printf(">>polygon %d\n",i);
				vector<Vector3d> vertices = polys[i].vertices;
				for(int j = 0; j < vertices.size(); j++) {
					cout << vertices[j][0] << "," <<  vertices[j][1] << "," <<  vertices[j][2] << endl;
				}
				*/
				if(isOverlapLine(lot,polys[i]))
					q.push_back(polys[i]);
			}
		}
		q.pop_front();
	}
}

double LotGenerator::getPerpendicularBisectorForLongestSide(Polygon &p, Vector3d& start, Vector3d& end) {
	//printf(":::getting bisector\n");
	Vector3d lStart, lEnd;
	double distance = getLongestEdge(p,lStart,lEnd);
	//cout << "longestEdge start: " << lStart[0] << " " << lStart[1] << " " << lStart[2] << ", l=" << distance << endl;
	//cout << "longestEdge end: " << lEnd[0] << " " << lEnd[1] << " " << lEnd[2] << ", l=" << distance << endl;

	double deviance = maxDeviance * ((rand() % 10) / 10.0);
	//cout << deviance << endl;
	Vector3d dirVector = lEnd - lStart;
	Vector3d midPt = lStart + (0.5 + deviance) * dirVector;
	//cout << "mid dirVector: " << dirVector[0] << " " << dirVector[1] << " " << dirVector[2] << endl;
	//cout << "mid point: " << midPt[0] << " " << midPt[1] << " " << midPt[2] << endl;

	Vector3d x1 = p.vertices[1] - p.vertices[0];
	Vector3d x2 = p.vertices[2] - p.vertices[0];
	Vector3d x = x1.cross(x2);
	Vector3d dirVector2 = x.cross(dirVector);
	//cout << "dir vector: " << dirVector2[0] << " " << dirVector2[1] << " " << dirVector2[2] << endl;

	start = midPt + -10*dirVector2;
	end = midPt + 10*dirVector2;
	//cout << "start vector: " << start[0] << " " << start[1] << " " << start[2] << endl;
	//cout << "end vector: " << end[0] << " " << end[1] << " " << end[2] << endl;

	return distance;
}

double LotGenerator::getLongestEdge(Polygon &lot, Vector3d &start, Vector3d &end) {
	//printf(":::getting longest edge\n");
	if(lot.vertices.size() <= 2) return -1;
	double maxlength = -1;
	//printf("---------------------------------\n");
	//printf("LotGenerator: num vertices: %i\n", lot.vertices.size());
	for(unsigned i = 0; i < lot.vertices.size() - 1; i++) {
		double length = getEdgeLength(lot.vertices[i],lot.vertices[i+1]);
		//printf("LotGenerator::getLongestEdge: length = %f\n", length);
		if(length > maxlength) {
			start = lot.vertices[i];
			end = lot.vertices [i+1];
			maxlength = length;
		}
	}
	if(maxlength == -1) {
		//printf("LotGenerator::getLongestEdge, edge not found\n");
	}
	//printf("LotGenerator::getLongestEdge: %f\n", maxlength);
	return maxlength;
}

//http://alienryderflex.com/polygon_area/
double LotGenerator::getArea(Polygon &p) {
	//printf(":::getting area\n");
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
void LotGenerator::splitPolygon(Polygon& p, Vector3d& slStart, Vector3d& slEnd, vector<Polygon>& results) {
	//printf(":::splitting polygon\n");
	vector<Vector3d> vertices = p.vertices;
	vector<Vector3d> points;
	vector<Vector3d> crossPts;

	//cout << "points" << endl;
	points.push_back(vertices[0]);
	//cout << "n: " << vertices[0][0] << "," << vertices[0][1] << "," << vertices[0][2] << endl;
	for(unsigned i = 0; i < vertices.size(); i++) {
		//printf(":::loop\n");
		
		/*
		printf(":::vertices: [\n");
		for(unsigned j = 0; j < vertices.size(); j++)
			cout << ":::   " << j << ". [" << vertices[j][0] << ", " << vertices[j][1] << ", " << vertices[j][2] << endl;
		printf(":::]\n");
		cout << ":::start = " << i << endl;
		cout << ":::end   = " << ((i+1)%vertices.size()) << endl;
		*/

		int startI = i;
		int endI = (i+1) % vertices.size();
		//if(endI < startI) printf("lol i might crash now\n");
		Vector3d& start = vertices[i];
		//printf(":::got start\n");
		Vector3d& end = vertices[(i + 1) % vertices.size()];
		//printf(":::got end\n");
		Vector3d intersection;
		bool intersects = intersection3d(start, end, slStart, slEnd, intersection);
		//printf(":::got intersection\n");
		if(intersects) {
			//printf("if 1\n");
			//cout << "i: " << intersection[0] << "," << intersection[1] << "," << intersection[2] << endl;
			if(end[0] == intersection[0] && end[1] == intersection[1] && end[2] == intersection[2]) {
				//printf("if 1, 1\n");
				//cout << "lol" << endl;
			} else {
				//printf("if 1, 0\n");
				crossPts.push_back(intersection);
				if(find(points.begin(),points.end(),intersection) == points.end()) {
					//cout << "in" << endl;
					points.push_back(intersection);
				}
			}
		}
		if((i + 1) % vertices.size() != 0) {
			//printf("if 2\n");
			//cout << "n: " << end[0] << "," << end[1] << "," << end[2] << endl;
			points.push_back(end);
		}
	}
	//printf(":::pq\n");

	typedef priority_queue<Vector3d, vector<Vector3d>, CompareDistanceFromPt> my_pq;
	my_pq pq(CompareDistanceFromPt(slStart[0], slStart[1], slStart[2]));

	//cout << "cuts" << endl;
	for(unsigned i = 0; i < crossPts.size(); i++) {
		pq.push(crossPts[i]);
	}
	//printf(":::pushed to pq\n");
	crossPts.clear();
	while(!pq.empty()) {
		crossPts.push_back(pq.top());
		//cout << "c: " << pq.top()[0] << "," << pq.top()[1] << "," << pq.top()[2] << endl;
		pq.pop();
	}
	//printf(":::pq crosspts pushpop done\n");
	//crossPtsIds now sorts crossPts by closest distance to start pt of split line
	//cout << "hi" << endl;
	vector<vector<Vector3d>> splitPolys;
	vector<Vector3d> done;
	//printf(":::initial split call\n");
	split(0,splitPolys,points,crossPts,done);
	//printf(":::split done");

	for(unsigned i = 0; i < splitPolys.size(); i++) {
		Polygon p;
		p.vertices = splitPolys[i];
		results.push_back(p);
	}

	/**
	for(int i = 0; i < results.size(); i++) {
		printf("polygon %d\n",i);
		vector<Vector3d> vertices = results[i].vertices;
		for(int j = 0; j < vertices.size(); j++) {
			cout << vertices[j][0] << "," <<  vertices[j][1] << "," <<  vertices[j][2] << endl;
		}
	}
	*/
}

//couldn't get right values from *crossBackItr, somehow pointing to wrong thing, hence var goal to store value
void LotGenerator::split(int startIndex, vector<vector<Vector3d>>& polys, vector<Vector3d>& points, 
	vector<Vector3d>& cross, vector<Vector3d>& done) {
	//cout << "::::::split points size: " << points.size() << endl;
	//printf("::::::split\n");
	vector<Vector3d> cut;
	int current = startIndex;
	vector<Vector3d>::iterator crossBackItr = cross.end();
	Vector3d goal;
	vector<Vector3d> result;
	result.push_back(points[current]);
	current = (current+1) % points.size();
	while(current != startIndex) {
		//printf("::::::   loop 1\n");
		Vector3d& currentPt = points[current];
		vector<Vector3d>::iterator itr = find(cross.begin(),cross.end(),currentPt);
		bool isCrossPt = itr != cross.end();
		if(crossBackItr == cross.end() && !isCrossPt) {
			//printf("::::::   if 1\n");
			result.push_back(currentPt);
		} else if(crossBackItr == cross.end() && isCrossPt) {
			//printf("::::::   elif 2\n");
			result.push_back(currentPt);
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
			//printf("::::::   elif 3\n");
			result.push_back(currentPt);
			crossBackItr = cross.end();
		}
		current = (current + 1) % points.size();
	}
	polys.push_back(result);
	//printf("::::::loop 1 done\n");
	for(unsigned i = 0; i < cut.size(); i++) {
		//printf("::::::   loop 2\n");
		if(find(done.begin(),done.end(),cut[i]) == done.end()) {
			int d = distance(points.begin(),find(points.begin(),points.end(),cut[i]));
			split(d,polys,points,cross,done);
		}
	}
	//printf("::::::split done\n");
}

//http://mathforum.org/library/drmath/view/62814.html
bool LotGenerator::intersection3d(Vector3d& x1, Vector3d& x2, 
	Vector3d& y1, Vector3d& y2, Vector3d& intersection) {
	Vector3d dirVectX = x2 - x1;
	Vector3d dirVectY = y2 - y1;
	Vector3d tmp1 = dirVectX.cross(dirVectY);
	Vector3d tmp2 = (y1 - x1).cross(dirVectY);
	double lSq1 = tmp1[0]*tmp1[0] + tmp1[1]*tmp1[1] + tmp1[2]*tmp1[2];
	double lSq2 = tmp2[0]*tmp2[0] + tmp2[1]*tmp2[1] + tmp2[2]*tmp2[2];
	double t = sqrt(lSq2/lSq1);
	bool signChange1 = (tmp1[0] > 0 && tmp2[0] < 0) || (tmp1[0] < 0 && tmp2[0] > 0);
	bool signChange2 = (tmp1[1] > 0 && tmp2[1] < 0) || (tmp1[1] < 0 && tmp2[0] > 0);
	bool signChange3 = (tmp1[2] > 0 && tmp2[2] < 0) || (tmp1[2] < 0 && tmp2[0] > 0);
	bool signChange = signChange1 || signChange2 || signChange3;
	if(signChange) {
		t = -t;
	}
	/**
	cout << dirVectX[0] << "," << dirVectX[1] << "," << dirVectX[2] << endl;
	cout << dirVectY[0] << "," << dirVectY[1] << "," << dirVectY[2] << endl;
	cout << tmp1[0] << "," << tmp1[1] << "," << tmp1[2] << endl;
	cout << tmp2[0] << "," << tmp2[1] << "," << tmp2[2] << endl;
	cout << lSq1 << endl;
	cout << lSq2 << endl;
	cout << t << endl;
	*/
	if(t < 0 || t > 1) {
		return false;
	} else {
		intersection = x1 + t*dirVectX;
		return true;
	}
}

bool LotGenerator::isOverlapLine(Polygon &lot, Vector3d& x, Vector3d& y) {
	for(unsigned i = 0; i < lot.vertices.size(); i++) {
		Vector3d s = lot.vertices[i];
		Vector3d e = lot.vertices[(i+1)%lot.vertices.size()];
		//cout << s[0] << " " << s[1] << " " << s[2] << "->" << e[0] << " " << e[1] << " " << e[2] << endl;
		if(isOverlapLine(x,y,s,e)) {
			return true;
		}
	}
	return false;
}

bool LotGenerator::isOverlapLine(Vector3d& x1, Vector3d& x2, Vector3d& y1, Vector3d& y2) {
	//Vector3d i;
	//bool first = intersection3d(x1,x1,y1,y2,i);
	//bool second = intersection3d(x2,x2,y1,y2,i);
	//cout << first << " " << second << endl;
	//return first && second;
	Vector3d r = x2 - x1;
	Vector3d q = y1 - x1;
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