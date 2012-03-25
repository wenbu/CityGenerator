#include "StdAfx.h"
#include "RoadSystem.h"

double RoadSystem::QUERY_DIST = 0.0125; //should be slightly more than half the max road length
double RoadSystem::SNAP_THRESHOLD = 0.0044; //0.0025
double RoadSystem::FINAL_SNAP_THRESHOLD = 0.00125;
double RoadSystem::ROAD_LENGTH_MAJOR = 0.0125; //0.03
double RoadSystem::ROAD_LENGTH_MINOR = 0.0125; //0.018
double RoadSystem::ROAD_LENGTH_HIGHWAY = 0.0125; //0.042
double RoadSystem::ROAD_LENGTH_JITTER = 0; //0.006
double RoadSystem::EVAL_DELAY_HIGHWAY = 1;
double RoadSystem::EVAL_DELAY_MAJOR = 10;
double RoadSystem::EVAL_DELAY_MINOR = 50;
double RoadSystem::CUTOFF_HIGHWAY_GROWTH = 100000;
double RoadSystem::NUM_HIGHWAY_PROBES = 10;
double RoadSystem::HIGHWAY_BRANCH_DELAY = 4;
double RoadSystem::ANGLE_HIGHWAY_PROBES = 0.8; // pi/4
double RoadSystem::HIGHWAY_PROBE_JITTER = 0.3;
double RoadSystem::WEIGHT_TERRAIN_CONTOUR = 0.4; //0.4
double RoadSystem::WEIGHT_TERRAIN_CONTOUR_JITTER = 0.2; //0.2
double RoadSystem::HIGHWAY_PROMOTION_PROBABILITY = 0.1;
double RoadSystem::GRID_DEVIATION_PROBABILITY = 0.1; //0.2    
double RoadSystem::GRID_DEVIATION_AMOUNT = 1.46; //1.46
double RoadSystem::SIDE_ROAD_PROBABILITY = 0.5;
double RoadSystem::COST_THRESHOLD = 0.5;
double RoadSystem::COST_THRESHOLD_JITTER = 0.2;
double RoadSystem::COST_INTERSECTION_OVERLOAD = 0.3;
double RoadSystem::COST_ROAD_SHORTENING = 0; //35
double RoadSystem::COST_ROAD_LENGTHENING = 10; //10
double RoadSystem::BONUS_INTERSECTION_UNDERLOAD = -10; //-0.3


/*-----------------------------------------------------
  EdgeKey
------------------------------------------------------*/

inline RoadSystem::EdgeKey::EdgeKey (int v0, int v1) {
	if(v0 < v1) {
		V[0] = v0;
		V[1] = v1;
	} else {
		V[0] = v1;
		V[1] = v0;
	}
}

inline bool RoadSystem::EdgeKey::operator < (const EdgeKey &key) const {
	if(V[1] < key.V[1]) return true;
	if(V[1] > key.V[1]) return false;
	return V[0] < key.V[0];
}

inline RoadSystem::EdgeKey::operator size_t () const {
	return V[0] | (V[1] << 16);
}

/*-----------------------------------------------------
 Vertex
------------------------------------------------------*/

RoadSystem::Vertex::Vertex(double x, double y, int i) {
	position = Vector3d(x, y, 0);
	index = i;
}

RoadSystem::Vertex::~Vertex() {
	adjacent.clear();
}

void RoadSystem::Vertex::insert(Vertex* adj) {
	adjacent.push_back(adj);
}

void RoadSystem::Vertex::remove(Vertex* adj) {
	int n = (int)adjacent.size();
	for(int i = 0; i < n; ++i) {
		if(adj == adjacent[i]) {
			--n;
			if(i < n)
				adjacent[i] = adjacent[n];
			adjacent.pop_back();
			return;
		}
	}
}

/*-----------------------------------------------------
 ProposedEdge constructors
------------------------------------------------------*/

RoadSystem::ProposedEdge::ProposedEdge(void) {
	t = 0;
	r.type = RT_GRID_MINOR;
	r.angle = 0;
	r.majorAxis = 0;
	r.iter = 1;
	r.branchDelay = 0;
	V[0] = -1; V[1] = -1;
	splitRoad = 0;
}

RoadSystem::ProposedEdge::ProposedEdge(int time, RoadType type,
									   double angle, double majorAxis) {
	t = time;
	r.type = type;
	r.angle = angle;
	r.majorAxis = majorAxis;
	r.iter = 1;
	r.branchDelay = 0;
	V[0] = -1; V[1] = -1;
	splitRoad = 0;
}

RoadSystem::ProposedEdge::ProposedEdge(int time, RoadType type,
									   double angle, double majorAxis,
									   double x0, double y0, double x1, double y1) {
	t = time;
	r.type = type;
	r.angle = angle;
	r.majorAxis = majorAxis;
	r.iter = 1;
	r.branchDelay = 0;
	V[0] = -1; V[1] = -1;
	newVertices[0] = Vector3d(x0, y0, 0);
	newVertices[1] = Vector3d(x1, y1, 0);
	splitRoad = 0;
}

RoadSystem::ProposedEdge::ProposedEdge(int time, RoadType type,
									   double angle, double majorAxis,
									   int v0, double x1, double y1) {
	t = time;
	r.type = type;
	r.angle = angle;
	r.majorAxis = majorAxis;
	r.iter = 1;
	r.branchDelay = 0;
	V[0] = v0; V[1] = -1;
	newVertices[1] = Vector3d(x1, y1, 0);
	splitRoad = 0;
}

RoadSystem::ProposedEdge::ProposedEdge(int time, RoadType type,
									   double angle, double majorAxis,
									   int v0, int v1) {
	t = time;
	r.type = type;
	r.angle = angle;
	r.majorAxis = majorAxis;
	r.iter = 1;
	r.branchDelay = 0;
	V[0] = v0; V[1] = v1;
	splitRoad = 0;
}

/*-----------------------------------------------------
 graph insertion/removal functions
------------------------------------------------------*/

bool RoadSystem::insertVertex(double x, double y, int i) {
	Vertices::iterator iter = vertices.find(i);
	if(iter != vertices.end()) {
		//vertex already exists
		return false;
	}

	Vertex* v = new Vertex(x, y, i);
	vertices[i] = v;
	qtree.insert(x, y, i);
	return true;
}

//don't think we'll need this one
//can't remove from qtree atm anyway
bool RoadSystem::removeVertex(int i) {
	Vertices::iterator iter = vertices.find(i);
	if(iter != vertices.end()) {
		Vertex* v = iter->second;
		if(v->adjacent.size() == 0) {
			vertices.erase(iter);
			delete v;
			return true;
		}
	}
	//vertex doesn't exist
	return false;
}

bool RoadSystem::insertEdge(int i0, int i1, roadAttr r) {
	Vertex* v0 = (Vertex*)getVertex(i0);
	if(!v0) {
		return false;
	}
	
	Vertex* v1 = (Vertex*)getVertex(i1);
	if(!v1) {
		return false;
	}

	EdgeKey e(i0, i1);
	map<EdgeKey, roadAttr>::iterator iter = edges.find(e);
	if(iter == edges.end()) {
		edges[e].type = r.type;
		edges[e].angle = r.angle;
		edges[e].majorAxis = r.majorAxis;
		edges[e].iter = r.iter;
		edges[e].branchDelay = r.branchDelay;
		v0->insert(v1);
		v1->insert(v0);
		return true;
	}
	return false;
}

bool RoadSystem::removeEdge(int i0, int i1) {
	Vertex* v0 = (Vertex*) getVertex(i0);
	if(!v0) return false;

	Vertex* v1 = (Vertex*) getVertex(i1);
	if(!v1) return false;

	EdgeKey e(i0, i1);
	map<EdgeKey, roadAttr>::iterator iter = edges.find(e);
	if(iter != edges.end()) {
		edges.erase(e);
		v0->remove(v1);
		v1->remove(v0);
		return true;
	}
	return false;
}

const RoadSystem::Vertex* RoadSystem::getVertex(int i) const {
	Vertices::const_iterator iter = vertices.find(i);
	return (iter != vertices.end() ? iter->second : 0);
}

//given qp, return vector edges containing EdgeKeys for connected edges
void RoadSystem::getRoads(vector<QuadPoint*> &qp, set<EdgeKey> &edges) const {
	for(unsigned i = 0; i < qp.size(); i++) {
		QuadPoint* q = qp[i];
		int i0 = q->index;
		Vertex* v0 = (Vertex*) getVertex(i0);
		for(unsigned j = 0; j < v0->adjacent.size(); j++) {
			int i1 = v0->adjacent[j]->index;
			edges.insert(EdgeKey(i0, i1));
		}
	}
}

/*-----------------------------------------------------
 various edge util functions
------------------------------------------------------*/

double RoadSystem::getLength(const ProposedEdge &e) const {
	// get intersections
	double x0, y0, x1, y1;
	if(e.V[0] < 0) {
		x0 = e.newVertices[0][0];
		y0 = e.newVertices[0][1];
	} else {
		Vertex* v0 = (Vertex*) getVertex(e.V[0]);
		x0 = v0->position[0];
		y0 = v0->position[1];
	}
	if(e.V[1] < 0) {
		x1 = e.newVertices[1][0];
		y1 = e.newVertices[1][1];
	} else {
		Vertex* v1 = (Vertex*) getVertex(e.V[1]);
		x1 = v1->position[0];
		y1 = v1->position[1];
	}

	return (Vector2d(x0, y0) - Vector2d(x1, y1)).norm();
}

//returns intersection (x,y) of the lines specified by
//v1,v2 and v3,v4; u specifies distance along first
//segment to intersection
//negative u => no intersection
void RoadSystem::intersectRoads(const Vertex* v0, const Vertex* v1,
								const Vertex* v2, const Vertex* v3,
								double &x, double &y, 
								double &u1, double &u2) const {
	Vector2d c1 = v0->position.head<2>();
	Vector2d c2 = v1->position.head<2>();
	Vector2d c3 = v2->position.head<2>();
	Vector2d c4 = v3->position.head<2>();

	Vector2d s43 = c4-c3;
	Vector2d s13 = c1-c3;
	Vector2d s21 = c2-c1;
	double numer1 = (s43[0]*s13[1]) - (s43[1]*s13[0]);
	double numer2 = (s21[0]*s13[1]) - (s21[1]*s13[0]);
	double denom  = (s43[1]*s21[0]) - (s43[0]*s21[1]);

	if(denom <= 1e-5 && denom >= -1e-5) {
		u1 = -1;
		u2 = -1;
		return;
	}

	u1 = numer1/denom;
	u2 = numer2/denom;

	x = c1[0]+u1*(c2[0]-c1[0]);
	y = c1[1]+u1*(c2[1]-c1[1]);
}

/*-----------------------------------------------------
 RoadSystem - constructor/destructor
------------------------------------------------------*/

RoadSystem::RoadSystem(void) {
	minCoord = Vector3d(-10, -10, 0);
	maxCoord = Vector3d(10, 10, 0);
	qtree = QuadNode(minCoord, maxCoord);
	nIterations = 0;
	gotPolys = false;
	axiomVerts = 0;
	is = ImageSampler("density.png", "elevation.png", "legality.png");
	srand(time(NULL));
}

RoadSystem::RoadSystem(ImageSampler &img, InputParser &p) {
	is = img;
	minCoord = Vector3d(-10, -10, 0);
	maxCoord = Vector3d(10, 10, 0);
	qtree = QuadNode(minCoord, maxCoord);
	nIterations = 0;
	gotPolys = false;
	axiomVerts = 0;
	srand(time(NULL));

	QUERY_DIST =						p.get(CityParam::QUERY_DIST);
	SNAP_THRESHOLD =					p.get(CityParam::SNAP_THRESHOLD);
	FINAL_SNAP_THRESHOLD =				p.get(CityParam::FINAL_SNAP_THRESHOLD);
	ROAD_LENGTH_MAJOR =					p.get(CityParam::ROAD_LENGTH_MAJOR);
	ROAD_LENGTH_MINOR =					p.get(CityParam::ROAD_LENGTH_MINOR);
	ROAD_LENGTH_HIGHWAY =				p.get(CityParam::ROAD_LENGTH_HIGHWAY);
	ROAD_LENGTH_JITTER =				p.get(CityParam::ROAD_LENGTH_JITTER);
	EVAL_DELAY_HIGHWAY =				p.get(CityParam::EVAL_DELAY_HIGHWAY);
	EVAL_DELAY_MAJOR =					p.get(CityParam::EVAL_DELAY_MAJOR);
	EVAL_DELAY_MINOR =					p.get(CityParam::EVAL_DELAY_MINOR);
	CUTOFF_HIGHWAY_GROWTH =				p.get(CityParam::CUTOFF_HIGHWAY_GROWTH);
	HIGHWAY_BRANCH_DELAY =				p.get(CityParam::HIGHWAY_BRANCH_DELAY);
	NUM_HIGHWAY_PROBES =				p.get(CityParam::NUM_HIGHWAY_PROBES);
	ANGLE_HIGHWAY_PROBES =				p.get(CityParam::ANGLE_HIGHWAY_PROBES);
	HIGHWAY_PROBE_JITTER =				p.get(CityParam::HIGHWAY_PROBE_JITTER);
	WEIGHT_TERRAIN_CONTOUR =			p.get(CityParam::WEIGHT_TERRAIN_CONTOUR);
	WEIGHT_TERRAIN_CONTOUR_JITTER =		p.get(CityParam::WEIGHT_TERRAIN_CONTOUR_JITTER);
	HIGHWAY_PROMOTION_PROBABILITY =		p.get(CityParam::HIGHWAY_PROMOTION_PROBABILITY);
	GRID_DEVIATION_PROBABILITY =		p.get(CityParam::GRID_DEVIATION_PROBABILITY);
	GRID_DEVIATION_AMOUNT =				p.get(CityParam::GRID_DEVIATION_AMOUNT);
	SIDE_ROAD_PROBABILITY =				p.get(CityParam::SIDE_ROAD_PROBABILITY);
	COST_THRESHOLD =					p.get(CityParam::COST_THRESHOLD);
	COST_THRESHOLD_JITTER =				p.get(CityParam::COST_THRESHOLD_JITTER);
	COST_INTERSECTION_OVERLOAD =		p.get(CityParam::COST_INTERSECTION_OVERLOAD);
	COST_ROAD_SHORTENING =				p.get(CityParam::COST_ROAD_SHORTENING);
	COST_ROAD_LENGTHENING =				p.get(CityParam::COST_ROAD_LENGTHENING);
	BONUS_INTERSECTION_UNDERLOAD =		p.get(CityParam::BONUS_INTERSECTION_UNDERLOAD);
}

RoadSystem::~RoadSystem(void) {}

/*-----------------------------------------------------
 RoadSystem - road generation
------------------------------------------------------*/

void RoadSystem::addAxiom(double x0, double y0, double x1, double y1) {
	//assuming one axiom segment for now
	insertVertex(x0, y0, 0);
	insertVertex(x1, y1, 1);
	
	ProposedEdge e(0, RT_AXIOM, 0, 0, axiomVerts, axiomVerts+1);
	roadQueue.push(e);
	axiomVerts += 2;
}

void RoadSystem::acceptProposal(ProposedEdge &pe, int &i0, int &i1) {
	//get intersections
	if(pe.V[0] >= 0) {
		i0 = pe.V[0];
	} else {
		i0 = vertices.size();
		insertVertex(pe.newVertices[0][0],
					 pe.newVertices[0][1],
					 i0);
	}
	if(pe.V[1] >= 0) {
		i1 = pe.V[1];
	} else {
		//double check quadtree for nearby snap targets
		//set up query box
		Vector3d qMin = Vector3d(pe.newVertices[1][0] - FINAL_SNAP_THRESHOLD,
								 pe.newVertices[1][1] - FINAL_SNAP_THRESHOLD,
								 0);
		Vector3d qMax = Vector3d(pe.newVertices[1][0] + FINAL_SNAP_THRESHOLD,
								 pe.newVertices[1][1] + FINAL_SNAP_THRESHOLD,
								 0);
		vector<QuadPoint*> targets;
		qtree.query(qMin, qMax, targets);
		if(targets.size() > 0 && !pe.splitRoad) {
			//find closest vertex and snap to it
			double dist = 1e20;
			double curDist;
			int closest = -1;
			Vector2d p0 = Vector2d(pe.newVertices[1][0],
									   pe.newVertices[1][1]);
			for(unsigned i = 0; i < targets.size(); i++) {
				Vertex* v = (Vertex*) getVertex(targets[i]->index);
				curDist = (Vector2d(v->position[0], v->position[1]) - p0).norm();
				if(curDist < dist) {
					closest = targets[i]->index;
					dist = curDist;
				}
			}
			i1 = closest;
		} else {
			i1 = vertices.size();
			insertVertex(pe.newVertices[1][0], 
						 pe.newVertices[1][1],
						 i1);
		}
	}
	insertEdge(i0, i1, pe.r);

	//check for split
	if(pe.splitRoad) {
		// need to remove edge [vi, vj] and
		// replace with [vi, vk] and [vk, vj]
		int vi = pe.splitVerts[0];
		int vj = pe.splitVerts[1];
		int vk = i1;

		EdgeKey e(vi, vj);
		bool huh1 = insertEdge(vi, vk, edges[e]);
		bool huh2 = insertEdge(vk, vj, edges[e]);
		removeEdge(vi, vj);
	}
}

void RoadSystem::generate(int iterations) {
 	gotPolys = false;
	for(int i = 0; i < iterations; i++) {
		if(!roadQueue.empty() && iterations > 0) {
			if(DEBUG) {
				debugQuads.clear();
				debugEdges.clear();
				debugVerts.clear();
			}
			ProposedEdge pr = roadQueue.top();
			roadQueue.pop();
			bool accepted = localConstraints(pr);
			if(accepted) {
				int i0, i1;
				pr.r.iter = nIterations+i;
				acceptProposal(pr, i0, i1);
				globalGoals(pr.t, i0, i1);
			}
		}
	}
	printf("Iterations: %i\n",nIterations += iterations);
}

bool RoadSystem::localConstraints(ProposedEdge &e) {
	double cost = 0;
	double origLen = getLength(e);

	//if either point is proposed, do clip/snap checks
	if(e.V[0] < 0 || e.V[1] < 0) {
		vector<QuadPoint*> nearV;
		vector<QuadPoint*> nearE;
		Vertex* vert0 = (Vertex*) getVertex(e.V[0]);
		Vector3d v0 = vert0->position;
		Vector3d v1 = e.newVertices[1];
		Vector3d unitEdge = (v1 - v0).normalized();
		Vector3d perpEdge = Vector3d(-unitEdge[1], unitEdge[0], 0);

		//edge query
		Vector3d qv = 2 * QUERY_DIST * unitEdge;
		Vector3d qvp = 2 * QUERY_DIST * perpEdge;
		Vector3d q0 = v0 - qvp;
		Vector3d q1 = v1 + qv - qvp;
		Vector3d q2 = v1 + qv + qvp;
		Vector3d q3 = v0 + qvp;
		qtree.query(q0, q1, q2, q3, nearE);
		if(DEBUG) {
			//edge pre-clip
			debugEdges.push_back((Vector6d() << v0, 0.3, 0.1, 0.1).finished());
			debugEdges.push_back((Vector6d() << v1, 0.3, 0.1, 0.1).finished());
		}
		clipRoad(e, nearE, cost);

		//get vert query box
		//assume second vert is proposed
		v1 = e.newVertices[1];
		unitEdge = (v1 - v0).normalized();
		perpEdge = Vector3d(-unitEdge[1], unitEdge[0], 0);
		qv = QUERY_DIST * unitEdge;
		qvp = QUERY_DIST * perpEdge;
		q0 = v1 - 0.5 * qv - 0.5 * qvp;
		q1 = v1 + qv - qvp;
		q2 = v1 + qv + qvp;
		q3 = v1 - 0.5 * qv + 0.5 * qvp;
		qtree.query(q0, q1, q2, q3, nearV);

		if(DEBUG) {
			debugQuads.push_back((Vector6d() << q0, 0.5, 0.25, 0).finished());
			debugQuads.push_back((Vector6d() << q1, 0.5, 0.25, 0).finished());
			debugQuads.push_back((Vector6d() << q2, 0.5, 0.25, 0).finished());
			debugQuads.push_back((Vector6d() << q3, 0.5, 0.25, 0).finished());
			for(unsigned i = 0; i < nearV.size(); i++) {
				debugVerts.push_back((Vector6d() << nearV[i]->position, 1, 0, 1).finished());
			}

			// edge post-clip, pre-snap
			debugEdges.push_back((Vector6d() << v0, 0.5, 0, 0).finished());
			debugEdges.push_back((Vector6d() << v1, 0.5, 0, 0).finished());
		}

		snapIntersections(e, nearV, cost);
		if(DEBUG) {
			// edge post-snap
			debugEdges.push_back((Vector6d() << v0, 0, 1, 0).finished());
			if(e.V[1] >= 0)
				v1 = ((Vertex*) getVertex(e.V[1]))->position;
			else
				v1 = e.newVertices[1];
			debugEdges.push_back((Vector6d() << v1, 0, 1, 0).finished());
		}
	}
	if(checkLegality(e, cost)) {
		double newLen = getLength(e);
		if(newLen < origLen)
			cost += (origLen - newLen) * COST_ROAD_SHORTENING;
		if(origLen < newLen)
			cost += (newLen - origLen) * COST_ROAD_LENGTHENING;
		return (cost < jitter(COST_THRESHOLD, COST_THRESHOLD_JITTER));
	}
	else return false;
}

void RoadSystem::globalGoals(int t, int i0, int i1) {
	Vertex* v0 = (Vertex*) getVertex(i0);
	Vertex* v1 = (Vertex*) getVertex(i1);
	if(v1->adjacent.size() > 1) return;

	EdgeKey e(i0, i1);
	double angle = edges[e].angle;
	double axis = edges[e].majorAxis;
	double branchDelay = edges[e].branchDelay;
	
	double j = 0;
	if(getRand() < GRID_DEVIATION_PROBABILITY) {
		j = jitter(0, GRID_DEVIATION_AMOUNT);
		j = j*j*j;
	}

	double angleDeltas[2] = {-1.57, 1.57};
	double x0 = v0->position[0];
	double y0 = v0->position[1];
	double x1 = v1->position[0];
	double y1 = v1->position[1];
	double z0 = elevation(x0, y0);
	double z1 = elevation(x1, y1);
	double rho = density(x1, y1);

	switch(edges[e].type) {
		case RT_AXIOM: {
			//axiom produces one highway from each end
			double a[2] = { angle+3.14, angle };
			int vi[2] = { i0, i1 };
			Vector3d positions[2] = { Vector3d(x0, y0, z0),
									  Vector3d(x1, y1, z1) };
			double delay = (t+EVAL_DELAY_HIGHWAY)/rho;
			for(unsigned i = 0; i < 2; i++) {
				double x, y, z;
				getGridPoint(a[i], positions[i], 0, 1, x, y, z);
				ProposedEdge ne(delay, RT_HIGHWAY, angle+a[i], angle+a[i],
								vi[i], x, y);
				ne.r.branchDelay = HIGHWAY_BRANCH_DELAY;
				roadQueue.push(ne);
			}

			//side roads
			double b[2] = {-1.57, 1.57};
			delay = (t+EVAL_DELAY_MAJOR)/rho;
			for(unsigned i = 0; i < 2; i++) {
				double x, y, z;
				getGridPoint(angle, Vector3d(x1, y1, z1), j, 
							 i == 0? 0 : 2, x, y, z);
				ProposedEdge ne = ProposedEdge(delay, RT_GRID_MAJOR,
											   angle+b[i]+j, axis+j, i1, x, y);
				roadQueue.push(ne);
			}
			break;
		}
		case RT_HIGHWAY: {
			//highway produces another highway segment from end
			//perpendicular gridmajor; these have small chance
			//to be promoted to highways
			double x, y, z;
			double fAngle = getHighwayPoint(angle, Vector3d(x1, y1, z1), x, y, z);
			double delay = (t+EVAL_DELAY_HIGHWAY)/rho;
			ProposedEdge newHighway(delay, RT_HIGHWAY, fAngle, fAngle,
									i1, x, y);
			if(branchDelay > 0)
				newHighway.r.branchDelay = branchDelay - 1;
			else
				newHighway.r.branchDelay = HIGHWAY_BRANCH_DELAY;
			roadQueue.push(newHighway);
	
			//side roads
			double a[2] = {-1.57, 1.57};
			for(unsigned i = 0; i < 2; i++) {
				getGridPoint(angle, Vector3d(x1, y1, z1), j, 
							 i == 0? 0 : 2, x, y, z);
				ProposedEdge ne;
				if(getRand() > HIGHWAY_PROMOTION_PROBABILITY) {
					//you lose the promotion roll
					delay = (t+EVAL_DELAY_MAJOR)/rho;
					ne = ProposedEdge(delay, RT_GRID_MAJOR,
									  angle+a[i]+j, axis+j, i1, x, y);
				} else if(nIterations < CUTOFF_HIGHWAY_GROWTH && branchDelay == 0) {
					delay = (t+EVAL_DELAY_HIGHWAY)/rho;
					ne = ProposedEdge(delay, RT_HIGHWAY,
									  angle+a[i]+j, angle+a[i]+j, i1, x, y);
					ne.r.branchDelay = HIGHWAY_BRANCH_DELAY;
				} else continue;
				roadQueue.push(ne);
			}
			break;
		}
		case RT_GRID_MAJOR: {
			double x, y, z;
			double a[3] = {-1.57, 0, 1.57};
			double delay;
			for(unsigned i = 0; i < 3; i++) {
				ProposedEdge ne;
				if(i == 1) {
					double fAngle = getHighwayPoint(angle, Vector3d(x1, y1, z1), x, y, z);
					delay = (t+EVAL_DELAY_MAJOR)/rho;
					ne = ProposedEdge(delay, RT_GRID_MAJOR, fAngle, fAngle,
									  i1, x, y);
				} else {
					getGridPoint(angle, Vector3d(x1, y1, z1), j, i, x, y, z);
					delay = (t+EVAL_DELAY_MINOR)/rho;
					ne = ProposedEdge(delay, RT_GRID_MINOR,
									  angle+a[i]+j, axis+j, i1, x, y);
				}
				roadQueue.push(ne);
			}
			break;
		}
		case RT_GRID_MINOR: {
			double x, y, z;
			double a[3] = {-1.57, 0, 1.57};
			double delay;
			for(unsigned i = 0; i < 3; i++) {
				ProposedEdge ne;
				if(i == 1) {
					getGridPoint(angle, Vector3d(x1, y1, z1), j, i, x, y, z);
					delay = (t+EVAL_DELAY_MINOR)/rho;
					ne = ProposedEdge(delay, RT_GRID_MINOR,
									  angle+a[i]+j, axis+j, i1, x, y);
				} else if(getRand() < SIDE_ROAD_PROBABILITY) {
					getGridPoint(angle, Vector3d(x1, y1, z1), j, i, x, y, z);
					delay = (t+EVAL_DELAY_MAJOR)/rho;
					ne = ProposedEdge(delay, RT_GRID_MAJOR,
									  angle+a[i]+j, axis+j, i1, x, y);
				} else continue;
				roadQueue.push(ne);
			}
			break;
		}
	}
}

bool RoadSystem::snapIntersections(ProposedEdge &e, 
								   vector<QuadPoint*> &nV, 
								   double &cost) {
	if(nV.size() == 0) return false;
	int i0 = e.V[0];
	if(e.V[1] >= 0) return false; //don't move existing verts
	double x1 = e.newVertices[1][0];
	double y1 = e.newVertices[1][1];
	Vector3d vpIsec(x1, y1, 0);

	//find nearest vertex
	double dist = 1e20;
	int best = -1;
	double x, y, curDist = 0;
	for(unsigned i = 0; i < nV.size(); i++) {
		x = nV[i]->position[0];
		y = nV[i]->position[1];
		Vector3d tIsec = Vector3d(x, y, 0);
		curDist = (tIsec - vpIsec).norm();
		if(curDist < dist && nV[i]->index != i0) {
			best = nV[i]->index;
			dist = curDist;
		}
	}

	if(best < 0) return false;
	Vertex* v2 = (Vertex*) getVertex(best);
	if(DEBUG) {
		debugVerts.push_back((Vector6d() << v2->position, 1, 0, 0).finished());
	}
	double x2 = v2->position[0];
	double y2 = v2->position[1];

	if(dist > SNAP_THRESHOLD) {
		return false;
	} else {
		//cost = distance moved
		cost += dist;
		//calculate cost/bonus for vertex valence
		if(v2->adjacent.size() >= 4) {
			int extra = v2->adjacent.size() - 3;
			cost += extra * COST_INTERSECTION_OVERLOAD;
		} else {
			int missing = 4 - v2->adjacent.size();
			cost += missing * BONUS_INTERSECTION_UNDERLOAD;
		}
		if(e.V[1] < 0)
			e.V[1] = best;
		else printf("HUH\n");
		e.splitRoad = false;
		return true;
	}
}

bool RoadSystem::clipRoad(ProposedEdge &e, vector<QuadPoint*> &qp,
						  double &cost) {
	//get edges
	set<EdgeKey> edges;
	getRoads(qp, edges);

	//get length
	Vertex* v0 = (Vertex*) getVertex(e.V[0]);
	// dummy since v1 isn't stored as a Vertex
	Vertex* v1 = new Vertex(e.newVertices[1][0], e.newVertices[1][1], -1);
	double x0 = v0->position[0];
	double y0 = v0->position[1];
	double x1 = v1->position[0];
	double y1 = v1->position[1];
	double prLen = (Vector2d(x1, y1) - Vector2d(x0, y0)).norm();

	Vector3d clipPoint, extensionPoint;
	EdgeKey clippingRoad;
	EdgeKey extendedRoad;
	double curU1_c = 1e20, curU2_c = 1e20;
	double curU1_e = 1e20, curU2_e = 1e20;
	bool clipped = false, extended = false;

	for(set<EdgeKey>::const_iterator iter = edges.begin();
		iter != edges.end(); ++iter) {

		EdgeKey edge = *iter;
		if(!areConnected(edge.V[0], edge.V[1]))
			continue;
		//get intersection
		Vertex* v2 = (Vertex*) getVertex(edge.V[0]);
		Vertex* v3 = (Vertex*) getVertex(edge.V[1]);

		double xi, yi, u1, u2;
		intersectRoads(v0, v1, v2, v3, xi, yi, u1, u2);

		//u1 < 0 -> intersection behind -> irrelevant
		if(u1 <= 0) continue;
		//u2 out of interval -> irrelevant
		if(0 >= u2 || u2 >= 1) continue;

		if(u1 < curU1_c && u1 < 1) {
			//potential clip target
			clipped = true;
			curU1_c = u1;
			curU2_c = u2;
			clippingRoad = EdgeKey(edge);
			clipPoint = Vector3d(xi, yi, 0);
		} else if(u1 < curU1_e && u1 > 1) {
			//potential extension target
			extended = true;
			curU1_e = u1;
			curU2_e = u2;
			extendedRoad = EdgeKey(edge);
			extensionPoint = Vector3d(xi, yi, 0);
		}
	}
	delete v1;

	if(!clipped && !extended) return 0;
	//got clip/extension point
	if(clipped) {
		e.newVertices[1] = clipPoint;
		e.splitVerts[0] = clippingRoad.V[0];
		e.splitVerts[1] = clippingRoad.V[1];
	} else if(extended) {
		e.newVertices[1] = extensionPoint;
		e.splitVerts[0] = extendedRoad.V[0];
		e.splitVerts[1] = extendedRoad.V[1];
	}
	e.splitRoad = true;
	double dist = (v0->position - e.newVertices[1]).norm();
	cost += dist;
	if(extended) cost += 10*dist;
	return true;
}

bool RoadSystem::checkLegality(ProposedEdge &e, double &cost) {
	//redundancy check if both endpoints already exist
	if(e.V[1] >= 0) {
		EdgeKey edge(e.V[0], e.V[1]);
		Edges::const_iterator iter = edges.find(edge);
		if(iter != edges.end()) {
			//edge already exists, reject
			return false;
		}
		return true;
	} else {
		//check against all other edges out of v0
		//get dot product, reject if too high
		//i.e. no slivers
		Vertex* vert0 = (Vertex*) getVertex(e.V[0]);
		Vector3d v0 = vert0->position;
		Vector3d v1;
		Vector3d unitEdge;
		if(e.V[1] < 0)
			v1 = e.newVertices[1];
		else
			v1 = ((Vertex*) getVertex(e.V[1]))->position;
		unitEdge = (v0 - v1).normalized();

		const vector<Vertex*> adj = vert0->adjacent;
		for(unsigned i = 0; i < adj.size(); i++) {
			//get edge
			Vertex* vert2 = adj[i];
			Vector3d v2 = vert2->position;
			Vector3d adjEdge = (v0 - v2).normalized();
			//cos(20 deg) ~= 0.94
			double dot = unitEdge.dot(adjEdge);
			if(dot > 0.94) return false;
		}
		
		//todo: have road follow contour of boundary
		return legal(e.newVertices[1][0], e.newVertices[1][1]);
	}
}

//returns new angle
double RoadSystem::getHighwayPoint(double srcAngle, Vector3d &src, 
								   double &x, double &y, double &z) const {
	double currMinScore = 1e20;
	double currMinAngle;
	double length = jitter(ROAD_LENGTH_HIGHWAY, ROAD_LENGTH_JITTER);
	//start sampling
	for(unsigned i = 0; i < NUM_HIGHWAY_PROBES; i++) {
		double probeAngle = ((ANGLE_HIGHWAY_PROBES*i)/NUM_HIGHWAY_PROBES) -
							(ANGLE_HIGHWAY_PROBES/2) + srcAngle;
		double px = src[0] + length*cos(probeAngle);
		double py = src[1] + length*sin(probeAngle);
		double pz = elevation(px, py);
		double deltaZ = abs(pz - src[2]);
		if(deltaZ < currMinScore) {
			currMinScore = deltaZ;
			currMinAngle = probeAngle;
		}
	}
	double terrainAngle = jitter(currMinAngle, HIGHWAY_PROBE_JITTER);
	double w = jitter(WEIGHT_TERRAIN_CONTOUR, WEIGHT_TERRAIN_CONTOUR_JITTER);
	double fAngle = w*terrainAngle + (1-w)*srcAngle;
	x = src[0] + length*cos(fAngle);
	y = src[1] + length*sin(fAngle);
	return fAngle;
}

// i in [0,2]
void RoadSystem::getGridPoint(double srcAngle, Vector3d &src, double j, int i, 
							  double &x, double &y, double &z) const {
	double dA[3] = {-1.57, 0, 1.57};
	double lJ = jitter(0, ROAD_LENGTH_JITTER);
	double l1 = ROAD_LENGTH_MAJOR + lJ;
	double l2 = ROAD_LENGTH_MINOR + lJ;
	
	x = src[0] + l1*cos(srcAngle+dA[i])*cos(j) - l2*sin(srcAngle+dA[i])*sin(j);
	y = src[1] + l1*cos(srcAngle+dA[i])*sin(j) + l2*sin(srcAngle+dA[i])*cos(j);
}

void RoadSystem::extractPolygons(vector<Polygon> &polygons) {
	if(!gotPolys)
		getPolygons();
	polygons = polys;
}

void RoadSystem::getPolygons() {
	if(!gotPolys) {
		vector<double> V; V.reserve(vertices.size() * 2);
		vector<double> E; E.reserve(edges.size() * 2);
		Vertices::const_iterator v_iter = vertices.begin();
		for(; v_iter != vertices.end(); ++v_iter) {
			V.push_back(v_iter->second->position[0]);
			V.push_back(v_iter->second->position[1]);
		}
		Edges::const_iterator e_iter = edges.begin();
		for(; e_iter != edges.end(); ++e_iter) {
			EdgeKey ek = e_iter->first;
			if(!areConnected(ek.V[0], ek.V[1]))
				continue;
			E.push_back(ek.V[0]);
			E.push_back(ek.V[1]);
		}

		//generate a new Graph with the same data
		//since minimal cycle basis extraction is destructive
		Graph G(V, E);
		G.getPolygons(polys);
		int origSize = polys.size();
		for(unsigned p = 0; p < polys.size(); p++) {
			if(polys[p].vertices.size() <= 2) {
				polys.erase(polys.begin()+p);
				p--;
			} else {
				//scale
				//todo: polygon skeleton?
				//accumulate polygon centroid
				Vector3d avg = Vector3d(0, 0, 0);
				for(unsigned i = 0; i < polys[p].vertices.size(); i++) {
					double x = polys[p].vertices[i][0];
					double y = polys[p].vertices[i][1];
					avg[0] += (x - avg[0])/(i+1);
					avg[1] += (y - avg[1])/(i+1);
				}

				double s = 0.7;
				for(unsigned i = 0; i < polys[p].vertices.size(); i++) {
					Vector3d newpoint = polys[p].vertices[i];
					newpoint -= avg;
					newpoint *= s;
					newpoint += avg;
					//newpoint[2] = elevation(newpoint[0], newpoint[1])/8;
					polys[p].vertices[i] = Vector3d(newpoint);

				}
			}
		}
		int size = polys.size();
		printf("got %i blocks\n", size);
		gotPolys = true;
	}
}

void RoadSystem::drawDebug() const {
	//draw quads
	//assume debugQuads.size%4 = 0
	for(unsigned i = 0; i < debugQuads.size(); i += 4) {
		Vector6d quads[4] = { debugQuads[i], 
							  debugQuads[i+1], 
							  debugQuads[i+2], 
							  debugQuads[i+3] };
		glBegin(GL_LINE_STRIP);
			for(unsigned j = 0; j < 4; j++) {
				glColor3f(quads[j][3], quads[j][4], quads[j][5]);
				glVertex3f(quads[j][0], quads[j][1], 0);
			}
			glColor3f(quads[0][3], quads[0][4], quads[0][5]);
			glVertex3f(quads[0][0], quads[0][1], 0);
		glEnd();
	}
	//draw edges
	//assume debugEdges has even size; consecutive verts = edge
	for(unsigned i = 0; i < debugEdges.size(); i += 2) {
		Vector6d pts[2] = { debugEdges[i],
							debugEdges[i+1] };
		glBegin(GL_LINE_STRIP);
			for(unsigned j = 0; j < 2; j++) {
				glColor3f(pts[j][3], pts[j][4], pts[j][5]);
				glVertex3f(pts[j][0], pts[j][1], 0);
			}
		glEnd();
	}
	//draw verts
	for(unsigned i = 0; i < debugVerts.size(); i++) {
		Vector6d v = debugVerts[i];
		if(v[3] == 1 && v[4] == 0 && v[5] == 0)
			glPointSize(5);
		else
			glPointSize(3);
		glBegin(GL_POINTS);
		glColor3f(v[3], v[4], v[5]);
		glVertex3f(v[0], v[1], 0);
		glEnd();
	}
}

void RoadSystem::draw(bool drawBlocks, bool drawQuadTree) const {
	double x0, y0, z0, x1, y1, z1;
	if(drawQuadTree)
		qtree.draw();
	if(DEBUG)
		drawDebug();
	for(Edges::const_iterator iter = edges.begin();
		iter != edges.end(); ++iter) {
		
		EdgeKey e = iter->first;
		roadAttr r = iter->second;
		Vertex* v0 = (Vertex*) getVertex(e.V[0]);
		Vertex* v1 = (Vertex*) getVertex(e.V[1]);
		x0 = v0->position[0];
		y0 = v0->position[1];
		x1 = v1->position[0];
		y1 = v1->position[1];
		z0 = elevation(x0,y0)/8;
		z1 = elevation(x1,y1)/8;
		glBegin(GL_LINE_STRIP);
			if(drawBlocks)
				glColor3f(0.7, 0.7, 0.7);
			else
				glColor3f(0.3, 0.3, 0.3);
			switch(r.type) {
			case RT_GRID_MINOR:
				glColor4f(0.5,0,0, 0.5);
				glVertex3f(x0, y0, z0);
				glColor4f(0.1,0.1,0.1, 0.5);
				glVertex3f(x1, y1, z1);
				break;
			case RT_GRID_MAJOR:
				glColor4f(0, 0.5, 0.5, 0.5);
				glVertex3f(x0, y0, z0);
				glColor4f(0.1, 0.1, 0.1, 0.5);
				glVertex3f(x1, y1, z1);
				break;
			case RT_HIGHWAY:
				glColor4f(1, 0.5, 0, 0.5);
				glVertex3f(x0, y0, z0);
				glColor4f(0.1, 0.1, 0.1, 0.5);
				glVertex3f(x1, y1, z1);
				break;
			case RT_AXIOM:
				glColor4f(1,0.5,0, 0.5);
				glVertex3f(x0, y0, z0);
				glColor4f(0.1, 0.1, 0.1, 0.5);
				glVertex3f(x1, y1, z1);
				break;
			}
		glEnd();
	}
}

void RoadSystem::printStats() const {
	printf("Vertices: %i\n", vertices.size());
	printf("Edges: %i\n", edges.size());
}

void RoadSystem::reset() {
	vertices.clear();
	edges.clear();
	roadQueue = priority_queue<ProposedEdge, vector<ProposedEdge>, ProposedEdgeComparator>();
	qtree.reset();
}

bool RoadSystem::areConnected(int i0, int i1) const {
	Vertex* v0 = (Vertex*) getVertex(i0);
	Vertex* v1 = (Vertex*) getVertex(i1);
	return areConnected(v0, v1);
}

bool RoadSystem::areConnected(Vertex* v0, Vertex* v1) const {
	bool b_v0 = false;
	bool b_v1 = false;
	bool b_e = false;
	//check verts list each other in adjacency lists
	vector<Vertex*> adj0 = v0->adjacent;
	vector<Vertex*> adj1 = v1->adjacent;
	for(unsigned i = 0; i < adj0.size(); i++) {
		if(adj0[i] == v1) {
			b_v0 = true;
			break;
		}
	}
	for(unsigned i = 0; i < adj1.size(); i++) {
		if(adj1[i] == v0) {
			b_v1 = true;
			break;
		}
	}

	//check edge exists in edges
	EdgeKey e = EdgeKey(v0->index, v1->index);
	b_e = edges.count(e);
	return b_v0 && b_v1 && b_e;
}

void RoadSystem::getVertices(Vector3d q0, Vector3d q1, Vector3d q2, Vector3d q3,
							 vector<const RoadSystem::Vertex*> &verts) {
	vector<QuadPoint*> qps;
	qtree.query(q0, q1, q2, q3, qps);
	for(unsigned i = 0; i < qps.size(); i++) {
		verts.push_back(getVertex(qps[i]->index));
	}
}

// O(n^2) rofl
void RoadSystem::getEdges(Vector3d q0, Vector3d q1, Vector3d q2, Vector3d q3,
						  set<const RoadSystem::EdgeKey> &edges) {
	vector<QuadPoint*> qps;
	qtree.query(q0, q1, q2, q3, qps);
	for(unsigned i = 0; i < qps.size(); i++) {
		const Vertex* v = (Vertex*) getVertex(qps[i]->index);
		vector<Vertex*> adj = v->adjacent;
		for(unsigned j = 0; j < adj.size(); j++) {
			const EdgeKey e = EdgeKey(v->index, adj[j]->index);
			edges.insert(e);
		}
	}
}

double RoadSystem::elevation(double x, double y) const {
	return ((double)is.getElevation(x, y))/255;
}

double RoadSystem::density(double x, double y) const {
	return (((double)is.getDensity(x, y))+1)/32;
}

bool RoadSystem::legal(double x, double y) const {
	bool valid = is.isValid(x, y);
	return valid;
}

/*-----------------------------------------------------
 RoadSystem - RNG functions
------------------------------------------------------*/

//returns a double in the range[c-r, c+r]
inline double RoadSystem::jitter(double c, double r) const {
	double x = 2*((double)rand()/RAND_MAX);
	x -= 1;
	return r*x + c;
}

//returns a double in the range [0,1]
inline double RoadSystem::getRand() const {
	return ((double)rand()/RAND_MAX);
}

/*------------------------------------------------------
 RoadSystem - output functions
 ------------------------------------------------------*/
void RoadSystem::dumpPolygons() const {
	ofstream outputFile;
	outputFile.open("blocks.obj");
	for(unsigned i = 0; i < polys.size(); i++) {
		for(unsigned j = 0; j < polys[i].vertices.size(); j++) {
			//maya uses y-up
			outputFile << "v ";
			outputFile << polys[i].vertices[j][0] << " ";
			outputFile << polys[i].vertices[j][2] << " ";
			outputFile << polys[i].vertices[j][1] << " ";
			outputFile << endl;
		}
	}
	int vertexCounter = 0;
	for(unsigned i = 0; i < polys.size(); i++) {
		outputFile << "f ";
		for(unsigned j = 0; j < polys[i].vertices.size(); j++)
			outputFile << ++vertexCounter << " ";
		outputFile << endl;
	}
	outputFile.close();
}

void RoadSystem::dumpRoads() const {
	ofstream outputFile;
	outputFile.open("roads.mel");
	Edges::const_iterator iter = iter = edges.begin();
	int scale = 25;
	for(; iter != edges.end(); ++iter) {
		EdgeKey e = iter->first;
		roadAttr r = iter->second;
		Vertex* v0 = (Vertex*) getVertex(e.V[0]);
		Vertex* v1 = (Vertex*) getVertex(e.V[1]);
		if(!areConnected(v0, v1))
			continue;
		Vector3d p0 = v0->position;
		Vector3d p1 = v1->position;

		//make the curve, assign to $c
		outputFile << "$c = `curve -d 1 -p ";
		outputFile << p0[0] << " " << p0[2] << " " << p0[1] << " -p ";
		outputFile << p1[0] << " " << p1[2] << " " << p1[1] << " ";
		outputFile << "-k 0 -k 1`;\n";

		//set keyframes
		outputFile << "setKeyframe -v 0 -t 0 -at visibility $c;\n";
		outputFile << "setKeyframe -v 1 -t " << r.iter/scale << " -at visibility $c;\n";
	}
}