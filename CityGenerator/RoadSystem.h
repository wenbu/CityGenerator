#pragma once
//#include "stdafx.h"
#include "Polygon.h"
#include "Intersection.h"
#include "Road.h"
#include "QuadNode.h"
#include "Graph.h"
#include "ImageSampler.h"
#include "InputParser.h"
#include <GL/glut.h>
#include <GL/glu.h>
#include <queue>
#include <map>
#include <time.h>
#include <iostream>
#include <fstream>

#define DEBUG false

using namespace std;
using namespace Eigen;

class RoadSystem {
	//all angles are in radians

	// quadtree query box (for intersection snaps) has side length
	// 2*QUERY_DIST
	static double QUERY_DIST;

	// will not snap intersections further than SNAP_THRESHOLD
	static double SNAP_THRESHOLD;

	static double FINAL_SNAP_THRESHOLD;

	// ideal road length, major axis
	static double ROAD_LENGTH_MAJOR;

	// ideal road length, minor axis
	static double ROAD_LENGTH_MINOR;

	static double ROAD_LENGTH_HIGHWAY;

	// jitter for road length
	static double ROAD_LENGTH_JITTER;

	static double EVAL_DELAY_HIGHWAY;
	static double EVAL_DELAY_MAJOR;
	static double EVAL_DELAY_MINOR;

	// cutoff for highway growth
	// after this many iterations, no more highways
	// are proposed
	static double CUTOFF_HIGHWAY_GROWTH;

	// need at least this many segments between branches
	static double HIGHWAY_BRANCH_DELAY;

	// number of probes for highway growth
	static double NUM_HIGHWAY_PROBES;

	// arc to probe in highway growth
	static double ANGLE_HIGHWAY_PROBES;

	// probe angle jitter
	static double HIGHWAY_PROBE_JITTER;

	// weight for highway contour following [0,1]
	static double WEIGHT_TERRAIN_CONTOUR;
	static double WEIGHT_TERRAIN_CONTOUR_JITTER;
	
	//weight for highway gradient following [0,1]
	static double WEIGHT_TERRAIN_GRADIENT;
	static double WEIGHT_TERRAIN_GRADIENT_JITTER;

	// probability that a grid major road will be promoted
	// to a highway
	static double HIGHWAY_PROMOTION_PROBABILITY;

	// probability that proposed roads will deviate from the
	// established pattern
	static double GRID_DEVIATION_PROBABILITY;

	// deviating roads will do so by [-GRID_DEVIATION_AMOUNT,
	// +GRID_DEVIATION_AMOUNT]^3 radians
	static double GRID_DEVIATION_AMOUNT;

	// probability that a side road will be proposed
	static double SIDE_ROAD_PROBABILITY;

	// cost threshold for road modification
	static double COST_THRESHOLD;

	// cost threshold jitter
	static double COST_THRESHOLD_JITTER;

	// cost for adding each additional road to an intersection
	// with more than 3 attached roads
	static double COST_INTERSECTION_OVERLOAD;

	// cost for shortening road
	// probably higher than for lengthening, since mostly we
	// just don't want stupidly short roads
	static double COST_ROAD_SHORTENING;

	// cost for lengthening road
	static double COST_ROAD_LENGTHENING;

	// bonus for adding road to an intersection
	// with less than 4 attached roads
	static double BONUS_INTERSECTION_UNDERLOAD;

	enum RoadType {
		RT_AXIOM,
		RT_GRID_MINOR,
		RT_GRID_MAJOR,
		RT_HIGHWAY
	};

	/**
	* when adding new roadAttr members
	* need to update:
	*  - ProposedRoad constructors
	*  - insertEdge
	*  - globalGoals
	**/
	struct roadAttr {
		RoadType type;
		double angle;
		double majorAxis;
		int iter; //iteration road was created
		int branchDelay; //don't want branching too often
						 //start at N, count down
	};

	class EdgeKey {
	public:
		EdgeKey(int v0 = -1, int v1 = -1);
		bool operator < (const EdgeKey &key) const;
		operator size_t () const;
		int V[2];
	};

	class Vertex {
	public:
		Vertex(double x, double y, int index);
		~Vertex();
		int index;
		Vector3d position;
		vector<Vertex*> adjacent;
		void insert(Vertex* adjacent);
		void remove(Vertex* adjacent);
	};

	class ProposedEdge {
	public:
		ProposedEdge(void);
		ProposedEdge(int t, RoadType type,
					 double angle, double majorAxis);
		ProposedEdge(int t, RoadType type,
					 double angle, double majorAxis,
					 double x0, double y0, double x1, double y1);
		ProposedEdge(int t, RoadType type,
					 double angle, double majorAxis,
					 int v0, double x1, double y1);
		ProposedEdge(int t, RoadType type,
					 double angle, double majorAxis,
					 int v0, int v1);
		int t; //evaluation priority
		int V[2]; //indices of vertices
		Vector3d newVertices[2]; //specification of nonexistent vertices
		roadAttr r;
		bool splitRoad;
		int splitVerts[2]; //indices of vertices defining
						   //road to split
	};

	class ProposedEdgeComparator {
	public:
		bool operator() (ProposedEdge &e1, ProposedEdge &e2) {
			return e1.t > e2.t;
		}
	};

	typedef map<int, Vertex*> Vertices;
	typedef map<EdgeKey, roadAttr> Edges;

	QuadNode qtree; //for spatial queries
	Vertices vertices;
	Edges edges;
	vector<Polygon> polys;
	vector<EdgeKey> deletedEdges;
	// vector6d: x, y, z, r, g, b
	typedef Matrix<double, 6, 1> Vector6d;
	vector<Vector6d> debugVerts;
	vector<Vector6d> debugEdges;
	vector<Vector6d> debugQuads;
	bool gotPolys;
	Vector3d minCoord, maxCoord;
	priority_queue<ProposedEdge, vector<ProposedEdge>, ProposedEdgeComparator> roadQueue;
	ImageSampler is;

	int axiomVerts;
	const Vertex* getVertex(int i) const;
	bool insertVertex(double x, double y, int i);
	bool removeVertex(int i);

	void getRoads(vector<QuadPoint*> &qp, set<EdgeKey> &roads) const;
	bool insertEdge(int i1, int i2, roadAttr r);
	bool removeEdge(int i1, int i2);

public:
	int nIterations;
	RoadSystem(void);
	~RoadSystem(void);
	RoadSystem(ImageSampler &img, InputParser &parser);
	void generate(int iterations);
	bool localConstraints(ProposedEdge &r);
	void globalGoals(int t, int i0, int i1);
	void addAxiom(double x1, double y1, double x2, double y2);
	void draw(bool drawBlocks = false, bool drawQuadTree = false) const;
	void printStats() const;
	void extractPolygons(vector<Polygon> &polys);
	void dumpPolygons() const;
	void dumpRoads() const;
	void getPolygons();
	void reset();
private:
	void acceptProposal(ProposedEdge &pe, int &i0, int &i1);
	
	// todo: need to re-query for verts again post-clip
	bool snapIntersections(ProposedEdge &e, vector<QuadPoint*> &nearIntersections, double &cost);
	bool clipRoad(ProposedEdge &e, vector<QuadPoint*> &nearRoads, double &cost);
	bool checkLegality(ProposedEdge &e, double &cost);
	
	void intersectRoads(const Vertex* v0, const Vertex* v1, 
						const Vertex* v2, const Vertex* v3, 
						double &x, double &y, 
						double &u1, double &u2) const;
	
	double getHighwayPoint(double srcAngle, Vector3d &src, double &x, double &y, double &z) const;
	double getSteepPoint(double srcAngle, Vector3d &src, double &x, double &y, double &z) const;

	//issue: getGridPoint is not aware of the graph structure. Ideas:
	// * Pass in v1, query edges for dot product, discard if too low (no slivers)
	// * Pass in v1, query quadtree for verts close by, discard if so (no tiny shit)
	// * Pass in v1, query quadtree for verts/edges in direction of road, discard if none
	//      (no roads to nowhere)

	//would it be better to throw this stuff in checkLegality instead?
	//todo: need to check these in snap/clip also
	void getGridPoint(double srcAngle, Vector3d &src, double rotJitter, int index, double &x, double &y, double &z) const;

	/************
	* util funcs
	************/
	//draw the edges contained in debugEdges
	void drawDebug() const;

	//query quadtree for verts/edges within the quad defined by [q0, q1, q2, q3]
	void getVertices(Vector3d q0, Vector3d q1, Vector3d q2, Vector3d q3, vector<const Vertex*> &verts);
	void getEdges(Vector3d q0, Vector3d q1, Vector3d q2, Vector3d q3, set<const EdgeKey> &edges);

	//query graph for connectivity
	bool areConnected(int i0, int i1) const;
	bool areConnected(Vertex* v0, Vertex* v1) const;

	double getLength(const ProposedEdge &e) const;

	//query image maps at (x, y)
	double elevation(double x, double y) const;
	double density(double x, double y) const;
	bool legal(double x, double y) const;

	//RNG stuff

	//returns number in range [c-r, c+r]
	double jitter(double center, double range) const;

	//returns number in range [0, 1]
	double getRand() const;
};

