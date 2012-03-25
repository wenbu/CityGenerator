/**---------------------
* All code in Graph.* more or less lifted from
* http://www.geometrictools.com
----------------------**/

#pragma once
#include "stdafx.h"
#include "Polygon.h"
#include "Intersection.h"
#include "Road.h"
#include "RoadSystem.h"
#include <vector>
#include <stack>
#include <set>
#include <map>

class Graph {
public:
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
		int index; //index into isecn vector in roadsystem!
		double x, y;
		vector<Vertex*> adjacent;
		void insert(Vertex* adjacent);
		void remove(Vertex* adjacent);
	};

	typedef map<int, Vertex*> Vertices;
	typedef map<EdgeKey, bool> Edges;

	Graph(const vector<const Intersection> &isecns, const vector<const Road> &roads);
	Graph(const vector<double> &V, const vector<double> &E);
	~Graph(void);

	const Vertices& getVertices() const;
	const Vertex* getVertex(int i) const;
	bool insertVertex(double x, double y, int i);
	bool removeVertex(int i);

	const Edges& getEdges() const;
	bool insertEdge(int i1, int i2);
	bool removeEdge(int i1, int i2);

	struct PrimitiveVertex {
		double x;
		double y;
		int index;
		PrimitiveVertex(double a, double b, int j) {
			x = a; y = b; index = j;
		}
	};

	enum PrimitiveType {
		PT_VERTEX,
		PT_FILAMENT,
		PT_CYCLE
	};

	class Primitive {
	public:
		Primitive(PrimitiveType type);
		PrimitiveType type;
		vector<PrimitiveVertex> seq;
	};

	void getPolygons(vector<Polygon> &polys);
	void extractPrimitives(vector<Primitive*> &primitives);

protected:
	class VertexPtr {
	public:
		VertexPtr(Vertex* vertex);
		inline operator Vertex* () { return mVertex; }
		bool operator < (const VertexPtr& vertexPtr) const;
	private:
		Vertex* mVertex;
	};

	void setCycleEdge(int v0, int v1);
	bool getCycleEdge(int v0, int v1) const;

	void extractVertex(Vertex* v0, set<VertexPtr> &heap, vector<Primitive*> &primitives);
	void extractFilament(Vertex* v0, Vertex* v1, set<VertexPtr> &heap, vector<Primitive*> &primitives);
	void extractPrimitive(Vertex* v0, set<VertexPtr> &heap, vector<Primitive*> &primitives);

	Vertex* getCW(Vertex* vprev, Vertex* vcurr);
	Vertex* getCCW(Vertex* vprex, Vertex* vcurr);

	Vertices mVertices;
	Edges mEdges;
};