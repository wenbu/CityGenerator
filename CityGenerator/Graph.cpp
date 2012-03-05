#pragma once
#include "stdafx.h"
#include "Graph.h"
#include "Intersection.h"
#include "Road.h"

Graph::Graph(const vector<const Intersection> &isecns, const vector<const Road> &roads) {
	for(unsigned i = 0; i < isecns.size(); i++)
		insertVertex(isecns[i][0], isecns[i][1], i);

	for(unsigned i = 0; i < roads.size(); i++) {
		int i1 = roads[i].getIntersection(0);
		int i2 = roads[i].getIntersection(1);
		/*
		if(i1 == i2) {
			printf("lol\n");
			continue;
		}
		*/
		insertEdge(i1, i2);
	}
}

Graph::Graph(const vector<double> &V, const vector<double> &E) {
	for(unsigned v = 0; v < V.size(); v += 2)
		insertVertex(V[v], V[v+1], v/2);

	for(unsigned e = 0; e < E.size(); e += 2)
		insertEdge(E[e], E[e+1]);
}

inline Graph::EdgeKey::EdgeKey (int v0, int v1) {
	if(v0 < v1) {
		V[0] = v0;
		V[1] = v1;
	} else {
		V[0] = v1;
		V[1] = v0;
	}
}

inline bool Graph::EdgeKey::operator < (const EdgeKey &key) const {
	if(V[1] < key.V[1])
		return true;
	if(V[1] > key.V[1])
		return false;
	return V[0] < key.V[0];
}

inline Graph::EdgeKey::operator size_t () const {
	return V[0] | (V[1] << 16);
}

Graph::~Graph() {
	Vertices::iterator iter = mVertices.begin();
	Vertices::iterator end = mVertices.end();
	for(; iter != end; ++iter) {
		Vertex* v = iter->second;
		delete v;
	}
}

const Graph::Vertices& Graph::getVertices() const {
	return mVertices;
}

const Graph::Vertex* Graph::getVertex(int i) const {
	Vertices::const_iterator iter = mVertices.find(i);
	return(iter != mVertices.end() ? iter->second : 0);
}

bool Graph::insertVertex(double x, double y, int i) {
	Vertices::iterator iter = mVertices.find(i);
	if(iter != mVertices.end())
		return false;

	Vertex* vertex = new Vertex(x, y, i);
	mVertices[i] = vertex;
	return true;
}

bool Graph::removeVertex(int i) {
	Vertices::iterator iter = mVertices.find(i);
	if(iter != mVertices.end()) {
		Vertex* v = iter->second;
		if(v->adjacent.size() == 0) {
			mVertices.erase(iter);
			delete v;
			return true;
		}
	}
	return false;
}

const Graph::Edges& Graph::getEdges() const {
	return mEdges;
}

bool Graph::insertEdge(int i0, int i1) {
	Vertex* v0 = (Vertex*)getVertex(i0);
	if(!v0)
		return false;

	Vertex* v1 = (Vertex*)getVertex(i1);
	if(!v1)
		return false;

	EdgeKey eKey(i0, i1);
	map<EdgeKey, bool>::iterator iter = mEdges.find(eKey);
	if(iter == mEdges.end()) {
		mEdges[eKey] = false;
		v0->insert(v1);
		v1->insert(v0);
		return true;
	}
	return false;
}

bool Graph::removeEdge(int i0, int i1) {
	Vertex* v0 = (Vertex*)getVertex(i0);
	if(!v0) return false;

	Vertex* v1 = (Vertex*)getVertex(i1);
	if(!v1) return false;

	EdgeKey eKey(i0, i1);
	map<EdgeKey, bool>::iterator iter = mEdges.find(eKey);
	if(iter != mEdges.end()) {
		mEdges.erase(iter);
		v0->remove(v1);
		v1->remove(v0);
		return true;
	}
	return false;
}

void Graph::getPolygons(vector<Polygon> &polys) {
	vector<Primitive*> prims;
	extractPrimitives(prims);
	//printf("got %i primitives\n", prims.size());
	for(unsigned i = 0; i < prims.size(); i++) {
		if(prims[i]->type == PT_CYCLE) {
			Polygon p;
			for(unsigned j = 0; j < prims[i]->seq.size(); j++) {
				Vector3d v = Vector3d(prims[i]->seq[j].x,
									  prims[i]->seq[j].y,
									  0);
				p.vertices.push_back(v);
			}
			polys.push_back(p);
		}
	}
}

void Graph::extractPrimitives(vector<Primitive*> &primitives) {
	set<VertexPtr> heap;
	Vertices::iterator iter = mVertices.begin();
	Vertices::iterator end = mVertices.end();
	for(; iter != end; ++iter)
		heap.insert(iter->second);

	while(!heap.empty()) {
		VertexPtr vptr = *heap.begin();
		Vertex* v0 = (Vertex*)vptr;

		if(v0->adjacent.size() == 0)
			extractVertex(v0, heap, primitives);
		else if(v0->adjacent.size() == 1)
			extractFilament(v0, v0->adjacent[0], heap, primitives);
		else
			extractPrimitive(v0, heap, primitives);
	}
}

void Graph::setCycleEdge(int i0, int i1) {
	EdgeKey eKey(i0, i1);
	Edges::iterator iter = mEdges.find(eKey);
	if(iter != mEdges.end())
		iter->second = true;
}

bool Graph::getCycleEdge(int i0, int i1) const {
	EdgeKey eKey(i0, i1);
	Edges::const_iterator iter = mEdges.find(eKey);
	if(iter != mEdges.end())
		return iter->second;
	return false;
}

void Graph::extractVertex(Vertex* v0, set<VertexPtr> &heap, vector<Primitive*> &primitives) {
	//printf("extractVertex %i\n", v0->index);
	Primitive* p = new Primitive(PT_VERTEX);

	p->seq.push_back(PrimitiveVertex(v0->x,v0->y,v0->index));
	heap.erase(v0);
	removeVertex(v0->index);
	primitives.push_back(p);
}

void Graph::extractFilament(Vertex* v0, Vertex* v1,
							set<VertexPtr> &heap, vector<Primitive*> &primitives) {
	//printf("extractFilament %i %i\n", v0->index, v1->index);
	if(getCycleEdge(v0->index, v1->index)) {
		if(v0->adjacent.size() >= 3) {
			removeEdge(v0->index, v1->index);
			v0 = v1;
			if(v0->adjacent.size() == 1)
				v1 = v0->adjacent[0];
		}
		while(v0->adjacent.size() == 1) {
			v1 = v0->adjacent[0];
			if(getCycleEdge(v0->index, v1->index)) {
				heap.erase(v0);
				removeEdge(v0->index, v1->index);
				removeVertex(v0->index);
				v0 = v1;
			}
			else break;
		}

		if(v0->adjacent.size() == 0) {
			heap.erase(v0);
			removeVertex(v0->index);
		}
	}

	else {
		Primitive* p = new Primitive(PT_FILAMENT);

		if(v0->adjacent.size() >= 3) {
			p->seq.push_back(PrimitiveVertex(v0->x,v0->y,v0->index));
			removeEdge(v0->index, v1->index);
			v0 = v1;
			if(v0->adjacent.size() == 1)
				v1 = v0->adjacent[0];
		}

		while(v0->adjacent.size() == 1) {
			v1 = v0->adjacent[0];
			p->seq.push_back(PrimitiveVertex(v0->x,v0->y,v0->index));
			heap.erase(v0);
			removeEdge(v0->index, v1->index);
			removeVertex(v0->index);
			v0 = v1;
		}

		p->seq.push_back(PrimitiveVertex(v0->x,v0->y,v0->index));
		if(v0->adjacent.size() == 0) {
			heap.erase(v0);
			removeVertex(v0->index);
		}
		primitives.push_back(p);
	}
}

void Graph::extractPrimitive(Vertex* v0, set<VertexPtr> &heap, vector<Primitive*> &primitives) {
	//printf("extractPrimitive %i\n", v0->index);
	set<Vertex*> visited;
	vector<PrimitiveVertex> sequence;
	sequence.push_back(PrimitiveVertex(v0->x,v0->y,v0->index));
	Vertex* v1 = getCW(0, v0);
	Vertex* vprev = v0;
	Vertex* vcurr = v1;

	while(vcurr && vcurr != v0 &&
		visited.find(vcurr) == visited.end()) {
			sequence.push_back(PrimitiveVertex(vcurr->x,vcurr->y,vcurr->index));
			visited.insert(vcurr);
			Vertex* vnext = getCCW(vprev, vcurr);
			vprev = vcurr;
			vcurr = vnext;
	}

	if(!vcurr) {
		extractFilament(vprev, vprev->adjacent[0], heap, primitives);
	}
	else if(vcurr == v0) {
		Primitive* p = new Primitive(PT_CYCLE);
		p->seq = sequence;
		primitives.push_back(p);

		int s = (int)sequence.size();
		for(int i0 = s-1, i1 = 0; i1 < s; i0 = i1++) {
			int iv0 = sequence[i0].index;
			int iv1 = sequence[i1].index;
			setCycleEdge(iv0, iv1);
		}

		removeEdge(v0->index, v1->index);

		if(v0->adjacent.size() == 1)
			extractFilament(v0, v0->adjacent[0], heap, primitives);
		if(v1->adjacent.size() == 1)
			extractFilament(v1, v1->adjacent[0], heap, primitives);
	} else {
		while(v0->adjacent.size() == 2) {
			if(v0->adjacent[0] != v1) {
				v1 = v0;
				v0 = v0->adjacent[0];
			} else {
				v1 = v0;
				v0 = v0->adjacent[1];
			}
		}
		extractFilament(v0, v1, heap, primitives);
	}
}

Graph::Vertex* Graph::getCW(Vertex* vprev, Vertex* vcurr) {
	Vertex* vnext = 0;
	Vector2d dcurr = vprev?
		Vector2d(vcurr->x-vprev->x,vcurr->y-vprev->y) :
		Vector2d(0, -1);
	Vector2d dnext;
	bool isConvex = false;

	for(int i = 0; i < (int)vcurr->adjacent.size(); ++i) {
		Vertex* vadj = vcurr->adjacent[i];

		if(vadj == vprev) continue;

		Vector2d dadj = Vector2d(vadj->x-vcurr->x,vadj->y-vcurr->y);

		if(!vnext) {
			vnext = vadj;
			dnext = dadj;
			isConvex = (dnext[0]*dcurr[1] - dnext[1]*dcurr[0] <= 0);
			continue;
		}

		if(isConvex) {
			if(dcurr[0]*dadj[1] - dcurr[1]*dadj[0] < 0 ||
			   dnext[0]*dadj[1] - dnext[1]*dadj[0] < 0) {
				   vnext = vadj;
				   dnext = dadj;
				   isConvex = (dnext[0]*dcurr[1] - dnext[1]*dcurr[0] <= 0);
			}
		} else {
			if(dcurr[0]*dadj[1] - dcurr[1]*dadj[0] < 0 &&
			   dnext[0]*dadj[1] - dnext[1]*dadj[0] < 0) {
				   vnext = vadj;
				   dnext = dadj;
				   isConvex = (dnext[0]*dcurr[1] - dnext[1]*dcurr[0] <= 0);
			}
		}
	}
	return vnext;
}

Graph::Vertex* Graph::getCCW(Vertex* vprev, Vertex* vcurr) {
	Vertex* vnext = 0;
	Vector2d dcurr = vprev?
		Vector2d(vcurr->x-vprev->x,vcurr->y-vprev->y) :
		Vector2d(0, -1);
	Vector2d dnext;
	bool isConvex = false;

	for(int i = 0; i < (int)vcurr->adjacent.size(); ++i) {
		Vertex* vadj = vcurr->adjacent[i];

		if(vadj == vprev) continue;

		Vector2d dadj = Vector2d(vadj->x-vcurr->x,vadj->y-vcurr->y);

		if(!vnext) {
			vnext = vadj;
			dnext = dadj;
			isConvex = (dnext[0]*dcurr[1] - dnext[1]*dcurr[0] <= 0);
			continue;
		}

		if(isConvex) {
			if(dcurr[0]*dadj[1] - dcurr[1]*dadj[0] > 0 &&
			   dnext[0]*dadj[1] - dnext[1]*dadj[0] > 0) {
				   vnext = vadj;
				   dnext = dadj;
				   isConvex = (dnext[0]*dcurr[1] - dnext[1]*dcurr[0] <= 0);
			}
		} else {
			if(dcurr[0]*dadj[1] - dcurr[1]*dadj[0] > 0 ||
			   dnext[0]*dadj[1] - dnext[1]*dadj[0] > 0) {
				   vnext = vadj;
				   dnext = dadj;
				   isConvex = (dnext[0]*dcurr[1] - dnext[1]*dcurr[0] <= 0);
			}
		}
	}
	return vnext;
}

Graph::Vertex::Vertex(double X, double Y, int I) {
	x = X;
	y = Y;
	index = I;
}

Graph::Vertex::~Vertex(void) {}

void Graph::Vertex::insert(Vertex* adj) {
	adjacent.push_back(adj);
}

void Graph::Vertex::remove(Vertex* adj) {
	int numA = (int)adjacent.size();
	for(int i = 0; i < numA; ++i) {
		if(adj == adjacent[i]) {
			--numA;
			if(i < numA)
				adjacent[i] = adjacent[numA];
			adjacent.pop_back();
			return;
		}
	}
}

Graph::Primitive::Primitive(PrimitiveType t) {
	type = t;
}

Graph::VertexPtr::VertexPtr(Vertex* v) {
	mVertex = v;
}

bool Graph::VertexPtr::operator < (const VertexPtr &other) const {
	if(mVertex->x < other.mVertex->x)
		return true;
	if(mVertex->x > other.mVertex->x)
		return false;
	return mVertex->y < other.mVertex->y;
}