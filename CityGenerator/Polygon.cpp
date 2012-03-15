#include "StdAfx.h"
#include "Polygon.h"


Polygon::Polygon(void)
{
}


Polygon::~Polygon(void)
{
}

void Polygon::insertCutPoints(vector<int> &indices, vector<Vector3d> &points) {
	printf("POLYGON:\n");
	for(unsigned i = 0, j = 0; i < indices.size(); i++) {
		while(j <= indices[i]) {
			printf("\t%i: [%f, %f]\n", j, vertices[j][0], vertices[j][1]);
			j++;
		}
		printf("*\t%i: [%f, %f]\n", indices[i], points[i][0], points[i][1]);
	}
	vector<Vector3d> newVerts;
	for(unsigned i = 0, j = 0; i < indices.size(); i++) {
		while(j <= indices[i]) {
			newVerts.push_back(vertices[j]);
			j++;
		}
		newVerts.push_back(points[i]);
	}
	vertices.clear();
	for(unsigned i = 0; i < newVerts.size(); i++)
		vertices.push_back(newVerts[i]);
	printf("result:\n");
	for(unsigned i = 0; i < vertices.size(); i++)
		printf("\t%i: [%f, %f]\n", i, vertices[i][0], vertices[i][1]);
	/*
	printf("\n------------------------------------\n");
	printf("POLYGON: inserting cut points\n");
	unsigned j = 0;
	for(unsigned i = 0; i < indices.size(); i++) {
		while(j <= indices[i]) {
			printf("\t%i: [%f, %f]\n", j, vertices[j][0], vertices[j][1]);
			j++;
		}
		printf("*\t%i:  [%f, %f]\n", indices[i], points[i][0], points[i][1]);
	}
	//assume indices.size == points.size; max(indices) <= vertices.size
	for(unsigned i = 0; i < indices.size(); i++) {
		int index = indices[i];
		Vector3d point = points[i];

		//want to insert point right after vertices[index];
		vertices.insert(vertices.begin()+index+1, point);

		//now we need to update indices: all entries > index need to be incremented
		//O(n^2) rofl
		for(unsigned j = 0; j < indices.size(); j++) {
			if(indices[j] > index)
				indices[j] += 1;
		}
		i += 1;
		printf("intermediate:\n");
		for(unsigned j = 0; j < vertices.size(); j++) {
			printf("\t%i: [%f, %f]\n", j, vertices[j][0], vertices[j][1]);
		}
	}
	printf("result:\n");
	for(unsigned i = 0; i < vertices.size(); i++) {
		printf("\t%i: [%f, %f]\n", i, vertices[i][0], vertices[i][1]);
	}
	*/
}
