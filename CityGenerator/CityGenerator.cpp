// CityGenerator.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "FreeImage.h"
#include "RoadSystem.h"
#include "LotGenerator.h"
#include "BuildingGenerator.h"
#include "InputParser.h"

#include <Windows.h>
#include <string.h>
#include <fstream>
#include <vector>
#include "targetver.h"
#include <GL/glut.h>
#include <GL/glu.h>
#include <iostream>
#include <time.h>
#include <math.h>

#ifdef _WIN32
static DWORD lastTime;
#else
static struct timeval lastTime;
#endif

using namespace std;
using namespace Eigen;

class Viewport {
public:
	int w, h; // width and height
};

enum CityStage {
	CS_ROADS,
	CS_LOTS,
	CS_BUILDINGS,
	CS_DONE
};

Viewport			viewport;
RoadSystem			roadsystem;
LotGenerator		lg;
BuildingGenerator	bg;
ImageSampler		is;
InputParser			parser;
CityStage			stage = CS_ROADS;

vector<Polygon>		blocks;
vector<Lot>			lots;
vector<Polygon>		buildings;

bool gotBlocks =	false;
bool gotLots =		false;
bool debug =		false;

double scale =			1;
double translateX =		0;
double translateY =		0;
double rotateX =		0;
double rotateY =		0;
double rotateZ =		0;
double wMult =			1;
double debugOffset =	0;

//****************************************************
// reshape viewport if the window is resized
//****************************************************
void myReshape(int w, int h) {
	viewport.w = w;
	viewport.h = h;

	glViewport(0,0,viewport.w,viewport.h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-1, 1, -1, 1, -100, 100);

	//------------------------------------------------------------
}

//****************************************************
// sets the window up
//****************************************************
void initScene(int &argc, char* argv[]){
	glClearColor(1,1,1, 0.0f);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_NORMALIZE);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glDepthMask(GL_TRUE);
	glDepthFunc(GL_LEQUAL);

	//set up road system
	roadsystem = RoadSystem(is, parser);
	roadsystem.addAxiom(0.5, 0.5, 0.5125, 0.5);
	printf("-------------------------------------\n");
	printf("Road Generation\n");
	printf("(Space) for one iteration, (x) for 100, (z) for 1000\n");
	printf("(i, k) to zoom, (arrows) to pan\n");
	printf("(d) to dump roads to roads.mel\n");
	printf("(p) to extract city blocks and proceed to lot subdivision\n");
}

void drawPolygon(Polygon &p, double zOffset = 0) {
	glBegin(GL_POLYGON);
		for(unsigned i = 0; i < p.vertices.size(); i++) {
			glVertex3f(p.vertices[i][0],
					   p.vertices[i][1],
					   p.vertices[i][2] + zOffset);
		}
	glEnd();
}

void drawWirePolygon(Polygon &p, double zOffset = 0) {
	for(unsigned i = 0; i < p.vertices.size(); i++) {
		glBegin(GL_LINE_STRIP);
		glVertex3f(p.vertices[i][0],
				   p.vertices[i][1],
				   p.vertices[i][2] + zOffset);
			if(i != p.vertices.size()-1) {
				glVertex3f(p.vertices[i+1][0],
						   p.vertices[i+1][1],
						   p.vertices[i+1][2] + zOffset);
			} else {
				glVertex3f(p.vertices[0][0],
						   p.vertices[0][1],
						   p.vertices[0][2] + zOffset);
			}
		glEnd();
	}
	if(debug) {
		glPointSize(3);
		glColor3f(0, 0, 0);
		for(unsigned i = 0; i < p.vertices.size(); i++) {
			glBegin(GL_POINTS);
			glVertex3f(p.vertices[i][0],
					   p.vertices[i][1],
					   p.vertices[i][2] + zOffset);
			glEnd();
		}
	}
}   

void drawStuff() {
	roadsystem.draw(gotBlocks);

	if(gotBlocks) {
		glColor3f(0.1, 0.1, 0.1);
		for(unsigned i = 0; i < blocks.size(); i++) {
			drawWirePolygon(blocks[i]);
		}
	}
	if(gotLots) {
		for(unsigned i = 0; i < lots.size(); i++) {
			glColor4f(0.16, 0.64, 0.32, 0.2);
			drawPolygon(lots[i].lot, i*debugOffset);
			glColor3f(0,0,0);
			drawWirePolygon(lots[i].lot, i*debugOffset);
		}
	}
	if(debug) {
		Polygon p;
		p.vertices.push_back(Vector3d(0.398803, 0.363752, 0));
		p.vertices.push_back(Vector3d(0.407103, 0.352854, 0));
		p.vertices.push_back(Vector3d(0.414602, 0.357363, 0));
		p.vertices.push_back(Vector3d(0.408383, 0.371349, 0));
		p.vertices.push_back(Vector3d(0.401586, 0.36584, 0));
		drawWirePolygon(p);
		lg.draw();
	}
}

void dumpLots() {
	ofstream outputFile;
	outputFile.open("lots.obj");
	for(unsigned i = 0; i < lots.size(); i++) {
		for(unsigned j = 0; j < lots[i].lot.vertices.size(); j++) {
			outputFile << "v ";
			outputFile << lots[i].lot.vertices[j][0] << " ";
			outputFile << lots[i].lot.vertices[j][2] << " ";
			outputFile << lots[i].lot.vertices[j][1] << " ";
			outputFile << endl;
		}
	}
	int vertexCounter = 0;
	for(unsigned i = 0; i < lots.size(); i++) {
		outputFile << "f ";
		for(unsigned j = 0; j < lots[i].lot.vertices.size(); j++) 
			outputFile << ++vertexCounter << " ";
		outputFile << endl;
	}
	outputFile.close();
}

//***************************************************
// function that does the actual drawing
//***************************************************
void myDisplay() {
	
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	glScalef(2*scale, 2*scale, 2*scale);
	glRotatef(rotateZ, 0, 0, 1);
	glRotatef(rotateY, 0, 1, 0);
	glRotatef(rotateX, 1, 0, 0);
	glTranslatef(translateX-0.5, translateY-0.5, 0);

	drawStuff();
	
	glFlush();
	glutSwapBuffers();
}

void myFrameMove() {
#ifdef _WIN32
	//Sleep(10);						//give ~10ms back to OS (so as not to waste the CPU)
#endif
	glutPostRedisplay();
}

void grow(int n) {
	roadsystem.generate(n);
	glutPostRedisplay();
	if(roadsystem.nIterations >= 20000) {
		roadsystem.reset();
		roadsystem.addAxiom(0.5, 0.5, 0.502, 0.5);
	}
	glutTimerFunc(50, grow, n);
}

void debugLots() {
	lots.clear();
	debugOffset = 0.005;
	vector<Polygon> debugBlocks;
	Polygon p;
	p.vertices.push_back(Vector3d(0.398803, 0.363752, 0));
	p.vertices.push_back(Vector3d(0.407103, 0.352854, 0));
	p.vertices.push_back(Vector3d(0.414602, 0.357363, 0));
	p.vertices.push_back(Vector3d(0.408383, 0.371349, 0));
	p.vertices.push_back(Vector3d(0.401586, 0.36584, 0));
	debugBlocks.push_back(p);
	lg = LotGenerator(debugBlocks,
									parser.get(LOT_EDGE_MAX_WIDTH)*wMult,
									0, //area
									parser.get(LOT_SPLIT_DEVIANCE),
									is);
	lg.getLots(lots);
	gotLots = true;
	debug = true;
	stage = CS_BUILDINGS;
}

void mySpecial(int key, int x, int y) {
	int mod = glutGetModifiers();
	if(key == GLUT_KEY_LEFT) {
		if(mod == GLUT_ACTIVE_SHIFT)
			rotateX += 1;
		else
			translateX += (0.1/scale);
	}
	if(key == GLUT_KEY_RIGHT) {
		if(mod == GLUT_ACTIVE_SHIFT)
			rotateX -= 1;
		else
			translateX -= (0.1/scale);
	}
	if(key == GLUT_KEY_UP) {
		if(mod == GLUT_ACTIVE_SHIFT)
			rotateY += 1;
		else
			translateY += (0.1/scale);
	}
	if(key == GLUT_KEY_DOWN) {
		if(mod == GLUT_ACTIVE_SHIFT)
			rotateY -= 1;
		else
			translateY -= (0.1/scale);
	}
}

void keyboard(unsigned char key, int x, int y) {
	if(key == ' ') {
		roadsystem.generate(1);
	}
	if(key == 'x') {
		roadsystem.generate(100);
	}
	if(key == 'z') {
		roadsystem.generate(1000);
	}
	if(key == 'q') {
		roadsystem.printStats();
	}
	if(key == 'i') {
		scale *= 1.1;
	}
	if(key == 'k') {
		scale /= 1.1;
	}
	/*
	if(key == 'u') {
		wMult += 0.05;
		debugLots();
	}
	if(key == 'j') {
		wMult -= 0.05;
		debugLots();
	}
	*/
	if(key == 'd') {
		if(stage >= CS_ROADS) {
			printf("dumping roads to roads.mel...");
			roadsystem.dumpRoads();
			printf("done!\n");
		}
	}
	if(key == 'f') {
		if(stage >= CS_LOTS) {
			printf("dumping blocks to blocks.obj...");
			roadsystem.dumpPolygons();
			printf("done!\n");
		}
	}
	if(key == 'g') {
		if(stage >= CS_BUILDINGS) {
			printf("dumping lots to lots.obj...");
			dumpLots();
			printf("done!\n");
		}
	}
	/*
	if(key == 'f') {
		printf("debugging\n");
		debugLots();
	}
	*/
	if(key == 'p') {
		switch(stage) {
		case(CS_ROADS):
			printf("extracting blocks...");
			roadsystem.getPolygons();
			blocks.clear();
			roadsystem.extractPolygons(blocks);
			gotBlocks = true;
			printf("done!\n");
			stage = CS_LOTS;
			printf("-------------------------------------\n");
			printf("Lot Subdivision\n");
			printf("(p) to generate lots\n");
			printf("(f) to dump blocks to blocks.obj\n");
			printf("(d) to dump roads to roads.mel\n");
			printf("The code likes to blow up at this stage :(\n");
			printf("If it does, just restart the program\n");
			break;
		case(CS_LOTS):
			lg = LotGenerator(blocks,
							  parser.get(LOT_EDGE_MAX_WIDTH),
							  parser.get(LOT_MIN_AREA),
							  parser.get(LOT_SPLIT_DEVIANCE),
							  is);
			lg.getLots(lots);
			gotLots = true;
			printf("Done subdividing blocks!\n");
			stage = CS_BUILDINGS;
			printf("-------------------------------------\n");
			printf("Building Generation\n");
			printf("(p) to generate buildings\n");
			printf("(g) to dump lots to lots.obj\n");
			printf("(f) to dump blocks to blocks.obj\n");
			printf("(d) to dump roads to roads.mel\n");
			printf("output will be written to outputFile.obj\n");
			break;
		case(CS_BUILDINGS):
			printf("generating buildings...");
			bg = BuildingGenerator(lots);
			bg.getBuildings(parser.get(BUILDING_MAX_HEIGHT), is);
			bg.generateObjFile();
			printf("done!\n");
			printf("(p) to quit\n");
			printf("(g) to dump lots to lots.obj\n");
			printf("(f) to dump blocks to blocks.obj\n");
			printf("(d) to dump roads to roads.mel\n");
			stage = CS_DONE;
			break;
		case(CS_DONE):
			exit(0);
			break;
		}
	}
	glutPostRedisplay();
}

int main(int argc, char* argv[])
{
	char* density = "";
	char* land = "";
	char* elevation = "";
	char* input = "";

	if(argc == 1) {
		cout << "Using default parameters" << endl;
		density = "density.png";
		land = "legality.png";
		elevation = "elevation.png";
		input = "input.txt";
	}

	else if(argc < 9) {
		cout << "Not enough args" << endl;
		exit(1);
	}
	else {
		for(int i = 1; i <= 8; i=i+2) {
			if(strcmp(argv[i],"-f") == 0) {
				input = argv[i+1];
			} else if(strcmp(argv[i],"-d") == 0) {
				density = argv[i+1];
			} else if(strcmp(argv[i],"-l") == 0) {
				land = argv[i+1];
			} else if(strcmp(argv[i],"-e") == 0) {
				elevation = argv[i+1];
			}
		}

		if(strcmp(input,"") == 0 || strcmp(density,"") == 0 || 
			strcmp(land,"") == 0 || strcmp(elevation,"") == 0) {
				cout << "Failed to initialize all args" << endl;
				exit(1);
		}
	}

	parser = InputParser(input);
	is = ImageSampler(density, elevation, land);

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	viewport.w = 600;
	viewport.h = 600;
	glutInitWindowSize(viewport.w, viewport.h);
  	glutInitWindowPosition(0, 0);
  	glutCreateWindow("City Generator");

	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  	initScene(argc, argv);
  
  	glutDisplayFunc(myDisplay);
  	glutReshapeFunc(myReshape);
  	glutIdleFunc(myFrameMove);
	glutKeyboardFunc(keyboard);
	glutSpecialFunc(mySpecial);
  	glutMainLoop();
	
	is.done();
	return 0;
}

