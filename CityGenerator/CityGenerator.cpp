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
//#include <stdlib.h>
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

double scale =		1;
double translateX =	0;
double translateY =	0;
double rotateX =	0;
double rotateY =	0;
double rotateZ =	0;

//****************************************************
// reshape viewport if the window is resized
//****************************************************
void myReshape(int w, int h) {
	viewport.w = w;
	viewport.h = h;

	glViewport(0,0,viewport.w,viewport.h);// sets the rectangle that will be the window
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();				// loading the identity matrix for the screen

	//----------- setting the projection -------------------------
	// glOrtho sets left, right, bottom, top, zNear, zFar of the chord system


	//glOrtho(-1, 1 + (w-400)/200.0 , -1 -(h-400)/200.0, 1, 1, -1); // resize type = add
	// glOrtho(-w/400.0, w/400.0, -h/400.0, h/400.0, 1, -1); // resize type = center

	glOrtho(-1, 1, -1, 1, -100, 100);	// resize type = stretch

	//------------------------------------------------------------
}

//****************************************************
// sets the window up
//****************************************************
void initScene(int &argc, char* argv[]){
	glClearColor(1, 1, 1, 0.0f); // Clear to black, fully transparent
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_NORMALIZE);
	glEnable(GL_BLEND);
	glDepthMask(GL_TRUE);
	glDepthFunc(GL_LEQUAL);

	//set up road system
	roadsystem = RoadSystem(is, parser);
	roadsystem.addAxiom(0.5, 0.5, 0.5125, 0.5); //coords?
	//roadsystem.addAxiom(0.5, 0.8, 0.5125, 0.8);
	//roadsystem.addAxiom(0.24, 0.35, 0.2525, 0.35);
	//roadsystem.addAxiom(0.76, 0.35, 0.7725, 0.35);
	printf("-------------------------------------\n");
	printf("Road Generation\n");
	printf("(Space) for one iteration, (x) for 100, (z) for 1000\n");
	printf("(i, k) to zoom, (arrows) to pan\n");
	printf("(d) to dump roads to MEL\n");
	printf("(p) to extract city blocks and proceed to lot subdivision\n");
}

void drawPolygon(Polygon &p) {
	glBegin(GL_POLYGON);
		for(unsigned i = 0; i < p.vertices.size(); i++) {
			glVertex3f(p.vertices[i][0],
					   p.vertices[i][1],
					   p.vertices[i][2]);
		}
	glEnd();
}

void drawWirePolygon(Polygon &p) {
	//printf("p: \n");
	for(unsigned i = 0; i < p.vertices.size(); i++) {
		//printf("    %f, %f\n", p.vertices[i][0], p.vertices[i][1]);
		glBegin(GL_LINE_STRIP);
		glVertex3f(p.vertices[i][0],
				   p.vertices[i][1],
				   p.vertices[i][2]);
			if(i != p.vertices.size()-1) {
				glVertex3f(p.vertices[i+1][0],
						   p.vertices[i+1][1],
						   p.vertices[i+1][2]);
			} else {
				glVertex3f(p.vertices[0][0],
						   p.vertices[0][1],
						   p.vertices[0][2]);
			}
		glEnd();
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
		//roadsystem.draw();
		//printf("   lots.size = %i\n", lots.size());
		//printf("lg.lots.size = %i\n", lg.lots.size());
		for(unsigned i = 0; i < lots.size(); i++) {
			//glColor3f(0, 1, 0);
			glColor3f(0.16, 0.64, 0.32);
			drawPolygon(lg.lots[i].lot);
			glColor3f(0,0,0);
			drawWirePolygon(lg.lots[i].lot);
		}
	}
}

void dumpLots() {
	if(stage == CS_LOTS && gotLots) {
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
}

//***************************************************
// function that does the actual drawing
//***************************************************
void myDisplay() {
	
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);				// clear the color buffer (sets everything to black)
	glMatrixMode(GL_MODELVIEW);					// indicate we are specifying camera transformations
	glLoadIdentity();							// make sure transformation is "zero'd"
	
	//----------------------- code to draw objects --------------------------
	glScalef(2*scale, 2*scale, 2*scale);
	glRotatef(rotateZ, 0, 0, 1);
	glRotatef(rotateY, 0, 1, 0);
	glRotatef(rotateX, 1, 0, 0);
	glTranslatef(translateX-0.5, translateY-0.5, 0);

	/*
	glBegin(GL_LINE_STRIP);
		glColor3f(1,0,0);
		glVertex3f(0,0,0);
		glVertex3f(10,0,0);
	glEnd();
	glBegin(GL_LINE_STRIP);
		glColor3f(0,1,0);
		glVertex3f(0,0,0);
		glVertex3f(0,10,0);
	glEnd();
	glBegin(GL_LINE_STRIP);
		glColor3f(0,0,1);
		glVertex3f(0,0,0);
		glVertex3f(0,0,10);
	glEnd();
	*/
	drawStuff();
	//-----------------------------------------------------------------------
	
	glFlush();
	glutSwapBuffers();					// swap buffers (we earlier set double buffer)
}

void myFrameMove() {
	//nothing here for now
#ifdef _WIN32
	//Sleep(10);						//give ~10ms back to OS (so as not to waste the CPU)
#endif
	glutPostRedisplay(); // forces glut to call the display function (myDisplay())
}

void grow(int n) {
	roadsystem.generate(n);
	glutPostRedisplay();
	if(roadsystem.nIterations >= 20000) {
		roadsystem.reset();
		roadsystem.addAxiom(0.5, 0.5, 0.502, 0.5);
		//roadsystem.addAxiom(0.5, 0.8, 0.502, 0.8);
		//roadsystem.addAxiom(0.24, 0.35, 0.242, 0.35);
		//roadsystem.addAxiom(0.76, 0.35, 0.762, 0.35);
	}
	glutTimerFunc(50, grow, n);
}

void mySpecial(int key, int x, int y) {
	int mod = glutGetModifiers();
	if(key == GLUT_KEY_LEFT) {
		if(mod == GLUT_ACTIVE_SHIFT)
			rotateX += 1;
		else
			translateX += 0.1;
	}
	if(key == GLUT_KEY_RIGHT) {
		if(mod == GLUT_ACTIVE_SHIFT)
			rotateX -= 1;
		else
			translateX -= 0.1;
	}
	if(key == GLUT_KEY_UP) {
		if(mod == GLUT_ACTIVE_SHIFT)
			rotateY += 1;
		else
			translateY += 0.1;
	}
	if(key == GLUT_KEY_DOWN) {
		if(mod == GLUT_ACTIVE_SHIFT)
			rotateY -= 1;
		else
			translateY -= 0.1;
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
	if(key == 'd') {
		if(stage == CS_ROADS) {
			roadsystem.dumpRoads();
		}
		if(stage == CS_LOTS) {
			roadsystem.dumpPolygons();
		}
	}
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
			printf("(d) to dump polygons\n");
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
			printf("output will be written to outputFile.obj\n");
			break;
		case(CS_BUILDINGS):
			printf("generating buildings...");
			bg = BuildingGenerator(lots);
			bg.getBuildings(parser.get(BUILDING_MAX_HEIGHT), is);
			bg.generateObjFile();
			printf("done!\n");
			printf("(p) to quit\n");
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

  	initScene(argc, argv);							// quick function to set up scene
  
  	glutDisplayFunc(myDisplay);				// function to run when its time to draw something
  	glutReshapeFunc(myReshape);				// function to run when the window gets resized
  	glutIdleFunc(myFrameMove);				// function to run when not handling any other task
	glutKeyboardFunc(keyboard);
	glutSpecialFunc(mySpecial);
  	glutMainLoop();							// infinite loop that will keep drawing and resizing and whatever else
	
	is.done();
	return 0;
}

