#include "StdAfx.h"
#include "BuildingGenerator.h"
#include "ImageSampler.h"



using namespace std;
BuildingGenerator::BuildingGenerator(void)
{
}


BuildingGenerator::~BuildingGenerator(void)
{
}

BuildingGenerator::BuildingGenerator(vector<Lot> &l)
{
	int s = l.size();
	int i;
	for (i = 0; i<s; i++) {
		lots.push_back(l[i]);
	}

}

void BuildingGenerator::generateObjFile() {
	ofstream outputFile;
	outputFile.open("outputFile.obj");
	int b = buildings.size();
	for (int k = 0; k<b; k++) {
		if (buildings[k].size() == 10 ) {	
						outputFile << "v ";
			outputFile << buildings[k][0].vertices[0].x() << " ";
			outputFile << buildings[k][0].vertices[0].z() << " ";
			outputFile << buildings[k][0].vertices[0].y() << endl;

			outputFile << "v ";
			outputFile << buildings[k][0].vertices[1].x() << " ";
			outputFile << buildings[k][0].vertices[1].z() << " ";
			outputFile << buildings[k][0].vertices[1].y() << endl;

			outputFile << "v ";
			outputFile << buildings[k][0].vertices[2].x() << " ";
			outputFile << buildings[k][0].vertices[2].z() << " ";
			outputFile << buildings[k][0].vertices[2].y() << endl;

			outputFile << "v ";
			outputFile << buildings[k][0].vertices[3].x() << " ";
			outputFile << buildings[k][0].vertices[3].z() << " ";
			outputFile << buildings[k][0].vertices[3].y() << endl;

			//second side
			outputFile << "v ";
			outputFile << buildings[k][1].vertices[1].x() << " ";
			outputFile << buildings[k][1].vertices[1].z() << " ";
			outputFile << buildings[k][1].vertices[1].y() << endl;

			outputFile << "v ";
			outputFile << buildings[k][1].vertices[2].x() << " ";
			outputFile << buildings[k][1].vertices[2].z() << " ";
			outputFile << buildings[k][1].vertices[2].y() << endl;

			//third side [building #][side #]
			outputFile << "v ";
			outputFile << buildings[k][2].vertices[1].x() << " ";
			outputFile << buildings[k][2].vertices[1].z() << " ";
			outputFile << buildings[k][2].vertices[1].y() << endl;

			outputFile << "v ";
			outputFile << buildings[k][2].vertices[2].x() << " ";
			outputFile << buildings[k][2].vertices[2].z() << " ";
			outputFile << buildings[k][2].vertices[2].y() << endl;

			//southside roofside
			outputFile << "v ";
			outputFile << buildings[k][5].vertices[0].x() << " ";
			outputFile << buildings[k][5].vertices[0].z() << " ";
			outputFile << buildings[k][5].vertices[0].y() << endl;

			outputFile << "v ";
			outputFile << buildings[k][5].vertices[1].x() << " ";
			outputFile << buildings[k][5].vertices[1].z() << " ";
			outputFile << buildings[k][5].vertices[1].y() << endl;

			outputFile << "v ";
			outputFile << buildings[k][5].vertices[2].x() << " ";
			outputFile << buildings[k][5].vertices[2].z() << " ";
			outputFile << buildings[k][5].vertices[2].y() << endl;

			outputFile << "v ";
			outputFile << buildings[k][5].vertices[3].x() << " ";
			outputFile << buildings[k][5].vertices[3].z() << " ";
			outputFile << buildings[k][5].vertices[3].y() << endl;
						
			//second roofside
			outputFile << "v ";
			outputFile << buildings[k][6].vertices[1].x() << " ";
			outputFile << buildings[k][6].vertices[1].z() << " ";
			outputFile << buildings[k][6].vertices[1].y() << endl;

			outputFile << "v ";
			outputFile << buildings[k][6].vertices[2].x() << " ";
			outputFile << buildings[k][6].vertices[2].z() << " ";
			outputFile << buildings[k][6].vertices[2].y() << endl;

			//third roofside [building #][side #]
			outputFile << "v ";
			outputFile << buildings[k][7].vertices[1].x() << " ";
			outputFile << buildings[k][7].vertices[1].z() << " ";
			outputFile << buildings[k][7].vertices[1].y() << endl;

			outputFile << "v ";
			outputFile << buildings[k][7].vertices[2].x() << " ";
			outputFile << buildings[k][7].vertices[2].z() << " ";
			outputFile << buildings[k][7].vertices[2].y() << endl;
		}
		if (buildings[k].size() == 9 ) {
			outputFile << "v ";
			outputFile << buildings[k][0].vertices[0].x() << " ";
			outputFile << buildings[k][0].vertices[0].z() << " ";
			outputFile << buildings[k][0].vertices[0].y() << endl;

			outputFile << "v ";
			outputFile << buildings[k][0].vertices[1].x() << " ";
			outputFile << buildings[k][0].vertices[1].z() << " ";
			outputFile << buildings[k][0].vertices[1].y() << endl;

			outputFile << "v ";
			outputFile << buildings[k][0].vertices[2].x() << " ";
			outputFile << buildings[k][0].vertices[2].z() << " ";
			outputFile << buildings[k][0].vertices[2].y() << endl;

			outputFile << "v ";
			outputFile << buildings[k][0].vertices[3].x() << " ";
			outputFile << buildings[k][0].vertices[3].z() << " ";
			outputFile << buildings[k][0].vertices[3].y() << endl;

			//second side
			outputFile << "v ";
			outputFile << buildings[k][1].vertices[1].x() << " ";
			outputFile << buildings[k][1].vertices[1].z() << " ";
			outputFile << buildings[k][1].vertices[1].y() << endl;

			outputFile << "v ";
			outputFile << buildings[k][1].vertices[2].x() << " ";
			outputFile << buildings[k][1].vertices[2].z() << " ";
			outputFile << buildings[k][1].vertices[2].y() << endl;

			//third side [building #][side #]
			outputFile << "v ";
			outputFile << buildings[k][2].vertices[1].x() << " ";
			outputFile << buildings[k][2].vertices[1].z() << " ";
			outputFile << buildings[k][2].vertices[1].y() << endl;

			outputFile << "v ";
			outputFile << buildings[k][2].vertices[2].x() << " ";
			outputFile << buildings[k][2].vertices[2].z() << " ";
			outputFile << buildings[k][2].vertices[2].y() << endl;

			//southside roofside
			outputFile << "v ";
			outputFile << buildings[k][4].vertices[0].x() << " ";
			outputFile << buildings[k][4].vertices[0].z() << " ";
			outputFile << buildings[k][4].vertices[0].y() << endl;

			outputFile << "v ";
			outputFile << buildings[k][4].vertices[1].x() << " ";
			outputFile << buildings[k][4].vertices[1].z() << " ";
			outputFile << buildings[k][4].vertices[1].y() << endl;

			outputFile << "v ";
			outputFile << buildings[k][4].vertices[2].x() << " ";
			outputFile << buildings[k][4].vertices[2].z() << " ";
			outputFile << buildings[k][4].vertices[2].y() << endl;

			outputFile << "v ";
			outputFile << buildings[k][4].vertices[3].x() << " ";
			outputFile << buildings[k][4].vertices[3].z() << " ";
			outputFile << buildings[k][4].vertices[3].y() << endl;
						
			//second roofside
			outputFile << "v ";
			outputFile << buildings[k][5].vertices[1].x() << " ";
			outputFile << buildings[k][5].vertices[1].z() << " ";
			outputFile << buildings[k][5].vertices[1].y() << endl;

			outputFile << "v ";
			outputFile << buildings[k][5].vertices[2].x() << " ";
			outputFile << buildings[k][5].vertices[2].z() << " ";
			outputFile << buildings[k][5].vertices[2].y() << endl;

			//third roofside [building #][side #]
			outputFile << "v ";
			outputFile << buildings[k][6].vertices[1].x() << " ";
			outputFile << buildings[k][6].vertices[1].z() << " ";
			outputFile << buildings[k][6].vertices[1].y() << endl;

			outputFile << "v ";
			outputFile << buildings[k][6].vertices[2].x() << " ";
			outputFile << buildings[k][6].vertices[2].z() << " ";
			outputFile << buildings[k][6].vertices[2].y() << endl;

		}
		
		if (buildings[k].size() == 5) {
			outputFile << "v ";
			outputFile << buildings[k][0].vertices[0].x() << " ";
			outputFile << buildings[k][0].vertices[0].z() << " ";
			outputFile << buildings[k][0].vertices[0].y() << endl;

			outputFile << "v ";
			outputFile << buildings[k][0].vertices[1].x() << " ";
			outputFile << buildings[k][0].vertices[1].z() << " ";
			outputFile << buildings[k][0].vertices[1].y() << endl;

			outputFile << "v ";
			outputFile << buildings[k][0].vertices[2].x() << " ";
			outputFile << buildings[k][0].vertices[2].z() << " ";
			outputFile << buildings[k][0].vertices[2].y() << endl;

			outputFile << "v ";
			outputFile << buildings[k][0].vertices[3].x() << " ";
			outputFile << buildings[k][0].vertices[3].z() << " ";
			outputFile << buildings[k][0].vertices[3].y() << endl;

			//second side
			outputFile << "v ";
			outputFile << buildings[k][1].vertices[1].x() << " ";
			outputFile << buildings[k][1].vertices[1].z() << " ";
			outputFile << buildings[k][1].vertices[1].y() << endl;

			outputFile << "v ";
			outputFile << buildings[k][1].vertices[2].x() << " ";
			outputFile << buildings[k][1].vertices[2].z() << " ";
			outputFile << buildings[k][1].vertices[2].y() << endl;

			//third side [building #][side #]
			outputFile << "v ";
			outputFile << buildings[k][2].vertices[1].x() << " ";
			outputFile << buildings[k][2].vertices[1].z() << " ";
			outputFile << buildings[k][2].vertices[1].y() << endl;

			outputFile << "v ";
			outputFile << buildings[k][2].vertices[2].x() << " ";
			outputFile << buildings[k][2].vertices[2].z() << " ";
			outputFile << buildings[k][2].vertices[2].y() << endl;
		}
		if (buildings[k].size() == 4) {
			outputFile << "v ";
			outputFile << buildings[k][0].vertices[0].x() << " ";
			outputFile << buildings[k][0].vertices[0].z() << " ";
			outputFile << buildings[k][0].vertices[0].y() << endl;

			outputFile << "v ";
			outputFile << buildings[k][0].vertices[1].x() << " ";
			outputFile << buildings[k][0].vertices[1].z() << " ";
			outputFile << buildings[k][0].vertices[1].y() << endl;

			outputFile << "v ";
			outputFile << buildings[k][0].vertices[2].x() << " ";
			outputFile << buildings[k][0].vertices[2].z() << " ";
			outputFile << buildings[k][0].vertices[2].y() << endl;

			outputFile << "v ";
			outputFile << buildings[k][0].vertices[3].x() << " ";
			outputFile << buildings[k][0].vertices[3].z() << " ";
			outputFile << buildings[k][0].vertices[3].y() << endl;

			//second side
			outputFile << "v ";
			outputFile << buildings[k][1].vertices[1].x() << " ";
			outputFile << buildings[k][1].vertices[1].z() << " ";
			outputFile << buildings[k][1].vertices[1].y() << endl;

			outputFile << "v ";
			outputFile << buildings[k][1].vertices[2].x() << " ";
			outputFile << buildings[k][1].vertices[2].z() << " ";
			outputFile << buildings[k][1].vertices[2].y() << endl;
		}
	}

	int counter = 0;
	for (int f = 0; f<b;  f++) {
		if (buildings[f].size() == 10) {
						outputFile << "f ";
			outputFile << counter+1 << " ";
			outputFile << counter+2 << " ";
			outputFile << counter+3 << " ";
			outputFile << counter+4 << endl;

			outputFile << "f ";
			outputFile << counter+2 << " ";
			outputFile << counter+5 << " ";
			outputFile << counter+6 << " ";
			outputFile << counter+3 << endl;

			outputFile << "f ";
			outputFile << counter+5 << " ";
			outputFile << counter+7 << " ";
			outputFile << counter+8 << " ";
			outputFile << counter+6 << endl;
		
			outputFile << "f ";
			outputFile << counter+7 << " ";
			outputFile << counter+1 << " ";
			outputFile << counter+4 << " ";
			outputFile << counter+8 << endl;

			outputFile << "f ";
			outputFile << counter+4 << " ";
			outputFile << counter+3 << " ";
			outputFile << counter+6 << " ";
			outputFile << counter+8 << endl;

			outputFile << "f ";
			outputFile << counter+9 << " ";
			outputFile << counter+10 << " ";
			outputFile << counter+11 << " ";
			outputFile << counter+12 << endl;

			outputFile << "f ";
			outputFile << counter+10 << " ";
			outputFile << counter+13 << " ";
			outputFile << counter+14 << " ";
			outputFile << counter+11 << endl;

			outputFile << "f ";
			outputFile << counter+13 << " ";
			outputFile << counter+15 << " ";
			outputFile << counter+16 << " ";
			outputFile << counter+14 << endl;

			outputFile << "f ";
			outputFile << counter+15 << " ";
			outputFile << counter+9 << " ";
			outputFile << counter+12 << " ";
			outputFile << counter+16 << endl;

			outputFile << "f ";
			outputFile << counter+12 << " ";
			outputFile << counter+11 << " ";
			outputFile << counter+14 << " ";
			outputFile << counter+16 << endl;
			counter+=16;

		}
		if (buildings[f].size() == 9) {
			outputFile << "f ";
			outputFile << counter+1 << " ";
			outputFile << counter+2 << " ";
			outputFile << counter+3 << " ";
			outputFile << counter+4 << endl;

			outputFile << "f ";
			outputFile << counter+2 << " ";
			outputFile << counter+5 << " ";
			outputFile << counter+6 << " ";
			outputFile << counter+3 << endl;

			outputFile << "f ";
			outputFile << counter+5 << " ";
			outputFile << counter+7 << " ";
			outputFile << counter+8 << " ";
			outputFile << counter+6 << endl;
		
			outputFile << "f ";
			outputFile << counter+7 << " ";
			outputFile << counter+1 << " ";
			outputFile << counter+4 << " ";
			outputFile << counter+8 << endl;

			outputFile << "f ";
			outputFile << counter+9 << " ";
			outputFile << counter+10 << " ";
			outputFile << counter+11 << " ";
			outputFile << counter+12 << endl;

			outputFile << "f ";
			outputFile << counter+10 << " ";
			outputFile << counter+13 << " ";
			outputFile << counter+14 << " ";
			outputFile << counter+11 << endl;

			outputFile << "f ";
			outputFile << counter+13 << " ";
			outputFile << counter+15 << " ";
			outputFile << counter+16 << " ";
			outputFile << counter+14 << endl;

			outputFile << "f ";
			outputFile << counter+15 << " ";
			outputFile << counter+9 << " ";
			outputFile << counter+12 << " ";
			outputFile << counter+16 << endl;

			outputFile << "f ";
			outputFile << counter+12 << " ";
			outputFile << counter+11 << " ";
			outputFile << counter+14 << " ";
			outputFile << counter+16 << endl;
			counter+=16;
		}
		if (buildings[f].size() == 5) {
			outputFile << "f ";
			outputFile << counter+1 << " ";
			outputFile << counter+2 << " ";
			outputFile << counter+3 << " ";
			outputFile << counter+4 << endl;

			outputFile << "f ";
			outputFile << counter+2 << " ";
			outputFile << counter+5 << " ";
			outputFile << counter+6 << " ";
			outputFile << counter+3 << endl;

			outputFile << "f ";
			outputFile << counter+5 << " ";
			outputFile << counter+7 << " ";
			outputFile << counter+8 << " ";
			outputFile << counter+6 << endl;
		
			outputFile << "f ";
			outputFile << counter+7 << " ";
			outputFile << counter+1 << " ";
			outputFile << counter+4 << " ";
			outputFile << counter+8 << endl;

			outputFile << "f ";
			outputFile << counter+4 << " ";
			outputFile << counter+3 << " ";
			outputFile << counter+6 << " ";
			outputFile << counter+8 << endl;
		
			counter += 8;
		}
		if (buildings[f].size() == 4) {
			outputFile << "f ";
			outputFile << counter+1 << " ";
			outputFile << counter+2 << " ";
			outputFile << counter+3 << " ";
			outputFile << counter+4 << endl;

			outputFile << "f ";
			outputFile << counter+2 << " ";
			outputFile << counter+5 << " ";
			outputFile << counter+6 << " ";
			outputFile << counter+3 << endl;

			outputFile << "f ";
			outputFile << counter+5 << " ";
			outputFile << counter+1 << " ";
			outputFile << counter+4 << " ";
			outputFile << counter+6 << endl;

			outputFile << "f ";
			outputFile << counter+4 << " ";
			outputFile << counter+3 << " ";
			outputFile << counter+6 << endl;

			counter += 6;
		}
	}
	outputFile.close();
}

void BuildingGenerator::getBuildings(double maxHeight, ImageSampler &is) {
	int numLots = lots.size(); //how many buildings must be made and added to input
	int power = 8;
	for (int i = 0; i<numLots; i++) {
		Polygon p = lots[i].lot;
		if (p.vertices.size() >= 4 && p.vertices.size() <= 6) {
			//printf("%i verts; trying to make KOOL BUILDIN\n", p.vertices.size());
			// index of lot polygon: 0 = ll(lower left) , 1 = lr, 2 = ur, 3 = ul

			vector<Polygon> toPush;

			Vector3d ll = p.vertices[2] + (p.vertices[0]-p.vertices[2])*0.80; //ur is opposite
			Vector3d lr = p.vertices[3] + (p.vertices[1]-p.vertices[3])*0.80; //ul is opposite
			Vector3d ur = p.vertices[0] + (p.vertices[2]-p.vertices[0])*0.80; //ll is opposite
			Vector3d ul = p.vertices[1] + (p.vertices[3]-p.vertices[1])*0.80; //lr is opposite

			Vector3d llf = p.vertices[2] + (p.vertices[0]-p.vertices[2])*0.85; //ur is opposite
			Vector3d lrf = p.vertices[3] + (p.vertices[1]-p.vertices[3])*0.85; //ul is opposite
			Vector3d urf = p.vertices[0] + (p.vertices[2]-p.vertices[0])*0.85; //ll is opposite
			Vector3d ulf = p.vertices[1] + (p.vertices[3]-p.vertices[1])*0.85; //lr is opposite

			Vector3d lln = p.vertices[2] + (p.vertices[0]-p.vertices[2])*0.70; //ur is opposite
			Vector3d lrn = p.vertices[3] + (p.vertices[1]-p.vertices[3])*0.70; //ul is opposite
			Vector3d urn = p.vertices[0] + (p.vertices[2]-p.vertices[0])*0.70; //ll is opposite
			Vector3d uln = p.vertices[1] + (p.vertices[3]-p.vertices[1])*0.70; //lr is opposite

			if (p.vertices.size() == 6) {
				ll = p.vertices[2] + (p.vertices[0]-p.vertices[2])*0.80; //ur is opposite
				lr = p.vertices[4] + (p.vertices[1]-p.vertices[4])*0.80; //ul is opposite
				ur = p.vertices[0] + (p.vertices[2]-p.vertices[0])*0.80; //ll is opposite
				ul = p.vertices[1] + (p.vertices[4]-p.vertices[1])*0.80; //lr is opposite

				llf = p.vertices[2] + (p.vertices[0]-p.vertices[2])*0.85; //ur is opposite
				lrf = p.vertices[4] + (p.vertices[1]-p.vertices[4])*0.85; //ul is opposite
				urf = p.vertices[0] + (p.vertices[2]-p.vertices[0])*0.85; //ll is opposite
				ulf = p.vertices[1] + (p.vertices[4]-p.vertices[1])*0.85; //lr is opposite

				lln = p.vertices[2] + (p.vertices[0]-p.vertices[2])*0.70; //ur is opposite
				lrn = p.vertices[4] + (p.vertices[1]-p.vertices[4])*0.70; //ul is opposite
				urn = p.vertices[0] + (p.vertices[2]-p.vertices[0])*0.70; //ll is opposite
				uln = p.vertices[1] + (p.vertices[4]-p.vertices[1])*0.70; //lr is opposite
			} 

			double height;
			if (is.isValid(p)==1) {
				//height = lr.z()+(is.getDensity(p)/255)*maxHeight;
				//printf("\tsample at %f, %f\n", p.vertices[0][0], p.vertices[0][1]);
				height = lr.z()+(is.getDensity(p.vertices[0][0], p.vertices[0][1])/255)*maxHeight;
				height = pow(height, power);
				height /= pow(maxHeight, power-1);
				height = jitter(0.75*height, 0.25*height);
				height += 0.1*maxHeight;
				//printf("\tgot height = %f\n", height);
			}
			else{
				height = 0;
			}

			if (height >= (0.85*maxHeight)) {
				//printf("\tbuilding high top\n");
				//building top facing in positive z direction
				//make southern side, then add to(push_back) buildings[i] 
				Polygon southSide;
				southSide.vertices.push_back(ll); southSide.vertices.push_back(lr); 
				southSide.vertices.push_back(Vector3d(lr.x(), lr.y(), height)); 
				southSide.vertices.push_back(Vector3d(ll.x(), ll.y(), height));
				toPush.push_back(southSide);

				//make eastern side, then add to toPush
				Polygon westSide;
				westSide.vertices.push_back(lr); westSide.vertices.push_back(ur); 
				westSide.vertices.push_back(Vector3d(ur.x(), ur.y(), height)); 
				westSide.vertices.push_back(Vector3d(lr.x(), lr.y(), height));
				toPush.push_back(westSide);

				//make northern side, then add to toPush
				Polygon northSide;
				northSide.vertices.push_back(ur); northSide.vertices.push_back(ul); 
				northSide.vertices.push_back(Vector3d(ul.x(), ul.y(), height)); 
				northSide.vertices.push_back(Vector3d(ur.x(), ur.y(), height));
				toPush.push_back(northSide);

				//make western side, then add to toPush
				Polygon eastSide;
				eastSide.vertices.push_back(ul); eastSide.vertices.push_back(ll); 
				eastSide.vertices.push_back(Vector3d(ll.x(), ll.y(), height)); 
				eastSide.vertices.push_back(Vector3d(ul.x(), ul.y(), height));
				toPush.push_back(eastSide);

				//make top side, then add to toPush
				Polygon topSide;
				topSide.vertices.push_back(Vector3d(ll.x(), ll.y(), height)); 
				topSide.vertices.push_back(Vector3d(lr.x(), lr.y(), height)); 
				topSide.vertices.push_back(Vector3d(ur.x(), ur.y(), height)); 
				topSide.vertices.push_back(Vector3d(ul.x(), ul.y(), height));
				toPush.push_back(topSide);

				//top level

				Polygon southSidetop;
				southSidetop.vertices.push_back(lln); southSidetop.vertices.push_back(lrn); 
				southSidetop.vertices.push_back(Vector3d(lrn.x(), lrn.y(), 1.7*height)); 
				southSidetop.vertices.push_back(Vector3d(lln.x(), lln.y(), 1.7*height));
				toPush.push_back(southSidetop);

				//make eastern side, then add to toPush
				Polygon westSidetop;
				westSidetop.vertices.push_back(lrn); westSidetop.vertices.push_back(urn); 
				westSidetop.vertices.push_back(Vector3d(urn.x(), urn.y(), 1.7*height)); 
				westSidetop.vertices.push_back(Vector3d(lrn.x(), lrn.y(), 1.7*height));
				toPush.push_back(westSidetop);

				//make northern side, then add to toPush
				Polygon northSidetop;
				northSidetop.vertices.push_back(urn); northSidetop.vertices.push_back(uln); 
				northSidetop.vertices.push_back(Vector3d(uln.x(), uln.y(), 1.7*height)); 
				northSidetop.vertices.push_back(Vector3d(urn.x(), urn.y(), 1.7*height));
				toPush.push_back(northSidetop);

				//make western side, then add to toPush
				Polygon eastSidetop;
				eastSidetop.vertices.push_back(uln); eastSidetop.vertices.push_back(lln); 
				eastSidetop.vertices.push_back(Vector3d(lln.x(), lln.y(), 1.7*height)); 
				eastSidetop.vertices.push_back(Vector3d(uln.x(), uln.y(), 1.7*height));
				toPush.push_back(eastSidetop);

				//make top side, then add to toPush
				Polygon topSidetop;
				topSidetop.vertices.push_back(Vector3d(lln.x(), lln.y(), 1.7*height)); 
				topSidetop.vertices.push_back(Vector3d(lrn.x(), lrn.y(), 1.7*height)); 
				topSidetop.vertices.push_back(Vector3d(urn.x(), urn.y(), 1.7*height)); 
				topSidetop.vertices.push_back(Vector3d(uln.x(), uln.y(), 1.7*height));
				toPush.push_back(topSidetop);

				if (is.isValid(p)==1) {
					buildings.push_back(toPush);
				}
			}

			if (height >= (0.33*maxHeight) && height < 0.85*maxHeight) {
				//printf("\tbuilding medium top\n");
				//building top facing in positive z direction
				//make southern side, then add to(push_back) buildings[i] 
				Polygon southSide;
				southSide.vertices.push_back(ll); southSide.vertices.push_back(lr); 
				southSide.vertices.push_back(Vector3d(lr.x(), lr.y(), height)); 
				southSide.vertices.push_back(Vector3d(ll.x(), ll.y(), height));
				toPush.push_back(southSide);

				//make eastern side, then add to toPush
				Polygon westSide;
				westSide.vertices.push_back(lr); westSide.vertices.push_back(ur); 
				westSide.vertices.push_back(Vector3d(ur.x(), ur.y(), height)); 
				westSide.vertices.push_back(Vector3d(lr.x(), lr.y(), height));
				toPush.push_back(westSide);

				//make northern side, then add to toPush
				Polygon northSide;
				northSide.vertices.push_back(ur); northSide.vertices.push_back(ul); 
				northSide.vertices.push_back(Vector3d(ul.x(), ul.y(), height)); 
				northSide.vertices.push_back(Vector3d(ur.x(), ur.y(), height));
				toPush.push_back(northSide);

				//make western side, then add to toPush
				Polygon eastSide;
				eastSide.vertices.push_back(ul); eastSide.vertices.push_back(ll); 
				eastSide.vertices.push_back(Vector3d(ll.x(), ll.y(), height)); 
				eastSide.vertices.push_back(Vector3d(ul.x(), ul.y(), height));
				toPush.push_back(eastSide);

				//make top side, then add to toPush
				Polygon topSide;
				topSide.vertices.push_back(Vector3d(ll.x(), ll.y(), height)); 
				topSide.vertices.push_back(Vector3d(lr.x(), lr.y(), height)); 
				topSide.vertices.push_back(Vector3d(ur.x(), ur.y(), height)); 
				topSide.vertices.push_back(Vector3d(ul.x(), ul.y(), height));
				toPush.push_back(topSide);

				if (is.isValid(p)==1) {
					buildings.push_back(toPush);
				}

			}

			if (height < 0.33*maxHeight) {
				//printf("\tbuilding low top\n");
				Polygon southSide;
				southSide.vertices.push_back(ll); southSide.vertices.push_back(lr); 
				southSide.vertices.push_back(Vector3d(lr.x(), lr.y(), height)); 
				southSide.vertices.push_back(Vector3d(ll.x(), ll.y(), height));
				toPush.push_back(southSide);

				//make eastern side, then add to toPush
				Polygon westSide;
				westSide.vertices.push_back(lr); westSide.vertices.push_back(ur); 
				westSide.vertices.push_back(Vector3d(ur.x(), ur.y(), height)); 
				westSide.vertices.push_back(Vector3d(lr.x(), lr.y(), height));
				toPush.push_back(westSide);

				//make northern side, then add to toPush
				Polygon northSide;
				northSide.vertices.push_back(ur); northSide.vertices.push_back(ul); 
				northSide.vertices.push_back(Vector3d(ul.x(), ul.y(), height)); 
				northSide.vertices.push_back(Vector3d(ur.x(), ur.y(), height));
				toPush.push_back(northSide);

				//make western side, then add to toPush
				Polygon eastSide;
				eastSide.vertices.push_back(ul); eastSide.vertices.push_back(ll); 
				eastSide.vertices.push_back(Vector3d(ll.x(), ll.y(), height)); 
				eastSide.vertices.push_back(Vector3d(ul.x(), ul.y(), height));
				toPush.push_back(eastSide);
				//
				Polygon southRoof;
				southRoof.vertices.push_back(Vector3d(llf.x(), llf.y(), height)); 
				southRoof.vertices.push_back(Vector3d(lrf.x(), lrf.y(), height)); 
				southRoof.vertices.push_back(Vector3d(lrn.x(), lrn.y(), 1.2*height)); 
				southRoof.vertices.push_back(Vector3d(lln.x(), lln.y(), 1.2*height));
				toPush.push_back(southRoof);

				Polygon eastRoof;
				eastRoof.vertices.push_back(Vector3d(lrf.x(), lrf.y(), height)); 
				eastRoof.vertices.push_back(Vector3d(urf.x(), urf.y(), height)); 
				eastRoof.vertices.push_back(Vector3d(urn.x(), urn.y(), 1.2*height)); 
				eastRoof.vertices.push_back(Vector3d(lln.x(), lln.y(), 1.2*height));
				toPush.push_back(eastRoof);

				Polygon northRoof;
				northRoof.vertices.push_back(Vector3d(urf.x(), urf.y(), height)); 
				northRoof.vertices.push_back(Vector3d(ulf.x(), ulf.y(), height)); 
				northRoof.vertices.push_back(Vector3d(uln.x(), uln.y(), 1.2*height)); 
				northRoof.vertices.push_back(Vector3d(urn.x(), urn.y(), 1.2*height));
				toPush.push_back(northRoof);

				Polygon westRoof;
				westRoof.vertices.push_back(Vector3d(ulf.x(), ulf.y(), height)); 
				westRoof.vertices.push_back(Vector3d(llf.x(), llf.y(), height)); 
				westRoof.vertices.push_back(Vector3d(lln.x(), lln.y(), 1.2*height)); 
				westRoof.vertices.push_back(Vector3d(uln.x(), uln.y(), 1.2*height));
				toPush.push_back(westRoof);

				Polygon topSide;
				topSide.vertices.push_back(Vector3d(lln.x(), lln.y(), 1.2*height)); 
				topSide.vertices.push_back(Vector3d(lrn.x(), lrn.y(), 1.2*height)); 
				topSide.vertices.push_back(Vector3d(urn.x(), urn.y(), 1.2*height)); 
				topSide.vertices.push_back(Vector3d(uln.x(), uln.y(), 1.2*height));
				toPush.push_back(topSide);
			}


			if (is.isValid(p)==1) {
				buildings.push_back(toPush);
			}

		}
		if (p.vertices.size() == 3) {
			//printf("%i verts; stuck with boring building\n", p.vertices.size());
			vector<Polygon> toPush;

			Vector3d a = (p.vertices[1] + (p.vertices[2]-p.vertices[1])*0.5) + (p.vertices[0]-((p.vertices[1] + (p.vertices[2]-p.vertices[1])*0.5)))*0.80;
			Vector3d b = (p.vertices[2] + (p.vertices[0]-p.vertices[2])*0.5) + (p.vertices[1]-((p.vertices[2] + (p.vertices[0]-p.vertices[2])*0.5)))*0.80;
			Vector3d c = (p.vertices[0] + (p.vertices[1]-p.vertices[0])*0.5) + (p.vertices[2]-((p.vertices[0] + (p.vertices[1]-p.vertices[0])*0.5)))*0.80;

			double height;
			if (is.isValid(p)==1) {
				//height =  b.z()+(is.getDensity(p)/255)*maxHeight;
				//printf("\tsample at %f, %f\n", p.vertices[0][0], p.vertices[0][1]);
				height =  b.z()+(is.getDensity(p.vertices[0][0], p.vertices[0][1])/255)*maxHeight;
				height = pow(height, power);
				height /= pow(maxHeight, power-1);
				height = jitter(0.75*height, 0.25*height);
				height += 0.1*maxHeight;
				//printf("\tgot height = %f\n", height);
			}else {
				height = 0;
			}

			Polygon southSide;
			southSide.vertices.push_back(a); southSide.vertices.push_back(b); 
			southSide.vertices.push_back(Vector3d(b.x(), b.y(), height)); 
			southSide.vertices.push_back(Vector3d(a.x(), a.y(), height));
			toPush.push_back(southSide);

			Polygon eastSide;
			eastSide.vertices.push_back(b); eastSide.vertices.push_back(c); 
			eastSide.vertices.push_back(Vector3d(c.x(), c.y(), height));
			eastSide.vertices.push_back(Vector3d(b.x(), b.y(), height));
			toPush.push_back(eastSide);

			Polygon westSide;
			westSide.vertices.push_back(c); westSide.vertices.push_back(a); 
			westSide.vertices.push_back(Vector3d(a.x(), a.y(), height)); 
			westSide.vertices.push_back(Vector3d(c.x(), c.y(), height));
			toPush.push_back(westSide);

			Polygon topSide;
			topSide.vertices.push_back(Vector3d(a.x(), a.y(), height)); 
			topSide.vertices.push_back(Vector3d(b.x(), b.y(), height)); 
			topSide.vertices.push_back(Vector3d(c.x(), c.y(), height)); 
			toPush.push_back(topSide);

			if (is.isValid(p)==1) {
				buildings.push_back(toPush);
			}
		}
	}
}

void BuildingGenerator::getPolygons(vector<Polygon> &polys) const {
	for(unsigned i = 0; i < buildings.size(); i++) {
		for(unsigned j = 0; j < buildings[i].size(); j++) {
			polys.push_back(buildings[i][j]);
		}
	}
}

//returns a double in the range[c-r, c+r]
inline double BuildingGenerator::jitter(double c, double r) const {
	double x = 2*((double)rand()/RAND_MAX);
	x -= 1;
	return r*x + c;
}