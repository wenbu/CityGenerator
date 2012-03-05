#include "StdAfx.h"
#include "InputParser.h"
#include <iostream>
#include <fstream>
using namespace std;

InputParser::InputParser(void) {
}

InputParser::~InputParser(void) {
}

InputParser::InputParser(char* fileName) {
	parse(fileName);
}

void InputParser::parse(char* fileName) {
	char line[1024];
	ifstream infile(fileName, ifstream::in);
	if(!infile) {
		cout << "Could not open given file " << fileName << endl;
		exit(1);
	}
	while(infile.good()) {
		infile.getline(line, 1023);
		if(!parseLine(line)) {
			cout << "Could not parse: " << line << endl;
		}
	}
	infile.close();
}

bool InputParser::parseLine(char* line) {
	string key;
	double value;
	if(string(line).empty())
		return true;
	stringstream ss(stringstream::in | stringstream::out);
	ss.str(line);
	ss >> key >> value;

	if(key.compare("QUERY_DIST") == 0) {
		params[QUERY_DIST] = value;
	} else if(key.compare("SNAP_THRESHOLD") == 0) {
		params[SNAP_THRESHOLD] = value;
	} else if(key.compare("FINAL_SNAP_THRESHOLD") == 0) {
		params[FINAL_SNAP_THRESHOLD] = value;
	} else if(key.compare("ROAD_LENGTH_MAJOR") == 0) {
		params[ROAD_LENGTH_MAJOR] = value;
	} else if(key.compare("ROAD_LENGTH_MINOR") == 0) {
		params[ROAD_LENGTH_MINOR] = value;
	} else if(key.compare("ROAD_LENGTH_HIGHWAY") == 0) {
		params[ROAD_LENGTH_HIGHWAY] = value;
	} else if(key.compare("ROAD_LENGTH_JITTER") == 0) {
		params[ROAD_LENGTH_JITTER] = value;
	} else if(key.compare("EVAL_DELAY_HIGHWAY") == 0) {
		params[EVAL_DELAY_HIGHWAY] = value;
	} else if(key.compare("EVAL_DELAY_MAJOR") == 0) {
		params[EVAL_DELAY_MAJOR] = value;
	} else if(key.compare("EVAL_DELAY_MINOR") == 0) {
		params[EVAL_DELAY_MINOR] = value;
	} else if(key.compare("CUTOFF_HIGHWAY_GROWTH") == 0) {
		params[CUTOFF_HIGHWAY_GROWTH] = value;
	} else if(key.compare("HIGHWAY_BRANCH_DELAY") == 0) {
		params[HIGHWAY_BRANCH_DELAY] = value;
	} else if(key.compare("NUM_HIGHWAY_PROBES") == 0) {
		params[NUM_HIGHWAY_PROBES] = value;
	} else if(key.compare("ANGLE_HIGHWAY_PROBES") == 0) {
		params[ANGLE_HIGHWAY_PROBES] = value;
	} else if(key.compare("HIGHWAY_PROBE_JITTER") == 0) {
		params[HIGHWAY_PROBE_JITTER] = value;
	} else if(key.compare("WEIGHT_TERRAIN_CONTOUR") == 0) {
		params[WEIGHT_TERRAIN_CONTOUR] = value;
	} else if(key.compare("WEIGHT_TERRAIN_CONTOUR_JITTER") == 0) {
		params[WEIGHT_TERRAIN_CONTOUR_JITTER] = value;
	} else if(key.compare("WEIGHT_TERRAIN_GRADIENT") == 0) {
		params[WEIGHT_TERRAIN_GRADIENT] = value;
	} else if(key.compare("WEIGHT_TERRAIN_GRADIENT_JITTER") == 0) {
		params[WEIGHT_TERRAIN_GRADIENT_JITTER] = value;
	} else if(key.compare("HIGHWAY_PROMOTION_PROBABILITY") == 0) {
		params[HIGHWAY_PROMOTION_PROBABILITY] = value;
	} else if(key.compare("GRID_DEVIATION_PROBABILITY") == 0) {
		params[GRID_DEVIATION_PROBABILITY] = value;
	} else if(key.compare("GRID_DEVIATION_AMOUNT") == 0) {
		params[GRID_DEVIATION_AMOUNT] = value;
	} else if(key.compare("SIDE_ROAD_PROBABILITY") == 0) {
		params[SIDE_ROAD_PROBABILITY] = value;
	} else if(key.compare("COST_THRESHOLD") == 0) {
		params[COST_THRESHOLD] = value;
	} else if(key.compare("COST_THRESHOLD_JITTER") == 0) {
		params[COST_THRESHOLD_JITTER] = value;
	} else if(key.compare("COST_INTERSECTION_OVERLOAD") == 0) {
		params[COST_INTERSECTION_OVERLOAD] = value;
	} else if(key.compare("COST_ROAD_SHORTENING") == 0) {
		params[COST_ROAD_SHORTENING] = value;
	} else if(key.compare("COST_ROAD_LENGTHENING") == 0) {
		params[COST_ROAD_LENGTHENING] = value;
	} else if(key.compare("BONUS_INTERSECTION_UNDERLOAD") == 0) {
		params[BONUS_INTERSECTION_UNDERLOAD] = value;
	} else if(key.compare("LOT_EDGE_MAX_WIDTH") == 0) {
		params[LOT_EDGE_MAX_WIDTH] = value;
	} else if(key.compare("LOT_MIN_AREA") == 0) {
		params[LOT_MIN_AREA] = value;
	} else if(key.compare("LOT_SPLIT_DEVIANCE") == 0) {
		params[LOT_SPLIT_DEVIANCE] = value;
	} else if(key.compare("BUILDING_MAX_HEIGHT") == 0) {
		params[BUILDING_MAX_HEIGHT] = value;
	} else {
		return false;
	}
	return true;
}

double InputParser::get(CityParam arg) {
	if(params.find(arg) != params.end()) {
		return params[arg];
	} else {
		return 0;
	}
}