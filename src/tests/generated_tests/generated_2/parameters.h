#ifndef DGGML_PARAMETERS_HPP
#define DGGML_PARAMETERS_HPP
#include <string>
#include "simdjson.h"
struct Parameters {
	double TOTAL_TIME = 10.0;
	double DELTA = 0.1;
	double MAXIMAL_REACTION_RADIUS = 0.4;
	int NUM_STEPS = TOTAL_TIME / DELTA;
	int MIN_DELTA_STEPS = 5;
	double DELTA_DELTA_T = DELTA / MIN_DELTA_STEPS;
	double DELTA_T_MIN = DELTA_DELTA_T;
	double NX_MIN = 8;
	double DX = 3 / NX_MIN;
	double CELL_NX = NX_MIN;
	double CELL_NY = NX_MIN;
	double CELL_NZ = NX_MIN;
	double CELL_DX = DX;
	double CELL_DY = DX;
	double CELL_DZ = DX;
	double varepsilon = 0.02;
};
#endif