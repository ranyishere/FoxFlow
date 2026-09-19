#ifndef DGGML_PARAMETERS_HPP
#define DGGML_PARAMETERS_HPP
#include <string>
#include "simdjson.h"
struct Parameters {
	double TOTAL_TIME = 400;
	double DELTA = 0.4;
	double MAXIMAL_REACTION_RADIUS = 0.1;
	double varepsilon = MAXIMAL_REACTION_RADIUS;
	int NUM_STEPS = TOTAL_TIME / DELTA;
	int MIN_DELTA_STEPS = 5;
	double DELTA_DELTA_T = DELTA / MIN_DELTA_STEPS;
	double DELTA_T_MIN = DELTA_DELTA_T;
	double CELL_NX = 10;
	double CELL_NY = 10;
	double CELL_DX = 1;
	double CELL_DY = 1;
	double rho_create = 0.026;
	double L_DIV = 0.075;
	double s_min = 0.005;
	double s_max = 0.01;
};
#endif