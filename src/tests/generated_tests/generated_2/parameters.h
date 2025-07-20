#ifndef DGGML_PARAMETERS_HPP
#define DGGML_PARAMETERS_HPP
#include <string>
#include "simdjson.h"
struct Parameters {
	int NUM_STEPS = 100;
	double DELTA = 0.4;
	double MAXIMAL_REACTION_RADIUS = 0.1;
	int MIN_DELTA_STEPS = 5;
	double DELTA_DELTA_T = DELTA / MIN_DELTA_STEPS;
	double DELTA_T_MIN = DELTA_DELTA_T;
	double CELL_NX = 10;
	double CELL_NY = 10;
	double CELL_DX = 1.66;
	double CELL_DY = 1.66;
	double total_time = 100;
};
#endif