#ifndef DGGML_PARAMETERS_HPP
#define DGGML_PARAMETERS_HPP
#include <string>
#include "simdjson.h"
struct Parameters {
	double TOTAL_TIME = 400;
	double DELTA = 0.1;
	double MAXIMAL_REACTION_RADIUS = 0.1;
	double varepsilon = MAXIMAL_REACTION_RADIUS;
	int NUM_STEPS = TOTAL_TIME / DELTA;
	int MIN_DELTA_STEPS = 5;
	double DELTA_DELTA_T = DELTA / MIN_DELTA_STEPS;
	double DELTA_T_MIN = DELTA_DELTA_T;
	double CELL_NX = 3;
	double CELL_NY = 3;
	double CELL_NZ = 3;
	double CELL_DX = 1;
	double CELL_DY = 1;
	double CELL_DZ = 1;
	double rho_create = 0.026;
	double rho_grow = 100.0;
	double rho_retract = 10.0;
	double L_DIV = 0.075;
	double L_min = 0.025;
	double s_min = 0.005;
	double s_max = 0.01;
	double buffer = CELL_NX / 10.0;
	double offset = (3.0 - buffer * 2) / 10.0;
	int boundary_pts = 21;
	double boundary_offset = 3.0 / 20.0;
	double creation_factor = 1.0;
	double creation_rate = 0.0026;
	double mt_min_segment_init = 0.005;
	double mt_max_segment_init = 0.01;
	double v_plus = 0.0615;
	double collision_distance = 5.0;
	double boundary_cic = 4000000;
	double destruction_factor = 0.0026;
	double gr_to_ret = 0.016;
	double ret_to_gr = 0.016;
	double s_col = 0.025;
	double rho_int_cic = 12000;
	double rho_grow_cic = 12000;
	double rho_retract_cic = 12000;
	double theta_cic = 40.0;
	double theta_cross = 40.0;
	double rho_cross = 200.0 * 100.0;
	double rho_zip_hit = 4000.0;
	double rho_zip_guard = 12000.0;
};
#endif