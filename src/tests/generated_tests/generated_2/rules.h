#ifndef DGGML_RULES_HPP
#define DGGML_RULES_HPP
#include "types.h"
#include "parameters.h"
#include "functions.h"
namespace FractureNetwork {
using GT = Microtubule::graph_type;
void start_rock_prop(DGGML::Grammar<Microtubule::graph_type> &gamma,
           Microtubule::graph_type &system_graph,
           Parameters &settings) {

GT start_rock_prop_lhs;
start_rock_prop_lhs.addNode({1, {Microtubule::StartType{} }});

GT start_rock_prop_rhs;
start_rock_prop_rhs.addNode({1, {Microtubule::RockStart{} }});

DGGML::WithRule<GT> start_rock_prop("start_rock_prop", start_rock_prop_lhs, start_rock_prop_rhs,
[&](auto &lhs, auto &m) {

return 
(
10000
)
;},
[&](auto &lhs, auto &rhs, auto &m1, auto &m2) {
double offset_distance = (0.1);
double lx = (0);
		rhs[m2[1]].position[0] = lx;
		std::get<Microtubule::RockStart>(rhs[m2[1]].data).fflow_a9c368[0] = lx;
double ly = (0);
		rhs[m2[1]].position[1] = ly;
		std::get<Microtubule::RockStart>(rhs[m2[1]].data).fflow_a9c368[1] = ly;
double lz = (0);
		rhs[m2[1]].position[2] = lz;
		std::get<Microtubule::RockStart>(rhs[m2[1]].data).fflow_a9c368[2] = lz;
double cx = (30);
		std::get<Microtubule::RockStart>(rhs[m2[1]].data).fflow_a13208 = cx;
double cy = (30);
		std::get<Microtubule::RockStart>(rhs[m2[1]].data).fflow_992eaa = cy;
double cz = (30);
		std::get<Microtubule::RockStart>(rhs[m2[1]].data).fflow_d20285 = cz;
}
);
gamma.addRule(start_rock_prop);
};
void propagate_rock_x(DGGML::Grammar<Microtubule::graph_type> &gamma,
           Microtubule::graph_type &system_graph,
           Parameters &settings) {

GT propagate_rock_x_lhs;
propagate_rock_x_lhs.addNode({1, {Microtubule::RockStart{} }});

GT propagate_rock_x_rhs;
propagate_rock_x_rhs.addNode({1, {Microtubule::RockStart{} }});

propagate_rock_x_rhs.addNode({2, {Microtubule::RockStart{} }});

propagate_rock_x_rhs.addEdge(1, 2);

DGGML::WithRule<GT> propagate_rock_x("propagate_rock_x", propagate_rock_x_lhs, propagate_rock_x_rhs,
[&](auto &lhs, auto &m) {

double cx0 = std::get<Microtubule::RockStart>(lhs[m[1]].data).fflow_a13208;

return 
(
1
 * 
( (cx0 < 300) ? 1.0 : 0.0)
)
;},
[&](auto &lhs, auto &rhs, auto &m1, auto &m2) {
double offset_distance = (0.1);
double x0 = std::get<Microtubule::RockStart>(lhs[m1[1]].data).fflow_a9c368[0];

double x0n = (x0);
		rhs[m2[1]].position[0] = x0n;
		std::get<Microtubule::RockStart>(rhs[m2[1]].data).fflow_a9c368[0] = x0n;
double y0 = std::get<Microtubule::RockStart>(lhs[m1[1]].data).fflow_a9c368[1];

double y0n = (y0);
		rhs[m2[1]].position[1] = y0n;
		std::get<Microtubule::RockStart>(rhs[m2[1]].data).fflow_a9c368[1] = y0n;
double z0 = std::get<Microtubule::RockStart>(lhs[m1[1]].data).fflow_a9c368[2];

double z0n = (z0);
		rhs[m2[1]].position[2] = z0n;
		std::get<Microtubule::RockStart>(rhs[m2[1]].data).fflow_a9c368[2] = z0n;
double cx0n = (0);
		std::get<Microtubule::RockStart>(rhs[m2[1]].data).fflow_a13208 = cx0n;
double cy0 = std::get<Microtubule::RockStart>(lhs[m1[1]].data).fflow_992eaa;

double cy0n = (cy0);
		std::get<Microtubule::RockStart>(rhs[m2[1]].data).fflow_992eaa = cy0n;
double cz0 = std::get<Microtubule::RockStart>(lhs[m1[1]].data).fflow_d20285;

double cz0n = (cz0);
		std::get<Microtubule::RockStart>(rhs[m2[1]].data).fflow_d20285 = cz0n;
double x1n = (x0 + offset_distance);
		rhs[m2[2]].position[0] = x1n;
		std::get<Microtubule::RockStart>(rhs[m2[2]].data).fflow_a9c368[0] = x1n;
double y1n = (y0);
		rhs[m2[2]].position[1] = y1n;
		std::get<Microtubule::RockStart>(rhs[m2[2]].data).fflow_a9c368[1] = y1n;
double z1n = (z0);
		rhs[m2[2]].position[2] = z1n;
		std::get<Microtubule::RockStart>(rhs[m2[2]].data).fflow_a9c368[2] = z1n;
double cx0 = std::get<Microtubule::RockStart>(lhs[m1[1]].data).fflow_a13208;

double cx1n = (cx0 - 1);
		std::get<Microtubule::RockStart>(rhs[m2[2]].data).fflow_a13208 = cx1n;
double cy1n = (cy0);
		std::get<Microtubule::RockStart>(rhs[m2[2]].data).fflow_992eaa = cy1n;
double cz1n = (cz0);
		std::get<Microtubule::RockStart>(rhs[m2[2]].data).fflow_d20285 = cz1n;
}
);
gamma.addRule(propagate_rock_x);
};
void propagate_rock_y(DGGML::Grammar<Microtubule::graph_type> &gamma,
           Microtubule::graph_type &system_graph,
           Parameters &settings) {

GT propagate_rock_y_lhs;
propagate_rock_y_lhs.addNode({1, {Microtubule::RockStart{} }});

GT propagate_rock_y_rhs;
propagate_rock_y_rhs.addNode({1, {Microtubule::RockStart{} }});

propagate_rock_y_rhs.addNode({2, {Microtubule::RockStart{} }});

propagate_rock_y_rhs.addEdge(1, 2);

DGGML::WithRule<GT> propagate_rock_y("propagate_rock_y", propagate_rock_y_lhs, propagate_rock_y_rhs,
[&](auto &lhs, auto &m) {

double cy0 = std::get<Microtubule::RockStart>(lhs[m[1]].data).fflow_992eaa;

return 
(
1
 * 
( (cy0 < 300) ? 1.0 : 0.0)
)
;},
[&](auto &lhs, auto &rhs, auto &m1, auto &m2) {
double offset_distance = (0.1);
double x0 = std::get<Microtubule::RockStart>(lhs[m1[1]].data).fflow_a9c368[0];

double x0n = (x0);
		rhs[m2[1]].position[0] = x0n;
		std::get<Microtubule::RockStart>(rhs[m2[1]].data).fflow_a9c368[0] = x0n;
double y0 = std::get<Microtubule::RockStart>(lhs[m1[1]].data).fflow_a9c368[1];

double y0n = (y0);
		rhs[m2[1]].position[1] = y0n;
		std::get<Microtubule::RockStart>(rhs[m2[1]].data).fflow_a9c368[1] = y0n;
double z0 = std::get<Microtubule::RockStart>(lhs[m1[1]].data).fflow_a9c368[2];

double z0n = (z0);
		rhs[m2[1]].position[2] = z0n;
		std::get<Microtubule::RockStart>(rhs[m2[1]].data).fflow_a9c368[2] = z0n;
double cx0 = std::get<Microtubule::RockStart>(lhs[m1[1]].data).fflow_a13208;

double cx0n = (cx0);
		std::get<Microtubule::RockStart>(rhs[m2[1]].data).fflow_a13208 = cx0n;
double cy0n = (0);
		std::get<Microtubule::RockStart>(rhs[m2[1]].data).fflow_992eaa = cy0n;
double cz0 = std::get<Microtubule::RockStart>(lhs[m1[1]].data).fflow_d20285;

double cz0n = (cz0);
		std::get<Microtubule::RockStart>(rhs[m2[1]].data).fflow_d20285 = cz0n;
double x1n = (x0);
		rhs[m2[2]].position[0] = x1n;
		std::get<Microtubule::RockStart>(rhs[m2[2]].data).fflow_a9c368[0] = x1n;
double y1n = (y0 + offset_distance);
		rhs[m2[2]].position[1] = y1n;
		std::get<Microtubule::RockStart>(rhs[m2[2]].data).fflow_a9c368[1] = y1n;
double z1n = (z0);
		rhs[m2[2]].position[2] = z1n;
		std::get<Microtubule::RockStart>(rhs[m2[2]].data).fflow_a9c368[2] = z1n;
double cx1n = (cx0);
		std::get<Microtubule::RockStart>(rhs[m2[2]].data).fflow_a13208 = cx1n;
double cy0 = std::get<Microtubule::RockStart>(lhs[m1[1]].data).fflow_992eaa;

double cy1n = (cy0 - 1);
		std::get<Microtubule::RockStart>(rhs[m2[2]].data).fflow_992eaa = cy1n;
double cz1n = (cz0);
		std::get<Microtubule::RockStart>(rhs[m2[2]].data).fflow_d20285 = cz1n;
}
);
gamma.addRule(propagate_rock_y);
};
void propagate_rock_z(DGGML::Grammar<Microtubule::graph_type> &gamma,
           Microtubule::graph_type &system_graph,
           Parameters &settings) {

GT propagate_rock_z_lhs;
propagate_rock_z_lhs.addNode({1, {Microtubule::RockStart{} }});

GT propagate_rock_z_rhs;
propagate_rock_z_rhs.addNode({1, {Microtubule::RockStart{} }});

propagate_rock_z_rhs.addNode({2, {Microtubule::RockStart{} }});

propagate_rock_z_rhs.addEdge(1, 2);

DGGML::WithRule<GT> propagate_rock_z("propagate_rock_z", propagate_rock_z_lhs, propagate_rock_z_rhs,
[&](auto &lhs, auto &m) {

double cy0 = std::get<Microtubule::RockStart>(lhs[m[1]].data).fflow_992eaa;

return 
(
1
 * 
( (cy0 < 300) ? 1.0 : 0.0)
)
;},
[&](auto &lhs, auto &rhs, auto &m1, auto &m2) {
double offset_distance = (0.1);
double x0 = std::get<Microtubule::RockStart>(lhs[m1[1]].data).fflow_a9c368[0];

double x0n = (x0);
		rhs[m2[1]].position[0] = x0n;
		std::get<Microtubule::RockStart>(rhs[m2[1]].data).fflow_a9c368[0] = x0n;
double y0 = std::get<Microtubule::RockStart>(lhs[m1[1]].data).fflow_a9c368[1];

double y0n = (y0);
		rhs[m2[1]].position[1] = y0n;
		std::get<Microtubule::RockStart>(rhs[m2[1]].data).fflow_a9c368[1] = y0n;
double z0 = std::get<Microtubule::RockStart>(lhs[m1[1]].data).fflow_a9c368[2];

double z0n = (z0);
		rhs[m2[1]].position[2] = z0n;
		std::get<Microtubule::RockStart>(rhs[m2[1]].data).fflow_a9c368[2] = z0n;
double cx0 = std::get<Microtubule::RockStart>(lhs[m1[1]].data).fflow_a13208;

double cx0n = (cx0);
		std::get<Microtubule::RockStart>(rhs[m2[1]].data).fflow_a13208 = cx0n;
double cy0 = std::get<Microtubule::RockStart>(lhs[m1[1]].data).fflow_992eaa;

double cy0n = (cy0);
		std::get<Microtubule::RockStart>(rhs[m2[1]].data).fflow_992eaa = cy0n;
double cz0n = (0);
		std::get<Microtubule::RockStart>(rhs[m2[1]].data).fflow_d20285 = cz0n;
double x1n = (x0);
		rhs[m2[2]].position[0] = x1n;
		std::get<Microtubule::RockStart>(rhs[m2[2]].data).fflow_a9c368[0] = x1n;
double y1n = (y0);
		rhs[m2[2]].position[1] = y1n;
		std::get<Microtubule::RockStart>(rhs[m2[2]].data).fflow_a9c368[1] = y1n;
double z1n = (z0 + offset_distance);
		rhs[m2[2]].position[2] = z1n;
		std::get<Microtubule::RockStart>(rhs[m2[2]].data).fflow_a9c368[2] = z1n;
double cx1n = (cx0);
		std::get<Microtubule::RockStart>(rhs[m2[2]].data).fflow_a13208 = cx1n;
double cy1n = (cy0);
		std::get<Microtubule::RockStart>(rhs[m2[2]].data).fflow_992eaa = cy1n;
double cz0 = std::get<Microtubule::RockStart>(lhs[m1[1]].data).fflow_d20285;

double cz1n = (cz0 - 1);
		std::get<Microtubule::RockStart>(rhs[m2[2]].data).fflow_d20285 = cz1n;
}
);
gamma.addRule(propagate_rock_z);
};
}
#endif