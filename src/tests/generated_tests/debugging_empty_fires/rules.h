#ifndef DGGML_RULES_HPP
#define DGGML_RULES_HPP
#include "types.h"
#include "parameters.h"
namespace Microtubule {
using GT = Microtubule::graph_type;
void start_node_dup(DGGML::Grammar<Microtubule::graph_type> &gamma,
           Microtubule::graph_type &system_graph,
           Parameters &settings) {

GT start_node_dup_lhs;
start_node_dup_lhs.addNode({1, {Microtubule::StartType{} }});

GT start_node_dup_rhs;
start_node_dup_rhs.addNode({1, {Microtubule::FractureSegment{} }});

start_node_dup_rhs.addNode({2, {Microtubule::FractureSegment{} }});

DGGML::WithRule<GT> start_node_dup("start_node_dup", start_node_dup_lhs, start_node_dup_rhs,
[&](auto &lhs, auto &m) {

return 
(
10
)
;},
[&](auto &lhs, auto &rhs, auto &m1, auto &m2) {
double start_xn = (3.0);

		rhs[m2[1]].position[0] = start_xn;
		std::get<Microtubule::FractureSegment>(rhs[m2[1]].data).fflow_bb7745[0] = start_xn;
        double start_yn = (3.0);
		rhs[m2[1]].position[1] = start_yn;

		std::get<Microtubule::FractureSegment>(rhs[m2[1]].data).fflow_bb7745[1] = start_yn;
        double pn = (1.0);
		std::get<Microtubule::FractureSegment>(rhs[m2[1]].data).fflow_87a940 = pn;
        double ux = (1.0);
		std::get<Microtubule::FractureSegment>(rhs[m2[1]].data).fflow_b62e5f[0] = ux;
        double uy = (1.0);
		std::get<Microtubule::FractureSegment>(rhs[m2[1]].data).fflow_b62e5f[1] = uy;
        double start_xn1 = (3);
		rhs[m2[2]].position[0] = start_xn1;
		std::get<Microtubule::FractureSegment>(rhs[m2[2]].data).fflow_bb7745[0] = start_xn1;

double start_yn1 = (7);
		rhs[m2[2]].position[1] = start_yn1;
		std::get<Microtubule::FractureSegment>(rhs[m2[2]].data).fflow_bb7745[1] = start_yn1;
double pn1 = (1.0);
		std::get<Microtubule::FractureSegment>(rhs[m2[2]].data).fflow_87a940 = pn1;
double ux1 = (1.0);
		std::get<Microtubule::FractureSegment>(rhs[m2[2]].data).fflow_b62e5f[0] = ux1;
double uy1 = (-1.0);
		std::get<Microtubule::FractureSegment>(rhs[m2[2]].data).fflow_b62e5f[1] = uy1;
}
);
gamma.addRule(start_node_dup);
};
void start_to_node_fracture(DGGML::Grammar<Microtubule::graph_type> &gamma,
           Microtubule::graph_type &system_graph,
           Parameters &settings) {

GT start_to_node_fracture_lhs;
start_to_node_fracture_lhs.addNode({1, {Microtubule::FractureSegment{} }});

GT start_to_node_fracture_rhs;
start_to_node_fracture_rhs.addNode({1, {Microtubule::FractureSegment{} }});

start_to_node_fracture_rhs.addNode({2, {Microtubule::FractureSegmentEnd{} }});

start_to_node_fracture_rhs.addEdge(1, 2);

DGGML::WithRule<GT> start_to_node_fracture("start_to_node_fracture", start_to_node_fracture_lhs, start_to_node_fracture_rhs,
[&](auto &lhs, auto &m) {

return 
(
DGGML::heaviside(10, 1)
)
;},
[&](auto &lhs, auto &rhs, auto &m1, auto &m2) {
double x = std::get<Microtubule::FractureSegment>(lhs[m1[1]].data).fflow_bb7745[0];

double xn0 = (x);
		rhs[m2[1]].position[0] = xn0;
		std::get<Microtubule::FractureSegment>(rhs[m2[1]].data).fflow_bb7745[0] = xn0;
double y = std::get<Microtubule::FractureSegment>(lhs[m1[1]].data).fflow_bb7745[1];

double yn0 = (y);
		rhs[m2[1]].position[1] = yn0;
		std::get<Microtubule::FractureSegment>(rhs[m2[1]].data).fflow_bb7745[1] = yn0;
double p = std::get<Microtubule::FractureSegment>(lhs[m1[1]].data).fflow_87a940;

double pn0 = (p);
		std::get<Microtubule::FractureSegment>(rhs[m2[1]].data).fflow_87a940 = pn0;
double ux = std::get<Microtubule::FractureSegment>(lhs[m1[1]].data).fflow_b62e5f[0];

double ux0 = (ux);
		std::get<Microtubule::FractureSegment>(rhs[m2[1]].data).fflow_b62e5f[0] = ux0;
double uy = std::get<Microtubule::FractureSegment>(lhs[m1[1]].data).fflow_b62e5f[1];

double uy0 = (uy);
		std::get<Microtubule::FractureSegment>(rhs[m2[1]].data).fflow_b62e5f[1] = uy0;
double x2 = ((xn0 + 1 * (settings.s_min + settings.s_min)) * ux0);
		rhs[m2[2]].position[0] = x2;
		std::get<Microtubule::FractureSegmentEnd>(rhs[m2[2]].data).fflow_f9bbc8[0] = x2;
double y2 = ((yn0 + 1 * (settings.s_min + settings.s_min)) * uy0);
		rhs[m2[2]].position[1] = y2;
		std::get<Microtubule::FractureSegmentEnd>(rhs[m2[2]].data).fflow_f9bbc8[1] = y2;
double p2 = (1.0);
		std::get<Microtubule::FractureSegmentEnd>(rhs[m2[2]].data).fflow_1e83fe = p2;
double ux2 = (ux0);
		std::get<Microtubule::FractureSegmentEnd>(rhs[m2[2]].data).fflow_f13ca0[0] = ux2;
double uy2 = (uy0);
		std::get<Microtubule::FractureSegmentEnd>(rhs[m2[2]].data).fflow_f13ca0[1] = uy2;
}
);
gamma.addRule(start_to_node_fracture);
};
void grow_fracture(DGGML::Grammar<Microtubule::graph_type> &gamma,
           Microtubule::graph_type &system_graph,
           Parameters &settings) {

GT grow_fracture_lhs;
grow_fracture_lhs.addNode({1, {Microtubule::FractureSegment{} }});

grow_fracture_lhs.addNode({2, {Microtubule::FractureSegmentEnd{} }});

grow_fracture_lhs.addEdge(1, 2);

GT grow_fracture_rhs;
grow_fracture_rhs.addNode({1, {Microtubule::FractureSegment{} }});

grow_fracture_rhs.addNode({2, {Microtubule::FractureSegment{} }});

grow_fracture_rhs.addEdge(1, 2);

grow_fracture_rhs.addNode({3, {Microtubule::FractureSegmentEnd{} }});

grow_fracture_rhs.addEdge(2, 3);

DGGML::WithRule<GT> grow_fracture("grow_fracture", grow_fracture_lhs, grow_fracture_rhs,
[&](auto &lhs, auto &m) {

double p2 = std::get<Microtubule::FractureSegmentEnd>(lhs[m[2]].data).fflow_1e83fe;

return 
(
1
 * 
DGGML::heaviside(p2, 1)
)
;},
[&](auto &lhs, auto &rhs, auto &m1, auto &m2) {
double x1 = std::get<Microtubule::FractureSegment>(lhs[m1[1]].data).fflow_bb7745[0];

double x1n = (x1);
		rhs[m2[1]].position[0] = x1n;
		std::get<Microtubule::FractureSegment>(rhs[m2[1]].data).fflow_bb7745[0] = x1n;
double y1 = std::get<Microtubule::FractureSegment>(lhs[m1[1]].data).fflow_bb7745[1];

double y1n = (y1);
		rhs[m2[1]].position[1] = y1n;
		std::get<Microtubule::FractureSegment>(rhs[m2[1]].data).fflow_bb7745[1] = y1n;
double p1 = std::get<Microtubule::FractureSegment>(lhs[m1[1]].data).fflow_87a940;

double p1n = (p1);
		std::get<Microtubule::FractureSegment>(rhs[m2[1]].data).fflow_87a940 = p1n;
double ux1 = std::get<Microtubule::FractureSegment>(lhs[m1[1]].data).fflow_b62e5f[0];

double ux1n = (ux1);
		std::get<Microtubule::FractureSegment>(rhs[m2[1]].data).fflow_b62e5f[0] = ux1n;
double uy1 = std::get<Microtubule::FractureSegment>(lhs[m1[1]].data).fflow_b62e5f[1];

double uy1n = (uy1);
		std::get<Microtubule::FractureSegment>(rhs[m2[1]].data).fflow_b62e5f[1] = uy1n;
double x2 = std::get<Microtubule::FractureSegmentEnd>(lhs[m1[2]].data).fflow_f9bbc8[0];

double x2n = (x2);
		rhs[m2[2]].position[0] = x2n;
		std::get<Microtubule::FractureSegment>(rhs[m2[2]].data).fflow_bb7745[0] = x2n;
double y2 = std::get<Microtubule::FractureSegmentEnd>(lhs[m1[2]].data).fflow_f9bbc8[1];

double y2n = (y2);
		rhs[m2[2]].position[1] = y2n;
		std::get<Microtubule::FractureSegment>(rhs[m2[2]].data).fflow_bb7745[1] = y2n;
double p2 = std::get<Microtubule::FractureSegmentEnd>(lhs[m1[2]].data).fflow_1e83fe;

double p2n = (p2);
		std::get<Microtubule::FractureSegment>(rhs[m2[2]].data).fflow_87a940 = p2n;
double ux2 = std::get<Microtubule::FractureSegmentEnd>(lhs[m1[2]].data).fflow_f13ca0[0];

double ux2n = (ux2);
		std::get<Microtubule::FractureSegment>(rhs[m2[2]].data).fflow_b62e5f[0] = ux2n;
double uy2 = std::get<Microtubule::FractureSegmentEnd>(lhs[m1[2]].data).fflow_f13ca0[1];

double uy2n = (uy2);
		std::get<Microtubule::FractureSegment>(rhs[m2[2]].data).fflow_b62e5f[1] = uy2n;
double x3n = ((x2 + settings.s_min) * ux2);
		rhs[m2[3]].position[0] = x3n;
		std::get<Microtubule::FractureSegmentEnd>(rhs[m2[3]].data).fflow_f9bbc8[0] = x3n;
double y3n = ((y2 + settings.s_min) * uy2);
		rhs[m2[3]].position[1] = y3n;
		std::get<Microtubule::FractureSegmentEnd>(rhs[m2[3]].data).fflow_f9bbc8[1] = y3n;
double p3n = (1.0);
		std::get<Microtubule::FractureSegmentEnd>(rhs[m2[3]].data).fflow_1e83fe = p3n;
double ux3 = (ux2);
		std::get<Microtubule::FractureSegmentEnd>(rhs[m2[3]].data).fflow_f13ca0[0] = ux3;
double uy3 = (uy2);
		std::get<Microtubule::FractureSegmentEnd>(rhs[m2[3]].data).fflow_f13ca0[1] = uy3;
}
);
gamma.addRule(grow_fracture);
};
void hit_boundary(DGGML::Grammar<Microtubule::graph_type> &gamma,
           Microtubule::graph_type &system_graph,
           Parameters &settings) {

GT hit_boundary_lhs;
hit_boundary_lhs.addNode({1, {Microtubule::FractureSegment{} }});

hit_boundary_lhs.addNode({2, {Microtubule::FractureSegmentEnd{} }});

hit_boundary_lhs.addEdge(1, 2);

hit_boundary_lhs.addNode({3, {Microtubule::Boundary{} }});

hit_boundary_lhs.addNode({4, {Microtubule::Boundary{} }});

hit_boundary_lhs.addEdge(3, 4);

GT hit_boundary_rhs;
hit_boundary_rhs.addNode({1, {Microtubule::FractureSegment{} }});

hit_boundary_rhs.addNode({2, {Microtubule::FractureSegmentEnd{} }});

hit_boundary_rhs.addEdge(1, 2);

DGGML::WithRule<GT> hit_boundary("hit_boundary", hit_boundary_lhs, hit_boundary_rhs,
[&](auto &lhs, auto &m) {

return 
(
0.0
)
;},
[&](auto &lhs, auto &rhs, auto &m1, auto &m2) {
double x0 = std::get<Microtubule::FractureSegment>(lhs[m1[1]].data).fflow_bb7745[0];

double x0n = (x0);
		rhs[m2[1]].position[0] = x0n;
		std::get<Microtubule::FractureSegment>(rhs[m2[1]].data).fflow_bb7745[0] = x0n;
double y0 = std::get<Microtubule::FractureSegment>(lhs[m1[1]].data).fflow_bb7745[1];

double y0n = (y0);
		rhs[m2[1]].position[1] = y0n;
		std::get<Microtubule::FractureSegment>(rhs[m2[1]].data).fflow_bb7745[1] = y0n;
double p0 = std::get<Microtubule::FractureSegment>(lhs[m1[1]].data).fflow_87a940;

double p0n = (p0);
		std::get<Microtubule::FractureSegment>(rhs[m2[1]].data).fflow_87a940 = p0n;
double x1 = std::get<Microtubule::FractureSegmentEnd>(lhs[m1[2]].data).fflow_f9bbc8[0];

double x1n = (x1);
		rhs[m2[2]].position[0] = x1n;
		std::get<Microtubule::FractureSegmentEnd>(rhs[m2[2]].data).fflow_f9bbc8[0] = x1n;
double y1 = std::get<Microtubule::FractureSegmentEnd>(lhs[m1[2]].data).fflow_f9bbc8[1];

double y1n = (y1);
		rhs[m2[2]].position[1] = y1n;
		std::get<Microtubule::FractureSegmentEnd>(rhs[m2[2]].data).fflow_f9bbc8[1] = y1n;
double p1n = (0.0);
		std::get<Microtubule::FractureSegmentEnd>(rhs[m2[2]].data).fflow_1e83fe = p1n;
}
);
gamma.addRule(hit_boundary);
};
}
#endif
