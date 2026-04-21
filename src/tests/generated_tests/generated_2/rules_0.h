#ifndef DGGML_RULES_STAGE_0_HPP
#define DGGML_RULES_STAGE_0_HPP
#include "types.h"
#include "parameters.h"
#include "functions.h"
namespace Dissolution {
using GT = Dissolution::graph_type;
void start_rock_prop(DGGML::Grammar<Dissolution::graph_type> &gamma,
           Dissolution::graph_type &system_graph,
           Parameters &settings) {

GT start_rock_prop_lhs;
start_rock_prop_lhs.addNode({1, {Dissolution::StartType{} }});

GT start_rock_prop_rhs;
start_rock_prop_rhs.addNode({1, {Dissolution::FluidSource{} }});

start_rock_prop_rhs.addNode({2, {Dissolution::RockStart{} }});

DGGML::WithRule<GT> start_rock_prop("start_rock_prop", start_rock_prop_lhs, start_rock_prop_rhs,
[&](auto &lhs, auto &m1) {

return 
(
10000.0
)
;},
[&](auto &lhs, auto &rhs, auto &m1, auto &m2) {
torch::Tensor fpos = (torch::tensor({1.35, 1.35, 0.0}, torch::kFloat64));
for (int _i = 0; _i < 3 && _i < fpos.numel(); _i++) { rhs[m2[1]].position[_i] = fpos[_i].template item<double>(); }
std::get<Dissolution::FluidSource>(rhs[m2[ 1 ]].data).Position = fpos;
torch::Tensor funit = (torch::tensor({0.0, 0.0, 1.0}, torch::kFloat64));
std::get<Dissolution::FluidSource>(rhs[m2[ 1 ]].data).Unit = funit;
torch::Tensor fc = (torch::tensor({10}, torch::kFloat64));
std::get<Dissolution::FluidSource>(rhs[m2[ 1 ]].data).FCount = fc;
double fP = (1.0);
std::get<Dissolution::FluidSource>(rhs[m2[ 1 ]].data).Pressure = fP;
double fC = (0.0);
std::get<Dissolution::FluidSource>(rhs[m2[ 1 ]].data).Concentration = fC;
torch::Tensor pos = (torch::tensor({0.0, 0.0, 0.0}, torch::kFloat64));
for (int _i = 0; _i < 3 && _i < pos.numel(); _i++) { rhs[m2[2]].position[_i] = pos[_i].template item<double>(); }
std::get<Dissolution::RockStart>(rhs[m2[ 2 ]].data).Position = pos;
torch::Tensor c_count = (torch::tensor({10, 10, 10}, torch::kFloat64));
std::get<Dissolution::RockStart>(rhs[m2[ 2 ]].data).Count = c_count;
torch::Tensor dir = (torch::tensor({0, 0, 0}, torch::kFloat64));
std::get<Dissolution::RockStart>(rhs[m2[ 2 ]].data).Dir = dir;
double rho = (1.0);
std::get<Dissolution::RockStart>(rhs[m2[ 2 ]].data).Density = rho;
}
);
gamma.addRule(start_rock_prop);
};
void propagate_rock_x(DGGML::Grammar<Dissolution::graph_type> &gamma,
           Dissolution::graph_type &system_graph,
           Parameters &settings) {

GT propagate_rock_x_lhs;
propagate_rock_x_lhs.addNode({1, {Dissolution::RockStart{} }});

GT propagate_rock_x_rhs;
propagate_rock_x_rhs.addNode({1, {Dissolution::RockStart{} }});

propagate_rock_x_rhs.addNode({2, {Dissolution::RockStart{} }});

DGGML::WithRule<GT> propagate_rock_x("propagate_rock_x", propagate_rock_x_lhs, propagate_rock_x_rhs,
[&](auto &lhs, auto &m1) {

torch::Tensor c_count0 = std::get<Dissolution::RockStart>(lhs[m1[1]].data).Count;

return 
(
10000.0
 * 
( (c_count0[0].template item<double>() > 0) ? 1.0 : 0.0)
)
;},
[&](auto &lhs, auto &rhs, auto &m1, auto &m2) {
double offset_distance = (0.3);
torch::Tensor pos0 = std::get<Dissolution::RockStart>(lhs[m1[1]].data).Position;

torch::Tensor pos0n = (pos0);
for (int _i = 0; _i < 3 && _i < pos0n.numel(); _i++) { rhs[m2[1]].position[_i] = pos0n[_i].template item<double>(); }
std::get<Dissolution::RockStart>(rhs[m2[ 1 ]].data).Position = pos0n;
double _arr_tmp_0 = 0;
torch::Tensor c_count0 = std::get<Dissolution::RockStart>(lhs[m1[1]].data).Count;

double _arr_tmp_1 = c_count0[1].template item<double>();
double _arr_tmp_2 = c_count0[2].template item<double>();
torch::Tensor c_count0n = (torch::tensor({_arr_tmp_0, _arr_tmp_1, _arr_tmp_2}, torch::kFloat64));
std::get<Dissolution::RockStart>(rhs[m2[ 1 ]].data).Count = c_count0n;
torch::Tensor dir0n = (torch::tensor({1, 0, 0}, torch::kFloat64));
std::get<Dissolution::RockStart>(rhs[m2[ 1 ]].data).Dir = dir0n;
double rho0 = std::get<Dissolution::RockStart>(lhs[m1[1]].data).Density;

double rho0n = (rho0);
std::get<Dissolution::RockStart>(rhs[m2[ 1 ]].data).Density = rho0n;
double _arr_tmp_3 = pos0[0].template item<double>() + offset_distance;
double _arr_tmp_4 = pos0[1].template item<double>();
double _arr_tmp_5 = pos0[2].template item<double>();
torch::Tensor pos1n = (torch::tensor({_arr_tmp_3, _arr_tmp_4, _arr_tmp_5}, torch::kFloat64));
for (int _i = 0; _i < 3 && _i < pos1n.numel(); _i++) { rhs[m2[2]].position[_i] = pos1n[_i].template item<double>(); }
std::get<Dissolution::RockStart>(rhs[m2[ 2 ]].data).Position = pos1n;
double _arr_tmp_6 = c_count0[0].template item<double>() - 1;
double _arr_tmp_7 = c_count0[1].template item<double>();
double _arr_tmp_8 = c_count0[2].template item<double>();
torch::Tensor c_count1n = (torch::tensor({_arr_tmp_6, _arr_tmp_7, _arr_tmp_8}, torch::kFloat64));
std::get<Dissolution::RockStart>(rhs[m2[ 2 ]].data).Count = c_count1n;
torch::Tensor dir1n = (torch::tensor({1, 0, 0}, torch::kFloat64));
std::get<Dissolution::RockStart>(rhs[m2[ 2 ]].data).Dir = dir1n;
double rho1n = (1.0);
std::get<Dissolution::RockStart>(rhs[m2[ 2 ]].data).Density = rho1n;
}
);
gamma.addRule(propagate_rock_x);
};
void propagate_rock_y(DGGML::Grammar<Dissolution::graph_type> &gamma,
           Dissolution::graph_type &system_graph,
           Parameters &settings) {

GT propagate_rock_y_lhs;
propagate_rock_y_lhs.addNode({1, {Dissolution::RockStart{} }});

GT propagate_rock_y_rhs;
propagate_rock_y_rhs.addNode({1, {Dissolution::RockStart{} }});

propagate_rock_y_rhs.addNode({2, {Dissolution::RockStart{} }});

DGGML::WithRule<GT> propagate_rock_y("propagate_rock_y", propagate_rock_y_lhs, propagate_rock_y_rhs,
[&](auto &lhs, auto &m1) {

torch::Tensor c_count0 = std::get<Dissolution::RockStart>(lhs[m1[1]].data).Count;

return 
(
10000.0
 * 
( (c_count0[0].template item<double>() == 0) ? 1.0 : 0.0)
 * 
( (c_count0[1].template item<double>() > 0) ? 1.0 : 0.0)
)
;},
[&](auto &lhs, auto &rhs, auto &m1, auto &m2) {
double offset_distance = (0.3);
torch::Tensor pos0 = std::get<Dissolution::RockStart>(lhs[m1[1]].data).Position;

torch::Tensor pos0n = (pos0);
for (int _i = 0; _i < 3 && _i < pos0n.numel(); _i++) { rhs[m2[1]].position[_i] = pos0n[_i].template item<double>(); }
std::get<Dissolution::RockStart>(rhs[m2[ 1 ]].data).Position = pos0n;
torch::Tensor c_count0 = std::get<Dissolution::RockStart>(lhs[m1[1]].data).Count;

double _arr_tmp_9 = c_count0[0].template item<double>();
double _arr_tmp_10 = 0;
double _arr_tmp_11 = c_count0[2].template item<double>();
torch::Tensor c_count0n = (torch::tensor({_arr_tmp_9, _arr_tmp_10, _arr_tmp_11}, torch::kFloat64));
std::get<Dissolution::RockStart>(rhs[m2[ 1 ]].data).Count = c_count0n;
torch::Tensor dir0n = (torch::tensor({0, 1, 0}, torch::kFloat64));
std::get<Dissolution::RockStart>(rhs[m2[ 1 ]].data).Dir = dir0n;
double rho0 = std::get<Dissolution::RockStart>(lhs[m1[1]].data).Density;

double rho0n = (rho0);
std::get<Dissolution::RockStart>(rhs[m2[ 1 ]].data).Density = rho0n;
double _arr_tmp_12 = pos0[0].template item<double>();
double _arr_tmp_13 = pos0[1].template item<double>() + offset_distance;
double _arr_tmp_14 = pos0[2].template item<double>();
torch::Tensor pos1n = (torch::tensor({_arr_tmp_12, _arr_tmp_13, _arr_tmp_14}, torch::kFloat64));
for (int _i = 0; _i < 3 && _i < pos1n.numel(); _i++) { rhs[m2[2]].position[_i] = pos1n[_i].template item<double>(); }
std::get<Dissolution::RockStart>(rhs[m2[ 2 ]].data).Position = pos1n;
double _arr_tmp_15 = c_count0[0].template item<double>();
double _arr_tmp_16 = c_count0[1].template item<double>() - 1;
double _arr_tmp_17 = c_count0[2].template item<double>();
torch::Tensor c_count1n = (torch::tensor({_arr_tmp_15, _arr_tmp_16, _arr_tmp_17}, torch::kFloat64));
std::get<Dissolution::RockStart>(rhs[m2[ 2 ]].data).Count = c_count1n;
torch::Tensor dir1n = (torch::tensor({0, 1, 0}, torch::kFloat64));
std::get<Dissolution::RockStart>(rhs[m2[ 2 ]].data).Dir = dir1n;
double rho1n = (1.0);
std::get<Dissolution::RockStart>(rhs[m2[ 2 ]].data).Density = rho1n;
}
);
gamma.addRule(propagate_rock_y);
};
void propagate_rock_z(DGGML::Grammar<Dissolution::graph_type> &gamma,
           Dissolution::graph_type &system_graph,
           Parameters &settings) {

GT propagate_rock_z_lhs;
propagate_rock_z_lhs.addNode({1, {Dissolution::RockStart{} }});

GT propagate_rock_z_rhs;
propagate_rock_z_rhs.addNode({1, {Dissolution::RockStart{} }});

propagate_rock_z_rhs.addNode({2, {Dissolution::RockStart{} }});

DGGML::WithRule<GT> propagate_rock_z("propagate_rock_z", propagate_rock_z_lhs, propagate_rock_z_rhs,
[&](auto &lhs, auto &m1) {

torch::Tensor c_count0 = std::get<Dissolution::RockStart>(lhs[m1[1]].data).Count;

return 
(
10000.0
 * 
( (c_count0[0].template item<double>() == 0) ? 1.0 : 0.0)
 * 
( (c_count0[1].template item<double>() == 0) ? 1.0 : 0.0)
 * 
( (c_count0[2].template item<double>() > 0) ? 1.0 : 0.0)
)
;},
[&](auto &lhs, auto &rhs, auto &m1, auto &m2) {
double offset_distance = (0.3);
torch::Tensor pos0 = std::get<Dissolution::RockStart>(lhs[m1[1]].data).Position;

torch::Tensor pos0n = (pos0);
for (int _i = 0; _i < 3 && _i < pos0n.numel(); _i++) { rhs[m2[1]].position[_i] = pos0n[_i].template item<double>(); }
std::get<Dissolution::RockStart>(rhs[m2[ 1 ]].data).Position = pos0n;
torch::Tensor c_count0 = std::get<Dissolution::RockStart>(lhs[m1[1]].data).Count;

double _arr_tmp_18 = c_count0[0].template item<double>();
double _arr_tmp_19 = c_count0[1].template item<double>();
double _arr_tmp_20 = 0;
torch::Tensor c_count0n = (torch::tensor({_arr_tmp_18, _arr_tmp_19, _arr_tmp_20}, torch::kFloat64));
std::get<Dissolution::RockStart>(rhs[m2[ 1 ]].data).Count = c_count0n;
torch::Tensor dir0n = (torch::tensor({0, 0, 1}, torch::kFloat64));
std::get<Dissolution::RockStart>(rhs[m2[ 1 ]].data).Dir = dir0n;
double rho0 = std::get<Dissolution::RockStart>(lhs[m1[1]].data).Density;

double rho0n = (rho0);
std::get<Dissolution::RockStart>(rhs[m2[ 1 ]].data).Density = rho0n;
double _arr_tmp_21 = pos0[0].template item<double>();
double _arr_tmp_22 = pos0[1].template item<double>();
double _arr_tmp_23 = pos0[2].template item<double>() + offset_distance;
torch::Tensor pos1n = (torch::tensor({_arr_tmp_21, _arr_tmp_22, _arr_tmp_23}, torch::kFloat64));
for (int _i = 0; _i < 3 && _i < pos1n.numel(); _i++) { rhs[m2[2]].position[_i] = pos1n[_i].template item<double>(); }
std::get<Dissolution::RockStart>(rhs[m2[ 2 ]].data).Position = pos1n;
double _arr_tmp_24 = c_count0[0].template item<double>();
double _arr_tmp_25 = c_count0[1].template item<double>();
double _arr_tmp_26 = c_count0[2].template item<double>() - 1;
torch::Tensor c_count1n = (torch::tensor({_arr_tmp_24, _arr_tmp_25, _arr_tmp_26}, torch::kFloat64));
std::get<Dissolution::RockStart>(rhs[m2[ 2 ]].data).Count = c_count1n;
torch::Tensor dir1n = (torch::tensor({0, 0, 1}, torch::kFloat64));
std::get<Dissolution::RockStart>(rhs[m2[ 2 ]].data).Dir = dir1n;
double rho1n = (1.0);
std::get<Dissolution::RockStart>(rhs[m2[ 2 ]].data).Density = rho1n;
}
);
gamma.addRule(propagate_rock_z);
};
void propagate_from_source(DGGML::Grammar<Dissolution::graph_type> &gamma,
           Dissolution::graph_type &system_graph,
           Parameters &settings) {

GT propagate_from_source_lhs;
propagate_from_source_lhs.addNode({1, {Dissolution::FluidSource{} }});

GT propagate_from_source_rhs;
propagate_from_source_rhs.addNode({1, {Dissolution::FluidSource{} }});

propagate_from_source_rhs.addNode({2, {Dissolution::Fluid{} }});

propagate_from_source_rhs.addEdge(1, 2);

DGGML::WithRule<GT> propagate_from_source("propagate_from_source", propagate_from_source_lhs, propagate_from_source_rhs,
[&](auto &lhs, auto &m1) {

torch::Tensor fc0 = std::get<Dissolution::FluidSource>(lhs[m1[1]].data).FCount;

return 
(
10000.0
 * 
( (fc0[0].template item<double>() > 1) ? 1.0 : 0.0)
)
;},
[&](auto &lhs, auto &rhs, auto &m1, auto &m2) {
double fluid_spacing = (0.3);
torch::Tensor fp0 = std::get<Dissolution::FluidSource>(lhs[m1[1]].data).Position;

torch::Tensor fp0n = (fp0);
for (int _i = 0; _i < 3 && _i < fp0n.numel(); _i++) { rhs[m2[1]].position[_i] = fp0n[_i].template item<double>(); }
std::get<Dissolution::FluidSource>(rhs[m2[ 1 ]].data).Position = fp0n;
torch::Tensor fu0 = std::get<Dissolution::FluidSource>(lhs[m1[1]].data).Unit;

torch::Tensor fu0n = (fu0);
std::get<Dissolution::FluidSource>(rhs[m2[ 1 ]].data).Unit = fu0n;
torch::Tensor fc0n = (torch::tensor({0}, torch::kFloat64));
std::get<Dissolution::FluidSource>(rhs[m2[ 1 ]].data).FCount = fc0n;
double fP0 = std::get<Dissolution::FluidSource>(lhs[m1[1]].data).Pressure;

double fP0n = (fP0);
std::get<Dissolution::FluidSource>(rhs[m2[ 1 ]].data).Pressure = fP0n;
double fC0 = std::get<Dissolution::FluidSource>(lhs[m1[1]].data).Concentration;

double fC0n = (fC0);
std::get<Dissolution::FluidSource>(rhs[m2[ 1 ]].data).Concentration = fC0n;
double _arr_tmp_27 = fp0[0].template item<double>() + fu0[0].template item<double>() * fluid_spacing;
double _arr_tmp_28 = fp0[1].template item<double>() + fu0[1].template item<double>() * fluid_spacing;
double _arr_tmp_29 = fp0[2].template item<double>() + fu0[2].template item<double>() * fluid_spacing;
torch::Tensor fp1n = (torch::tensor({_arr_tmp_27, _arr_tmp_28, _arr_tmp_29}, torch::kFloat64));
for (int _i = 0; _i < 3 && _i < fp1n.numel(); _i++) { rhs[m2[2]].position[_i] = fp1n[_i].template item<double>(); }
std::get<Dissolution::Fluid>(rhs[m2[ 2 ]].data).Position = fp1n;
torch::Tensor fu1n = (fu0);
std::get<Dissolution::Fluid>(rhs[m2[ 2 ]].data).Unit = fu1n;
torch::Tensor fc0 = std::get<Dissolution::FluidSource>(lhs[m1[1]].data).FCount;

double _arr_tmp_30 = fc0[0].template item<double>() - 1;
torch::Tensor fc1n = (torch::tensor({_arr_tmp_30}, torch::kFloat64));
std::get<Dissolution::Fluid>(rhs[m2[ 2 ]].data).FCount = fc1n;
double fP1n = (0.0);
std::get<Dissolution::Fluid>(rhs[m2[ 2 ]].data).Pressure = fP1n;
double fC1n = (0.0);
std::get<Dissolution::Fluid>(rhs[m2[ 2 ]].data).Concentration = fC1n;
}
);
gamma.addRule(propagate_from_source);
};
void propagate_source_to_sink(DGGML::Grammar<Dissolution::graph_type> &gamma,
           Dissolution::graph_type &system_graph,
           Parameters &settings) {

GT propagate_source_to_sink_lhs;
propagate_source_to_sink_lhs.addNode({1, {Dissolution::FluidSource{} }});

GT propagate_source_to_sink_rhs;
propagate_source_to_sink_rhs.addNode({1, {Dissolution::FluidSource{} }});

propagate_source_to_sink_rhs.addNode({2, {Dissolution::FluidSink{} }});

propagate_source_to_sink_rhs.addEdge(1, 2);

DGGML::WithRule<GT> propagate_source_to_sink("propagate_source_to_sink", propagate_source_to_sink_lhs, propagate_source_to_sink_rhs,
[&](auto &lhs, auto &m1) {

torch::Tensor fc0 = std::get<Dissolution::FluidSource>(lhs[m1[1]].data).FCount;

return 
(
10000.0
 * 
( (fc0[0].template item<double>() == 1) ? 1.0 : 0.0)
)
;},
[&](auto &lhs, auto &rhs, auto &m1, auto &m2) {
double fluid_spacing = (0.3);
torch::Tensor fp0 = std::get<Dissolution::FluidSource>(lhs[m1[1]].data).Position;

torch::Tensor fp0n = (fp0);
for (int _i = 0; _i < 3 && _i < fp0n.numel(); _i++) { rhs[m2[1]].position[_i] = fp0n[_i].template item<double>(); }
std::get<Dissolution::FluidSource>(rhs[m2[ 1 ]].data).Position = fp0n;
torch::Tensor fu0 = std::get<Dissolution::FluidSource>(lhs[m1[1]].data).Unit;

torch::Tensor fu0n = (fu0);
std::get<Dissolution::FluidSource>(rhs[m2[ 1 ]].data).Unit = fu0n;
torch::Tensor fc0n = (torch::tensor({0}, torch::kFloat64));
std::get<Dissolution::FluidSource>(rhs[m2[ 1 ]].data).FCount = fc0n;
double fP0 = std::get<Dissolution::FluidSource>(lhs[m1[1]].data).Pressure;

double fP0n = (fP0);
std::get<Dissolution::FluidSource>(rhs[m2[ 1 ]].data).Pressure = fP0n;
double fC0 = std::get<Dissolution::FluidSource>(lhs[m1[1]].data).Concentration;

double fC0n = (fC0);
std::get<Dissolution::FluidSource>(rhs[m2[ 1 ]].data).Concentration = fC0n;
double _arr_tmp_31 = fp0[0].template item<double>() + fu0[0].template item<double>() * fluid_spacing;
double _arr_tmp_32 = fp0[1].template item<double>() + fu0[1].template item<double>() * fluid_spacing;
double _arr_tmp_33 = fp0[2].template item<double>() + fu0[2].template item<double>() * fluid_spacing;
torch::Tensor fp1n = (torch::tensor({_arr_tmp_31, _arr_tmp_32, _arr_tmp_33}, torch::kFloat64));
for (int _i = 0; _i < 3 && _i < fp1n.numel(); _i++) { rhs[m2[2]].position[_i] = fp1n[_i].template item<double>(); }
std::get<Dissolution::FluidSink>(rhs[m2[ 2 ]].data).Position = fp1n;
torch::Tensor fu1n = (fu0);
std::get<Dissolution::FluidSink>(rhs[m2[ 2 ]].data).Unit = fu1n;
torch::Tensor fc1n = (torch::tensor({0}, torch::kFloat64));
std::get<Dissolution::FluidSink>(rhs[m2[ 2 ]].data).FCount = fc1n;
double fP1n = (0.0);
std::get<Dissolution::FluidSink>(rhs[m2[ 2 ]].data).Pressure = fP1n;
double fC1n = (0.0);
std::get<Dissolution::FluidSink>(rhs[m2[ 2 ]].data).Concentration = fC1n;
}
);
gamma.addRule(propagate_source_to_sink);
};
void propagate_initial_fluid(DGGML::Grammar<Dissolution::graph_type> &gamma,
           Dissolution::graph_type &system_graph,
           Parameters &settings) {

GT propagate_initial_fluid_lhs;
propagate_initial_fluid_lhs.addNode({1, {Dissolution::Fluid{} }});

GT propagate_initial_fluid_rhs;
propagate_initial_fluid_rhs.addNode({1, {Dissolution::Fluid{} }});

propagate_initial_fluid_rhs.addNode({2, {Dissolution::Fluid{} }});

propagate_initial_fluid_rhs.addEdge(1, 2);

DGGML::WithRule<GT> propagate_initial_fluid("propagate_initial_fluid", propagate_initial_fluid_lhs, propagate_initial_fluid_rhs,
[&](auto &lhs, auto &m1) {

torch::Tensor fc0 = std::get<Dissolution::Fluid>(lhs[m1[1]].data).FCount;

return 
(
10000.0
 * 
( (fc0[0].template item<double>() > 1) ? 1.0 : 0.0)
)
;},
[&](auto &lhs, auto &rhs, auto &m1, auto &m2) {
double fluid_spacing = (0.3);
torch::Tensor fp0 = std::get<Dissolution::Fluid>(lhs[m1[1]].data).Position;

torch::Tensor fp0n = (fp0);
for (int _i = 0; _i < 3 && _i < fp0n.numel(); _i++) { rhs[m2[1]].position[_i] = fp0n[_i].template item<double>(); }
std::get<Dissolution::Fluid>(rhs[m2[ 1 ]].data).Position = fp0n;
torch::Tensor fu0 = std::get<Dissolution::Fluid>(lhs[m1[1]].data).Unit;

torch::Tensor fu0n = (fu0);
std::get<Dissolution::Fluid>(rhs[m2[ 1 ]].data).Unit = fu0n;
torch::Tensor fc0n = (torch::tensor({0}, torch::kFloat64));
std::get<Dissolution::Fluid>(rhs[m2[ 1 ]].data).FCount = fc0n;
double fP0 = std::get<Dissolution::Fluid>(lhs[m1[1]].data).Pressure;

double fP0n = (fP0);
std::get<Dissolution::Fluid>(rhs[m2[ 1 ]].data).Pressure = fP0n;
double fC0 = std::get<Dissolution::Fluid>(lhs[m1[1]].data).Concentration;

double fC0n = (fC0);
std::get<Dissolution::Fluid>(rhs[m2[ 1 ]].data).Concentration = fC0n;
double _arr_tmp_34 = fp0[0].template item<double>() + fu0[0].template item<double>() * fluid_spacing;
double _arr_tmp_35 = fp0[1].template item<double>() + fu0[1].template item<double>() * fluid_spacing;
double _arr_tmp_36 = fp0[2].template item<double>() + fu0[2].template item<double>() * fluid_spacing;
torch::Tensor fp1n = (torch::tensor({_arr_tmp_34, _arr_tmp_35, _arr_tmp_36}, torch::kFloat64));
for (int _i = 0; _i < 3 && _i < fp1n.numel(); _i++) { rhs[m2[2]].position[_i] = fp1n[_i].template item<double>(); }
std::get<Dissolution::Fluid>(rhs[m2[ 2 ]].data).Position = fp1n;
torch::Tensor fu1n = (fu0);
std::get<Dissolution::Fluid>(rhs[m2[ 2 ]].data).Unit = fu1n;
torch::Tensor fc0 = std::get<Dissolution::Fluid>(lhs[m1[1]].data).FCount;

double _arr_tmp_37 = fc0[0].template item<double>() - 1;
torch::Tensor fc1n = (torch::tensor({_arr_tmp_37}, torch::kFloat64));
std::get<Dissolution::Fluid>(rhs[m2[ 2 ]].data).FCount = fc1n;
double fP1n = (0.0);
std::get<Dissolution::Fluid>(rhs[m2[ 2 ]].data).Pressure = fP1n;
double fC1n = (0.0);
std::get<Dissolution::Fluid>(rhs[m2[ 2 ]].data).Concentration = fC1n;
}
);
gamma.addRule(propagate_initial_fluid);
};
void propagate_final_fluid(DGGML::Grammar<Dissolution::graph_type> &gamma,
           Dissolution::graph_type &system_graph,
           Parameters &settings) {

GT propagate_final_fluid_lhs;
propagate_final_fluid_lhs.addNode({1, {Dissolution::Fluid{} }});

GT propagate_final_fluid_rhs;
propagate_final_fluid_rhs.addNode({1, {Dissolution::Fluid{} }});

propagate_final_fluid_rhs.addNode({2, {Dissolution::FluidSink{} }});

propagate_final_fluid_rhs.addEdge(1, 2);

DGGML::WithRule<GT> propagate_final_fluid("propagate_final_fluid", propagate_final_fluid_lhs, propagate_final_fluid_rhs,
[&](auto &lhs, auto &m1) {

torch::Tensor fc0 = std::get<Dissolution::Fluid>(lhs[m1[1]].data).FCount;

return 
(
10000.0
 * 
( (fc0[0].template item<double>() == 1) ? 1.0 : 0.0)
)
;},
[&](auto &lhs, auto &rhs, auto &m1, auto &m2) {
double fluid_spacing = (0.3);
torch::Tensor fp0 = std::get<Dissolution::Fluid>(lhs[m1[1]].data).Position;

torch::Tensor fp0n = (fp0);
for (int _i = 0; _i < 3 && _i < fp0n.numel(); _i++) { rhs[m2[1]].position[_i] = fp0n[_i].template item<double>(); }
std::get<Dissolution::Fluid>(rhs[m2[ 1 ]].data).Position = fp0n;
torch::Tensor fu0 = std::get<Dissolution::Fluid>(lhs[m1[1]].data).Unit;

torch::Tensor fu0n = (fu0);
std::get<Dissolution::Fluid>(rhs[m2[ 1 ]].data).Unit = fu0n;
torch::Tensor fc0n = (torch::tensor({0}, torch::kFloat64));
std::get<Dissolution::Fluid>(rhs[m2[ 1 ]].data).FCount = fc0n;
double fP0 = std::get<Dissolution::Fluid>(lhs[m1[1]].data).Pressure;

double fP0n = (fP0);
std::get<Dissolution::Fluid>(rhs[m2[ 1 ]].data).Pressure = fP0n;
double fC0 = std::get<Dissolution::Fluid>(lhs[m1[1]].data).Concentration;

double fC0n = (fC0);
std::get<Dissolution::Fluid>(rhs[m2[ 1 ]].data).Concentration = fC0n;
double _arr_tmp_38 = fp0[0].template item<double>() + fu0[0].template item<double>() * fluid_spacing;
double _arr_tmp_39 = fp0[1].template item<double>() + fu0[1].template item<double>() * fluid_spacing;
double _arr_tmp_40 = fp0[2].template item<double>() + fu0[2].template item<double>() * fluid_spacing;
torch::Tensor fp1n = (torch::tensor({_arr_tmp_38, _arr_tmp_39, _arr_tmp_40}, torch::kFloat64));
for (int _i = 0; _i < 3 && _i < fp1n.numel(); _i++) { rhs[m2[2]].position[_i] = fp1n[_i].template item<double>(); }
std::get<Dissolution::FluidSink>(rhs[m2[ 2 ]].data).Position = fp1n;
torch::Tensor fu1n = (fu0);
std::get<Dissolution::FluidSink>(rhs[m2[ 2 ]].data).Unit = fu1n;
torch::Tensor fc1n = (torch::tensor({0}, torch::kFloat64));
std::get<Dissolution::FluidSink>(rhs[m2[ 2 ]].data).FCount = fc1n;
double fP1n = (0.0);
std::get<Dissolution::FluidSink>(rhs[m2[ 2 ]].data).Pressure = fP1n;
double fC1n = (0.0);
std::get<Dissolution::FluidSink>(rhs[m2[ 2 ]].data).Concentration = fC1n;
}
);
gamma.addRule(propagate_final_fluid);
};
}
#endif