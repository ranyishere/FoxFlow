#ifndef DGGML_RULES_STAGE_1_HPP
#define DGGML_RULES_STAGE_1_HPP
#include "types.h"
#include "parameters.h"
#include "functions.h"
namespace Dissolution {
using GT = Dissolution::graph_type;
void source_flow(DGGML::Grammar<Dissolution::graph_type> &gamma,
           Dissolution::graph_type &system_graph,
           Parameters &settings) {

GT source_flow_lhs;
source_flow_lhs.addNode({1, {Dissolution::FluidSource{} }});

source_flow_lhs.addNode({2, {Dissolution::Fluid{} }});

source_flow_lhs.addEdge(1, 2);

GT source_flow_rhs;
source_flow_rhs.addNode({1, {Dissolution::FluidSource{} }});

source_flow_rhs.addNode({2, {Dissolution::Fluid{} }});

source_flow_rhs.addEdge(1, 2);

DGGML::SolvingRule<GT> source_flow("source_flow", source_flow_lhs, source_flow_lhs,
2,
[](auto &lhs, auto &m1, auto &varset) {
auto &node_2_4 = std::get<Dissolution::Fluid>(lhs[m1[2]].data).Pressure;
varset.insert(&node_2_4);
auto &node_1_4 = std::get<Dissolution::FluidSource>(lhs[m1[1]].data).Pressure;
varset.insert(&node_1_4);
},
[&](auto &lhs, auto &m1, auto y, auto ydot, auto &varmap) {
auto &ix_P1 = std::get<Dissolution::Fluid>(lhs[m1[2]].data).Pressure;
auto &ix_P0 = std::get<Dissolution::FluidSource>(lhs[m1[1]].data).Pressure;
NV_Ith_S(ydot, varmap[&ix_P1]) += 11.11 * (NV_Ith_S(y, varmap.at(&ix_P0)) - NV_Ith_S(y, varmap.at(&ix_P1)));
NV_Ith_S(ydot, varmap[&ix_P0]) +=  0.0 ;
// ── Symbolic ODE system ──
// d(dP1)/dt += 11.11 * (P0 - P1)
// d(dP0)/dt += 0.0
{ static bool _sym_dumped = false;
  if (!_sym_dumped && std::getenv("ODE_DUMP")) { _sym_dumped = true;
    std::cout << "  ── Symbolic ODE ──" << std::endl;
    std::cout << "    d(dP1)/dt += 11.11 * (P0 - P1)" << std::endl;
    std::cout << "    d(dP0)/dt += 0.0" << std::endl;
  }
}
}
);
gamma.addRule(source_flow);
};
void fluid_flow(DGGML::Grammar<Dissolution::graph_type> &gamma,
           Dissolution::graph_type &system_graph,
           Parameters &settings) {

GT fluid_flow_lhs;
fluid_flow_lhs.addNode({1, {Dissolution::Fluid{} }});

fluid_flow_lhs.addNode({2, {Dissolution::Fluid{} }});

fluid_flow_lhs.addEdge(1, 2);

GT fluid_flow_rhs;
fluid_flow_rhs.addNode({1, {Dissolution::Fluid{} }});

fluid_flow_rhs.addNode({2, {Dissolution::Fluid{} }});

fluid_flow_rhs.addEdge(1, 2);

DGGML::SolvingRule<GT> fluid_flow("fluid_flow", fluid_flow_lhs, fluid_flow_lhs,
2,
[](auto &lhs, auto &m1, auto &varset) {
auto &node_2_4 = std::get<Dissolution::Fluid>(lhs[m1[2]].data).Pressure;
varset.insert(&node_2_4);
auto &node_1_4 = std::get<Dissolution::Fluid>(lhs[m1[1]].data).Pressure;
varset.insert(&node_1_4);
},
[&](auto &lhs, auto &m1, auto y, auto ydot, auto &varmap) {
auto &ix_P1 = std::get<Dissolution::Fluid>(lhs[m1[2]].data).Pressure;
auto &ix_P0 = std::get<Dissolution::Fluid>(lhs[m1[1]].data).Pressure;
NV_Ith_S(ydot, varmap[&ix_P1]) += 11.11 * (NV_Ith_S(y, varmap.at(&ix_P0)) - NV_Ith_S(y, varmap.at(&ix_P1)));
NV_Ith_S(ydot, varmap[&ix_P0]) +=  0.0 ;
// ── Symbolic ODE system ──
// d(dP1)/dt += 11.11 * (P0 - P1)
// d(dP0)/dt += 0.0
{ static bool _sym_dumped = false;
  if (!_sym_dumped && std::getenv("ODE_DUMP")) { _sym_dumped = true;
    std::cout << "  ── Symbolic ODE ──" << std::endl;
    std::cout << "    d(dP1)/dt += 11.11 * (P0 - P1)" << std::endl;
    std::cout << "    d(dP0)/dt += 0.0" << std::endl;
  }
}
}
);
gamma.addRule(fluid_flow);
};
void sink_flow(DGGML::Grammar<Dissolution::graph_type> &gamma,
           Dissolution::graph_type &system_graph,
           Parameters &settings) {

GT sink_flow_lhs;
sink_flow_lhs.addNode({1, {Dissolution::Fluid{} }});

sink_flow_lhs.addNode({2, {Dissolution::FluidSink{} }});

sink_flow_lhs.addEdge(1, 2);

GT sink_flow_rhs;
sink_flow_rhs.addNode({1, {Dissolution::Fluid{} }});

sink_flow_rhs.addNode({2, {Dissolution::FluidSink{} }});

sink_flow_rhs.addEdge(1, 2);

DGGML::SolvingRule<GT> sink_flow("sink_flow", sink_flow_lhs, sink_flow_lhs,
2,
[](auto &lhs, auto &m1, auto &varset) {
auto &node_1_4 = std::get<Dissolution::Fluid>(lhs[m1[1]].data).Pressure;
varset.insert(&node_1_4);
auto &node_2_4 = std::get<Dissolution::FluidSink>(lhs[m1[2]].data).Pressure;
varset.insert(&node_2_4);
},
[&](auto &lhs, auto &m1, auto y, auto ydot, auto &varmap) {
auto &ix_P1 = std::get<Dissolution::FluidSink>(lhs[m1[2]].data).Pressure;
auto &ix_P0 = std::get<Dissolution::Fluid>(lhs[m1[1]].data).Pressure;
NV_Ith_S(ydot, varmap[&ix_P0]) += 11.11 * (NV_Ith_S(y, varmap.at(&ix_P1)) - NV_Ith_S(y, varmap.at(&ix_P0)));
NV_Ith_S(ydot, varmap[&ix_P1]) +=  0.0 ;
// ── Symbolic ODE system ──
// d(dP0)/dt += 11.11 * (P1 - P0)
// d(dP1)/dt += 0.0
{ static bool _sym_dumped = false;
  if (!_sym_dumped && std::getenv("ODE_DUMP")) { _sym_dumped = true;
    std::cout << "  ── Symbolic ODE ──" << std::endl;
    std::cout << "    d(dP0)/dt += 11.11 * (P1 - P0)" << std::endl;
    std::cout << "    d(dP1)/dt += 0.0" << std::endl;
  }
}
}
);
gamma.addRule(sink_flow);
};
void erode_rock(DGGML::Grammar<Dissolution::graph_type> &gamma,
           Dissolution::graph_type &system_graph,
           Parameters &settings) {

GT erode_rock_lhs;
erode_rock_lhs.addNode({1, {Dissolution::Fluid{} }});

erode_rock_lhs.addNode({2, {Dissolution::Fluid{} }});

erode_rock_lhs.addEdge(1, 2);

erode_rock_lhs.addNode({3, {Dissolution::RockStart{} }});

GT erode_rock_rhs;
erode_rock_rhs.addNode({1, {Dissolution::Fluid{} }});

erode_rock_rhs.addNode({2, {Dissolution::Fluid{} }});

erode_rock_rhs.addEdge(1, 2);

erode_rock_rhs.addNode({3, {Dissolution::RockStart{} }});

DGGML::WithRule<GT> erode_rock("erode_rock", erode_rock_lhs, erode_rock_rhs,
[&](auto &lhs, auto &m1) {

double P0 = std::get<Dissolution::Fluid>(lhs[m1[1]].data).Pressure;

double rho = std::get<Dissolution::RockStart>(lhs[m1[3]].data).Density;

torch::Tensor fp0 = std::get<Dissolution::Fluid>(lhs[m1[1]].data).Position;

torch::Tensor rpos = std::get<Dissolution::RockStart>(lhs[m1[3]].data).Position;

return 
(
10.0
 * 
P0
 * 
( (rho > 0.1) ? 1.0 : 0.0)
 * 
( (HELP::distance(fp0[0].template item<double>(), fp0[1].template item<double>(), fp0[2].template item<double>(), rpos[0].template item<double>(), rpos[1].template item<double>(), rpos[2].template item<double>()) <= 0.35) ? 1.0 : 0.0)
)
;},
[&](auto &lhs, auto &rhs, auto &m1, auto &m2) {
torch::Tensor fp0 = std::get<Dissolution::Fluid>(lhs[m1[1]].data).Position;

torch::Tensor fp0n = (fp0);
for (int _i = 0; _i < 3 && _i < fp0n.numel(); _i++) { rhs[m2[1]].position[_i] = fp0n[_i].template item<double>(); }
std::get<Dissolution::Fluid>(rhs[m2[ 1 ]].data).Position = fp0n;
torch::Tensor fu0 = std::get<Dissolution::Fluid>(lhs[m1[1]].data).Unit;

torch::Tensor fu0n = (fu0);
std::get<Dissolution::Fluid>(rhs[m2[ 1 ]].data).Unit = fu0n;
torch::Tensor fc0 = std::get<Dissolution::Fluid>(lhs[m1[1]].data).FCount;

torch::Tensor fc0n = (fc0);
std::get<Dissolution::Fluid>(rhs[m2[ 1 ]].data).FCount = fc0n;
double P0 = std::get<Dissolution::Fluid>(lhs[m1[1]].data).Pressure;

double P0n = (P0);
std::get<Dissolution::Fluid>(rhs[m2[ 1 ]].data).Pressure = P0n;
double C0 = std::get<Dissolution::Fluid>(lhs[m1[1]].data).Concentration;

double C0n = (C0);
std::get<Dissolution::Fluid>(rhs[m2[ 1 ]].data).Concentration = C0n;
torch::Tensor fp1 = std::get<Dissolution::Fluid>(lhs[m1[2]].data).Position;

torch::Tensor fp1n = (fp1);
for (int _i = 0; _i < 3 && _i < fp1n.numel(); _i++) { rhs[m2[2]].position[_i] = fp1n[_i].template item<double>(); }
std::get<Dissolution::Fluid>(rhs[m2[ 2 ]].data).Position = fp1n;
torch::Tensor fu1 = std::get<Dissolution::Fluid>(lhs[m1[2]].data).Unit;

torch::Tensor fu1n = (fu1);
std::get<Dissolution::Fluid>(rhs[m2[ 2 ]].data).Unit = fu1n;
torch::Tensor fc1 = std::get<Dissolution::Fluid>(lhs[m1[2]].data).FCount;

torch::Tensor fc1n = (fc1);
std::get<Dissolution::Fluid>(rhs[m2[ 2 ]].data).FCount = fc1n;
double P1 = std::get<Dissolution::Fluid>(lhs[m1[2]].data).Pressure;

double P1n = (P1);
std::get<Dissolution::Fluid>(rhs[m2[ 2 ]].data).Pressure = P1n;
double C1 = std::get<Dissolution::Fluid>(lhs[m1[2]].data).Concentration;

double C1n = (C1);
std::get<Dissolution::Fluid>(rhs[m2[ 2 ]].data).Concentration = C1n;
torch::Tensor rpos = std::get<Dissolution::RockStart>(lhs[m1[3]].data).Position;

torch::Tensor rposn = (rpos);
for (int _i = 0; _i < 3 && _i < rposn.numel(); _i++) { rhs[m2[3]].position[_i] = rposn[_i].template item<double>(); }
std::get<Dissolution::RockStart>(rhs[m2[ 3 ]].data).Position = rposn;
torch::Tensor rc = std::get<Dissolution::RockStart>(lhs[m1[3]].data).Count;

torch::Tensor rcn = (rc);
std::get<Dissolution::RockStart>(rhs[m2[ 3 ]].data).Count = rcn;
torch::Tensor rdir = std::get<Dissolution::RockStart>(lhs[m1[3]].data).Dir;

torch::Tensor rdirn = (rdir);
std::get<Dissolution::RockStart>(rhs[m2[ 3 ]].data).Dir = rdirn;
double rho = std::get<Dissolution::RockStart>(lhs[m1[3]].data).Density;

double rhon = (rho - 0.1);
std::get<Dissolution::RockStart>(rhs[m2[ 3 ]].data).Density = rhon;
}
);
gamma.addRule(erode_rock);
};
}
#endif