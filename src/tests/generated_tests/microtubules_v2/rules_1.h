#ifndef DGGML_RULES_STAGE_1_HPP
#define DGGML_RULES_STAGE_1_HPP
#include "types.h"
#include "parameters.h"
#include "functions.h"
namespace Microtubule {
using GT = Microtubule::graph_type;
void stochastic_mt_growth(DGGML::Grammar<Microtubule::graph_type> &gamma,
           Microtubule::graph_type &system_graph,
           Parameters &settings) {

GT stochastic_mt_growth_lhs;
stochastic_mt_growth_lhs.addNode({1, {Microtubule::Intermediate{} }});

stochastic_mt_growth_lhs.addNode({2, {Microtubule::Positive{} }});

stochastic_mt_growth_lhs.addEdge(1, 2);

GT stochastic_mt_growth_rhs;
stochastic_mt_growth_rhs.addNode({1, {Microtubule::Intermediate{} }});

stochastic_mt_growth_rhs.addNode({3, {Microtubule::Intermediate{} }});

stochastic_mt_growth_rhs.addEdge(1, 3);

stochastic_mt_growth_rhs.addNode({2, {Microtubule::Positive{} }});

stochastic_mt_growth_rhs.addEdge(3, 2);

DGGML::WithRule<GT> stochastic_mt_growth("stochastic_mt_growth", stochastic_mt_growth_lhs, stochastic_mt_growth_rhs,
[&](auto &lhs, auto &m1) {

torch::Tensor im_pos = torch::tensor({lhs[m1[1]].position[0], lhs[m1[1]].position[1], lhs[m1[1]].position[2]}, torch::kFloat64);

auto DIV_LENGTH =  settings.DIV_LENGTH;

torch::Tensor pos_pos = torch::tensor({lhs[m1[2]].position[0], lhs[m1[2]].position[1], lhs[m1[2]].position[2]}, torch::kFloat64);

return 
(
settings.WITH_GROWTH_RATE_FACTOR
 * 
DGGML::heaviside(HELP::distance(im_pos[0].template item<double>(), im_pos[1].template item<double>(), pos_pos[0].template item<double>(), pos_pos[1].template item<double>()), DIV_LENGTH)
)
;},
[&](auto &lhs, auto &rhs, auto &m1, auto &m2) {
std::copy(std::begin(lhs[m1[1]].position), std::end(lhs[m1[1]].position), std::begin(rhs[m2[1]].position));

std::get<Microtubule::Intermediate>(rhs[m2[ 1 ]].data).Direction = std::get<Microtubule::Intermediate>(lhs[m1[ 1 ]].data).Direction.clone();

std::copy(std::begin(lhs[m1[2]].position), std::end(lhs[m1[2]].position), std::begin(rhs[m2[2]].position));

std::get<Microtubule::Positive>(rhs[m2[ 2 ]].data).Direction = std::get<Microtubule::Positive>(lhs[m1[ 2 ]].data).Direction.clone();

torch::Tensor pos_pos = torch::tensor({rhs[m2[2]].position[0], rhs[m2[2]].position[1], rhs[m2[2]].position[2]}, torch::kFloat64);

torch::Tensor im_pos = torch::tensor({rhs[m2[1]].position[0], rhs[m2[1]].position[1], rhs[m2[1]].position[2]}, torch::kFloat64);

rhs[m2[3]].position[(0)] = static_cast<double>((pos_pos[0].template item<double>() - (pos_pos[0].template item<double>() - im_pos[0].template item<double>()) / 100.0));
rhs[m2[3]].position[(1)] = static_cast<double>((pos_pos[1].template item<double>() - (pos_pos[1].template item<double>() - im_pos[1].template item<double>()) / 100.0));
rhs[m2[3]].position[(2)] = static_cast<double>((pos_pos[2].template item<double>() - (pos_pos[2].template item<double>() - im_pos[2].template item<double>()) / 100.0));
torch::Tensor im_dir = std::get<Microtubule::Intermediate>(rhs[m2[1]].data).Direction.clone();

std::get<Microtubule::Intermediate>(rhs[m2[ 3 ]].data).Direction = std::get<Microtubule::Intermediate>(rhs[m2[ 3 ]].data).Direction.clone();
std::get<Microtubule::Intermediate>(rhs[m2[ 3 ]].data).Direction[(0)] = (im_dir[0].template item<double>());
std::get<Microtubule::Intermediate>(rhs[m2[ 3 ]].data).Direction[(1)] = (im_dir[1].template item<double>());
std::get<Microtubule::Intermediate>(rhs[m2[ 3 ]].data).Direction[(2)] = (im_dir[2].template item<double>());
std::random_device random_device;

std::mt19937 random_engine(random_device());

double angle = (std::uniform_real_distribution<double>(-settings.WOBBLE_ANGLE * 0.0174533, settings.WOBBLE_ANGLE * 0.0174533)(random_engine));
torch::Tensor pos_dir = std::get<Microtubule::Positive>(rhs[m2[2]].data).Direction.clone();

std::get<Microtubule::Positive>(rhs[m2[ 2 ]].data).Direction = std::get<Microtubule::Positive>(rhs[m2[ 2 ]].data).Direction.clone();
std::get<Microtubule::Positive>(rhs[m2[ 2 ]].data).Direction[(0)] = (pos_dir[0].template item<double>());
std::get<Microtubule::Positive>(rhs[m2[ 2 ]].data).Direction[(1)] = (pos_dir[1].template item<double>());
}
);
gamma.addRule(stochastic_mt_growth);
};
void ode_mt_growth(DGGML::Grammar<Microtubule::graph_type> &gamma,
           Microtubule::graph_type &system_graph,
           Parameters &settings) {

GT ode_mt_growth_lhs;
ode_mt_growth_lhs.addNode({1, {Microtubule::Intermediate{} }});

ode_mt_growth_lhs.addNode({2, {Microtubule::Positive{} }});

ode_mt_growth_lhs.addEdge(1, 2);

GT ode_mt_growth_rhs;
ode_mt_growth_rhs.addNode({1, {Microtubule::Intermediate{} }});

ode_mt_growth_rhs.addNode({2, {Microtubule::Positive{} }});

ode_mt_growth_rhs.addEdge(1, 2);

DGGML::SolvingRule<GT> ode_mt_growth("ode_mt_growth", ode_mt_growth_lhs, ode_mt_growth_lhs,
3,
[](auto &lhs, auto &m1, auto &varset) {
double* tensor_ptr_2_1 = &lhs[m1[2]].position[0];
varset.insert(&tensor_ptr_2_1[0]);
varset.insert(&tensor_ptr_2_1[1]);
varset.insert(&tensor_ptr_2_1[2]);
},
[&](auto &lhs, auto &m1, auto y, auto ydot, auto &varmap) {
double* tensor_ptr_2_1 = &lhs[m1[2]].position[0];
torch::Tensor pos_dir = std::get<Microtubule::Positive>(lhs[m1[2]].data).Direction.clone();

NV_Ith_S(ydot, varmap[&tensor_ptr_2_1[0]]) += settings.V_PLUS *  pos_dir [ 0 ].template item<double>();
NV_Ith_S(ydot, varmap[&tensor_ptr_2_1[1]]) += settings.V_PLUS *  pos_dir [ 1 ].template item<double>();
NV_Ith_S(ydot, varmap[&tensor_ptr_2_1[2]]) += settings.V_PLUS *  pos_dir [ 2 ].template item<double>();
// ── Symbolic ODE system ──
// d(dpos[0])/dt += V_PLUS * pos_dir[0]
// d(dpos[1])/dt += V_PLUS * pos_dir[1]
// d(dpos[2])/dt += V_PLUS * pos_dir[2]
{ static bool _sym_dumped = false;
  if (!_sym_dumped && std::getenv("ODE_DUMP")) { _sym_dumped = true;
    std::cout << "  ── Symbolic ODE ──" << std::endl;
    std::cout << "    d(dpos[0])/dt += V_PLUS * pos_dir[0]" << std::endl;
    std::cout << "    d(dpos[1])/dt += V_PLUS * pos_dir[1]" << std::endl;
    std::cout << "    d(dpos[2])/dt += V_PLUS * pos_dir[2]" << std::endl;
  }
}
}
);
gamma.addRule(ode_mt_growth);
};
void mt_stochastic_retraction(DGGML::Grammar<Microtubule::graph_type> &gamma,
           Microtubule::graph_type &system_graph,
           Parameters &settings) {

GT mt_stochastic_retraction_lhs;
mt_stochastic_retraction_lhs.addNode({1, {Microtubule::Negative{} }});

mt_stochastic_retraction_lhs.addNode({2, {Microtubule::Intermediate{} }});

mt_stochastic_retraction_lhs.addEdge(1, 2);

mt_stochastic_retraction_lhs.addNode({3, {Microtubule::Intermediate{} }});

mt_stochastic_retraction_lhs.addEdge(2, 3);

GT mt_stochastic_retraction_rhs;
mt_stochastic_retraction_rhs.addNode({1, {Microtubule::Negative{} }});

mt_stochastic_retraction_rhs.addNode({3, {Microtubule::Intermediate{} }});

mt_stochastic_retraction_rhs.addEdge(1, 3);

DGGML::WithRule<GT> mt_stochastic_retraction("mt_stochastic_retraction", mt_stochastic_retraction_lhs, mt_stochastic_retraction_rhs,
[&](auto &lhs, auto &m1) {

auto DIV_LENGTH_RETRACT =  settings.DIV_LENGTH_RETRACT;

torch::Tensor im1_pos = torch::tensor({lhs[m1[2]].position[0], lhs[m1[2]].position[1], lhs[m1[2]].position[2]}, torch::kFloat64);

torch::Tensor neg_pos = torch::tensor({lhs[m1[1]].position[0], lhs[m1[1]].position[1], lhs[m1[1]].position[2]}, torch::kFloat64);

return 
(
settings.WITH_RETRACTION_RATE_FACTOR
 * 
DGGML::heaviside(DIV_LENGTH_RETRACT, HELP::distance(neg_pos[0].template item<double>(), neg_pos[1].template item<double>(), im1_pos[0].template item<double>(), im1_pos[1].template item<double>()))
)
;},
[&](auto &lhs, auto &rhs, auto &m1, auto &m2) {
std::copy(std::begin(lhs[m1[1]].position), std::end(lhs[m1[1]].position), std::begin(rhs[m2[1]].position));

std::get<Microtubule::Negative>(rhs[m2[ 1 ]].data).Direction = std::get<Microtubule::Negative>(lhs[m1[ 1 ]].data).Direction.clone();

std::copy(std::begin(lhs[m1[3]].position), std::end(lhs[m1[3]].position), std::begin(rhs[m2[3]].position));

torch::Tensor im2_pos = torch::tensor({lhs[m1[3]].position[0], lhs[m1[3]].position[1], lhs[m1[3]].position[2]}, torch::kFloat64);

torch::Tensor neg_pos = torch::tensor({lhs[m1[1]].position[0], lhs[m1[1]].position[1], lhs[m1[1]].position[2]}, torch::kFloat64);

double dx = (im2_pos[0].template item<double>() - neg_pos[0].template item<double>());
double dy = (im2_pos[1].template item<double>() - neg_pos[0].template item<double>());
double dz = (im2_pos[2].template item<double>() - neg_pos[2].template item<double>());
double len = (sqrt((((dx * dx) + (dy * dy)) + (dz * dz))));
std::get<Microtubule::Negative>(rhs[m2[ 1 ]].data).Direction = std::get<Microtubule::Negative>(rhs[m2[ 1 ]].data).Direction.clone();
std::get<Microtubule::Negative>(rhs[m2[ 1 ]].data).Direction[(0)] = (dx / len);
std::get<Microtubule::Negative>(rhs[m2[ 1 ]].data).Direction[(1)] = (dy / len);
std::get<Microtubule::Negative>(rhs[m2[ 1 ]].data).Direction[(2)] = (dz / len);
}
);
gamma.addRule(mt_stochastic_retraction);
};
void mt_ode_retraction(DGGML::Grammar<Microtubule::graph_type> &gamma,
           Microtubule::graph_type &system_graph,
           Parameters &settings) {

GT mt_ode_retraction_lhs;
mt_ode_retraction_lhs.addNode({1, {Microtubule::Negative{} }});

mt_ode_retraction_lhs.addNode({2, {Microtubule::Intermediate{} }});

mt_ode_retraction_lhs.addEdge(1, 2);

GT mt_ode_retraction_rhs;
mt_ode_retraction_rhs.addNode({1, {Microtubule::Negative{} }});

mt_ode_retraction_rhs.addNode({2, {Microtubule::Intermediate{} }});

mt_ode_retraction_rhs.addEdge(1, 2);

DGGML::SolvingRule<GT> mt_ode_retraction("mt_ode_retraction", mt_ode_retraction_lhs, mt_ode_retraction_lhs,
3,
[](auto &lhs, auto &m1, auto &varset) {
double* tensor_ptr_1_1 = &lhs[m1[1]].position[0];
varset.insert(&tensor_ptr_1_1[0]);
varset.insert(&tensor_ptr_1_1[1]);
varset.insert(&tensor_ptr_1_1[2]);
},
[&](auto &lhs, auto &m1, auto y, auto ydot, auto &varmap) {
double* tensor_ptr_1_1 = &lhs[m1[1]].position[0];
torch::Tensor neg_dir = std::get<Microtubule::Negative>(lhs[m1[1]].data).Direction.clone();

torch::Tensor im_pos = torch::tensor({lhs[m1[2]].position[0], lhs[m1[2]].position[1], lhs[m1[2]].position[2]}, torch::kFloat64);

auto DIV_LENGTH_RETRACT =  settings.DIV_LENGTH_RETRACT;

torch::Tensor neg_pos = torch::tensor({lhs[m1[1]].position[0], lhs[m1[1]].position[1], lhs[m1[1]].position[2]}, torch::kFloat64);

NV_Ith_S(ydot, varmap[&tensor_ptr_1_1[0]]) +=  -  1  * settings.V_MINUS *  neg_dir [ 0 ].template item<double>() * DGGML::heaviside(HELP::distance(neg_pos[0].template item<double>(), neg_pos[1].template item<double>(), im_pos[0].template item<double>(), im_pos[1].template item<double>()), (DIV_LENGTH_RETRACT / 2.0));
NV_Ith_S(ydot, varmap[&tensor_ptr_1_1[1]]) +=  -  1  * settings.V_MINUS *  neg_dir [ 1 ].template item<double>() * DGGML::heaviside(HELP::distance(neg_pos[0].template item<double>(), neg_pos[1].template item<double>(), im_pos[0].template item<double>(), im_pos[1].template item<double>()), (DIV_LENGTH_RETRACT / 2.0));
NV_Ith_S(ydot, varmap[&tensor_ptr_1_1[2]]) +=  -  1  * settings.V_MINUS *  neg_dir [ 2 ].template item<double>() * DGGML::heaviside(HELP::distance(neg_pos[0].template item<double>(), neg_pos[1].template item<double>(), im_pos[0].template item<double>(), im_pos[1].template item<double>()), (DIV_LENGTH_RETRACT / 2.0));
// ── Symbolic ODE system ──
// d(dneg[0])/dt += -1 * V_MINUS * neg_dir[0] * heaviside(HELP::distance(neg_pos[0], neg_pos[1], im_pos[0], im_pos[1]), DIV_LENGTH_RETRACT / 2.0)
// d(dneg[1])/dt += -1 * V_MINUS * neg_dir[1] * heaviside(HELP::distance(neg_pos[0], neg_pos[1], im_pos[0], im_pos[1]), DIV_LENGTH_RETRACT / 2.0)
// d(dneg[2])/dt += -1 * V_MINUS * neg_dir[2] * heaviside(HELP::distance(neg_pos[0], neg_pos[1], im_pos[0], im_pos[1]), DIV_LENGTH_RETRACT / 2.0)
{ static bool _sym_dumped = false;
  if (!_sym_dumped && std::getenv("ODE_DUMP")) { _sym_dumped = true;
    std::cout << "  ── Symbolic ODE ──" << std::endl;
    std::cout << "    d(dneg[0])/dt += -1 * V_MINUS * neg_dir[0] * heaviside(HELP::distance(neg_pos[0], neg_pos[1], im_pos[0], im_pos[1]), DIV_LENGTH_RETRACT / 2.0)" << std::endl;
    std::cout << "    d(dneg[1])/dt += -1 * V_MINUS * neg_dir[1] * heaviside(HELP::distance(neg_pos[0], neg_pos[1], im_pos[0], im_pos[1]), DIV_LENGTH_RETRACT / 2.0)" << std::endl;
    std::cout << "    d(dneg[2])/dt += -1 * V_MINUS * neg_dir[2] * heaviside(HELP::distance(neg_pos[0], neg_pos[1], im_pos[0], im_pos[1]), DIV_LENGTH_RETRACT / 2.0)" << std::endl;
  }
}
}
);
gamma.addRule(mt_ode_retraction);
};
void boundary_catastrophe1(DGGML::Grammar<Microtubule::graph_type> &gamma,
           Microtubule::graph_type &system_graph,
           Parameters &settings) {

GT boundary_catastrophe1_lhs;
boundary_catastrophe1_lhs.addNode({1, {Microtubule::Intermediate{} }});

boundary_catastrophe1_lhs.addNode({2, {Microtubule::Positive{} }});

boundary_catastrophe1_lhs.addEdge(1, 2);

boundary_catastrophe1_lhs.addNode({3, {Microtubule::CellBoundary{} }});

boundary_catastrophe1_lhs.addNode({4, {Microtubule::CellBoundary{} }});

boundary_catastrophe1_lhs.addEdge(3, 4);

GT boundary_catastrophe1_rhs;
boundary_catastrophe1_rhs.addNode({1, {Microtubule::Intermediate{} }});

boundary_catastrophe1_rhs.addNode({2, {Microtubule::Negative{} }});

boundary_catastrophe1_rhs.addEdge(1, 2);

boundary_catastrophe1_rhs.addNode({3, {Microtubule::CellBoundary{} }});

boundary_catastrophe1_rhs.addNode({4, {Microtubule::CellBoundary{} }});

boundary_catastrophe1_rhs.addEdge(3, 4);

DGGML::WithRule<GT> boundary_catastrophe1("boundary_catastrophe1", boundary_catastrophe1_lhs, boundary_catastrophe1_rhs,
[&](auto &lhs, auto &m1) {

auto COLLISION_DISTANCE_BOUNDARY =  settings.COLLISION_DISTANCE_BOUNDARY;

torch::Tensor b0_pos = torch::tensor({lhs[m1[3]].position[0], lhs[m1[3]].position[1], lhs[m1[3]].position[2]}, torch::kFloat64);

torch::Tensor b1_pos = torch::tensor({lhs[m1[4]].position[0], lhs[m1[4]].position[1], lhs[m1[4]].position[2]}, torch::kFloat64);

torch::Tensor pos_pos = torch::tensor({lhs[m1[2]].position[0], lhs[m1[2]].position[1], lhs[m1[2]].position[2]}, torch::kFloat64);

return 
(
( (DGGML::distanceToLineSegment(b0_pos[0].template item<double>(), b0_pos[1].template item<double>(), b1_pos[0].template item<double>(), b1_pos[1].template item<double>(), pos_pos[0].template item<double>(), pos_pos[1].template item<double>()) <= COLLISION_DISTANCE_BOUNDARY) ? 1.0 : 0.0)
 * 
settings.STANDARD_BOUNDARY_CATASTROPHE_RATE
)
;},
[&](auto &lhs, auto &rhs, auto &m1, auto &m2) {
std::copy(std::begin(lhs[m1[1]].position), std::end(lhs[m1[1]].position), std::begin(rhs[m2[1]].position));

std::copy(std::begin(lhs[m1[3]].position), std::end(lhs[m1[3]].position), std::begin(rhs[m2[3]].position));

std::copy(std::begin(lhs[m1[4]].position), std::end(lhs[m1[4]].position), std::begin(rhs[m2[4]].position));

torch::Tensor pos_pos = torch::tensor({lhs[m1[2]].position[0], lhs[m1[2]].position[1], lhs[m1[2]].position[2]}, torch::kFloat64);

rhs[m2[2]].position[(0)] = static_cast<double>((pos_pos[0].template item<double>()));
rhs[m2[2]].position[(1)] = static_cast<double>((pos_pos[1].template item<double>()));
rhs[m2[2]].position[(2)] = static_cast<double>((pos_pos[2].template item<double>()));
torch::Tensor pos_dir = std::get<Microtubule::Positive>(lhs[m1[2]].data).Direction.clone();

std::get<Microtubule::Negative>(rhs[m2[ 2 ]].data).Direction = std::get<Microtubule::Negative>(rhs[m2[ 2 ]].data).Direction.clone();
std::get<Microtubule::Negative>(rhs[m2[ 2 ]].data).Direction[(0)] = (pos_dir[0].template item<double>());
std::get<Microtubule::Negative>(rhs[m2[ 2 ]].data).Direction[(1)] = (pos_dir[1].template item<double>());
std::get<Microtubule::Negative>(rhs[m2[ 2 ]].data).Direction[(2)] = (pos_dir[2].template item<double>());
}
);
gamma.addRule(boundary_catastrophe1);
};
void boundary_clamp(DGGML::Grammar<Microtubule::graph_type> &gamma,
           Microtubule::graph_type &system_graph,
           Parameters &settings) {

GT boundary_clamp_lhs;
boundary_clamp_lhs.addNode({1, {Microtubule::Intermediate{} }});

boundary_clamp_lhs.addNode({2, {Microtubule::Positive{} }});

boundary_clamp_lhs.addEdge(1, 2);

GT boundary_clamp_rhs;
boundary_clamp_rhs.addNode({1, {Microtubule::Intermediate{} }});

boundary_clamp_rhs.addNode({2, {Microtubule::Negative{} }});

boundary_clamp_rhs.addEdge(1, 2);

DGGML::WithRule<GT> boundary_clamp("boundary_clamp", boundary_clamp_lhs, boundary_clamp_rhs,
[&](auto &lhs, auto &m1) {

auto DOMAIN_MAX =  settings.DOMAIN_MAX;

torch::Tensor pos_pos = torch::tensor({lhs[m1[2]].position[0], lhs[m1[2]].position[1], lhs[m1[2]].position[2]}, torch::kFloat64);

auto DOMAIN_MIN =  settings.DOMAIN_MIN;

return 
(
(
DGGML::heaviside(pos_pos[0].template item<double>(), DOMAIN_MAX)
 + 
DGGML::heaviside(DOMAIN_MIN, pos_pos[0].template item<double>())
 + 
DGGML::heaviside(pos_pos[1].template item<double>(), DOMAIN_MAX)
 + 
DGGML::heaviside(DOMAIN_MIN, pos_pos[1].template item<double>())
)
 * 
settings.BOUNDARY_CLAMP_RATE
)
;},
[&](auto &lhs, auto &rhs, auto &m1, auto &m2) {
std::copy(std::begin(lhs[m1[1]].position), std::end(lhs[m1[1]].position), std::begin(rhs[m2[1]].position));

torch::Tensor pos_pos = torch::tensor({lhs[m1[2]].position[0], lhs[m1[2]].position[1], lhs[m1[2]].position[2]}, torch::kFloat64);

rhs[m2[2]].position[(0)] = static_cast<double>((pos_pos[0].template item<double>()));
rhs[m2[2]].position[(1)] = static_cast<double>((pos_pos[1].template item<double>()));
rhs[m2[2]].position[(2)] = static_cast<double>((pos_pos[2].template item<double>()));
torch::Tensor pos_dir = std::get<Microtubule::Positive>(lhs[m1[2]].data).Direction.clone();

std::get<Microtubule::Negative>(rhs[m2[ 2 ]].data).Direction = std::get<Microtubule::Negative>(rhs[m2[ 2 ]].data).Direction.clone();
std::get<Microtubule::Negative>(rhs[m2[ 2 ]].data).Direction[(0)] = (pos_dir[0].template item<double>());
std::get<Microtubule::Negative>(rhs[m2[ 2 ]].data).Direction[(1)] = (pos_dir[1].template item<double>());
std::get<Microtubule::Negative>(rhs[m2[ 2 ]].data).Direction[(2)] = (pos_dir[2].template item<double>());
}
);
gamma.addRule(boundary_clamp);
};
void catastrophe2_case1(DGGML::Grammar<Microtubule::graph_type> &gamma,
           Microtubule::graph_type &system_graph,
           Parameters &settings) {

GT catastrophe2_case1_lhs;
catastrophe2_case1_lhs.addNode({1, {Microtubule::Intermediate{} }});

catastrophe2_case1_lhs.addNode({2, {Microtubule::Positive{} }});

catastrophe2_case1_lhs.addEdge(1, 2);

catastrophe2_case1_lhs.addNode({3, {Microtubule::Intermediate{} }});

catastrophe2_case1_lhs.addNode({4, {Microtubule::Intermediate{} }});

catastrophe2_case1_lhs.addEdge(3, 4);

GT catastrophe2_case1_rhs;
catastrophe2_case1_rhs.addNode({1, {Microtubule::Intermediate{} }});

catastrophe2_case1_rhs.addNode({2, {Microtubule::Negative{} }});

catastrophe2_case1_rhs.addEdge(1, 2);

catastrophe2_case1_rhs.addNode({3, {Microtubule::Intermediate{} }});

catastrophe2_case1_rhs.addNode({4, {Microtubule::Intermediate{} }});

catastrophe2_case1_rhs.addEdge(3, 4);

DGGML::WithRule<GT> catastrophe2_case1("catastrophe2_case1", catastrophe2_case1_lhs, catastrophe2_case1_rhs,
[&](auto &lhs, auto &m1) {

auto COLLISION_DISTANCE =  settings.COLLISION_DISTANCE;

torch::Tensor im3_pos = torch::tensor({lhs[m1[4]].position[0], lhs[m1[4]].position[1], lhs[m1[4]].position[2]}, torch::kFloat64);

torch::Tensor pos_pos = torch::tensor({lhs[m1[2]].position[0], lhs[m1[2]].position[1], lhs[m1[2]].position[2]}, torch::kFloat64);

torch::Tensor im2_pos = torch::tensor({lhs[m1[3]].position[0], lhs[m1[3]].position[1], lhs[m1[3]].position[2]}, torch::kFloat64);

return 
(
( (DGGML::distanceToLineSegment(im2_pos[0].template item<double>(), im2_pos[1].template item<double>(), im3_pos[0].template item<double>(), im3_pos[1].template item<double>(), pos_pos[0].template item<double>(), pos_pos[1].template item<double>()) <= COLLISION_DISTANCE) ? 1.0 : 0.0)
 * 
settings.INTERMEDIATE_CIC_RATE
)
;},
[&](auto &lhs, auto &rhs, auto &m1, auto &m2) {
std::copy(std::begin(lhs[m1[1]].position), std::end(lhs[m1[1]].position), std::begin(rhs[m2[1]].position));

std::copy(std::begin(lhs[m1[3]].position), std::end(lhs[m1[3]].position), std::begin(rhs[m2[3]].position));

std::copy(std::begin(lhs[m1[4]].position), std::end(lhs[m1[4]].position), std::begin(rhs[m2[4]].position));

torch::Tensor pos_pos = torch::tensor({lhs[m1[2]].position[0], lhs[m1[2]].position[1], lhs[m1[2]].position[2]}, torch::kFloat64);

rhs[m2[2]].position[(0)] = static_cast<double>((pos_pos[0].template item<double>()));
rhs[m2[2]].position[(1)] = static_cast<double>((pos_pos[1].template item<double>()));
rhs[m2[2]].position[(2)] = static_cast<double>((pos_pos[2].template item<double>()));
torch::Tensor pos_dir = std::get<Microtubule::Positive>(lhs[m1[2]].data).Direction.clone();

std::get<Microtubule::Negative>(rhs[m2[ 2 ]].data).Direction = std::get<Microtubule::Negative>(rhs[m2[ 2 ]].data).Direction.clone();
std::get<Microtubule::Negative>(rhs[m2[ 2 ]].data).Direction[(0)] = (-pos_dir[0].template item<double>());
std::get<Microtubule::Negative>(rhs[m2[ 2 ]].data).Direction[(1)] = (-pos_dir[1].template item<double>());
std::get<Microtubule::Negative>(rhs[m2[ 2 ]].data).Direction[(2)] = (-pos_dir[2].template item<double>());
}
);
gamma.addRule(catastrophe2_case1);
};
void creation_case1(DGGML::Grammar<Microtubule::graph_type> &gamma,
           Microtubule::graph_type &system_graph,
           Parameters &settings) {

GT creation_case1_lhs;
creation_case1_lhs.addNode({1, {Microtubule::Nucleator{} }});

GT creation_case1_rhs;
creation_case1_rhs.addNode({1, {Microtubule::Nucleator{} }});

creation_case1_rhs.addNode({2, {Microtubule::Negative{} }});

creation_case1_rhs.addNode({3, {Microtubule::Intermediate{} }});

creation_case1_rhs.addEdge(2, 3);

creation_case1_rhs.addNode({4, {Microtubule::Positive{} }});

creation_case1_rhs.addEdge(3, 4);

DGGML::WithRule<GT> creation_case1("creation_case1", creation_case1_lhs, creation_case1_rhs,
[&](auto &lhs, auto &m1) {

return 
(
settings.CREATION_FACTOR
 * 
settings.CREATION_RATE
)
;},
[&](auto &lhs, auto &rhs, auto &m1, auto &m2) {
std::copy(std::begin(lhs[m1[1]].position), std::end(lhs[m1[1]].position), std::begin(rhs[m2[1]].position));

std::random_device random_device;

std::mt19937 random_engine(random_device());

double seg_len = (std::uniform_real_distribution<double>(settings.MT_MIN_SEGMENT_INIT, settings.MT_MAX_SEGMENT_INIT)(random_engine));
double theta = (std::uniform_real_distribution<double>(0.0, 6.2831853)(random_engine));
torch::Tensor nuc_pos = torch::tensor({rhs[m2[1]].position[0], rhs[m2[1]].position[1], rhs[m2[1]].position[2]}, torch::kFloat64);

rhs[m2[3]].position[(0)] = static_cast<double>((nuc_pos[0].template item<double>()));
rhs[m2[3]].position[(1)] = static_cast<double>((nuc_pos[1].template item<double>()));
rhs[m2[3]].position[(2)] = static_cast<double>((0.0));
torch::Tensor pos_pos = torch::tensor({rhs[m2[4]].position[0], rhs[m2[4]].position[1], rhs[m2[4]].position[2]}, torch::kFloat64);

torch::Tensor neg_pos = torch::tensor({rhs[m2[2]].position[0], rhs[m2[2]].position[1], rhs[m2[2]].position[2]}, torch::kFloat64);

torch::Tensor im_dir_vector = (HELP::calculate_unit_vector(neg_pos[0].template item<double>(), neg_pos[1].template item<double>(), pos_pos[0].template item<double>(), pos_pos[1].template item<double>()));
rhs[m2[4]].position[(0)] = static_cast<double>((nuc_pos[0].template item<double>() + seg_len * sin(theta)));
pos_pos = torch::tensor({rhs[m2[4]].position[0], rhs[m2[4]].position[1], rhs[m2[4]].position[2]}, torch::kFloat64);
rhs[m2[4]].position[(1)] = static_cast<double>((nuc_pos[1].template item<double>() + seg_len * cos(theta)));
pos_pos = torch::tensor({rhs[m2[4]].position[0], rhs[m2[4]].position[1], rhs[m2[4]].position[2]}, torch::kFloat64);
rhs[m2[4]].position[(2)] = static_cast<double>((0.0));
pos_pos = torch::tensor({rhs[m2[4]].position[0], rhs[m2[4]].position[1], rhs[m2[4]].position[2]}, torch::kFloat64);
torch::Tensor im_pos = torch::tensor({rhs[m2[3]].position[0], rhs[m2[3]].position[1], rhs[m2[3]].position[2]}, torch::kFloat64);

torch::Tensor pos_dir_vector = (HELP::calculate_unit_vector(im_pos[0].template item<double>(), im_pos[1].template item<double>(), pos_pos[0].template item<double>(), pos_pos[1].template item<double>()));
std::get<Microtubule::Positive>(rhs[m2[ 4 ]].data).Direction = std::get<Microtubule::Positive>(rhs[m2[ 4 ]].data).Direction.clone();
std::get<Microtubule::Positive>(rhs[m2[ 4 ]].data).Direction[(0)] = (pos_dir_vector[0].template item<double>());
std::get<Microtubule::Positive>(rhs[m2[ 4 ]].data).Direction[(1)] = (pos_dir_vector[1].template item<double>());
std::get<Microtubule::Positive>(rhs[m2[ 4 ]].data).Direction[(2)] = (0.0);
rhs[m2[2]].position[(0)] = static_cast<double>((nuc_pos[0].template item<double>() - seg_len * sin(theta)));
neg_pos = torch::tensor({rhs[m2[2]].position[0], rhs[m2[2]].position[1], rhs[m2[2]].position[2]}, torch::kFloat64);
rhs[m2[2]].position[(1)] = static_cast<double>((nuc_pos[1].template item<double>() - seg_len * cos(theta)));
neg_pos = torch::tensor({rhs[m2[2]].position[0], rhs[m2[2]].position[1], rhs[m2[2]].position[2]}, torch::kFloat64);
rhs[m2[2]].position[(2)] = static_cast<double>((0.0));
neg_pos = torch::tensor({rhs[m2[2]].position[0], rhs[m2[2]].position[1], rhs[m2[2]].position[2]}, torch::kFloat64);
torch::Tensor neg_dir_vector = (HELP::calculate_unit_vector(im_pos[0].template item<double>(), im_pos[1].template item<double>(), neg_pos[0].template item<double>(), neg_pos[1].template item<double>()));
std::get<Microtubule::Negative>(rhs[m2[ 2 ]].data).Direction = std::get<Microtubule::Negative>(rhs[m2[ 2 ]].data).Direction.clone();
std::get<Microtubule::Negative>(rhs[m2[ 2 ]].data).Direction[(0)] = (-neg_dir_vector[0].template item<double>());
std::get<Microtubule::Negative>(rhs[m2[ 2 ]].data).Direction[(1)] = (-neg_dir_vector[1].template item<double>());
std::get<Microtubule::Negative>(rhs[m2[ 2 ]].data).Direction[(2)] = (0.0);
}
);
gamma.addRule(creation_case1);
};
void destruction_case2(DGGML::Grammar<Microtubule::graph_type> &gamma,
           Microtubule::graph_type &system_graph,
           Parameters &settings) {

GT destruction_case2_lhs;
destruction_case2_lhs.addNode({1, {Microtubule::Negative{} }});

destruction_case2_lhs.addNode({2, {Microtubule::Intermediate{} }});

destruction_case2_lhs.addEdge(1, 2);

destruction_case2_lhs.addNode({3, {Microtubule::Positive{} }});

destruction_case2_lhs.addEdge(2, 3);

destruction_case2_lhs.addNode({4, {Microtubule::Intermediate{} }});

GT destruction_case2_rhs;
destruction_case2_rhs.addNode({4, {Microtubule::Intermediate{} }});

DGGML::WithRule<GT> destruction_case2("destruction_case2", destruction_case2_lhs, destruction_case2_rhs,
[&](auto &lhs, auto &m1) {

auto COLLISION_DISTANCE =  settings.COLLISION_DISTANCE;

torch::Tensor im1_pos = torch::tensor({lhs[m1[2]].position[0], lhs[m1[2]].position[1], lhs[m1[2]].position[2]}, torch::kFloat64);

torch::Tensor im2_pos = torch::tensor({lhs[m1[4]].position[0], lhs[m1[4]].position[1], lhs[m1[4]].position[2]}, torch::kFloat64);

return 
(
( (HELP::distance(im1_pos[0].template item<double>(), im1_pos[1].template item<double>(), im2_pos[0].template item<double>(), im2_pos[1].template item<double>()) < COLLISION_DISTANCE) ? 1.0 : 0.0)
 * 
settings.MT_DESTRUCTION_RATE
)
;},
[&](auto &lhs, auto &rhs, auto &m1, auto &m2) {
std::copy(std::begin(lhs[m1[4]].position), std::end(lhs[m1[4]].position), std::begin(rhs[m2[4]].position));

}
);
gamma.addRule(destruction_case2);
};
}
#endif