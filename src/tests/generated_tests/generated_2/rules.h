#ifndef DGGML_RULES_HPP
#define DGGML_RULES_HPP
#include "types.h"
#include "parameters.h"
namespace particle_rules {
using GT = Particles::graph_type;
void start_to_node(DGGML::Grammar<Particles::graph_type> &gamma,
           Particles::graph_type &system_graph,
           Parameters &settings) {

GT start_to_node_lhs;
start_to_node_lhs.addNode({1, {Particles::StartType{} }});

GT start_to_node_rhs;
start_to_node_rhs.addNode({1, {Particles::ParticleNodeCreator{} }});

DGGML::WithRule<GT> start_to_node("start_to_node", start_to_node_lhs, start_to_node_rhs,
                                 [&](auto &lhs, auto &m) {

return DGGML::heaviside(10,1);},
[&](auto &lhs, auto &rhs, auto &m1, auto &m2) {
float start_x = std::get<Particles::StartType>(lhs[m1[1]].data).start_location[0];

float start_y = std::get<Particles::StartType>(lhs[m1[1]].data).start_location[1];

Particles::ParticleNodeCreator rhs_node_1 = std::get<Particles::ParticleNodeCreator>(rhs[m2[1]].data);
		std::get<Particles::ParticleNodeCreator>(rhs[m2[1]].data).fflow_1c7d3f = 9;
		rhs[m2[1]].position[0] = start_x+1;
		std::get<Particles::ParticleNodeCreator>(rhs[m2[1]].data).fflow_8df356[0] = start_x+1;
		rhs[m2[1]].position[1] = start_y+1;
		std::get<Particles::ParticleNodeCreator>(rhs[m2[1]].data).fflow_8df356[1] = start_y+1;
}
);
gamma.addRule(start_to_node);
};
void create_particle(DGGML::Grammar<Particles::graph_type> &gamma,
           Particles::graph_type &system_graph,
           Parameters &settings) {

GT create_particle_lhs;
create_particle_lhs.addNode({1, {Particles::ParticleNodeCreator{} }});

GT create_particle_rhs;
create_particle_rhs.addNode({1, {Particles::ParticleNodeCreator{} }});

create_particle_rhs.addNode({2, {Particles::ParticleNode{} }});

create_particle_rhs.addNode({3, {Particles::ParticleNode{} }});

create_particle_rhs.addEdge(2, 3);

create_particle_rhs.addNode({4, {Particles::ParticleNode{} }});

create_particle_rhs.addEdge(3, 4);

create_particle_rhs.addNode({5, {Particles::ParticleNode{} }});

create_particle_rhs.addEdge(3, 5);

DGGML::WithRule<GT> create_particle("create_particle", create_particle_lhs, create_particle_rhs,
                                 [&](auto &lhs, auto &m) {

int counter = std::get<Particles::ParticleNodeCreator>(lhs[m[1]].data).fflow_1c7d3f;

return DGGML::heaviside(counter,1);},
[&](auto &lhs, auto &rhs, auto &m1, auto &m2) {
int counter = std::get<Particles::ParticleNodeCreator>(lhs[m1[1]].data).fflow_1c7d3f;

double x = std::get<Particles::ParticleNodeCreator>(lhs[m1[1]].data).fflow_8df356[0];

double y = std::get<Particles::ParticleNodeCreator>(lhs[m1[1]].data).fflow_8df356[1];

Particles::ParticleNodeCreator rhs_node_1 = std::get<Particles::ParticleNodeCreator>(rhs[m2[1]].data);
		std::get<Particles::ParticleNodeCreator>(rhs[m2[1]].data).fflow_1c7d3f = counter-1;
		rhs[m2[1]].position[0] = x-0.2;
		std::get<Particles::ParticleNodeCreator>(rhs[m2[1]].data).fflow_8df356[0] = x-0.2;
		rhs[m2[1]].position[1] = y-0.2;
		std::get<Particles::ParticleNodeCreator>(rhs[m2[1]].data).fflow_8df356[1] = y-0.2;
Particles::ParticleNode rhs_node_2 = std::get<Particles::ParticleNode>(rhs[m2[2]].data);
		rhs[m2[2]].position[0] = x+0.3;
		std::get<Particles::ParticleNode>(rhs[m2[2]].data).fflow_3717b6[0] = x+0.3;
		rhs[m2[2]].position[1] = y+0.3;
		std::get<Particles::ParticleNode>(rhs[m2[2]].data).fflow_3717b6[1] = y+0.3;
Particles::ParticleNode rhs_node_5 = std::get<Particles::ParticleNode>(rhs[m2[5]].data);
		rhs[m2[5]].position[0] = x+0.6;
		std::get<Particles::ParticleNode>(rhs[m2[5]].data).fflow_3717b6[0] = x+0.6;
		rhs[m2[5]].position[1] = y+0.6;
		std::get<Particles::ParticleNode>(rhs[m2[5]].data).fflow_3717b6[1] = y+0.6;
Particles::ParticleNode rhs_node_4 = std::get<Particles::ParticleNode>(rhs[m2[4]].data);
		rhs[m2[4]].position[0] = x+0.6;
		std::get<Particles::ParticleNode>(rhs[m2[4]].data).fflow_3717b6[0] = x+0.6;
		rhs[m2[4]].position[1] = y+0.1;
		std::get<Particles::ParticleNode>(rhs[m2[4]].data).fflow_3717b6[1] = y+0.1;
}
);
gamma.addRule(create_particle);
};
}
#endif