#ifndef DGGML_RULES_HPP
#define DGGML_RULES_HPP
#include "types.h"
#include "parameters.h"
namespace particle_rules {
using GT = Particles::graph_type;
void transform(DGGML::Grammar<Particles::graph_type> &gamma,
           Particles::graph_type &system_graph,
           Parameters &settings) {

GT transform_lhs;
transform_lhs.addNode({1, {Particles::StartType{} }});

GT transform_rhs;
transform_rhs.addNode({1, {Particles::ParabolaParticle{} }});

DGGML::WithRule<GT> transform("transform", transform_lhs, transform_rhs,
[&](auto &lhs, auto &m) {

return DGGML::heaviside(10,1);},
[&](auto &lhs, auto &rhs, auto &m1, auto &m2) {
float x = std::get<Particles::StartType>(lhs[m1[1]].data).start_location[0];

float y = std::get<Particles::StartType>(lhs[m1[1]].data).start_location[1];

    Particles::ParabolaParticle rhs_node_1 = std::get<Particles::ParabolaParticle>(rhs[m2[1]].data);
        rhs[m2[1]].position[0] = x;
        std::get<Particles::ParabolaParticle>(rhs[m2[1]].data).fflow_bde751[0] = x;
        rhs[m2[1]].position[1] = y;
        std::get<Particles::ParabolaParticle>(rhs[m2[1]].data).fflow_bde751[1] = y;
        std::get<Particles::ParabolaParticle>(rhs[m2[1]].data).fflow_2121c5[0] = 0.5;
        std::get<Particles::ParabolaParticle>(rhs[m2[1]].data).fflow_2121c5[1] = 0.5;
    }

);
    gamma.addRule(transform);
};
void move_node(DGGML::Grammar<Particles::graph_type> &gamma,
           Particles::graph_type &system_graph,
           Parameters &settings) {

    GT move_node_lhs;
    move_node_lhs.addNode({1, {Particles::ParabolaParticle{} }});

    GT move_node_rhs;
    move_node_rhs.addNode({1, {Particles::ParabolaParticle{} }});

    DGGML::SolvingRule<GT> move_node("move_node", move_node_lhs, move_node_lhs,

4,
    [](auto &lhs, auto &m1, auto &varset) {

        varset.insert(&lhs[m1[1]].position[0]);
        varset.insert(&lhs[m1[1]].position[1]);

        auto &data0 = std::get<Particles::ParabolaParticle>(
                lhs[m1[1]].data
            ).fflow_2121c5[0];

        auto &data1 = std::get<Particles::ParabolaParticle>(
                lhs[m1[1]].data
            ).fflow_2121c5[1];

        varset.insert(&data0);
        varset.insert(&data1);

    },
[&](auto &lhs, auto &m1, auto y, auto ydot, auto &varmap) {

    auto &oof = std::get<Particles::ParabolaParticle>(lhs[m1[1]].data).fflow_2121c5[0];
    auto &oof1 = std::get<Particles::ParabolaParticle>(lhs[m1[1]].data).fflow_2121c5[1];

    // auto ix_v = varmap.find(&oof)->second;
    // auto ix_w = varmap.find(&oof1)->second;

    auto ix_v = varmap.at(&oof);
    auto ix_w = varmap.at(&oof1);
    // std::cout << "ix_v: " << ix_v << std::endl;

    auto ix_x = varmap.at(&lhs[m1[1]].position[0]);
    auto ix_y = varmap.at(&lhs[m1[1]].position[1]);

    // Moves in a circle
    NV_Ith_S(ydot, varmap[&lhs[m1[1]].position[0]]) +=  0.1 * NV_Ith_S(y, ix_v) * (
        NV_Ith_S(y, ix_y)  - 0.5);

    NV_Ith_S(ydot, varmap[&lhs[m1[1]].position[1]]) +=  0.1 * - NV_Ith_S(y, ix_w) * ( 
        NV_Ith_S(y, ix_x) - 0.5);

    NV_Ith_S(ydot, varmap[&oof]) +=  0;
    NV_Ith_S(ydot, varmap[&oof1]) +=  0.6;

    }
);
    gamma.addRule(move_node);
};
}
#endif
