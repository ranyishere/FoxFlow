#ifndef DGGML_RULES_HPP
#define DGGML_RULES_HPP
#include "types.h"
#include "parameters.h"

namespace particle_rules {

    using GT = Particles::graph_type;

    void change_to_particle(DGGML::Grammar<Particles::graph_type> &gamma, 
                Particles::graph_type &system_graph, Parameters &settings) {
                GT change_to_particle_lhs;
                change_to_particle_lhs.addNode({1, {Particles::ParticleNodeCreator{} }});

                // GT change_to_particle_rhs;
                GT change_to_particle_rhs;
                change_to_particle_rhs.addNode({1, {Particles::ParticleNodeCreator{} }});

                // change_to_particle_rhs.addNode({2, {Particles::Particle{} }});

                DGGML::WithRule<GT> rule(
                    "change_to_particle", change_to_particle_lhs, change_to_particle_rhs,
                [&](auto &lhs, auto &m) {
                    auto &node = lhs.findNode(m[1])->second.getData();
                    auto oof = std::get<Particles::ParticleNodeCreator>(node.data).a4c6400;
		    std::cout << "change_to_particle: " << oof << std::endl;
		    // return 9999999;
		    return 1;
                    // return DGGML::heaviside(oof,  0);
                },
                [&](auto &lhs, auto &rhs, auto &m1, auto &m2) {

                    int counter = std::get<Particles::ParticleNodeCreator>(lhs[m1[1]].data).a4c6400;
                    int new_counter = counter-1;
                    int x = counter+1;
                    int y = counter+1;

                    Particles::ParticleNodeCreator tmp = std::get<Particles::ParticleNodeCreator>(rhs[m2[1]].data);
                    tmp.a4c6400 = new_counter;
                    tmp.ad89e53[0] = x;
                    tmp.ad89e53[1] = y;

                    rhs[m2[1]].position[0] = lhs[m1[1]].position[0] + x*0.3;
                    rhs[m2[1]].position[1] = lhs[m1[1]].position[1] + y*0.3;

                    /*
                    std::get<Particles::Particle>(rhs[m2[2]].data).a35817d[0] = x+1;
                    std::get<Particles::Particle>(rhs[m2[2]].data).a35817d[1] = y+1;
                    */

            });

        gamma.addRule(rule);
    };
};

#endif
