#ifndef DGGML_Particles_TYPES_HPP
#define DGGML_Particles_TYPES_HPP
#include "YAGL_Graph.hpp" 
#include "YAGL_Node.hpp" 
#include "SpatialData3D.hpp" 
#include <variant>


namespace Particles {

	struct Type {};

    /**
     * @brief A class representing a boundary node in the graph.
     * Default Boundary Type
     * */
    struct Boundary : Type {
    };

	struct ParticleNodeCreator : Type {

		int a4c6400;
		float ad89e53[2];

		std::variant<int, float*> operator[](std::size_t index) {
			if (index == 0) return a4c6400;
			if (index == 1) return ad89e53;
			throw std::out_of_range("Index out of bounds");
		};

	};

	struct Particle : Type {

		float a35817d[2];

		std::variant<int, float*> operator[](std::size_t index) {
			if (index == 0) return a35817d;
			throw std::out_of_range("Index out of bounds");
		};

};
	using graph_type = YAGL::Graph<std::size_t,
        SpatialNode3D<Type, Boundary, ParticleNodeCreator,
        Particle>
    >;
};
#endif
