#ifndef DGGML_Particles_TYPES_HPP
#define DGGML_Particles_TYPES_HPP
#include "YAGL_Graph.hpp" 
#include "YAGL_Node.hpp" 
#include "SpatialData3D.hpp" 
namespace Particles {
	struct Type {};
	struct Boundary {};
	struct StartType {
            float start_location[2];

       };
	struct ParticleNodeCreator : Type {
		float fflow_8df356[2];
		int fflow_1c7d3f;
		void* operator[](std::size_t index) const {
			if (index == 0) return (void*)&fflow_8df356;
			if (index == 1) return (void*)&fflow_1c7d3f;
			throw std::out_of_range("Index out of bounds");
		};

};	struct ParticleNode : Type {
		float fflow_3717b6[2];
		void* operator[](std::size_t index) const {
			if (index == 0) return (void*)&fflow_3717b6;
			throw std::out_of_range("Index out of bounds");
		};

};	using graph_type = YAGL::Graph<std::size_t,	SpatialNode3D<StartType,Boundary,ParticleNodeCreator,ParticleNode>>;
};
#endif