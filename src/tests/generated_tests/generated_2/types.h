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
		float fflow_a55d52[2];
		int fflow_658e14;
		void* operator[](std::size_t index) const {
			if (index == 0) return (void*)&fflow_a55d52;
			if (index == 1) return (void*)&fflow_658e14;
			throw std::out_of_range("Index out of bounds");
		};

};	struct ParticleNode : Type {
		float fflow_bc174a[2];
		void* operator[](std::size_t index) const {
			if (index == 0) return (void*)&fflow_bc174a;
			throw std::out_of_range("Index out of bounds");
		};

};	using graph_type = YAGL::Graph<std::size_t,	SpatialNode3D<StartType,Boundary,ParticleNodeCreator,ParticleNode>>;
};
#endif