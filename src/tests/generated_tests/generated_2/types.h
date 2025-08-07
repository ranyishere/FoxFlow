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
		float fflow_261168[2];
		int fflow_936c91;
		void* operator[](std::size_t index) const {
			if (index == 0) return (void*)&fflow_261168;
			if (index == 1) return (void*)&fflow_936c91;
			throw std::out_of_range("Index out of bounds");
		};

};	struct ParticleNode : Type {
		float fflow_87fbe5[2];
		void* operator[](std::size_t index) const {
			if (index == 0) return (void*)&fflow_87fbe5;
			throw std::out_of_range("Index out of bounds");
		};

};	using graph_type = YAGL::Graph<std::size_t,	SpatialNode3D<StartType,Boundary,ParticleNodeCreator,ParticleNode>>;
};
#endif