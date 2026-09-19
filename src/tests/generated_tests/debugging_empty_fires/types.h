#ifndef DGGML_Microtubule_TYPES_HPP
#define DGGML_Microtubule_TYPES_HPP
#include "YAGL_Graph.hpp" 
#include "YAGL_Node.hpp" 
#include "SpatialData3D.hpp" 
namespace Microtubule {
	struct Type {};
	struct Boundary {
            float boundary_location[2];

        };
	struct StartType {
            float start_location[2];

       };
	struct FractureSegment : Type {
		double fflow_bb7745[2];
		double fflow_87a940;
		double fflow_b62e5f[2];
		void* operator[](std::size_t index) const {
			if (index == 0) return (void*)&fflow_bb7745;
			if (index == 1) return (void*)&fflow_87a940;
			if (index == 2) return (void*)&fflow_b62e5f;
			throw std::out_of_range("Index out of bounds");
		};

};	struct FractureSegmentEnd : Type {
		double fflow_f9bbc8[2];
		double fflow_1e83fe;
		double fflow_f13ca0[2];
		void* operator[](std::size_t index) const {
			if (index == 0) return (void*)&fflow_f9bbc8;
			if (index == 1) return (void*)&fflow_1e83fe;
			if (index == 2) return (void*)&fflow_f13ca0;
			throw std::out_of_range("Index out of bounds");
		};

};	using graph_type = YAGL::Graph<std::size_t,	SpatialNode3D<StartType,Boundary,FractureSegment,FractureSegmentEnd>>;
};
#endif