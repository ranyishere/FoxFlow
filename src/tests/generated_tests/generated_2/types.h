#ifndef DGGML_Dissolution_TYPES_HPP
#define DGGML_Dissolution_TYPES_HPP
#include "YAGL_Graph.hpp" 
#include "YAGL_Node.hpp" 
#include "SpatialData3D.hpp" 
#include "torch/torch.h"
namespace Dissolution {
struct Type {
	template <class Archive>
	void serialize(Archive& archive) {
	}
};
struct StartType {
	torch::Tensor start_location = torch::zeros(3, torch::kFloat64);
	template <class Archive>
	void serialize(Archive& archive) {
		archive(start_location);
	}
};
struct Boundary {
	torch::Tensor boundary_location = torch::zeros(3, torch::kFloat64);
	template <class Archive>
	void serialize(Archive& archive) {
		archive(boundary_location);
	}
};
	struct Fluid : Type {
		torch::Tensor Position = torch::zeros({3}, torch::kFloat64);
		torch::Tensor Unit = torch::zeros({3}, torch::kFloat64);
		torch::Tensor FCount = torch::zeros({1}, torch::kInt64);
		double Pressure;
		double Concentration;

	template <class Archive>
	void serialize(Archive& archive) {
		archive(Position);
		archive(Unit);
		archive(FCount);
		archive(Pressure);
		archive(Concentration);
	}

};	struct FluidSource : Type {
		torch::Tensor Position = torch::zeros({3}, torch::kFloat64);
		torch::Tensor Unit = torch::zeros({3}, torch::kFloat64);
		torch::Tensor FCount = torch::zeros({1}, torch::kInt64);
		double Pressure;
		double Concentration;

	template <class Archive>
	void serialize(Archive& archive) {
		archive(Position);
		archive(Unit);
		archive(FCount);
		archive(Pressure);
		archive(Concentration);
	}

};	struct FluidSink : Type {
		torch::Tensor Position = torch::zeros({3}, torch::kFloat64);
		torch::Tensor Unit = torch::zeros({3}, torch::kFloat64);
		torch::Tensor FCount = torch::zeros({1}, torch::kInt64);
		double Pressure;
		double Concentration;

	template <class Archive>
	void serialize(Archive& archive) {
		archive(Position);
		archive(Unit);
		archive(FCount);
		archive(Pressure);
		archive(Concentration);
	}

};	struct RockStart : Type {
		torch::Tensor Position = torch::zeros({3}, torch::kFloat64);
		torch::Tensor Count = torch::zeros({3}, torch::kInt64);
		torch::Tensor Dir = torch::zeros({3}, torch::kInt64);
		double Density;

	template <class Archive>
	void serialize(Archive& archive) {
		archive(Position);
		archive(Count);
		archive(Dir);
		archive(Density);
	}

};	using graph_type = YAGL::Graph<std::size_t,	SpatialNode3D<StartType,Boundary,Fluid,FluidSource,FluidSink,RockStart>>;
};
#endif