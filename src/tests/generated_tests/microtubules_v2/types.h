#ifndef DGGML_Microtubule_TYPES_HPP
#define DGGML_Microtubule_TYPES_HPP
#include "YAGL_Graph.hpp" 
#include "YAGL_Node.hpp" 
#include "SpatialData3D.hpp" 
#include "torch/torch.h"
namespace Microtubule {
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
	struct Intermediate : Type {
		torch::Tensor Direction = torch::zeros({3}, torch::kFloat64);

	template <class Archive>
	void serialize(Archive& archive) {
		archive(Direction);
	}

};	struct Positive : Type {
		torch::Tensor Direction = torch::zeros({3}, torch::kFloat64);

	template <class Archive>
	void serialize(Archive& archive) {
		archive(Direction);
	}

};	struct Negative : Type {
		torch::Tensor Direction = torch::zeros({3}, torch::kFloat64);

	template <class Archive>
	void serialize(Archive& archive) {
		archive(Direction);
	}

};	struct Nucleator : Type {
		torch::Tensor Direction = torch::zeros({2}, torch::kInt64);

	template <class Archive>
	void serialize(Archive& archive) {
		archive(Direction);
	}

};	struct CellBoundary : Type {
		torch::Tensor Direction = torch::zeros({3}, torch::kFloat64);
		torch::Tensor Count = torch::zeros({2}, torch::kInt64);

	template <class Archive>
	void serialize(Archive& archive) {
		archive(Direction);
		archive(Count);
	}

};	struct Zipper : Type {
		torch::Tensor Direction = torch::zeros({3}, torch::kFloat64);

	template <class Archive>
	void serialize(Archive& archive) {
		archive(Direction);
	}

};	struct Junction : Type {
		torch::Tensor Direction = torch::zeros({3}, torch::kFloat64);

	template <class Archive>
	void serialize(Archive& archive) {
		archive(Direction);
	}

};	using graph_type = YAGL::Graph<std::size_t,	SpatialNode3D<StartType,Boundary,Intermediate,Positive,Negative,Nucleator,CellBoundary,Zipper,Junction>>;
};
#endif