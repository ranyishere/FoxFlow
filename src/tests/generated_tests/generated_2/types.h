#ifndef DGGML_NeuralNetwork_TYPES_HPP
#define DGGML_NeuralNetwork_TYPES_HPP
#include "YAGL_Graph.hpp" 
#include "YAGL_Node.hpp" 
#include "SpatialData3D.hpp" 
#include "torch/torch.h"
namespace NeuralNetwork {
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
	struct InputLayer : Type {
		torch::Tensor Position = torch::zeros({3}, torch::kFloat64);
		torch::Tensor ImageBatch = torch::zeros({64, 1, 28, 28}, torch::kFloat64);
		int InputID;

	template <class Archive>
	void serialize(Archive& archive) {
		archive(Position);
		archive(ImageBatch);
		archive(InputID);
	}

};	struct Layer : Type {
		torch::Tensor Position = torch::zeros({3}, torch::kFloat64);
		torch::Tensor Weights = torch::zeros({100}, torch::kFloat64);
		int LayerID;

	template <class Archive>
	void serialize(Archive& archive) {
		archive(Position);
		archive(Weights);
		archive(LayerID);
	}

};	struct OutputLayer : Type {
		torch::Tensor Position = torch::zeros({3}, torch::kFloat64);
		torch::Tensor DigitClass = torch::zeros({10}, torch::kFloat64);
		int OutputID;

	template <class Archive>
	void serialize(Archive& archive) {
		archive(Position);
		archive(DigitClass);
		archive(OutputID);
	}

};	using graph_type = YAGL::Graph<std::size_t,	SpatialNode3D<StartType,Boundary,InputLayer,Layer,OutputLayer>>;
};
#endif