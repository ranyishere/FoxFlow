#ifndef DGGML_RULES_STAGE_0_HPP
#define DGGML_RULES_STAGE_0_HPP
#include "types.h"
#include "parameters.h"
#include "functions.h"
namespace NeuralNetwork {
using GT = NeuralNetwork::graph_type;
void build_nn(DGGML::Grammar<NeuralNetwork::graph_type> &gamma,
           NeuralNetwork::graph_type &system_graph,
           Parameters &settings) {

GT build_nn_lhs;
build_nn_lhs.addNode({1, {NeuralNetwork::StartType{} }});

GT build_nn_rhs;
build_nn_rhs.addNode({1, {NeuralNetwork::InputLayer{} }});

DGGML::WithRule<GT> build_nn("build_nn", build_nn_lhs, build_nn_rhs,
[&](auto &lhs, auto &m1) {

return 
(
1.0
)
;},
[&](auto &lhs, auto &rhs, auto &m1, auto &m2) {
int input_id = (0);
rhs[m2[ 1 ]].position[2] = input_id;
std::get<NeuralNetwork::InputLayer>(rhs[m2[ 1 ]].data).InputID = input_id;
torch::Tensor pos = (torch::tensor({1.0, 1.0, 1.0}, torch::kFloat64));
for (int _i = 0; _i < 3 && _i < pos.numel(); _i++) { rhs[m2[1]].position[_i] = pos[_i].template item<double>(); }
std::get<NeuralNetwork::InputLayer>(rhs[m2[ 1 ]].data).Position = pos;
torch::Tensor some_vec = (torch::tensor({1.0, 2.0, 3.0, 4.0}, torch::kFloat64));
torch::Tensor W = (torch::zeros({4, 4}, torch::kFloat64));
torch::Tensor out = (torch::mv(W, some_vec));
torch::Tensor W2 = ((W + torch::rand({4, 4}, torch::kFloat64)));
torch::Tensor W3 = (torch::mm(W, (W2).t()));
torch::Tensor A = (torch::rand({2, 3}, torch::kFloat64));
torch::Tensor B = (torch::rand({3, 4}, torch::kFloat64));
torch::Tensor C = (torch::einsum("ij,jk->ik", {A, B}));
torch::Tensor A3 = (torch::rand({8, 2, 3}, torch::kFloat64));
torch::Tensor B3 = (torch::rand({8, 3, 4}, torch::kFloat64));
torch::Tensor C3 = (torch::einsum("bij,bjk->bik", {A3, B3}));
torch::Tensor E = (A3.permute({1, 0, 2}));
torch::Tensor image_batch = (settings.images.index({torch::indexing::Slice(0, 64)}));
std::get<NeuralNetwork::InputLayer>(rhs[m2[ 1 ]].data).ImageBatch = image_batch;
torch::Tensor W0 = (torch::rand({4, 4}, torch::kFloat64));
torch::Tensor v = (torch::rand({4}, torch::kFloat64));
torch::Tensor dloss_dW = ([&]() {
    W0.requires_grad_(true);
    auto _fflow_ad_tmp = (torch::einsum("i->", {torch::mv(W0, v)}));
    _fflow_ad_tmp.backward();
    return W0.grad();
})();
}
);
gamma.addRule(build_nn);
};
void add_layer(DGGML::Grammar<NeuralNetwork::graph_type> &gamma,
           NeuralNetwork::graph_type &system_graph,
           Parameters &settings) {

GT add_layer_lhs;
add_layer_lhs.addNode({1, {NeuralNetwork::InputLayer{} }});

GT add_layer_rhs;
add_layer_rhs.addNode({1, {NeuralNetwork::InputLayer{} }});

add_layer_rhs.addNode({2, {NeuralNetwork::Layer{} }});

add_layer_rhs.addEdge(1, 2);

DGGML::WithRule<GT> add_layer("add_layer", add_layer_lhs, add_layer_rhs,
[&](auto &lhs, auto &m1) {

int input_id = std::get<NeuralNetwork::InputLayer>(lhs[m1[1]].data).InputID;

return 
(
1.0
 * 
( (input_id == 0) ? 1.0 : 0.0)
)
;},
[&](auto &lhs, auto &rhs, auto &m1, auto &m2) {
std::get<NeuralNetwork::InputLayer>(rhs[m2[ 1 ]].data).Position = std::get<NeuralNetwork::InputLayer>(lhs[m1[ 1 ]].data).Position;

std::copy(std::begin(lhs[m1[1]].position), std::end(lhs[m1[1]].position), std::begin(rhs[m2[1]].position));

std::get<NeuralNetwork::InputLayer>(rhs[m2[ 1 ]].data).ImageBatch = std::get<NeuralNetwork::InputLayer>(lhs[m1[ 1 ]].data).ImageBatch;

int input_id = std::get<NeuralNetwork::InputLayer>(lhs[m1[1]].data).InputID;

int new_input_id = (input_id + 1);
rhs[m2[ 1 ]].position[2] = new_input_id;
std::get<NeuralNetwork::InputLayer>(rhs[m2[ 1 ]].data).InputID = new_input_id;
torch::Tensor weights = (torch::rand({100}, torch::kFloat64));
std::get<NeuralNetwork::Layer>(rhs[m2[ 2 ]].data).Weights = weights;
int layer_id = (0);
rhs[m2[ 2 ]].position[2] = layer_id;
std::get<NeuralNetwork::Layer>(rhs[m2[ 2 ]].data).LayerID = layer_id;
torch::Tensor layer_pos = (torch::tensor({1.0, 1.0, 1.5}, torch::kFloat64));
for (int _i = 0; _i < 3 && _i < layer_pos.numel(); _i++) { rhs[m2[2]].position[_i] = layer_pos[_i].template item<double>(); }
std::get<NeuralNetwork::Layer>(rhs[m2[ 2 ]].data).Position = layer_pos;
}
);
gamma.addRule(add_layer);
};
void add_output_layer(DGGML::Grammar<NeuralNetwork::graph_type> &gamma,
           NeuralNetwork::graph_type &system_graph,
           Parameters &settings) {

GT add_output_layer_lhs;
add_output_layer_lhs.addNode({1, {NeuralNetwork::InputLayer{} }});

add_output_layer_lhs.addNode({2, {NeuralNetwork::Layer{} }});

add_output_layer_lhs.addEdge(1, 2);

GT add_output_layer_rhs;
add_output_layer_rhs.addNode({1, {NeuralNetwork::InputLayer{} }});

add_output_layer_rhs.addNode({2, {NeuralNetwork::Layer{} }});

add_output_layer_rhs.addEdge(1, 2);

add_output_layer_rhs.addNode({3, {NeuralNetwork::OutputLayer{} }});

add_output_layer_rhs.addEdge(2, 3);

DGGML::WithRule<GT> add_output_layer("add_output_layer", add_output_layer_lhs, add_output_layer_rhs,
[&](auto &lhs, auto &m1) {

int layer_id = std::get<NeuralNetwork::Layer>(lhs[m1[2]].data).LayerID;

return 
(
1.0
 * 
( (layer_id == 0) ? 1.0 : 0.0)
)
;},
[&](auto &lhs, auto &rhs, auto &m1, auto &m2) {
std::get<NeuralNetwork::InputLayer>(rhs[m2[ 1 ]].data).Position = std::get<NeuralNetwork::InputLayer>(lhs[m1[ 1 ]].data).Position;

std::copy(std::begin(lhs[m1[1]].position), std::end(lhs[m1[1]].position), std::begin(rhs[m2[1]].position));

std::get<NeuralNetwork::InputLayer>(rhs[m2[ 1 ]].data).ImageBatch = std::get<NeuralNetwork::InputLayer>(lhs[m1[ 1 ]].data).ImageBatch;

std::get<NeuralNetwork::InputLayer>(rhs[m2[ 1 ]].data).InputID = std::get<NeuralNetwork::InputLayer>(lhs[m1[ 1 ]].data).InputID;

rhs[m2[ 1 ]].position[2] = lhs[m1[ 1 ]].position[2];

std::get<NeuralNetwork::Layer>(rhs[m2[ 2 ]].data).Position = std::get<NeuralNetwork::Layer>(lhs[m1[ 2 ]].data).Position;

std::copy(std::begin(lhs[m1[2]].position), std::end(lhs[m1[2]].position), std::begin(rhs[m2[2]].position));

std::get<NeuralNetwork::Layer>(rhs[m2[ 2 ]].data).Weights = std::get<NeuralNetwork::Layer>(lhs[m1[ 2 ]].data).Weights;

int layer_id = std::get<NeuralNetwork::Layer>(lhs[m1[2]].data).LayerID;

int new_layer_id = (layer_id + 1);
rhs[m2[ 2 ]].position[2] = new_layer_id;
std::get<NeuralNetwork::Layer>(rhs[m2[ 2 ]].data).LayerID = new_layer_id;
int output_id = (0);
rhs[m2[ 3 ]].position[2] = output_id;
std::get<NeuralNetwork::OutputLayer>(rhs[m2[ 3 ]].data).OutputID = output_id;
torch::Tensor output_pos = (torch::tensor({1.0, 1.0, 2.0}, torch::kFloat64));
for (int _i = 0; _i < 3 && _i < output_pos.numel(); _i++) { rhs[m2[3]].position[_i] = output_pos[_i].template item<double>(); }
std::get<NeuralNetwork::OutputLayer>(rhs[m2[ 3 ]].data).Position = output_pos;
torch::Tensor digit_class = (torch::zeros({10}, torch::kFloat64));
std::get<NeuralNetwork::OutputLayer>(rhs[m2[ 3 ]].data).DigitClass = digit_class;
}
);
gamma.addRule(add_output_layer);
};
}
#endif