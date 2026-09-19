#ifndef DGGML_RULES_STAGE_1_HPP
#define DGGML_RULES_STAGE_1_HPP
#include "types.h"
#include "parameters.h"
#include "functions.h"
namespace NeuralNetwork {
using GT = NeuralNetwork::graph_type;
void backprop(DGGML::Grammar<NeuralNetwork::graph_type> &gamma,
           NeuralNetwork::graph_type &system_graph,
           Parameters &settings) {

GT backprop_lhs;
backprop_lhs.addNode({1, {NeuralNetwork::InputLayer{} }});

backprop_lhs.addNode({2, {NeuralNetwork::Layer{} }});

backprop_lhs.addEdge(1, 2);

backprop_lhs.addNode({3, {NeuralNetwork::OutputLayer{} }});

backprop_lhs.addEdge(2, 3);

GT backprop_rhs;
backprop_rhs.addNode({1, {NeuralNetwork::InputLayer{} }});

backprop_rhs.addNode({2, {NeuralNetwork::Layer{} }});

backprop_rhs.addEdge(1, 2);

backprop_rhs.addNode({3, {NeuralNetwork::OutputLayer{} }});

backprop_rhs.addEdge(2, 3);

DGGML::SolvingRule<GT> backprop("backprop", backprop_lhs, backprop_lhs,
103,
[](auto &lhs, auto &m1, auto &varset) {
auto &tensor_ref_2_2 = std::get<NeuralNetwork::Layer>(lhs[m1[2]].data).Weights;
double* tensor_ptr_2_2 = tensor_ref_2_2.template data_ptr<double>();
varset.insert(&tensor_ptr_2_2[0]);
varset.insert(&tensor_ptr_2_2[1]);
varset.insert(&tensor_ptr_2_2[2]);
varset.insert(&tensor_ptr_2_2[3]);
varset.insert(&tensor_ptr_2_2[4]);
varset.insert(&tensor_ptr_2_2[5]);
varset.insert(&tensor_ptr_2_2[6]);
varset.insert(&tensor_ptr_2_2[7]);
varset.insert(&tensor_ptr_2_2[8]);
varset.insert(&tensor_ptr_2_2[9]);
varset.insert(&tensor_ptr_2_2[10]);
varset.insert(&tensor_ptr_2_2[11]);
varset.insert(&tensor_ptr_2_2[12]);
varset.insert(&tensor_ptr_2_2[13]);
varset.insert(&tensor_ptr_2_2[14]);
varset.insert(&tensor_ptr_2_2[15]);
varset.insert(&tensor_ptr_2_2[16]);
varset.insert(&tensor_ptr_2_2[17]);
varset.insert(&tensor_ptr_2_2[18]);
varset.insert(&tensor_ptr_2_2[19]);
varset.insert(&tensor_ptr_2_2[20]);
varset.insert(&tensor_ptr_2_2[21]);
varset.insert(&tensor_ptr_2_2[22]);
varset.insert(&tensor_ptr_2_2[23]);
varset.insert(&tensor_ptr_2_2[24]);
varset.insert(&tensor_ptr_2_2[25]);
varset.insert(&tensor_ptr_2_2[26]);
varset.insert(&tensor_ptr_2_2[27]);
varset.insert(&tensor_ptr_2_2[28]);
varset.insert(&tensor_ptr_2_2[29]);
varset.insert(&tensor_ptr_2_2[30]);
varset.insert(&tensor_ptr_2_2[31]);
varset.insert(&tensor_ptr_2_2[32]);
varset.insert(&tensor_ptr_2_2[33]);
varset.insert(&tensor_ptr_2_2[34]);
varset.insert(&tensor_ptr_2_2[35]);
varset.insert(&tensor_ptr_2_2[36]);
varset.insert(&tensor_ptr_2_2[37]);
varset.insert(&tensor_ptr_2_2[38]);
varset.insert(&tensor_ptr_2_2[39]);
varset.insert(&tensor_ptr_2_2[40]);
varset.insert(&tensor_ptr_2_2[41]);
varset.insert(&tensor_ptr_2_2[42]);
varset.insert(&tensor_ptr_2_2[43]);
varset.insert(&tensor_ptr_2_2[44]);
varset.insert(&tensor_ptr_2_2[45]);
varset.insert(&tensor_ptr_2_2[46]);
varset.insert(&tensor_ptr_2_2[47]);
varset.insert(&tensor_ptr_2_2[48]);
varset.insert(&tensor_ptr_2_2[49]);
varset.insert(&tensor_ptr_2_2[50]);
varset.insert(&tensor_ptr_2_2[51]);
varset.insert(&tensor_ptr_2_2[52]);
varset.insert(&tensor_ptr_2_2[53]);
varset.insert(&tensor_ptr_2_2[54]);
varset.insert(&tensor_ptr_2_2[55]);
varset.insert(&tensor_ptr_2_2[56]);
varset.insert(&tensor_ptr_2_2[57]);
varset.insert(&tensor_ptr_2_2[58]);
varset.insert(&tensor_ptr_2_2[59]);
varset.insert(&tensor_ptr_2_2[60]);
varset.insert(&tensor_ptr_2_2[61]);
varset.insert(&tensor_ptr_2_2[62]);
varset.insert(&tensor_ptr_2_2[63]);
varset.insert(&tensor_ptr_2_2[64]);
varset.insert(&tensor_ptr_2_2[65]);
varset.insert(&tensor_ptr_2_2[66]);
varset.insert(&tensor_ptr_2_2[67]);
varset.insert(&tensor_ptr_2_2[68]);
varset.insert(&tensor_ptr_2_2[69]);
varset.insert(&tensor_ptr_2_2[70]);
varset.insert(&tensor_ptr_2_2[71]);
varset.insert(&tensor_ptr_2_2[72]);
varset.insert(&tensor_ptr_2_2[73]);
varset.insert(&tensor_ptr_2_2[74]);
varset.insert(&tensor_ptr_2_2[75]);
varset.insert(&tensor_ptr_2_2[76]);
varset.insert(&tensor_ptr_2_2[77]);
varset.insert(&tensor_ptr_2_2[78]);
varset.insert(&tensor_ptr_2_2[79]);
varset.insert(&tensor_ptr_2_2[80]);
varset.insert(&tensor_ptr_2_2[81]);
varset.insert(&tensor_ptr_2_2[82]);
varset.insert(&tensor_ptr_2_2[83]);
varset.insert(&tensor_ptr_2_2[84]);
varset.insert(&tensor_ptr_2_2[85]);
varset.insert(&tensor_ptr_2_2[86]);
varset.insert(&tensor_ptr_2_2[87]);
varset.insert(&tensor_ptr_2_2[88]);
varset.insert(&tensor_ptr_2_2[89]);
varset.insert(&tensor_ptr_2_2[90]);
varset.insert(&tensor_ptr_2_2[91]);
varset.insert(&tensor_ptr_2_2[92]);
varset.insert(&tensor_ptr_2_2[93]);
varset.insert(&tensor_ptr_2_2[94]);
varset.insert(&tensor_ptr_2_2[95]);
varset.insert(&tensor_ptr_2_2[96]);
varset.insert(&tensor_ptr_2_2[97]);
varset.insert(&tensor_ptr_2_2[98]);
varset.insert(&tensor_ptr_2_2[99]);
},
[&](auto &lhs, auto &m1, auto y, auto ydot, auto &varmap) {
double* tensor_ptr_2_2 = std::get<NeuralNetwork::Layer>(lhs[m1[2]].data).Weights.template data_ptr<double>();
{  // Neural ODE: dw
const int _base_dw = varmap.at(&tensor_ptr_2_2[0]);
// Zero-copy views into SUNDIALS y / ydot memory
torch::Tensor _y_dw = torch::from_blob(
    N_VGetArrayPointer(y) + _base_dw, {100}, torch::kFloat64);
torch::Tensor _ydot_dw = torch::from_blob(
    N_VGetArrayPointer(ydot) + _base_dw, {100}, torch::kFloat64);
auto& Weights = _y_dw;
torch::Tensor weights = std::get<NeuralNetwork::Layer>(lhs[m1[2]].data).Weights;

torch::Tensor _rhs_dw = HELP::activation(100, weights, weights);
_ydot_dw += _rhs_dw;
}  // end Neural ODE: dw
// ── Symbolic ODE system ──
// d(dw[0..99])/dt += HELP::activation(100, weights, weights)
{ static bool _sym_dumped = false;
  if (!_sym_dumped && std::getenv("ODE_DUMP")) { _sym_dumped = true;
    std::cout << "  ── Symbolic ODE ──" << std::endl;
    std::cout << "    d(dw[0..99])/dt += HELP::activation(100, weights, weights)" << std::endl;
  }
}
}
);
gamma.addRule(backprop);
};
}
#endif