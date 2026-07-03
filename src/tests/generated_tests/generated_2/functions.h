#ifndef DGGML_FUNCTIONS_HELP_HPP
#define DGGML_FUNCTIONS_HELP_HPP
#include<cmath>
#include <torch/script.h>
 namespace HELP {

torch::Tensor activation(int M, torch::Tensor x, torch::Tensor w){return torch::einsum("i->", (x * w));}
;


}
#endif