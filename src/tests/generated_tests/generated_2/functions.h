#ifndef DGGML_FUNCTIONS_HPP
#define DGGML_FUNCTIONS_HPP
#include<cmath>
 namespace FractureNetwork {

double distance(double x0, double y0, double x1, double y1){double res = (sqrt((std::pow((x1 - x0), 2) + std::pow((y1 - y0), 2))));return res;}
;


}
#endif