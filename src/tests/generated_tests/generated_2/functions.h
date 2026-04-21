#ifndef DGGML_FUNCTIONS_HELP_HPP
#define DGGML_FUNCTIONS_HELP_HPP
#include<cmath>
#include <torch/script.h>
 namespace HELP {

double distance(double x0, double y0, double z0, double x1, double y1, double z1){double res = (sqrt(((std::pow((x1 - x0), 2) + std::pow((y1 - y0), 2)) + std::pow((z1 - z0), 2))));return res;}
;


double line_intersection_x(double x1, double y1, double ux1, double uy1, double x2, double y2, double x3, double y3){double varepsilon = (0.001);double dx = (x3 - x2);double dy = (y3 - y2);double D = (dx * uy1 - dy * ux1);double res = (((x2 - x1) * (-dy) - (-dx) * (y2 - y1)) / D);return res;}
;


double line_intersection_y(double x1, double y1, double ux1, double uy1, double x2, double y2, double x3, double y3){double varepsilon = (0.001);double dx = (x3 - x2);double dy = (y3 - y2);double D = (dx * uy1 - dy * ux1);double res = ((ux1 * (y2 - y1) - uy1 * (x2 - x1)) / D);return res;}
;


double rotate_x(double theta, double x){return x * cos(theta) + x * sin(theta);}
;


double rotate_y(double theta, double x){return -x * sin(theta) + x * cos(theta);}
;


double minimum_distance_2d(double x1, double y1, double x2, double y2, double x3, double y3){double num = (abs((((x3 - x2) * (y1 - y2)) - ((y3 - y2) * (x1 - x2)))));double denom = (sqrt((std::pow((x3 - x2), 2) + std::pow((y3 - y2), 2))));return (num) / (denom + 0.0000001);}
;


double calculate_angle(double ux_1, double uy_1, double x_1, double y_1, double ux_2, double uy_2){double num = ((ux_1 * ux_2) + (uy_1 * uy_2));double denom = (sqrt((std::pow(ux_1, 2) + std::pow(uy_1, 2))) * sqrt((std::pow(ux_2, 2) + std::pow(uy_2, 2))));double deg_conversion = ((180.0) / (3.1415926535));return acos((num / denom)) * deg_conversion;}
;


double darcy_velocity(double P0, double P1){double res = (3.33 * (P0 - P1));return res;}
;


}
#endif