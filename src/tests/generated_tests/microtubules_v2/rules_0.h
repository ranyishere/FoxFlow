#ifndef DGGML_RULES_STAGE_0_HPP
#define DGGML_RULES_STAGE_0_HPP
#include "types.h"
#include "parameters.h"
#include "functions.h"
namespace Microtubule {
using GT = Microtubule::graph_type;
void start_to_node(DGGML::Grammar<Microtubule::graph_type> &gamma,
           Microtubule::graph_type &system_graph,
           Parameters &settings) {

GT start_to_node_lhs;
start_to_node_lhs.addNode({1, {Microtubule::StartType{} }});

GT start_to_node_rhs;
start_to_node_rhs.addNode({1, {Microtubule::Nucleator{} }});

start_to_node_rhs.addNode({2, {Microtubule::CellBoundary{} }});

DGGML::WithRule<GT> start_to_node("start_to_node", start_to_node_lhs, start_to_node_rhs,
[&](auto &lhs, auto &m1) {

return 
(
DGGML::heaviside(10, 1)
)
;},
[&](auto &lhs, auto &rhs, auto &m1, auto &m2) {
rhs[m2[1]].position[(0)] = static_cast<double>((0 + settings.buffer));
rhs[m2[1]].position[(1)] = static_cast<double>((0 + settings.buffer));
rhs[m2[1]].position[(2)] = static_cast<double>((0));
double _arr_tmp_0 = settings.nuc_count_n;
double _arr_tmp_1 = settings.nuc_count_n;
torch::Tensor nuc_count = (torch::tensor({_arr_tmp_0, _arr_tmp_1}, torch::kFloat64));
std::get<Microtubule::Nucleator>(rhs[m2[ 1 ]].data).Direction = (nuc_count).clone();
double _arr_tmp_2 = settings.boundary_margin;
double _arr_tmp_3 = settings.boundary_margin;
double _arr_tmp_4 = 0.0;
torch::Tensor b_pos = (torch::tensor({_arr_tmp_2, _arr_tmp_3, _arr_tmp_4}, torch::kFloat64));
for (int _i = 0; _i < 3 && _i < b_pos.numel(); _i++) { rhs[m2[2]].position[_i] = b_pos[_i].template item<double>(); }
torch::Tensor b_unit = (torch::tensor({1.0, 1.0, 1.0}, torch::kFloat64));
std::get<Microtubule::CellBoundary>(rhs[m2[ 2 ]].data).Direction = (b_unit).clone();
torch::Tensor b_count = (torch::tensor({1, 1}, torch::kFloat64));
std::get<Microtubule::CellBoundary>(rhs[m2[ 2 ]].data).Count = (b_count).clone();
}
);
gamma.addRule(start_to_node);
};
void create_nucleator_grid_x(DGGML::Grammar<Microtubule::graph_type> &gamma,
           Microtubule::graph_type &system_graph,
           Parameters &settings) {

GT create_nucleator_grid_x_lhs;
create_nucleator_grid_x_lhs.addNode({1, {Microtubule::Nucleator{} }});

GT create_nucleator_grid_x_rhs;
create_nucleator_grid_x_rhs.addNode({1, {Microtubule::Nucleator{} }});

create_nucleator_grid_x_rhs.addNode({2, {Microtubule::Nucleator{} }});

DGGML::WithRule<GT> create_nucleator_grid_x("create_nucleator_grid_x", create_nucleator_grid_x_lhs, create_nucleator_grid_x_rhs,
[&](auto &lhs, auto &m1) {

torch::Tensor nuc_count = std::get<Microtubule::Nucleator>(lhs[m1[1]].data).Direction.clone();

return 
(
( (0 != nuc_count[0].template item<double>()) ? 1.0 : 0.0)
)
;},
[&](auto &lhs, auto &rhs, auto &m1, auto &m2) {
std::copy(std::begin(lhs[m1[1]].position), std::end(lhs[m1[1]].position), std::begin(rhs[m2[1]].position));

std::get<Microtubule::Nucleator>(rhs[m2[ 1 ]].data).Direction = std::get<Microtubule::Nucleator>(lhs[m1[ 1 ]].data).Direction.clone();

torch::Tensor nuc_pos = torch::tensor({rhs[m2[1]].position[0], rhs[m2[1]].position[1], rhs[m2[1]].position[2]}, torch::kFloat64);

double _arr_tmp_5 = nuc_pos[0].template item<double>() + settings.offset;
double _arr_tmp_6 = nuc_pos[1].template item<double>();
double _arr_tmp_7 = 0;
torch::Tensor nuc_pos_2 = (torch::tensor({_arr_tmp_5, _arr_tmp_6, _arr_tmp_7}, torch::kFloat64));
for (int _i = 0; _i < 3 && _i < nuc_pos_2.numel(); _i++) { rhs[m2[2]].position[_i] = nuc_pos_2[_i].template item<double>(); }
torch::Tensor nuc_count = std::get<Microtubule::Nucleator>(rhs[m2[1]].data).Direction.clone();

double _arr_tmp_8 = nuc_count[0].template item<double>() - 1;
double _arr_tmp_9 = nuc_count[1].template item<double>();
torch::Tensor nuc_count_2 = (torch::tensor({_arr_tmp_8, _arr_tmp_9}, torch::kFloat64));
std::get<Microtubule::Nucleator>(rhs[m2[ 2 ]].data).Direction = (nuc_count_2).clone();
std::get<Microtubule::Nucleator>(rhs[m2[ 1 ]].data).Direction = std::get<Microtubule::Nucleator>(rhs[m2[ 1 ]].data).Direction.clone();
std::get<Microtubule::Nucleator>(rhs[m2[ 1 ]].data).Direction[(0)] = (0);
}
);
gamma.addRule(create_nucleator_grid_x);
};
void create_nucleator_grid_y(DGGML::Grammar<Microtubule::graph_type> &gamma,
           Microtubule::graph_type &system_graph,
           Parameters &settings) {

GT create_nucleator_grid_y_lhs;
create_nucleator_grid_y_lhs.addNode({1, {Microtubule::Nucleator{} }});

GT create_nucleator_grid_y_rhs;
create_nucleator_grid_y_rhs.addNode({1, {Microtubule::Nucleator{} }});

create_nucleator_grid_y_rhs.addNode({2, {Microtubule::Nucleator{} }});

DGGML::WithRule<GT> create_nucleator_grid_y("create_nucleator_grid_y", create_nucleator_grid_y_lhs, create_nucleator_grid_y_rhs,
[&](auto &lhs, auto &m1) {

torch::Tensor nuc_count = std::get<Microtubule::Nucleator>(lhs[m1[1]].data).Direction.clone();

return 
(
( (nuc_count[0].template item<double>() == 0) ? 1.0 : 0.0)
 * 
( (0 != nuc_count[1].template item<double>()) ? 1.0 : 0.0)
)
;},
[&](auto &lhs, auto &rhs, auto &m1, auto &m2) {
std::copy(std::begin(lhs[m1[1]].position), std::end(lhs[m1[1]].position), std::begin(rhs[m2[1]].position));

std::get<Microtubule::Nucleator>(rhs[m2[ 1 ]].data).Direction = std::get<Microtubule::Nucleator>(lhs[m1[ 1 ]].data).Direction.clone();

torch::Tensor nuc_pos = torch::tensor({rhs[m2[1]].position[0], rhs[m2[1]].position[1], rhs[m2[1]].position[2]}, torch::kFloat64);

double _arr_tmp_10 = nuc_pos[0].template item<double>();
double _arr_tmp_11 = nuc_pos[1].template item<double>() + settings.offset;
double _arr_tmp_12 = 0;
torch::Tensor nuc_pos_2 = (torch::tensor({_arr_tmp_10, _arr_tmp_11, _arr_tmp_12}, torch::kFloat64));
for (int _i = 0; _i < 3 && _i < nuc_pos_2.numel(); _i++) { rhs[m2[2]].position[_i] = nuc_pos_2[_i].template item<double>(); }
torch::Tensor nuc_count = std::get<Microtubule::Nucleator>(rhs[m2[1]].data).Direction.clone();

double _arr_tmp_13 = nuc_count[0].template item<double>();
double _arr_tmp_14 = nuc_count[1].template item<double>() - 1;
torch::Tensor nuc_count_2 = (torch::tensor({_arr_tmp_13, _arr_tmp_14}, torch::kFloat64));
std::get<Microtubule::Nucleator>(rhs[m2[ 2 ]].data).Direction = (nuc_count_2).clone();
std::get<Microtubule::Nucleator>(rhs[m2[ 1 ]].data).Direction = std::get<Microtubule::Nucleator>(rhs[m2[ 1 ]].data).Direction.clone();
std::get<Microtubule::Nucleator>(rhs[m2[ 1 ]].data).Direction[(1)] = (0);
}
);
gamma.addRule(create_nucleator_grid_y);
};
void make_boundary_bottom(DGGML::Grammar<Microtubule::graph_type> &gamma,
           Microtubule::graph_type &system_graph,
           Parameters &settings) {

GT make_boundary_bottom_lhs;
make_boundary_bottom_lhs.addNode({1, {Microtubule::CellBoundary{} }});

GT make_boundary_bottom_rhs;
make_boundary_bottom_rhs.addNode({1, {Microtubule::CellBoundary{} }});

make_boundary_bottom_rhs.addNode({2, {Microtubule::CellBoundary{} }});

make_boundary_bottom_rhs.addEdge(1, 2);

DGGML::WithRule<GT> make_boundary_bottom("make_boundary_bottom", make_boundary_bottom_lhs, make_boundary_bottom_rhs,
[&](auto &lhs, auto &m1) {

torch::Tensor cb_count = std::get<Microtubule::CellBoundary>(lhs[m1[1]].data).Count.clone();

auto boundary_pts =  settings.boundary_pts;

torch::Tensor cb_unit = std::get<Microtubule::CellBoundary>(lhs[m1[1]].data).Direction.clone();

return 
(
( (cb_count[0].template item<double>() != 0) ? 1.0 : 0.0)
 * 
( (cb_count[0].template item<double>() != boundary_pts) ? 1.0 : 0.0)
 * 
( (cb_unit[0].template item<double>() == 1) ? 1.0 : 0.0)
 * 
( (cb_count[0].template item<double>() <= boundary_pts) ? 1.0 : 0.0)
 * 
settings.BOUNDARY_BUILD_SPEED
)
;},
[&](auto &lhs, auto &rhs, auto &m1, auto &m2) {
std::copy(std::begin(lhs[m1[1]].position), std::end(lhs[m1[1]].position), std::begin(rhs[m2[1]].position));

std::get<Microtubule::CellBoundary>(rhs[m2[ 1 ]].data).Direction = std::get<Microtubule::CellBoundary>(lhs[m1[ 1 ]].data).Direction.clone();

std::get<Microtubule::CellBoundary>(rhs[m2[ 1 ]].data).Count = std::get<Microtubule::CellBoundary>(lhs[m1[ 1 ]].data).Count.clone();

torch::Tensor cb_pos = torch::tensor({rhs[m2[1]].position[0], rhs[m2[1]].position[1], rhs[m2[1]].position[2]}, torch::kFloat64);

double _arr_tmp_15 = cb_pos[0].template item<double>() + settings.boundary_offset;
double _arr_tmp_16 = cb_pos[1].template item<double>();
double _arr_tmp_17 = cb_pos[2].template item<double>();
torch::Tensor cb_pos_1 = (torch::tensor({_arr_tmp_15, _arr_tmp_16, _arr_tmp_17}, torch::kFloat64));
for (int _i = 0; _i < 3 && _i < cb_pos_1.numel(); _i++) { rhs[m2[2]].position[_i] = cb_pos_1[_i].template item<double>(); }
torch::Tensor cb_unit = std::get<Microtubule::CellBoundary>(rhs[m2[1]].data).Direction.clone();

double _arr_tmp_18 = cb_unit[0].template item<double>();
double _arr_tmp_19 = cb_unit[1].template item<double>();
double _arr_tmp_20 = cb_unit[2].template item<double>();
torch::Tensor cb_unit_1 = (torch::tensor({_arr_tmp_18, _arr_tmp_19, _arr_tmp_20}, torch::kFloat64));
std::get<Microtubule::CellBoundary>(rhs[m2[ 2 ]].data).Direction = (cb_unit_1).clone();
torch::Tensor cb_count = std::get<Microtubule::CellBoundary>(rhs[m2[1]].data).Count.clone();

double _arr_tmp_21 = cb_count[0].template item<double>() + 1;
double _arr_tmp_22 = cb_count[1].template item<double>();
torch::Tensor cb_count_1 = (torch::tensor({_arr_tmp_21, _arr_tmp_22}, torch::kFloat64));
std::get<Microtubule::CellBoundary>(rhs[m2[ 2 ]].data).Count = (cb_count_1).clone();
std::get<Microtubule::CellBoundary>(rhs[m2[ 1 ]].data).Count = std::get<Microtubule::CellBoundary>(rhs[m2[ 1 ]].data).Count.clone();
std::get<Microtubule::CellBoundary>(rhs[m2[ 1 ]].data).Count[(0)] = (0);
}
);
gamma.addRule(make_boundary_bottom);
};
void make_boundary_right(DGGML::Grammar<Microtubule::graph_type> &gamma,
           Microtubule::graph_type &system_graph,
           Parameters &settings) {

GT make_boundary_right_lhs;
make_boundary_right_lhs.addNode({1, {Microtubule::CellBoundary{} }});

GT make_boundary_right_rhs;
make_boundary_right_rhs.addNode({1, {Microtubule::CellBoundary{} }});

make_boundary_right_rhs.addNode({2, {Microtubule::CellBoundary{} }});

make_boundary_right_rhs.addEdge(1, 2);

DGGML::WithRule<GT> make_boundary_right("make_boundary_right", make_boundary_right_lhs, make_boundary_right_rhs,
[&](auto &lhs, auto &m1) {

auto boundary_pts =  settings.boundary_pts;

torch::Tensor cb_count = std::get<Microtubule::CellBoundary>(lhs[m1[1]].data).Count.clone();

torch::Tensor cb_unit = std::get<Microtubule::CellBoundary>(lhs[m1[1]].data).Direction.clone();

return 
(
( (cb_count[0].template item<double>() >= boundary_pts) ? 1.0 : 0.0)
 * 
( (cb_count[1].template item<double>() != 0) ? 1.0 : 0.0)
 * 
( (cb_count[1].template item<double>() != boundary_pts) ? 1.0 : 0.0)
 * 
( (cb_unit[1].template item<double>() == 1) ? 1.0 : 0.0)
 * 
settings.BOUNDARY_BUILD_SPEED
)
;},
[&](auto &lhs, auto &rhs, auto &m1, auto &m2) {
std::copy(std::begin(lhs[m1[1]].position), std::end(lhs[m1[1]].position), std::begin(rhs[m2[1]].position));

std::get<Microtubule::CellBoundary>(rhs[m2[ 1 ]].data).Direction = std::get<Microtubule::CellBoundary>(lhs[m1[ 1 ]].data).Direction.clone();

std::get<Microtubule::CellBoundary>(rhs[m2[ 1 ]].data).Count = std::get<Microtubule::CellBoundary>(lhs[m1[ 1 ]].data).Count.clone();

torch::Tensor cb_pos = torch::tensor({rhs[m2[1]].position[0], rhs[m2[1]].position[1], rhs[m2[1]].position[2]}, torch::kFloat64);

double _arr_tmp_23 = cb_pos[0].template item<double>();
double _arr_tmp_24 = cb_pos[1].template item<double>() + settings.boundary_offset;
double _arr_tmp_25 = cb_pos[2].template item<double>();
torch::Tensor cb_pos_1 = (torch::tensor({_arr_tmp_23, _arr_tmp_24, _arr_tmp_25}, torch::kFloat64));
for (int _i = 0; _i < 3 && _i < cb_pos_1.numel(); _i++) { rhs[m2[2]].position[_i] = cb_pos_1[_i].template item<double>(); }
double _arr_tmp_26 = -1;
torch::Tensor cb_unit = std::get<Microtubule::CellBoundary>(rhs[m2[1]].data).Direction.clone();

double _arr_tmp_27 = cb_unit[1].template item<double>();
double _arr_tmp_28 = cb_unit[2].template item<double>();
torch::Tensor cb_unit_1 = (torch::tensor({_arr_tmp_26, _arr_tmp_27, _arr_tmp_28}, torch::kFloat64));
std::get<Microtubule::CellBoundary>(rhs[m2[ 2 ]].data).Direction = (cb_unit_1).clone();
torch::Tensor cb_count = std::get<Microtubule::CellBoundary>(rhs[m2[1]].data).Count.clone();

double _arr_tmp_29 = cb_count[0].template item<double>();
double _arr_tmp_30 = cb_count[1].template item<double>() + 1;
torch::Tensor cb_count_1 = (torch::tensor({_arr_tmp_29, _arr_tmp_30}, torch::kFloat64));
std::get<Microtubule::CellBoundary>(rhs[m2[ 2 ]].data).Count = (cb_count_1).clone();
std::get<Microtubule::CellBoundary>(rhs[m2[ 1 ]].data).Count = std::get<Microtubule::CellBoundary>(rhs[m2[ 1 ]].data).Count.clone();
std::get<Microtubule::CellBoundary>(rhs[m2[ 1 ]].data).Count[(1)] = (0);
}
);
gamma.addRule(make_boundary_right);
};
void make_boundary_top(DGGML::Grammar<Microtubule::graph_type> &gamma,
           Microtubule::graph_type &system_graph,
           Parameters &settings) {

GT make_boundary_top_lhs;
make_boundary_top_lhs.addNode({1, {Microtubule::CellBoundary{} }});

GT make_boundary_top_rhs;
make_boundary_top_rhs.addNode({1, {Microtubule::CellBoundary{} }});

make_boundary_top_rhs.addNode({2, {Microtubule::CellBoundary{} }});

make_boundary_top_rhs.addEdge(1, 2);

DGGML::WithRule<GT> make_boundary_top("make_boundary_top", make_boundary_top_lhs, make_boundary_top_rhs,
[&](auto &lhs, auto &m1) {

auto boundary_pts =  settings.boundary_pts;

torch::Tensor cb_count = std::get<Microtubule::CellBoundary>(lhs[m1[1]].data).Count.clone();

torch::Tensor cb_unit = std::get<Microtubule::CellBoundary>(lhs[m1[1]].data).Direction.clone();

return 
(
( (cb_count[0].template item<double>() <= ((2 * boundary_pts) - 2)) ? 1.0 : 0.0)
 * 
( (cb_count[0].template item<double>() >= boundary_pts) ? 1.0 : 0.0)
 * 
( (cb_count[1].template item<double>() >= boundary_pts) ? 1.0 : 0.0)
 * 
( (cb_unit[0].template item<double>() == (-1)) ? 1.0 : 0.0)
 * 
( (cb_unit[1].template item<double>() == 1) ? 1.0 : 0.0)
 * 
settings.BOUNDARY_BUILD_SPEED
)
;},
[&](auto &lhs, auto &rhs, auto &m1, auto &m2) {
std::copy(std::begin(lhs[m1[1]].position), std::end(lhs[m1[1]].position), std::begin(rhs[m2[1]].position));

std::get<Microtubule::CellBoundary>(rhs[m2[ 1 ]].data).Direction = std::get<Microtubule::CellBoundary>(lhs[m1[ 1 ]].data).Direction.clone();

std::get<Microtubule::CellBoundary>(rhs[m2[ 1 ]].data).Count = std::get<Microtubule::CellBoundary>(lhs[m1[ 1 ]].data).Count.clone();

torch::Tensor cb_pos = torch::tensor({rhs[m2[1]].position[0], rhs[m2[1]].position[1], rhs[m2[1]].position[2]}, torch::kFloat64);

double _arr_tmp_31 = cb_pos[0].template item<double>() - settings.boundary_offset;
double _arr_tmp_32 = cb_pos[1].template item<double>();
double _arr_tmp_33 = cb_pos[2].template item<double>();
torch::Tensor cb_pos_1 = (torch::tensor({_arr_tmp_31, _arr_tmp_32, _arr_tmp_33}, torch::kFloat64));
for (int _i = 0; _i < 3 && _i < cb_pos_1.numel(); _i++) { rhs[m2[2]].position[_i] = cb_pos_1[_i].template item<double>(); }
double _arr_tmp_34 = -1;
torch::Tensor cb_unit = std::get<Microtubule::CellBoundary>(rhs[m2[1]].data).Direction.clone();

double _arr_tmp_35 = cb_unit[1].template item<double>();
double _arr_tmp_36 = cb_unit[2].template item<double>();
torch::Tensor cb_unit_1 = (torch::tensor({_arr_tmp_34, _arr_tmp_35, _arr_tmp_36}, torch::kFloat64));
std::get<Microtubule::CellBoundary>(rhs[m2[ 2 ]].data).Direction = (cb_unit_1).clone();
torch::Tensor cb_count = std::get<Microtubule::CellBoundary>(rhs[m2[1]].data).Count.clone();

double _arr_tmp_37 = cb_count[0].template item<double>() + 1;
double _arr_tmp_38 = cb_count[1].template item<double>();
torch::Tensor cb_count_1 = (torch::tensor({_arr_tmp_37, _arr_tmp_38}, torch::kFloat64));
std::get<Microtubule::CellBoundary>(rhs[m2[ 2 ]].data).Count = (cb_count_1).clone();
std::get<Microtubule::CellBoundary>(rhs[m2[ 1 ]].data).Direction = std::get<Microtubule::CellBoundary>(rhs[m2[ 1 ]].data).Direction.clone();
std::get<Microtubule::CellBoundary>(rhs[m2[ 1 ]].data).Direction[(0)] = (1);
}
);
gamma.addRule(make_boundary_top);
};
void make_boundary_left(DGGML::Grammar<Microtubule::graph_type> &gamma,
           Microtubule::graph_type &system_graph,
           Parameters &settings) {

GT make_boundary_left_lhs;
make_boundary_left_lhs.addNode({1, {Microtubule::CellBoundary{} }});

GT make_boundary_left_rhs;
make_boundary_left_rhs.addNode({1, {Microtubule::CellBoundary{} }});

make_boundary_left_rhs.addNode({2, {Microtubule::CellBoundary{} }});

make_boundary_left_rhs.addEdge(1, 2);

DGGML::WithRule<GT> make_boundary_left("make_boundary_left", make_boundary_left_lhs, make_boundary_left_rhs,
[&](auto &lhs, auto &m1) {

auto boundary_pts =  settings.boundary_pts;

torch::Tensor cb_count = std::get<Microtubule::CellBoundary>(lhs[m1[1]].data).Count.clone();

return 
(
( (cb_count[0].template item<double>() > ((2 * boundary_pts) - 2)) ? 1.0 : 0.0)
 * 
( (cb_count[1].template item<double>() >= boundary_pts) ? 1.0 : 0.0)
 * 
( (cb_count[1].template item<double>() <= ((2 * boundary_pts) - 2)) ? 1.0 : 0.0)
 * 
settings.BOUNDARY_BUILD_SPEED
)
;},
[&](auto &lhs, auto &rhs, auto &m1, auto &m2) {
std::copy(std::begin(lhs[m1[1]].position), std::end(lhs[m1[1]].position), std::begin(rhs[m2[1]].position));

std::get<Microtubule::CellBoundary>(rhs[m2[ 1 ]].data).Direction = std::get<Microtubule::CellBoundary>(lhs[m1[ 1 ]].data).Direction.clone();

std::get<Microtubule::CellBoundary>(rhs[m2[ 1 ]].data).Count = std::get<Microtubule::CellBoundary>(lhs[m1[ 1 ]].data).Count.clone();

torch::Tensor cb_pos = torch::tensor({rhs[m2[1]].position[0], rhs[m2[1]].position[1], rhs[m2[1]].position[2]}, torch::kFloat64);

double _arr_tmp_39 = cb_pos[0].template item<double>();
double _arr_tmp_40 = cb_pos[1].template item<double>() - settings.boundary_offset;
double _arr_tmp_41 = cb_pos[2].template item<double>();
torch::Tensor cb_pos_1 = (torch::tensor({_arr_tmp_39, _arr_tmp_40, _arr_tmp_41}, torch::kFloat64));
for (int _i = 0; _i < 3 && _i < cb_pos_1.numel(); _i++) { rhs[m2[2]].position[_i] = cb_pos_1[_i].template item<double>(); }
double _arr_tmp_42 = -1;
double _arr_tmp_43 = -1;
torch::Tensor cb_unit = std::get<Microtubule::CellBoundary>(rhs[m2[1]].data).Direction.clone();

double _arr_tmp_44 = cb_unit[2].template item<double>();
torch::Tensor cb_unit_1 = (torch::tensor({_arr_tmp_42, _arr_tmp_43, _arr_tmp_44}, torch::kFloat64));
std::get<Microtubule::CellBoundary>(rhs[m2[ 2 ]].data).Direction = (cb_unit_1).clone();
torch::Tensor cb_count = std::get<Microtubule::CellBoundary>(rhs[m2[1]].data).Count.clone();

double _arr_tmp_45 = cb_count[0].template item<double>();
double _arr_tmp_46 = cb_count[1].template item<double>() + 1;
torch::Tensor cb_count_1 = (torch::tensor({_arr_tmp_45, _arr_tmp_46}, torch::kFloat64));
std::get<Microtubule::CellBoundary>(rhs[m2[ 2 ]].data).Count = (cb_count_1).clone();
std::get<Microtubule::CellBoundary>(rhs[m2[ 1 ]].data).Count = std::get<Microtubule::CellBoundary>(rhs[m2[ 1 ]].data).Count.clone();
std::get<Microtubule::CellBoundary>(rhs[m2[ 1 ]].data).Count[(1)] = (0);
}
);
gamma.addRule(make_boundary_left);
};
}
#endif