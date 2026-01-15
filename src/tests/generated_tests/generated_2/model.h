#ifndef DGGML_MODELS_HPP
#define DGGML_MODELS_HPP
#include <fstream>
#include "DGGML.h"
#include "rules.h"
#include "simdjson.h"
#include "ExpandedComplex2D.hpp"
#include "ExpandedComplex3D.hpp"
#include "YAGL_Algorithms.hpp"
#include "parameters.h"
#include <cereal/archives/binary.hpp>
#include <cereal/types/variant.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/unordered_set.hpp>
namespace cereal {
        // This function tells cereal how to handle SpatialNode3D without
        // modifying the DGGML source code.
        template <class Archive, typename... Ts>
        void serialize(Archive& ar, SpatialNode3D<Ts...>& node) {
            // Map the internal members of the class to the archive
            ar(node.position);
            ar(node.data);
        }

    }
    namespace Microtubule {
using graph_grammar_t = DGGML::Grammar<Microtubule::graph_type>;
class Model : public DGGML::Model3D<graph_grammar_t> {
	public:
	Parameters settings;
template<typename GraphType, typename CplexType, typename ParamType, typename GenType>
    void add_default_type(
                    GraphType &graph,
                    CplexType &cplex,
                    ParamType &settings, GenType &gen
                    ) {

            using node_type = typename GraphType::node_type;
            const double min_x = 0.0;
            const double min_y = 0.0;
        const double min_z = 0.0;

            // const double max_x = 4.98;
            const double max_x = settings.CELL_NX-1;
            // const double max_y = 4.98;
            const double max_y = settings.CELL_NY-1;

            const std::size_t max_nx = 16;
            const std::size_t max_ny = 16;
        const std::size_t max_nz = 16;

            // DGGML::CartesianGrid2D &reaction_grid = cplex.reaction_grid;
            DGGML::CartesianGrid3D &reaction_grid = cplex.reaction_grid;
            std::size_t max_particles = reaction_grid.totalNumCells();

            //step 1 create an ordered number of selectable cells
            std::vector<std::size_t> selected_cells;
            for (int i = max_particles/2; i < max_particles; i++) {
                selected_cells.push_back(i);
            }

        auto centerCardinal = reaction_grid.cardinalCellIndex(
                        reaction_grid._nx/2, reaction_grid._ny/2, reaction_grid._nz/2
                    );

            double center_x, center_y, center_z;
            // reaction_grid.cardinalCellToPoint(center_x, center_y, centerCardinal);
            reaction_grid.cardinalCellToPoint(center_x, center_y, center_z, centerCardinal);
            Microtubule::graph_type tg;
            Microtubule::StartType tmp = Microtubule::StartType{};

            // Initial placement of particle node creator
            tmp.start_location[0] = center_x; // X
            tmp.start_location[1] = center_y; // Y
            tmp.start_location[2] = center_z; // Y
                                                       //
            node_type oof_node = {gen.get_key(),//i*segments,
                                {tmp, 
                                 center_x, center_y, center_z}
                        };

            graph.addNode(oof_node);
            };template <typename GraphType>
auto get_type_attributes(GraphType &system_graph) {
std::vector<std::pair<std::string, std::vector<double>>> type_attributes;
double nan_value = std::numeric_limits<double>::quiet_NaN();
const std::array<std::string, 48> col_names = {
"StartType.2", 
"StartType.3", 
"StartType.1", 
"Boundary.2", 
"Boundary.3", 
"Boundary.1", 
"FractureSegment.5", 
"FractureSegment.4", 
"FractureSegment.6", 
"FractureSegment.7", 
"FractureSegment.2", 
"FractureSegment.3", 
"FractureSegment.1", 
"FractureSegmentEnd.5", 
"FractureSegmentEnd.4", 
"FractureSegmentEnd.6", 
"FractureSegmentEnd.7", 
"FractureSegmentEnd.2", 
"FractureSegmentEnd.3", 
"FractureSegmentEnd.1", 
"Junction.5", 
"Junction.4", 
"Junction.6", 
"Junction.7", 
"Junction.2", 
"Junction.3", 
"Junction.1", 
"PressureSegment.4", 
"PressureSegment.2", 
"PressureSegment.3", 
"PressureSegment.1", 
"PressureSource.5", 
"PressureSource.4", 
"PressureSource.6", 
"PressureSource.7", 
"PressureSource.2", 
"PressureSource.3", 
"PressureSource.1", 
"RockStart.5", 
"RockStart.4", 
"RockStart.6", 
"RockStart.2", 
"RockStart.3", 
"RockStart.1", 
"Rock.4", 
"Rock.2", 
"Rock.3", 
"Rock.1", 
};

std::array<std::vector<double>, 48> columns;
std::size_t N = system_graph.numNodes();
for (auto &col : columns) {
        col.resize(N, nan_value);
        }
std::unordered_map<std::string, std::size_t> col_idx;
for (std::size_t i = 0; i < col_names.size(); ++i) {
        col_idx[col_names[i]] = i;
        }
    auto set = [&](const std::string &col_name, std::size_t row, double value) {
        columns[col_idx.at(col_name)][row] = value;
    };

std::size_t row = 0;
for (auto it = system_graph.node_list_begin(); it != system_graph.node_list_end(); ++it, ++row) {
auto &n = it->second.getData();
std::visit([&](auto &alt) {
using T = std::decay_t<decltype(alt)>;
if constexpr (std::is_same_v<T, Microtubule::PressureSource>) {
set("PressureSource.6", row, alt.fflow_297c72[1]);
set("PressureSource.3", row, alt.fflow_583de4[2]);
set("PressureSource.2", row, alt.fflow_583de4[1]);
set("PressureSource.4", row, alt.fflow_f25522);
set("PressureSource.7", row, alt.fflow_297c72[2]);
set("PressureSource.5", row, alt.fflow_297c72[0]);
set("PressureSource.1", row, alt.fflow_583de4[0]);
}
if constexpr (std::is_same_v<T, Microtubule::FractureSegmentEnd>) {
set("FractureSegmentEnd.3", row, alt.fflow_8227fa[2]);
set("FractureSegmentEnd.2", row, alt.fflow_8227fa[1]);
set("FractureSegmentEnd.5", row, alt.fflow_30aaee[0]);
set("FractureSegmentEnd.4", row, alt.fflow_e957c5);
set("FractureSegmentEnd.6", row, alt.fflow_30aaee[1]);
set("FractureSegmentEnd.7", row, alt.fflow_30aaee[2]);
set("FractureSegmentEnd.1", row, alt.fflow_8227fa[0]);
}
if constexpr (std::is_same_v<T, Microtubule::Junction>) {
set("Junction.7", row, alt.fflow_f00b58[2]);
set("Junction.1", row, alt.fflow_390a7c[0]);
set("Junction.6", row, alt.fflow_f00b58[1]);
set("Junction.5", row, alt.fflow_f00b58[0]);
set("Junction.4", row, alt.fflow_b78c6a);
set("Junction.3", row, alt.fflow_390a7c[2]);
set("Junction.2", row, alt.fflow_390a7c[1]);
}
if constexpr (std::is_same_v<T, Microtubule::RockStart>) {
set("RockStart.3", row, alt.fflow_a9c368[2]);
set("RockStart.4", row, alt.fflow_a13208);
set("RockStart.5", row, alt.fflow_992eaa);
set("RockStart.2", row, alt.fflow_a9c368[1]);
set("RockStart.6", row, alt.fflow_d20285);
set("RockStart.1", row, alt.fflow_a9c368[0]);
}
if constexpr (std::is_same_v<T, Microtubule::StartType>) {
set("StartType.1", row, alt.start_location[0]);
set("StartType.2", row, alt.start_location[1]);
set("StartType.3", row, alt.start_location[2]);
}
if constexpr (std::is_same_v<T, Microtubule::FractureSegment>) {
set("FractureSegment.6", row, alt.fflow_6e2a05[1]);
set("FractureSegment.7", row, alt.fflow_6e2a05[2]);
set("FractureSegment.1", row, alt.fflow_6b3cfc[0]);
set("FractureSegment.5", row, alt.fflow_6e2a05[0]);
set("FractureSegment.4", row, alt.fflow_122ce5);
set("FractureSegment.2", row, alt.fflow_6b3cfc[1]);
set("FractureSegment.3", row, alt.fflow_6b3cfc[2]);
}
if constexpr (std::is_same_v<T, Microtubule::PressureSegment>) {
set("PressureSegment.1", row, alt.fflow_26d4e5[0]);
set("PressureSegment.4", row, alt.fflow_37b7e1);
set("PressureSegment.2", row, alt.fflow_26d4e5[1]);
set("PressureSegment.3", row, alt.fflow_26d4e5[2]);
}
if constexpr (std::is_same_v<T, Microtubule::Rock>) {
set("Rock.4", row, alt.fflow_9ca0e1);
set("Rock.1", row, alt.fflow_e73169[0]);
set("Rock.2", row, alt.fflow_e73169[1]);
set("Rock.3", row, alt.fflow_e73169[2]);
}
if constexpr (std::is_same_v<T, Microtubule::Boundary>) {
set("Boundary.2", row, alt.boundary_location[1]);
set("Boundary.3", row, alt.boundary_location[2]);
set("Boundary.1", row, alt.boundary_location[0]);
}
}, n.data);
}
for (std::size_t i = 0; i < col_names.size(); ++i) {
        type_attributes.emplace_back(col_names[i], std::move(columns[i]));
    }
return type_attributes;}void initialize() override {
int geoplex_size = 3;
         int cell_nx = geoplex_size;
         int cell_ny = geoplex_size;
         int cell_nz = geoplex_size;

         double cell_dx = 1.0;
         double cell_dy = 1.0;
         double cell_dz = 1.0;

         double maximal_rx_radius = 0.05;
        geoplex2D.init(
                        cell_nx,
                        cell_ny,
                        cell_nz,
                        cell_dx,
                        cell_dy,
                        cell_dz,
                        false,
                        maximal_rx_radius
                );
        	FractureNetwork::propagate_rock_y(gamma, this->system_graph, settings);
	FractureNetwork::propagate_rock_x(gamma, this->system_graph, settings);
	FractureNetwork::propagate_rock_z(gamma, this->system_graph, settings);
	FractureNetwork::start_rock_prop(gamma, this->system_graph, settings);

		this->add_default_type(this->system_graph,
                    geoplex2D,
                    settings,
                    this->gen);}
void save_graph(DGGML::Grammar<Microtubule::graph_type> &grammar,
            Microtubule::graph_type &system_graph,
            Parameters &settings, std::string filename = "simulation_state.bin"
        ) {

        std::ofstream os(filename, std::ios::binary);
        cereal::BinaryOutputArchive archive(os);

        size_t num_nodes = system_graph.numNodes();
        archive(num_nodes);

        auto check = system_graph.getNodeSetRef();

        int count_nodes = 0;
        // Loop through system graph
        for (auto i=system_graph.node_list_begin(); i != system_graph.node_list_end(); ++i) {
            count_nodes++;
            auto& key =  i->first;
            auto& node_data =  i->second.getData();

            archive(key, node_data);
            archive(system_graph.out_neighbors(i->second));
            archive(system_graph.in_neighbors(i->second));

        }

    };
    void checkpoint(std::size_t step) override {std::string results_dir_name = "my_results";if(step == 0)
                {
                    // Create the local save directory
                    std::filesystem::remove_all(results_dir_name);
                    std::filesystem::create_directory(results_dir_name);

                    DGGML::GridFileWriter3D grid_writer;
                    grid_writer.save({geoplex2D.reaction_grid,geoplex2D.dim_label},
                                     results_dir_name+"/expanded_cell_complex");

                    std::string title = results_dir_name+"/simulation_step_";
                    DGGML::VtkFileWriterComplete<graph_type> vtk_writer;
                    auto attr = this->get_type_attributes(this->system_graph);
                    vtk_writer.set_extra_point_data(attr);

                    vtk_writer.save(system_graph, title+std::to_string(step));
                    collect(step);
                }
                if( step != 0 & step % 1 == 0)
                {
                    std::string title = results_dir_name+"/simulation_step_";
                    DGGML::VtkFileWriterComplete<graph_type> vtk_writer;

                    auto attr = this->get_type_attributes(this->system_graph);
                    vtk_writer.set_extra_point_data(attr);

                    vtk_writer.save(system_graph, title+std::to_string(step));
                    collect(step);}
 }};
};
#endif