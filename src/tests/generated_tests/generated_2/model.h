#ifndef DGGML_MODELS_HPP
#define DGGML_MODELS_HPP
#include <fstream>
#include "DGGML.h"
#include "rules_0.h"
#include "rules_1.h"
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
    #include <cereal/archives/binary.hpp>
    #include <cereal/types/vector.hpp>
    #include <torch/torch.h>

    namespace cereal {
        // SAVE function
        template<class Archive>
        void save(Archive& ar, const torch::Tensor& tensor) {
            // 1. Save metadata (dimensions and type)
            std::vector<int64_t> size = tensor.sizes().vec();
            int64_t type = static_cast<int64_t>(tensor.scalar_type());
            ar(size, type);

            // 2. Save raw data
            auto contiguous_tensor = tensor.contiguous();
            size_t bytes = contiguous_tensor.nbytes();
            ar(binary_data(contiguous_tensor.data_ptr(), bytes));
        }

        // LOAD function
        template<class Archive>
        void load(Archive& ar, torch::Tensor& tensor) {
            // 1. Load metadata
            std::vector<int64_t> size;
            int64_t type_int;
            ar(size, type_int);

            // 2. Prepare tensor and load raw data
            tensor = torch::empty(size, torch::TensorOptions().dtype(static_cast<torch::ScalarType>(type_int)));
            ar(binary_data(tensor.data_ptr(), tensor.nbytes()));
        }
    }
namespace Dissolution {
using graph_grammar_t = DGGML::Grammar<Dissolution::graph_type>;
class Model_0 : public DGGML::Model3D<graph_grammar_t> {
	public:
	Parameters settings;
bool save_system_graph = true;
bool load_initial_state = false;
std::string initial_state_filename = "simulation_state.bin";
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

            const double max_x = settings.CELL_NX-1;
            const double max_y = settings.CELL_NY-1;

            const std::size_t max_nx = 16;
            const std::size_t max_ny = 16;
            const std::size_t max_nz = 16;

            DGGML::CartesianGrid3D &reaction_grid = cplex.reaction_grid;
            std::size_t max_particles = reaction_grid.totalNumCells();

            std::vector<std::size_t> selected_cells;
            for (int i = max_particles/2; i < max_particles; i++) {
                selected_cells.push_back(i);
            }

            auto centerCardinal = reaction_grid.cardinalCellIndex(
                        reaction_grid._nx/2, reaction_grid._ny/2, reaction_grid._nz/2
                    );

            double center_x, center_y, center_z;
            reaction_grid.cardinalCellToPoint(center_x, center_y, center_z, centerCardinal);
            Dissolution::graph_type tg;
            Dissolution::StartType tmp = Dissolution::StartType{};

            tmp.start_location[0] = center_x;
            tmp.start_location[1] = center_y;
            tmp.start_location[2] = center_z;

            node_type oof_node = {gen.get_key(),
                                {tmp, 
                                 center_x, center_y, center_z}
                        };

            graph.addNode(oof_node);
            };template <typename GraphType>
auto get_type_attributes(GraphType &system_graph) {
std::vector<std::pair<std::string, std::vector<double>>> type_attributes;
double nan_value = std::numeric_limits<double>::quiet_NaN();
std::vector<std::string> col_names;
std::size_t N = system_graph.numNodes();
// Register column names from type definitions
col_names.push_back("FluidSource.Concentration");
col_names.push_back("FluidSource.Pressure");
col_names.push_back("FluidSource.Unit.0");
col_names.push_back("FluidSource.Unit.1");
col_names.push_back("FluidSource.Unit.2");
col_names.push_back("FluidSource.FCount.0");
col_names.push_back("FluidSource.Position.0");
col_names.push_back("FluidSource.Position.1");
col_names.push_back("FluidSource.Position.2");
col_names.push_back("FluidSink.Concentration");
col_names.push_back("FluidSink.Pressure");
col_names.push_back("FluidSink.Unit.0");
col_names.push_back("FluidSink.Unit.1");
col_names.push_back("FluidSink.Unit.2");
col_names.push_back("FluidSink.FCount.0");
col_names.push_back("FluidSink.Position.0");
col_names.push_back("FluidSink.Position.1");
col_names.push_back("FluidSink.Position.2");
col_names.push_back("RockStart.Density");
col_names.push_back("RockStart.Count.0");
col_names.push_back("RockStart.Count.1");
col_names.push_back("RockStart.Count.2");
col_names.push_back("RockStart.Dir.0");
col_names.push_back("RockStart.Dir.1");
col_names.push_back("RockStart.Dir.2");
col_names.push_back("RockStart.Position.0");
col_names.push_back("RockStart.Position.1");
col_names.push_back("RockStart.Position.2");
col_names.push_back("StartType.start_location.0");
col_names.push_back("StartType.start_location.1");
col_names.push_back("StartType.start_location.2");
col_names.push_back("Boundary.boundary_location.0");
col_names.push_back("Boundary.boundary_location.1");
col_names.push_back("Boundary.boundary_location.2");
col_names.push_back("Fluid.Concentration");
col_names.push_back("Fluid.Pressure");
col_names.push_back("Fluid.Unit.0");
col_names.push_back("Fluid.Unit.1");
col_names.push_back("Fluid.Unit.2");
col_names.push_back("Fluid.FCount.0");
col_names.push_back("Fluid.Position.0");
col_names.push_back("Fluid.Position.1");
col_names.push_back("Fluid.Position.2");
std::vector<std::vector<double>> columns(col_names.size());
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
if constexpr (std::is_same_v<T, Dissolution::FluidSource>) {
set("FluidSource.Concentration", row, alt.Concentration);
set("FluidSource.Pressure", row, alt.Pressure);
{
auto flat_tensor = alt.Unit.contiguous().view(-1);
for (int64_t j = 0; j < flat_tensor.size(0); ++j) {
std::string col_key = "FluidSource.Unit." + std::to_string(j);
set(col_key, row, flat_tensor[j].template item<double>());
}
}
{
auto flat_tensor = alt.FCount.contiguous().view(-1);
for (int64_t j = 0; j < flat_tensor.size(0); ++j) {
std::string col_key = "FluidSource.FCount." + std::to_string(j);
set(col_key, row, flat_tensor[j].template item<double>());
}
}
{
auto flat_tensor = alt.Position.contiguous().view(-1);
for (int64_t j = 0; j < flat_tensor.size(0); ++j) {
std::string col_key = "FluidSource.Position." + std::to_string(j);
set(col_key, row, flat_tensor[j].template item<double>());
}
}
}
if constexpr (std::is_same_v<T, Dissolution::FluidSink>) {
set("FluidSink.Concentration", row, alt.Concentration);
set("FluidSink.Pressure", row, alt.Pressure);
{
auto flat_tensor = alt.Unit.contiguous().view(-1);
for (int64_t j = 0; j < flat_tensor.size(0); ++j) {
std::string col_key = "FluidSink.Unit." + std::to_string(j);
set(col_key, row, flat_tensor[j].template item<double>());
}
}
{
auto flat_tensor = alt.FCount.contiguous().view(-1);
for (int64_t j = 0; j < flat_tensor.size(0); ++j) {
std::string col_key = "FluidSink.FCount." + std::to_string(j);
set(col_key, row, flat_tensor[j].template item<double>());
}
}
{
auto flat_tensor = alt.Position.contiguous().view(-1);
for (int64_t j = 0; j < flat_tensor.size(0); ++j) {
std::string col_key = "FluidSink.Position." + std::to_string(j);
set(col_key, row, flat_tensor[j].template item<double>());
}
}
}
if constexpr (std::is_same_v<T, Dissolution::RockStart>) {
set("RockStart.Density", row, alt.Density);
{
auto flat_tensor = alt.Count.contiguous().view(-1);
for (int64_t j = 0; j < flat_tensor.size(0); ++j) {
std::string col_key = "RockStart.Count." + std::to_string(j);
set(col_key, row, flat_tensor[j].template item<double>());
}
}
{
auto flat_tensor = alt.Dir.contiguous().view(-1);
for (int64_t j = 0; j < flat_tensor.size(0); ++j) {
std::string col_key = "RockStart.Dir." + std::to_string(j);
set(col_key, row, flat_tensor[j].template item<double>());
}
}
{
auto flat_tensor = alt.Position.contiguous().view(-1);
for (int64_t j = 0; j < flat_tensor.size(0); ++j) {
std::string col_key = "RockStart.Position." + std::to_string(j);
set(col_key, row, flat_tensor[j].template item<double>());
}
}
}
if constexpr (std::is_same_v<T, Dissolution::StartType>) {
{
auto flat_tensor = alt.start_location.contiguous().view(-1);
for (int64_t j = 0; j < flat_tensor.size(0); ++j) {
std::string col_key = "StartType.start_location." + std::to_string(j);
set(col_key, row, flat_tensor[j].template item<double>());
}
}
}
if constexpr (std::is_same_v<T, Dissolution::Boundary>) {
{
auto flat_tensor = alt.boundary_location.contiguous().view(-1);
for (int64_t j = 0; j < flat_tensor.size(0); ++j) {
std::string col_key = "Boundary.boundary_location." + std::to_string(j);
set(col_key, row, flat_tensor[j].template item<double>());
}
}
}
if constexpr (std::is_same_v<T, Dissolution::Fluid>) {
set("Fluid.Concentration", row, alt.Concentration);
set("Fluid.Pressure", row, alt.Pressure);
{
auto flat_tensor = alt.Unit.contiguous().view(-1);
for (int64_t j = 0; j < flat_tensor.size(0); ++j) {
std::string col_key = "Fluid.Unit." + std::to_string(j);
set(col_key, row, flat_tensor[j].template item<double>());
}
}
{
auto flat_tensor = alt.FCount.contiguous().view(-1);
for (int64_t j = 0; j < flat_tensor.size(0); ++j) {
std::string col_key = "Fluid.FCount." + std::to_string(j);
set(col_key, row, flat_tensor[j].template item<double>());
}
}
{
auto flat_tensor = alt.Position.contiguous().view(-1);
for (int64_t j = 0; j < flat_tensor.size(0); ++j) {
std::string col_key = "Fluid.Position." + std::to_string(j);
set(col_key, row, flat_tensor[j].template item<double>());
}
}
}
}, n.data);
}
for (std::size_t i = 0; i < col_names.size(); ++i) {
        type_attributes.emplace_back(col_names[i], std::move(columns[i]));
    }
return type_attributes;}
void load_graph(DGGML::Grammar<Dissolution::graph_type> &grammar,
               Dissolution::graph_type &system_graph,
               Parameters &settings, DGGML::KeyGenerator<key_type> &gen,
               std::string filename = "simulation_state.bin"
               ) {

    std::ifstream is(filename, std::ios::binary);
    if (!is) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    // Create the input archive
    cereal::BinaryInputArchive archive(is);

    size_t num_nodes;
    archive(num_nodes);
    std::cout << "Number of nodes to load: " << num_nodes << std::endl;

    // Check to see if key already exists in the graph
    std::unordered_map<key_type, key_type> existing_keys;

    std::unordered_map<key_type, std::unordered_set<key_type>> out_edges_map;
    std::unordered_map<key_type, std::unordered_set<key_type>> in_edges_map;

    for (unsigned int i=0; i < num_nodes; i++) {

        key_type old_key;

        SpatialNode3D<Dissolution::StartType, Dissolution::Boundary, Dissolution::Fluid, Dissolution::FluidSource, Dissolution::FluidSink, Dissolution::RockStart> node_data;
        std::unordered_set<key_type> out_neighbors;
        std::unordered_set<key_type> in_neighbors;

        archive(old_key, node_data);
        archive(out_neighbors);
        archive(in_neighbors);

        key_type new_node_key;
        new_node_key = gen.get_key();
        existing_keys[old_key] = new_node_key;

        out_edges_map[new_node_key] = out_neighbors;
        in_edges_map[new_node_key] = in_neighbors;

        std::visit([&](auto &alt) {
using T = std::decay_t<decltype(alt)>;
if constexpr (std::is_same_v<T, Dissolution::StartType>) {
Dissolution::StartType copy = alt;
system_graph.addNode({new_node_key, {copy,
        node_data.position[0], node_data.position[1], node_data.position[2]} });
} else if constexpr (std::is_same_v<T, Dissolution::Boundary>) {
Dissolution::Boundary copy = alt;
system_graph.addNode({new_node_key, {copy,
        node_data.position[0], node_data.position[1], node_data.position[2]} });
} else if constexpr (std::is_same_v<T, Dissolution::Fluid>) {
Dissolution::Fluid copy = alt;
system_graph.addNode({new_node_key, {copy,
        node_data.position[0], node_data.position[1], node_data.position[2]} });
} else if constexpr (std::is_same_v<T, Dissolution::FluidSource>) {
Dissolution::FluidSource copy = alt;
system_graph.addNode({new_node_key, {copy,
        node_data.position[0], node_data.position[1], node_data.position[2]} });
} else if constexpr (std::is_same_v<T, Dissolution::FluidSink>) {
Dissolution::FluidSink copy = alt;
system_graph.addNode({new_node_key, {copy,
        node_data.position[0], node_data.position[1], node_data.position[2]} });
} else if constexpr (std::is_same_v<T, Dissolution::RockStart>) {
Dissolution::RockStart copy = alt;
system_graph.addNode({new_node_key, {copy,
        node_data.position[0], node_data.position[1], node_data.position[2]} });
} else {
std::cout << "  Type: Unknown\n";exit(0);
}
}, node_data.data);

    }

    // Adding outgoing edges
    for(const auto& [old_key, out_neighbors] : out_edges_map) {
        for (const auto& old_neighbor_key : out_neighbors) {
            key_type new_neighbor_key = existing_keys[old_neighbor_key];
            system_graph.addEdge(old_key, new_neighbor_key);
        }
    }

    // Add incoming edges
    for(const auto& [old_key, in_neighbors] : in_edges_map) {
        for (const auto& old_neighbor_key : in_neighbors) {
            key_type new_neighbor_key = existing_keys[old_neighbor_key];
            system_graph.addEdge(new_neighbor_key, old_key);
        }
    }

};
void initialize() override {

             // Per-stage simulation time override
             settings.TOTAL_TIME = 3;
             settings.NUM_STEPS = static_cast<int>(settings.TOTAL_TIME / settings.DELTA);
            int geoplex_size = settings.CELL_NX;
             int cell_nx = geoplex_size;
             int cell_ny = geoplex_size;
             int cell_nz = geoplex_size;

             double cell_dx = settings.CELL_DX;
             double cell_dy = settings.CELL_DY;
             double cell_dz = settings.CELL_DZ;

             double geoplex_epsilon = settings.varepsilon;
            geoplex2D.init(
                            cell_nx,
                            cell_ny,
                            cell_nz,
                            cell_dx,
                            cell_dy,
                            cell_dz,
                            false,
                            geoplex_epsilon
                    );
            	Dissolution::propagate_rock_y(gamma, this->system_graph, settings);
	Dissolution::propagate_from_source(gamma, this->system_graph, settings);
	Dissolution::propagate_source_to_sink(gamma, this->system_graph, settings);
	Dissolution::propagate_rock_x(gamma, this->system_graph, settings);
	Dissolution::propagate_rock_z(gamma, this->system_graph, settings);
	Dissolution::propagate_initial_fluid(gamma, this->system_graph, settings);
	Dissolution::propagate_final_fluid(gamma, this->system_graph, settings);
	Dissolution::start_rock_prop(gamma, this->system_graph, settings);

		 if (load_initial_state == false) {this->add_default_type(this->system_graph,
                        geoplex2D,
                        settings,
                        this->gen);} else {load_graph(gamma, this->system_graph, settings, this->gen, 
                        this->initial_state_filename);}
}
void save_graph(DGGML::Grammar<Dissolution::graph_type> &grammar,
            Dissolution::graph_type &system_graph,
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
 if (save_system_graph) {save_graph(gamma, this->system_graph, settings, results_dir_name+"/simulation_state_"+std::to_string(step)+".bin");
save_graph(gamma, this->system_graph, settings, results_dir_name+"/simulation_state_latest.bin");}
 }
}; // end Model_0

class Model_1 : public DGGML::Model3D<graph_grammar_t> {
	public:
	Parameters settings;
bool save_system_graph = true;
bool load_initial_state = true;
std::string initial_state_filename = "my_results/simulation_state_latest.bin";
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

            const double max_x = settings.CELL_NX-1;
            const double max_y = settings.CELL_NY-1;

            const std::size_t max_nx = 16;
            const std::size_t max_ny = 16;
            const std::size_t max_nz = 16;

            DGGML::CartesianGrid3D &reaction_grid = cplex.reaction_grid;
            std::size_t max_particles = reaction_grid.totalNumCells();

            std::vector<std::size_t> selected_cells;
            for (int i = max_particles/2; i < max_particles; i++) {
                selected_cells.push_back(i);
            }

            auto centerCardinal = reaction_grid.cardinalCellIndex(
                        reaction_grid._nx/2, reaction_grid._ny/2, reaction_grid._nz/2
                    );

            double center_x, center_y, center_z;
            reaction_grid.cardinalCellToPoint(center_x, center_y, center_z, centerCardinal);
            Dissolution::graph_type tg;
            Dissolution::StartType tmp = Dissolution::StartType{};

            tmp.start_location[0] = center_x;
            tmp.start_location[1] = center_y;
            tmp.start_location[2] = center_z;

            node_type oof_node = {gen.get_key(),
                                {tmp, 
                                 center_x, center_y, center_z}
                        };

            graph.addNode(oof_node);
            };template <typename GraphType>
auto get_type_attributes(GraphType &system_graph) {
std::vector<std::pair<std::string, std::vector<double>>> type_attributes;
double nan_value = std::numeric_limits<double>::quiet_NaN();
std::vector<std::string> col_names;
std::size_t N = system_graph.numNodes();
// Register column names from type definitions
col_names.push_back("FluidSource.Concentration");
col_names.push_back("FluidSource.Pressure");
col_names.push_back("FluidSource.Unit.0");
col_names.push_back("FluidSource.Unit.1");
col_names.push_back("FluidSource.Unit.2");
col_names.push_back("FluidSource.FCount.0");
col_names.push_back("FluidSource.Position.0");
col_names.push_back("FluidSource.Position.1");
col_names.push_back("FluidSource.Position.2");
col_names.push_back("FluidSink.Concentration");
col_names.push_back("FluidSink.Pressure");
col_names.push_back("FluidSink.Unit.0");
col_names.push_back("FluidSink.Unit.1");
col_names.push_back("FluidSink.Unit.2");
col_names.push_back("FluidSink.FCount.0");
col_names.push_back("FluidSink.Position.0");
col_names.push_back("FluidSink.Position.1");
col_names.push_back("FluidSink.Position.2");
col_names.push_back("RockStart.Density");
col_names.push_back("RockStart.Count.0");
col_names.push_back("RockStart.Count.1");
col_names.push_back("RockStart.Count.2");
col_names.push_back("RockStart.Dir.0");
col_names.push_back("RockStart.Dir.1");
col_names.push_back("RockStart.Dir.2");
col_names.push_back("RockStart.Position.0");
col_names.push_back("RockStart.Position.1");
col_names.push_back("RockStart.Position.2");
col_names.push_back("StartType.start_location.0");
col_names.push_back("StartType.start_location.1");
col_names.push_back("StartType.start_location.2");
col_names.push_back("Boundary.boundary_location.0");
col_names.push_back("Boundary.boundary_location.1");
col_names.push_back("Boundary.boundary_location.2");
col_names.push_back("Fluid.Concentration");
col_names.push_back("Fluid.Pressure");
col_names.push_back("Fluid.Unit.0");
col_names.push_back("Fluid.Unit.1");
col_names.push_back("Fluid.Unit.2");
col_names.push_back("Fluid.FCount.0");
col_names.push_back("Fluid.Position.0");
col_names.push_back("Fluid.Position.1");
col_names.push_back("Fluid.Position.2");
std::vector<std::vector<double>> columns(col_names.size());
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
if constexpr (std::is_same_v<T, Dissolution::FluidSource>) {
set("FluidSource.Concentration", row, alt.Concentration);
set("FluidSource.Pressure", row, alt.Pressure);
{
auto flat_tensor = alt.Unit.contiguous().view(-1);
for (int64_t j = 0; j < flat_tensor.size(0); ++j) {
std::string col_key = "FluidSource.Unit." + std::to_string(j);
set(col_key, row, flat_tensor[j].template item<double>());
}
}
{
auto flat_tensor = alt.FCount.contiguous().view(-1);
for (int64_t j = 0; j < flat_tensor.size(0); ++j) {
std::string col_key = "FluidSource.FCount." + std::to_string(j);
set(col_key, row, flat_tensor[j].template item<double>());
}
}
{
auto flat_tensor = alt.Position.contiguous().view(-1);
for (int64_t j = 0; j < flat_tensor.size(0); ++j) {
std::string col_key = "FluidSource.Position." + std::to_string(j);
set(col_key, row, flat_tensor[j].template item<double>());
}
}
}
if constexpr (std::is_same_v<T, Dissolution::FluidSink>) {
set("FluidSink.Concentration", row, alt.Concentration);
set("FluidSink.Pressure", row, alt.Pressure);
{
auto flat_tensor = alt.Unit.contiguous().view(-1);
for (int64_t j = 0; j < flat_tensor.size(0); ++j) {
std::string col_key = "FluidSink.Unit." + std::to_string(j);
set(col_key, row, flat_tensor[j].template item<double>());
}
}
{
auto flat_tensor = alt.FCount.contiguous().view(-1);
for (int64_t j = 0; j < flat_tensor.size(0); ++j) {
std::string col_key = "FluidSink.FCount." + std::to_string(j);
set(col_key, row, flat_tensor[j].template item<double>());
}
}
{
auto flat_tensor = alt.Position.contiguous().view(-1);
for (int64_t j = 0; j < flat_tensor.size(0); ++j) {
std::string col_key = "FluidSink.Position." + std::to_string(j);
set(col_key, row, flat_tensor[j].template item<double>());
}
}
}
if constexpr (std::is_same_v<T, Dissolution::RockStart>) {
set("RockStart.Density", row, alt.Density);
{
auto flat_tensor = alt.Count.contiguous().view(-1);
for (int64_t j = 0; j < flat_tensor.size(0); ++j) {
std::string col_key = "RockStart.Count." + std::to_string(j);
set(col_key, row, flat_tensor[j].template item<double>());
}
}
{
auto flat_tensor = alt.Dir.contiguous().view(-1);
for (int64_t j = 0; j < flat_tensor.size(0); ++j) {
std::string col_key = "RockStart.Dir." + std::to_string(j);
set(col_key, row, flat_tensor[j].template item<double>());
}
}
{
auto flat_tensor = alt.Position.contiguous().view(-1);
for (int64_t j = 0; j < flat_tensor.size(0); ++j) {
std::string col_key = "RockStart.Position." + std::to_string(j);
set(col_key, row, flat_tensor[j].template item<double>());
}
}
}
if constexpr (std::is_same_v<T, Dissolution::StartType>) {
{
auto flat_tensor = alt.start_location.contiguous().view(-1);
for (int64_t j = 0; j < flat_tensor.size(0); ++j) {
std::string col_key = "StartType.start_location." + std::to_string(j);
set(col_key, row, flat_tensor[j].template item<double>());
}
}
}
if constexpr (std::is_same_v<T, Dissolution::Boundary>) {
{
auto flat_tensor = alt.boundary_location.contiguous().view(-1);
for (int64_t j = 0; j < flat_tensor.size(0); ++j) {
std::string col_key = "Boundary.boundary_location." + std::to_string(j);
set(col_key, row, flat_tensor[j].template item<double>());
}
}
}
if constexpr (std::is_same_v<T, Dissolution::Fluid>) {
set("Fluid.Concentration", row, alt.Concentration);
set("Fluid.Pressure", row, alt.Pressure);
{
auto flat_tensor = alt.Unit.contiguous().view(-1);
for (int64_t j = 0; j < flat_tensor.size(0); ++j) {
std::string col_key = "Fluid.Unit." + std::to_string(j);
set(col_key, row, flat_tensor[j].template item<double>());
}
}
{
auto flat_tensor = alt.FCount.contiguous().view(-1);
for (int64_t j = 0; j < flat_tensor.size(0); ++j) {
std::string col_key = "Fluid.FCount." + std::to_string(j);
set(col_key, row, flat_tensor[j].template item<double>());
}
}
{
auto flat_tensor = alt.Position.contiguous().view(-1);
for (int64_t j = 0; j < flat_tensor.size(0); ++j) {
std::string col_key = "Fluid.Position." + std::to_string(j);
set(col_key, row, flat_tensor[j].template item<double>());
}
}
}
}, n.data);
}
for (std::size_t i = 0; i < col_names.size(); ++i) {
        type_attributes.emplace_back(col_names[i], std::move(columns[i]));
    }
return type_attributes;}
void load_graph(DGGML::Grammar<Dissolution::graph_type> &grammar,
               Dissolution::graph_type &system_graph,
               Parameters &settings, DGGML::KeyGenerator<key_type> &gen,
               std::string filename = "simulation_state.bin"
               ) {

    std::ifstream is(filename, std::ios::binary);
    if (!is) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    // Create the input archive
    cereal::BinaryInputArchive archive(is);

    size_t num_nodes;
    archive(num_nodes);
    std::cout << "Number of nodes to load: " << num_nodes << std::endl;

    // Check to see if key already exists in the graph
    std::unordered_map<key_type, key_type> existing_keys;

    std::unordered_map<key_type, std::unordered_set<key_type>> out_edges_map;
    std::unordered_map<key_type, std::unordered_set<key_type>> in_edges_map;

    for (unsigned int i=0; i < num_nodes; i++) {

        key_type old_key;

        SpatialNode3D<Dissolution::StartType, Dissolution::Boundary, Dissolution::Fluid, Dissolution::FluidSource, Dissolution::FluidSink, Dissolution::RockStart> node_data;
        std::unordered_set<key_type> out_neighbors;
        std::unordered_set<key_type> in_neighbors;

        archive(old_key, node_data);
        archive(out_neighbors);
        archive(in_neighbors);

        key_type new_node_key;
        new_node_key = gen.get_key();
        existing_keys[old_key] = new_node_key;

        out_edges_map[new_node_key] = out_neighbors;
        in_edges_map[new_node_key] = in_neighbors;

        std::visit([&](auto &alt) {
using T = std::decay_t<decltype(alt)>;
if constexpr (std::is_same_v<T, Dissolution::StartType>) {
Dissolution::StartType copy = alt;
system_graph.addNode({new_node_key, {copy,
        node_data.position[0], node_data.position[1], node_data.position[2]} });
} else if constexpr (std::is_same_v<T, Dissolution::Boundary>) {
Dissolution::Boundary copy = alt;
system_graph.addNode({new_node_key, {copy,
        node_data.position[0], node_data.position[1], node_data.position[2]} });
} else if constexpr (std::is_same_v<T, Dissolution::Fluid>) {
Dissolution::Fluid copy = alt;
system_graph.addNode({new_node_key, {copy,
        node_data.position[0], node_data.position[1], node_data.position[2]} });
} else if constexpr (std::is_same_v<T, Dissolution::FluidSource>) {
Dissolution::FluidSource copy = alt;
system_graph.addNode({new_node_key, {copy,
        node_data.position[0], node_data.position[1], node_data.position[2]} });
} else if constexpr (std::is_same_v<T, Dissolution::FluidSink>) {
Dissolution::FluidSink copy = alt;
system_graph.addNode({new_node_key, {copy,
        node_data.position[0], node_data.position[1], node_data.position[2]} });
} else if constexpr (std::is_same_v<T, Dissolution::RockStart>) {
Dissolution::RockStart copy = alt;
system_graph.addNode({new_node_key, {copy,
        node_data.position[0], node_data.position[1], node_data.position[2]} });
} else {
std::cout << "  Type: Unknown\n";exit(0);
}
}, node_data.data);

    }

    // Adding outgoing edges
    for(const auto& [old_key, out_neighbors] : out_edges_map) {
        for (const auto& old_neighbor_key : out_neighbors) {
            key_type new_neighbor_key = existing_keys[old_neighbor_key];
            system_graph.addEdge(old_key, new_neighbor_key);
        }
    }

    // Add incoming edges
    for(const auto& [old_key, in_neighbors] : in_edges_map) {
        for (const auto& old_neighbor_key : in_neighbors) {
            key_type new_neighbor_key = existing_keys[old_neighbor_key];
            system_graph.addEdge(new_neighbor_key, old_key);
        }
    }

};
void initialize() override {

             // Per-stage simulation time override
             settings.TOTAL_TIME = 10;
             settings.NUM_STEPS = static_cast<int>(settings.TOTAL_TIME / settings.DELTA);
            int geoplex_size = settings.CELL_NX;
             int cell_nx = geoplex_size;
             int cell_ny = geoplex_size;
             int cell_nz = geoplex_size;

             double cell_dx = settings.CELL_DX;
             double cell_dy = settings.CELL_DY;
             double cell_dz = settings.CELL_DZ;

             double geoplex_epsilon = settings.varepsilon;
            geoplex2D.init(
                            cell_nx,
                            cell_ny,
                            cell_nz,
                            cell_dx,
                            cell_dy,
                            cell_dz,
                            false,
                            geoplex_epsilon
                    );
            	Dissolution::fluid_flow(gamma, this->system_graph, settings);
	Dissolution::source_flow(gamma, this->system_graph, settings);
	Dissolution::erode_rock(gamma, this->system_graph, settings);
	Dissolution::sink_flow(gamma, this->system_graph, settings);

		 if (load_initial_state == false) {this->add_default_type(this->system_graph,
                        geoplex2D,
                        settings,
                        this->gen);} else {load_graph(gamma, this->system_graph, settings, this->gen, 
                        this->initial_state_filename);}
}
void save_graph(DGGML::Grammar<Dissolution::graph_type> &grammar,
            Dissolution::graph_type &system_graph,
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
 if (save_system_graph) {save_graph(gamma, this->system_graph, settings, results_dir_name+"/simulation_state_"+std::to_string(step)+".bin");
save_graph(gamma, this->system_graph, settings, results_dir_name+"/simulation_state_latest.bin");}
 }
}; // end Model_1

};
#endif
