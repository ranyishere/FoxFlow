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
namespace Microtubule {
using graph_grammar_t = DGGML::Grammar<Microtubule::graph_type>;
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
            Microtubule::graph_type tg;
            Microtubule::StartType tmp = Microtubule::StartType{};

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
// Register column names
std::vector<std::vector<double>> columns(col_names.size());
for (auto &col : columns) { col.resize(N, nan_value); }
std::unordered_map<std::string, std::size_t> col_idx;
for (std::size_t i = 0; i < col_names.size(); ++i) { col_idx[col_names[i]] = i; }
    auto set = [&](const std::string &col_name, std::size_t row, double value) {
        columns[col_idx.at(col_name)][row] = value;
    };

std::size_t row = 0;
for (auto it = system_graph.node_list_begin(); it != system_graph.node_list_end(); ++it, ++row) {
auto &n = it->second.getData();
std::visit([&](auto &alt) {
using T = std::decay_t<decltype(alt)>;
if constexpr (std::is_same_v<T, Microtubule::CellBoundary>) {
}
if constexpr (std::is_same_v<T, Microtubule::Nucleator>) {
}
if constexpr (std::is_same_v<T, Microtubule::Zipper>) {
}
if constexpr (std::is_same_v<T, Microtubule::Intermediate>) {
}
if constexpr (std::is_same_v<T, Microtubule::Negative>) {
}
if constexpr (std::is_same_v<T, Microtubule::StartType>) {
}
if constexpr (std::is_same_v<T, Microtubule::Junction>) {
}
if constexpr (std::is_same_v<T, Microtubule::Positive>) {
}
if constexpr (std::is_same_v<T, Microtubule::Boundary>) {
}
}, n.data);
}
for (std::size_t i = 0; i < col_names.size(); ++i) {
        type_attributes.emplace_back(col_names[i], std::move(columns[i]));
    }
return type_attributes;}
void load_graph(DGGML::Grammar<Microtubule::graph_type> &grammar,
               Microtubule::graph_type &system_graph,
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

        SpatialNode3D<Microtubule::StartType, Microtubule::Boundary, Microtubule::Intermediate, Microtubule::Positive, Microtubule::Negative, Microtubule::Nucleator, Microtubule::CellBoundary, Microtubule::Zipper, Microtubule::Junction> node_data;
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
if constexpr (std::is_same_v<T, Microtubule::StartType>) {
Microtubule::StartType copy = alt;
system_graph.addNode({new_node_key, {copy,
        node_data.position[0], node_data.position[1], node_data.position[2]} });
} else if constexpr (std::is_same_v<T, Microtubule::Boundary>) {
Microtubule::Boundary copy = alt;
system_graph.addNode({new_node_key, {copy,
        node_data.position[0], node_data.position[1], node_data.position[2]} });
} else if constexpr (std::is_same_v<T, Microtubule::Intermediate>) {
Microtubule::Intermediate copy = alt;
system_graph.addNode({new_node_key, {copy,
        node_data.position[0], node_data.position[1], node_data.position[2]} });
} else if constexpr (std::is_same_v<T, Microtubule::Positive>) {
Microtubule::Positive copy = alt;
system_graph.addNode({new_node_key, {copy,
        node_data.position[0], node_data.position[1], node_data.position[2]} });
} else if constexpr (std::is_same_v<T, Microtubule::Negative>) {
Microtubule::Negative copy = alt;
system_graph.addNode({new_node_key, {copy,
        node_data.position[0], node_data.position[1], node_data.position[2]} });
} else if constexpr (std::is_same_v<T, Microtubule::Nucleator>) {
Microtubule::Nucleator copy = alt;
system_graph.addNode({new_node_key, {copy,
        node_data.position[0], node_data.position[1], node_data.position[2]} });
} else if constexpr (std::is_same_v<T, Microtubule::CellBoundary>) {
Microtubule::CellBoundary copy = alt;
system_graph.addNode({new_node_key, {copy,
        node_data.position[0], node_data.position[1], node_data.position[2]} });
} else if constexpr (std::is_same_v<T, Microtubule::Zipper>) {
Microtubule::Zipper copy = alt;
system_graph.addNode({new_node_key, {copy,
        node_data.position[0], node_data.position[1], node_data.position[2]} });
} else if constexpr (std::is_same_v<T, Microtubule::Junction>) {
Microtubule::Junction copy = alt;
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
             settings.TOTAL_TIME = 40;
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
            	Microtubule::start_to_node(gamma, this->system_graph, settings);
	Microtubule::create_nucleator_grid_x(gamma, this->system_graph, settings);
	Microtubule::make_boundary_top(gamma, this->system_graph, settings);
	Microtubule::make_boundary_bottom(gamma, this->system_graph, settings);
	Microtubule::make_boundary_right(gamma, this->system_graph, settings);
	Microtubule::create_nucleator_grid_y(gamma, this->system_graph, settings);
	Microtubule::make_boundary_left(gamma, this->system_graph, settings);

		 if (load_initial_state == false) {this->add_default_type(this->system_graph,
                        geoplex2D,
                        settings,
                        this->gen);} else {load_graph(gamma, this->system_graph, settings, this->gen, 
                        this->initial_state_filename);}
}
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
                        { std::ofstream ts(results_dir_name+"/timesteps.csv", std::ios::trunc);
                          ts << "step,time\n" << step << "," << (static_cast<double>(step) * settings.DELTA) << "\n"; }
                    }
                    if( step != 0 & step % 1 == 0)
                    {
                        std::string title = results_dir_name+"/simulation_step_";
                        DGGML::VtkFileWriterComplete<graph_type> vtk_writer;

                        auto attr = this->get_type_attributes(this->system_graph);
                        vtk_writer.set_extra_point_data(attr);

                        vtk_writer.save(system_graph, title+std::to_string(step));
                        collect(step);
                        { std::ofstream ts(results_dir_name+"/timesteps.csv", std::ios::app);
                          ts << step << "," << (static_cast<double>(step) * settings.DELTA) << "\n"; }}
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
            Microtubule::graph_type tg;
            Microtubule::StartType tmp = Microtubule::StartType{};

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
// Register column names
std::vector<std::vector<double>> columns(col_names.size());
for (auto &col : columns) { col.resize(N, nan_value); }
std::unordered_map<std::string, std::size_t> col_idx;
for (std::size_t i = 0; i < col_names.size(); ++i) { col_idx[col_names[i]] = i; }
    auto set = [&](const std::string &col_name, std::size_t row, double value) {
        columns[col_idx.at(col_name)][row] = value;
    };

std::size_t row = 0;
for (auto it = system_graph.node_list_begin(); it != system_graph.node_list_end(); ++it, ++row) {
auto &n = it->second.getData();
std::visit([&](auto &alt) {
using T = std::decay_t<decltype(alt)>;
if constexpr (std::is_same_v<T, Microtubule::CellBoundary>) {
}
if constexpr (std::is_same_v<T, Microtubule::Nucleator>) {
}
if constexpr (std::is_same_v<T, Microtubule::Zipper>) {
}
if constexpr (std::is_same_v<T, Microtubule::Intermediate>) {
}
if constexpr (std::is_same_v<T, Microtubule::Negative>) {
}
if constexpr (std::is_same_v<T, Microtubule::StartType>) {
}
if constexpr (std::is_same_v<T, Microtubule::Junction>) {
}
if constexpr (std::is_same_v<T, Microtubule::Positive>) {
}
if constexpr (std::is_same_v<T, Microtubule::Boundary>) {
}
}, n.data);
}
for (std::size_t i = 0; i < col_names.size(); ++i) {
        type_attributes.emplace_back(col_names[i], std::move(columns[i]));
    }
return type_attributes;}
void load_graph(DGGML::Grammar<Microtubule::graph_type> &grammar,
               Microtubule::graph_type &system_graph,
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

        SpatialNode3D<Microtubule::StartType, Microtubule::Boundary, Microtubule::Intermediate, Microtubule::Positive, Microtubule::Negative, Microtubule::Nucleator, Microtubule::CellBoundary, Microtubule::Zipper, Microtubule::Junction> node_data;
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
if constexpr (std::is_same_v<T, Microtubule::StartType>) {
Microtubule::StartType copy = alt;
system_graph.addNode({new_node_key, {copy,
        node_data.position[0], node_data.position[1], node_data.position[2]} });
} else if constexpr (std::is_same_v<T, Microtubule::Boundary>) {
Microtubule::Boundary copy = alt;
system_graph.addNode({new_node_key, {copy,
        node_data.position[0], node_data.position[1], node_data.position[2]} });
} else if constexpr (std::is_same_v<T, Microtubule::Intermediate>) {
Microtubule::Intermediate copy = alt;
system_graph.addNode({new_node_key, {copy,
        node_data.position[0], node_data.position[1], node_data.position[2]} });
} else if constexpr (std::is_same_v<T, Microtubule::Positive>) {
Microtubule::Positive copy = alt;
system_graph.addNode({new_node_key, {copy,
        node_data.position[0], node_data.position[1], node_data.position[2]} });
} else if constexpr (std::is_same_v<T, Microtubule::Negative>) {
Microtubule::Negative copy = alt;
system_graph.addNode({new_node_key, {copy,
        node_data.position[0], node_data.position[1], node_data.position[2]} });
} else if constexpr (std::is_same_v<T, Microtubule::Nucleator>) {
Microtubule::Nucleator copy = alt;
system_graph.addNode({new_node_key, {copy,
        node_data.position[0], node_data.position[1], node_data.position[2]} });
} else if constexpr (std::is_same_v<T, Microtubule::CellBoundary>) {
Microtubule::CellBoundary copy = alt;
system_graph.addNode({new_node_key, {copy,
        node_data.position[0], node_data.position[1], node_data.position[2]} });
} else if constexpr (std::is_same_v<T, Microtubule::Zipper>) {
Microtubule::Zipper copy = alt;
system_graph.addNode({new_node_key, {copy,
        node_data.position[0], node_data.position[1], node_data.position[2]} });
} else if constexpr (std::is_same_v<T, Microtubule::Junction>) {
Microtubule::Junction copy = alt;
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
            	Microtubule::ode_mt_growth(gamma, this->system_graph, settings);
	Microtubule::creation_case1(gamma, this->system_graph, settings);
	Microtubule::boundary_clamp(gamma, this->system_graph, settings);
	Microtubule::stochastic_mt_growth(gamma, this->system_graph, settings);
	Microtubule::mt_stochastic_retraction(gamma, this->system_graph, settings);
	Microtubule::catastrophe2_case1(gamma, this->system_graph, settings);
	Microtubule::destruction_case2(gamma, this->system_graph, settings);
	Microtubule::mt_ode_retraction(gamma, this->system_graph, settings);
	Microtubule::boundary_catastrophe1(gamma, this->system_graph, settings);

		 if (load_initial_state == false) {this->add_default_type(this->system_graph,
                        geoplex2D,
                        settings,
                        this->gen);} else {load_graph(gamma, this->system_graph, settings, this->gen, 
                        this->initial_state_filename);}
}
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
                        { std::ofstream ts(results_dir_name+"/timesteps.csv", std::ios::trunc);
                          ts << "step,time\n" << step << "," << (static_cast<double>(step) * settings.DELTA) << "\n"; }
                    }
                    if( step != 0 & step % 1 == 0)
                    {
                        std::string title = results_dir_name+"/simulation_step_";
                        DGGML::VtkFileWriterComplete<graph_type> vtk_writer;

                        auto attr = this->get_type_attributes(this->system_graph);
                        vtk_writer.set_extra_point_data(attr);

                        vtk_writer.save(system_graph, title+std::to_string(step));
                        collect(step);
                        { std::ofstream ts(results_dir_name+"/timesteps.csv", std::ios::app);
                          ts << step << "," << (static_cast<double>(step) * settings.DELTA) << "\n"; }}
 if (save_system_graph) {save_graph(gamma, this->system_graph, settings, results_dir_name+"/simulation_state_"+std::to_string(step)+".bin");
save_graph(gamma, this->system_graph, settings, results_dir_name+"/simulation_state_latest.bin");}
 }
}; // end Model_1

};
#endif
