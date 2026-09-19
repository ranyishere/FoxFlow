#ifndef DGGML_MODELS_HPP
#define DGGML_MODELS_HPP
#include <fstream>
#include "DGGML.h"
#include "rules.h"
#include "simdjson.h"
#include "ExpandedComplex2D.hpp"
#include "YAGL_Algorithms.hpp"
#include "parameters.h"
#include <unistd.h> // For sleep() and usleep()


namespace Microtubule {

    /*
    // Assuming cplex2D is your cell complex
auto& cplex_graph = cplex2D.getGraph();

for (auto& pair : graph.getNodeSetRef()) {
    auto& node = pair.second;
    double x = node.getData().position[0];
    double y = node.getData().position[1];

    int ic, jc;
    grid.locatePoint(x, y, ic, jc);

    // Map (ic, jc) to cell complex key (depends on your implementation)
    std::size_t cell_key = grid.cardinalCellIndex(ic, jc);

    // Now try to find the complex node
    auto cell_node_it = cplex_graph.findNode(cell_key);
    if (cell_node_it != cplex_graph.node_list_end()) {
        bool is_interior = cell_node_it->second.getData().interior;
        std::cout << "Node at (" << x << "," << y << ") is "
                  << (is_interior ? "INTERIOR" : "EXTERIOR") << std::endl;
    } else {
        std::cout << "Node at (" << x << "," << y << ") is NOT IN THE CELL COMPLEX" << std::endl;
    }
}
*/

    // Suppose: grid is your CartesianGrid2D, graph is YAGL::Graph
void print_node_grid_assignments(DGGML::CartesianGrid2D& grid, Microtubule::graph_type& graph, auto geoplex2D) {

    auto& geoplex_graph = geoplex2D.getGraph();

    for (auto& pair : graph.getNodeSetRef()) {
        auto key = pair.first;
        
        auto& node = pair.second;


        auto& spatial = node.getData();


        double x = spatial.position[0]; // or spatial.position[0], whichever is correct
        double y = spatial.position[1];

        // double x = node.position[0];
        // double y = node.position[1];

        auto rx_grid_max = grid._nx * grid._dx;
        auto ry_grid_max = grid._ny * grid._dy;

        int ic, jc;


        grid.locatePoint(x, y, ic, jc);
        auto cardinal = grid.cardinalCellIndex(ic, jc);

        auto type = spatial.type;

        // Check if it is interior
        auto cell_node_it = geoplex_graph.findNode(cardinal);
        bool interior = false;

        if (cell_node_it != geoplex_graph.node_list_end()) {
            interior = cell_node_it->second.getData().interior;
            std::cout << "found: "<< interior << " cardinal: " << cardinal << " key: "<< key << "x: "<< x << " y: "<< y << std::endl;
        } else {
            std::cout << "Node " << key << " at (" << x << "," << y << ") maps to cell_key " << cardinal << " which is NOT IN THE CELL COMPLEX" << std::endl;
            // exit(0);
        }

        // auto interior = node.getData().interior;

        // bool in_bounds = (ic >= 0 && ic < rx_grid_max) && (jc >= 0 && jc < ry_grid_max);
        /*
        bool in_bounds = (ic >= 0 && x < rx_grid_max) && (jc >= 0 && y < ry_grid_max);
        if (!in_bounds) {
            std::cout << "rx_grid_max: " << rx_grid_max << ", ry_grid_max: " << ry_grid_max << "\n";
            std::cout << "Node " << key << " at (" << x << "," << y << ") is OUT OF GRID! (ic=" << ic << ", jc=" << jc << ")\n";
            exit(0);
        } else {
            std::cout << "type " << type << " Node " << key << " at (" << x << "," << y << ") is binned to (ic=" << ic << ", jc=" << jc << ")\n";
        }
        */

        // sleep(10);
    }

}

    /*
    void print_node_grid_assignments(DGGML::CartesianGrid2D& reaction_grid, Microtubule::graph_type& graph) {
        for (const auto& [key, node] : graph.numNodes()) {
            double x = node.position[0];
            double y = node.position[1];
            int ic, jc;
            bool assigned = reaction_grid.locatePoint(x, y, ic, jc); // Adjust signature as needed
            if (!assigned) {
                std::cout << "Node " << key << " at (" << x << "," << y << ") NOT assigned to any grid cell!\n";
            } else {
                std::cout << "Node " << key << " at (" << x << "," << y << ") is in grid cell (" << ic << "," << jc << ")\n";
            }
        }
    }
    */

using graph_grammar_t = DGGML::Grammar<Microtubule::graph_type>;
class Model : public DGGML::Model<graph_grammar_t> {
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
            // const double max_x = 4.98;
            const double max_x = settings.CELL_NX-1;
            // const double max_y = 4.98;
            const double max_y = settings.CELL_NY-1;
            const std::size_t max_nx = 16;
            const std::size_t max_ny = 16;

            DGGML::CartesianGrid2D &reaction_grid = cplex.reaction_grid;
            std::size_t max_particles = reaction_grid.totalNumCells();

            //step 1 create an ordered number of selectable cells
            std::vector<std::size_t> selected_cells;
            for (int i = max_particles/2; i < max_particles; i++) {
                selected_cells.push_back(i);
            }

        auto centerCardinal = reaction_grid.cardinalCellIndex(
                        reaction_grid._nx/2, reaction_grid._ny/2
                    );

            double center_x, center_y;
            reaction_grid.cardinalCellToPoint(center_x, center_y, centerCardinal);
            Microtubule::graph_type tg;
            Microtubule::StartType tmp = Microtubule::StartType{};

            // Initial placement of particle node creator
            tmp.start_location[0] = center_x; // X
            tmp.start_location[1] = center_y; // Y
                                                       //
            node_type oof_node = {gen.get_key(),//i*segments,
                                {tmp, 
                                 center_x, center_y, 0.0}
                        };

            graph.addNode(oof_node);
            };
template<typename GraphType, typename CplexType, typename ParamType, typename GenType>
void add_boundary(
                    DGGML::CartesianGrid2D &reaction_grid,
                    GraphType &graph,
                    CplexType &cplex,
                    ParamType &settings, GenType &gen
            ) {

            using node_type = typename GraphType::node_type;
            using key_type = typename GraphType::key_type;
            key_type prev_key; //only the first step doesn't have a prev
            key_type first_key; //needed to complete the loop around the boundary
            // Creating Boundary
            for (auto i = 0; i < reaction_grid._nx; i++) {
                auto cardinal = reaction_grid.cardinalCellIndex(i, 0);
                double px, py;
                reaction_grid.cardinalCellToPoint(px, py, cardinal);
                key_type curr_key = gen.get_key();
                node_type node_n = {curr_key, {Microtubule::Boundary{}, px, py, 0.0}};
                graph.addNode(node_n);
                //connect to previous node or its the first
                if (i >= 1) graph.addEdge(prev_key, curr_key);
                else first_key = curr_key;
                prev_key = curr_key;
            }

            // Note: the previous gets carried over!
            // right side interior, bottom to top
            // Creates boundary
            for (auto j = 1; j < reaction_grid._ny - 1; j++) {
                auto cardinal = reaction_grid.cardinalCellIndex(reaction_grid._nx - 1, j);
                double px, py;
                reaction_grid.cardinalCellToPoint(px, py, cardinal);
                key_type curr_key = gen.get_key();
                node_type node_n = {curr_key, {Microtubule::Boundary{}, px, py, 0.0}};
                graph.addNode(node_n);
                //connect to previous node
                graph.addEdge(prev_key, curr_key);
                prev_key = curr_key;
            }

            // Note: the previous gets carried over!
            // top, right to left
            for (auto i = reaction_grid._nx - 1; i >= 0; i--) {
                auto cardinal = reaction_grid.cardinalCellIndex(i, reaction_grid._ny - 1);
                double px, py;
                reaction_grid.cardinalCellToPoint(px, py, cardinal);
                key_type curr_key = gen.get_key();
                node_type node_n = {curr_key, {Microtubule::Boundary{}, px, py, 0.0}};
                graph.addNode(node_n);
                //connect to previous node
                graph.addEdge(prev_key, curr_key);
                prev_key = curr_key;
            }
            //note: the previous gets carried over!
            //left side interior, bottom to top
            // Creating boundary
            for (auto j = reaction_grid._ny - 2; j > 0; j--) {
                auto cardinal = reaction_grid.cardinalCellIndex(0, j);
                double px, py;
                reaction_grid.cardinalCellToPoint(px, py, cardinal);
                key_type curr_key = gen.get_key();
                node_type node_n = {curr_key, {Microtubule::Boundary{}, px, py, 0.0}};
                graph.addNode(node_n);
                //connect to previous node
                graph.addEdge(prev_key, curr_key);
                prev_key = curr_key;
            }

            //complete the loop with the first
            graph.addEdge(prev_key, first_key);

        };

    void initialize() override {
        geoplex2D.init(
                        settings.CELL_NX,
                        settings.CELL_NY,
                        settings.CELL_DX,
                        settings.CELL_DY,
                        settings.MAXIMAL_REACTION_RADIUS
                );
    Microtubule::hit_boundary(gamma, this->system_graph, settings);
    Microtubule::start_node_dup(gamma, this->system_graph, settings);
    Microtubule::start_to_node_fracture(gamma, this->system_graph, settings);
    Microtubule::grow_fracture(gamma, this->system_graph, settings);

        this->add_default_type(this->system_graph,
                    geoplex2D,
                    settings,
                    this->gen);
        this->add_boundary(
                        geoplex2D.reaction_grid,
                        this->system_graph,
                        geoplex2D,
                        settings,
                        this->gen);}
        void checkpoint(std::size_t step) override {std::string results_dir_name = "my_results";if(step == 0)
                {
                    // Create the local save directory
                    std::filesystem::remove_all(results_dir_name);
                    std::filesystem::create_directory(results_dir_name);

                        DGGML::GridFileWriter grid_writer;
                    grid_writer.save({geoplex2D.reaction_grid,geoplex2D.dim_label},
                                     results_dir_name+"/expanded_cell_complex");

                    std::string title = results_dir_name+"/simulation_step_";
                    DGGML::VtkFileWriter<graph_type> vtk_writer;
                    vtk_writer.save(system_graph, title+std::to_string(step));
                    collect(step);
                }

                print_node_grid_assignments(geoplex2D.reaction_grid, this->system_graph, geoplex2D);
                if( step != 0 & step % 5 == 0)
                {

                    // sleep(10); // Pause for 1 second to ensure distinct timestamps if needed

                    DGGML::GridFileWriter grid_writer;
                    grid_writer.save({geoplex2D.reaction_grid,geoplex2D.dim_label},
                                     results_dir_name+"/expanded_cell_complex");

                    std::string title = results_dir_name+"/simulation_step_";
                    DGGML::VtkFileWriter<graph_type> vtk_writer;
                    vtk_writer.save(system_graph, title+std::to_string(step));
                    collect(step);
                }
                }};
};
#endif
