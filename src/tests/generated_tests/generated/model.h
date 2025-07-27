#ifndef DGGML_MODELS_HPP
#define DGGML_MODELS_HPP
#include <fstream>
#include "DGGML.h"
#include "rules.h"
#include "simdjson.h"
#include "ExpandedComplex2D.hpp"
#include "YAGL_Algorithms.hpp"
#include "parameters.h"


namespace Particles {

    using graph_grammar_t = DGGML::Grammar<Particles::graph_type>;
    class Model : public DGGML::Model<graph_grammar_t> {

        public:

            Parameters settings;

            // template<typename GraphType, typename CplexType, typename ParamType, typename GenType>
            /**
             *  Add a boundary to the graph
             *  @param graph The graph to which the boundary is added
             *  @param system_graph The system graph
             *  @param settings The settings for the simulation
             *  @param gen The random number generator
             * */
            /*
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
                node_type node_n = {curr_key, {Particles::Boundary{}, px, py, 0.0}};
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
                node_type node_n = {curr_key, {Particles::Boundary{}, px, py, 0.0}};
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
                node_type node_n = {curr_key, {Particles::Boundary{}, px, py, 0.0}};
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
                node_type node_n = {curr_key, {Particles::Boundary{}, px, py, 0.0}};
                graph.addNode(node_n);
                //connect to previous node
                graph.addEdge(prev_key, curr_key);
                prev_key = curr_key;
            }
            //complete the loop with the first
            graph.addEdge(prev_key, first_key);

            auto tmp = 0.075;

            // TODO: How do you determine the max and the min??
            // int max_particles = 10;
            // int min_particles = 1;

            const double min_x = 0.0;
            const double min_y = 0.0;
            // const double max_x = 4.98;
            const double max_x = settings.CELL_NX-1;
            // const double max_y = 4.98;
            const double max_y = settings.CELL_NY-1;
            const std::size_t max_nx = 16;
            const std::size_t max_ny = 16;

            std::size_t max_particles = reaction_grid.totalNumCells();
            //step 1 create an ordered number of selectable cells
            std::vector<std::size_t> selected_cells;
            for (int i = max_particles/2; i < max_particles; i++) {
                selected_cells.push_back(i);
            }

            std::size_t segments = 3;
            // int num_mt = 100;
            int num_mt = 50;

            for (auto i = 0; i < num_mt; i++) {
                double x_c, y_c;

                x_c = i*0.01;
                y_c = i*0.01;

                reaction_grid.cardinalCellToPoint(x_c, y_c, selected_cells[i]);

                Particles::graph_type tg;
                Particles::ParticleNodeCreator tmp = Particles::ParticleNodeCreator{};

                tmp.a4c6400 = 10;

            // Initial placement of particle node creator
                tmp.ad89e53[0] = max_x/2 + 1; // X
                tmp.ad89e53[1] = max_y/2 + i+1; // Y

                node_type oof_node = {gen.get_key(),//i*segments,
                                    {tmp, 
                                     x_c, y_c, 0.0}
                            };

                graph.addNode(oof_node);
            };

        };
    */

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

            // Fetching Reaction Grid
            std::size_t segments = 3;
            // int num_mt = 100;
            int num_mt = 50;

            for (auto i = 0; i < num_mt; i++) {
                double x_c, y_c;

                x_c = i*0.01;
                y_c = i*0.01;

                reaction_grid.cardinalCellToPoint(x_c, y_c, selected_cells[i]);

                Particles::graph_type tg;
                Particles::ParticleNodeCreator tmp = Particles::ParticleNodeCreator{};

                tmp.a4c6400 = 10;

                // Initial placement of particle node creator
                tmp.ad89e53[0] = max_x/2 + 1; // X
                tmp.ad89e53[1] = max_y/2 + i+1; // Y

                node_type oof_node = {gen.get_key(),//i*segments,
                                    {tmp, 
                                     x_c, y_c, 0.0}
                            };

                graph.addNode(oof_node);
            };

                /*
                this->add_boundary(
                        reaction_grid,
                        graph,
                        cplex,
                        // system_graph,
                        settings, gen);
                */

            };

            void initialize() override {
                
                // Default Rule
                // Initialize the model
                particle_rules::change_to_particle(
                            gamma, this->system_graph, settings
                );

                geoplex2D.init(
                        settings.CELL_NX,
                        settings.CELL_NY,
                        settings.CELL_DX,
                        settings.CELL_DY,
                        settings.MAXIMAL_REACTION_RADIUS
                );

                this->add_default_type(this->system_graph,
                    geoplex2D,
                    settings,
                    this->gen);

                // std::string filename = settings.json;
                //change_to_particle();
            }

            void checkpoint(std::size_t step) override {

                // auto results_dir_name = settings.RESULTS_DIR;
                std::string results_dir_name = "my_results";

                if(step == 0)
                {
                    std::cout << "Running the checkpoint for step " << step << "\n";
                    // Create the local save directory
                    std::cout << "Cleaning up old results folder if it exists and creating a new one\n";
                    std::filesystem::remove_all(results_dir_name);
                    std::filesystem::create_directory(results_dir_name);

                    // Save expanded cell complex graph
                    //DGGML::VtkFileWriter<typename DGGML::ExpandedComplex2D<>::types::graph_type> writer;
                    //writer.save(model->geoplex2D.getGraph(), results_dir_name+"/factory_geoplex");
                    DGGML::GridFileWriter grid_writer;
                    grid_writer.save({geoplex2D.reaction_grid,geoplex2D.dim_label},
                                     results_dir_name+"/expanded_cell_complex");

                    std::string title = results_dir_name+"/simulation_step_";
                    std::cout << "Saving the initial state of the system graph\n";
                    DGGML::VtkFileWriter<graph_type> vtk_writer;
                    vtk_writer.save(system_graph, title+std::to_string(step));
                    collect(step);
                }
                if( step != 0 & step % 5 == 0)
                {
                    std::cout << "Running the checkpoint for step " << step << "\n";
                    // Create the local save directory
                    // std::filesystem::remove_all(results_dir_name);
                    // std::filesystem::create_directory(results_dir_name);

                    DGGML::GridFileWriter grid_writer;
                    grid_writer.save({geoplex2D.reaction_grid,geoplex2D.dim_label},
                                     results_dir_name+"/expanded_cell_complex");

                    std::string title = results_dir_name+"/simulation_step_";
                    std::cout << "Saving the initial state of the system graph\n";
                    DGGML::VtkFileWriter<graph_type> vtk_writer;
                    std::cout << "title: " << title + std::to_string(step) << std::endl;
                    vtk_writer.save(system_graph, title+std::to_string(step));
                    collect(step);
                }

            }

            void set_parameters(const simdjson::ondemand::document &interface) {
                // settings.set_parameters(interface);
                // Set parameters from the JSON object
                // Example: obj["param_name"].get_to(param_name);
            };

    };
}
#endif
