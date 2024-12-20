using Mustache

"""
Generates Simulator
"""

main_tmp = mt"""
#include <iostream>
#include <chrono>

#include "DggFactory.hpp"
#include "{{header_name}}.h"
#include "simdjson.h"

int main(int argc, char* argv[])
{

    // Check if filename argument is provided
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <filename>" << std::endl;
        return 1;
    }
    std::cout << "Running {{sim_name}} Dynamic Graph Grammar Simulator\n";

    std::string filename = argv[1];
    // The user defines the parser, loads, and parses configuration file.
    simdjson::ondemand::parser parser;
    simdjson::padded_string json = simdjson::padded_string::load(filename);
    simdjson::ondemand::document settings_file = parser.iterate(json);

    auto start = std::chrono::high_resolution_clock::now();
    DGGML::SimulatorInterface<CMA::cmaModel> cma_simulation;
    CMA::cmaModel experiment1;
    experiment1.set_parameters(settings_file);

    //return 0;
    cma_simulation.setModel(experiment1);
    cma_simulation.simulate();
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration =
            std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
    std::cout << "\n\nSimulation took " << duration.count() / 1000.0 << " seconds\n";

}
"""

function generate_main(header_name, sim_name)
    sim_name = "Microtuble";
    header_name = "model";

    data = Dict(
             "sim_name" => sim_name,
             "header_name" => header_name
         );

    Mustache.render(main_tmp, data)
end
