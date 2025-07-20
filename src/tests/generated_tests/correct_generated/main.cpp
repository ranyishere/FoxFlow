#include <iostream>
#include <chrono>

#include "DggFactory.hpp"

#include "model.h"
#include "simdjson.h"

int main(int argc, char **argv) {

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <json_file>" << std::endl;
        return 1;
    }

    // Load the model from the JSON file
    std::string file_name = argv[1];

    // simdjson::ondemand::parser parser;
    // simdjson::padded_string json = simdjson::padded_string::load(file_name);
    // simdjson::ondemand::document doc = parser.iterate(json);

    DGGML::SimulatorInterface<Particles::Model> model_simulator;

    Particles::Model experiment;
    // experiment.set_parameters(doc);

    model_simulator.setModel(experiment);
    model_simulator.simulate();

    std::cout << "========== end ============" << std::endl;
    // model_simulator.getModel().print();
    // model_simulator.getModel().printParticles();
    model_simulator.simRunner->model->gamma.print();
    // model_simulator.model->gamma.print();
    std::cout << "========== end ============" << std::endl;

    /*
	std::cout << "Running simulation" << std::endl;
	std::string filename = settings.json;
	std::string json = simdjson::get_json(filename);
	DGGML::Model model(json);
	model.run();
    */
	return 0;
}
