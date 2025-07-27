#include <iostream>
#include <chrono>
#include "DggFactory.hpp"
#include "model.h"
#include "simdjson.h"
int main(int argc, char **argv) {
	 if (argc != 2) {
		 std::cerr << "Usage: " << argv[0] << " <json_file>" << std::endl;
		 return 1;
		}
	 std::string filename = argv[1];
	 DGGML::SimulatorInterface<Particles::Model> model_simulator;
	 Particles::Model current_model;
	 model_simulator.setModel(current_model);
	 model_simulator.simulate();
	return 0;
}
