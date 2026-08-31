#include <iostream>
#include <chrono>
#include "DggFactory.hpp"
#include "model.h"
#include "functions.h"
#include "simdjson.h"
int main(int argc, char **argv) {
	 if (argc != 2) {
		 std::cerr << "Usage: " << argv[0] << " <json_file>" << std::endl;
		 return 1;
		}
	 std::string filename = argv[1];

	 // ===== Stage 0 =====
	 {
		 DGGML::SimulatorInterface<Microtubule::Model_0> simulator_0;
		 Microtubule::Model_0 model_0;
		 simulator_0.setModel(model_0);
		 simulator_0.simulate();
	 }

	 // ===== Stage 1 =====
	 {
		 DGGML::SimulatorInterface<Microtubule::Model_1> simulator_1;
		 Microtubule::Model_1 model_1;
		 simulator_1.setModel(model_1);
		 simulator_1.simulate();
	 }
	return 0;
}
