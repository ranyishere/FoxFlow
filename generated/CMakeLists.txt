configure_file("settings.json" "settings.json" COPYONLY)
# include_directories(PlantModel)

# add_executable(main main.cpp ../../thirdparty/simdjson/simdjson.cpp)
add_executable(main main.cpp ../include/DGGML/thirdparty/simdjson/simdjson.cpp)
# add_executable(configuration_file_generator configuration_file_generator.cpp ../../thirdparty/simdjson/simdjson.cpp)

# target_link_libraries(mt_dgg_simulator PRIVATE vtu11::vtu11)
# target_link_libraries(mt_dgg_simulator PUBLIC SUNDIALS::arkode SUNDIALS::nvecserial)
target_link_libraries(main PRIVATE vtu11::vtu11)
target_link_libraries(main PUBLIC SUNDIALS::arkode SUNDIALS::nvecserial)
