using Mustache

include("./templates/sim_main.jl")
include("./templates/sim_model.jl")
include("./templates/sim_params.jl")
include("./templates/sim_rules.jl")
include("./templates/sim_settings.jl")
include("./templates/sim_types.jl")

include("./templates/sim_utils.jl")


function generate_im(grammars, symbol_table)
    """
    Generates C++ Code
    """

    header_name = "cmaModel";
    sim_name = "Microtuble";

    data = Dict();
    main_file_data = generate_main(header_name, sim_name);

    # println("main: ", main_file_data)
    model_file_data = generate_model(data);
    parameters_file_data = generate_parameters(data);
    types_file_data = generate_types(data);
    rules_file_data = generate_rules(data);
    settings_file_data = generate_settings(data);
    utils_file_data = generate_utils(data);

    Dict(
         "main" => main_file_data,
         "model" => model_file_data,
         "parameters" => parameters_file_data,
         "types" => types_file_data,
         "rules" => rules_file_data,
         "settings" => settings_file_data,
         "utils" => utils_file_data
        )
end
