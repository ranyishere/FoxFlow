# include("../main_parser_op.jl")
# include("../ast_nodes.jl")
# include("ir_builder.jl")
# include("ir.jl")
# include("utils.jl")

include("ir_rule_generation.jl")
include("ir_parameter_section.jl")
include("ir_types_section.jl")
include("ir_models_section.jl")
include("ir_functions_section.jl")
include("ir_simulations_section.jl")

import .IRRuleGeneration: ir_rules_section!
import .IRBuildUtils: emit, build, IRBuilder, build_sameline
import ..AstNodes: IntegerNode, FloatNode, IdentifierNode, BinaryOpNode, GroupNode, CallNode, IndexAccessNode, ArrayLiteralNode, NamedParameterNode

import ..Tokens: IntegerToken, FloatToken, PositionToken, LiteralToken, ErrorToken, OperatorToken
import .IRUtils: get_value, convert_type_name, write_file
using OrderedCollections

# include("ir_types_section.jl")
# FIXME: Spaces should not be relevant.
# FIXME: need to fix that the rules are updating the wrong nodes.
#    seems like node counts are wrong.
using UUIDs

rules_table = Dict()
propensity_table = Dict()
type_namespace = ""
rule_namespace = ""

# Maps ix to the attribute
parameter_table = Dict()

function ir_value(ast)
    if ast == nothing
        return ""
    elseif isa(ast, IntegerNode)
        return ast.token.position.value
    elseif isa(ast, FloatNode)
        return ast.token.position.value
    elseif isa(ast, IdentifierNode)
        return get_value(ast)
    elseif isa(ast, BinaryOpNode)
        left_value = ir_value(ast.lhs)
        right_value = ir_value(ast.rhs)
        op = ast.expression.position.value
        return "$(left_value) $(op) $(right_value)"
    elseif isa(ast, GroupNode)
        inner_value = ir_value(ast.expression)
        return "($(inner_value))"
    elseif isa(ast, CallNode)
        func_name = get_value(ast.func)
        args = [ir_value(arg) for arg in ast.args]
        return "$(func_name)($(join(args, ", ")))"
    else
        # For unknown node types, attempt to get the token value if possible
        try
            return ast.token.position.value
        catch
            return string(ast)
        end
    end
end

function ir_fixed_list_define(ast)
    # Handles defining a fixed list.
    list_size = ast.token[1].position.value
    list_type = ast.token[2].position.value
    return list_size, list_type
end


function ir_parameter(ast)

    params = []
    param_names = []
    param_types = []

    # (is a list, size, type)
    is_list_params = []

    if isempty(ast.token)
        return nothing, nothing, nothing, nothing
    else

        # TODO: Finish popping this.
        # should be generating the parameters
        # you need to lookahead and grab the parameter related
        # to the identifier.
        # =============================

        copy_ast = deepcopy(ast.token)
        while length(copy_ast) > 0

            each_token = popfirst!(copy_ast)

	    println("each_token: ", each_token)

            if !(each_token isa NamedParameterNode)
                throw(ErrorException(
                    "Type attributes must be named. Found unnamed attribute: '$(get_value(each_token))'.\n" *
                    "  Expected:  MyType << name : FixedList<<3, Float>> >> : Type\n" *
                    "  Got:       MyType << FixedList<<3, Float>> >> : Type\n" *
                    "  Please add a name before the ':' for each attribute."
                ))
            end

            name = get_value(each_token.name)
            type_name = get_value(each_token.parameter.name)

            println("=========> type_name: ", type_name, " name: ", name, " each_token: ", each_token, "<=========")

            # if isa(type_name, IntegerNode)
            if type_name == "Integer"
                # #int_name = "fflow_"*string(UUIDs.uuid4())[1:6]
                int_name = name

                push!(params,"\t\tint $(int_name);\n")
                push!(param_names,"$(int_name)")
                push!(param_types, "int")
                push!(is_list_params, (false, nothing, nothing))

            # elseif isa(type_name, FloatNode)
            elseif type_name == "Float"

                # float_name = "fflow_"*string(UUIDs.uuid4())[1:6]
                float_name = name

                push!(params,"\t\tdouble $(float_name);\n")
                push!(param_names,"$(float_name)")
                push!(param_types, "double")

                push!(is_list_params, (false, nothing, nothing))

            elseif type_name ==  "FixedList"

                # TODO: Must have a parameter node afterwards and parse it
                # Assuming the that Fixed list is of type double

                # name = get_value(each_token.name)
                # type_name = get_value(each_token.parameter.name)

                if type_name == "FixedList"
                    # A fixed list has two parameters
                    # the first one is the size
                    # the second one is the type

                    # fl_params = popfirst!(copy_ast)
                    # param_name = get_value(fl_params.name)
                    # list_params = fl_params.parameter.parameter
                    # check_params = fl_params.parameter.name

                    list_params = each_token.parameter.parameter.token
                    println("list_params: ", list_params)

                    list_size = []
                    list_type = nothing
                    while length(list_params) > 0
                        cur_par = popfirst!(list_params)
                        val = cur_par.token.position.value
                        if val == "Integer" || val == "Float" 
                            list_type = cur_par
                            break
                        else
                            push!(list_size, cur_par.token.position.value)
                        end
                    end

                    println(
                        "list_size: ", list_size, " param_name: ", list_type, " name: ", name
                    )

                    # list_size = list_params.token[1]
                    # loop until you hit the list type
                    # list_type = list_params.token[2]

                    # (Is a list, size of list, list_type)
                    push!(is_list_params, (true, list_size, list_type))

                    # list_name ="fflow_"* string(UUIDs.uuid4())[1:6]
                    list_name = name

                    # dim = list_size.token.position.value
                    # join the sizes by comma
                    dim = join(list_size, ", ")

                    list_type_str = list_type.token.position.value

                    convert_type = convert_type_name(list_type_str, true)

                    # println("convert_type: ", convert_type)
                    # exit(0)

                    # push!(params, "\t\tdouble $(list_name)[$(dim)];\n")
                    ir = "\t\ttorch::Tensor $(list_name) = torch::zeros({$(dim)}, $convert_type);\n"

                    # println("ir: ", ir)
                    # exit(0)
                    # push!(params, "\t\tdouble $(list_name)[$(dim)];\n")
                    push!(params, ir)
                    push!(param_names,"$(list_name)")
                    push!(param_types, list_type)

                else
                    println("Each Token: ", each_token)
                    println("Unsupported parameter type for named parameter: $(type_name)")
                    throw(ErrorException("Only FixedList is supported as a list type for now."))
                    # TODO: Handle if the Identifer is another FoxFlow type
                    # Make sure it also has a parameter node
                end
            else
                throw(ErrorException("Unsupported parameter type: $(each_token)"))

            end

        end
        return params, param_names, param_types, is_list_params
    end
    params, param_names, param_types, is_list_params
end

function ir_type_class(ast)

    type_class_name = get_value(ast.name)
    type_class_param, param_names, param_types, is_list = ir_parameter(ast.parameter)

    return type_class_name
end

function ir_struct_like_array(param_names)
    variant_func = [
        # "\t\tstd::variant<int, float*> operator[](std::size_t index) const {\n",
        "\t\tvoid* operator[](std::size_t index) const {\n"
    ]

    for (ix, name) in enumerate(param_names)
        push!(variant_func,"\t\t\tif (index == $(ix-1)) return (void*)&$(name);\n")
    end

    variant_func = [variant_func;["\t\t\tthrow std::out_of_range(\"Index out of bounds\");\n","\t\t};\n"]]

    return variant_func
end


function ir_type_instance(ast, define_type=false, symbol_tables=nothing)
    """
    IR Type Instance
    """

    type_name = get_value(ast.name)

    type_parameters = []

    # Track for each type, what parameters does it belong to?
    if define_type == true
        symbol_tables[type_name] = Dict()
    end

    # Associate each parameter with the type and parameter position.
    type_parameters, param_names, param_types, is_list = ir_parameter(ast.parameter)

    if define_type == true
        param_count = 1
        for (ix, name) in enumerate(param_names)

            if is_list[ix][1] == true

                # Store tensor as a single entry with type "torch::Tensor"
                # instead of unrolling into individual elements.
                symbol_tables[type_name][param_count] = ("torch::Tensor", name, is_list[ix])
                param_count += 1

            else
                # It's not a list
                symbol_tables[type_name][param_count] = (param_types[ix], name, is_list[ix])
                param_count += 1
            end

        end
    end

    type_class_name = ir_type_class(ast.type)
    begin_struct = []

    # TODO: Handle identifier names
    if ast.value == nothing
        begin_struct = ["\tstruct $(type_name) : $(type_class_name) {\n"]
        begin_struct = [begin_struct;type_parameters]

        operator_like_array = ir_struct_like_array(param_names)
        begin_struct = [begin_struct;operator_like_array]

    elseif isa(ast.value, LiteralToken)

        ir_type = convert_type_name(type_class_name)
        begin_struct = ["\t"*ir_type*" $(type_name) "*"= $(ast.value.position.value);"]
    elseif isa(ast.value, IntegerNode)

        ir_type = convert_type_name(get_value(ast.type.name))
        begin_struct = ["\t"*ir_type*" $(type_name) "*"= $(ast.value.token.position.value);"]
    elseif isa(ast.value, BinaryOpNode)
        # It's a binary operation
        expression = ir_value(ast.value)
        ir_type = convert_type_name(type_class_name)
        begin_struct = ["\t"*ir_type*" $(type_name) "*"= $(expression);"]

    elseif isa(ast.value, IdentifierNode)
        ir_type = convert_type_name(type_class_name)
        begin_struct = ["\t"*ir_type*" $(type_name) "*"= $(ast.value.token.position.value);"]
    end

    push!(begin_struct,"\n")
    join(begin_struct), param_names
end


function create_spatial_node_3d(instance_names)

    code_line = ["\tSpatialNode3D<StartType,Boundary,",]
    for (ix, name) in enumerate(instance_names)
        if ix != length(instance_names)
            push!(code_line, name*",")
        else
            push!(code_line, name)
        end
    end

    code_line = [code_line; ">"]
    join(code_line)
end

function link_parameter_to_node(param, node_loc, params_var)
    """
    params_var[param_name] = (param_count, node_loc)

    Node location represents in your rule, either on the lhs or rhs, where it appears in order.

    param_count represents the parameter index as it appears in the rule.
    """

    # Count params and match them to nodes
    # if the param count goes beyond 
    # the number of attributes it has for that node
    # we should jump to the next node

    param_count = 1
    for tok in param.token


        # Need to check if the parameter belongs to an edge node
        if tok isa ParameterNode
            node_loc += 1
            link_parameter_to_node(tok, node_loc, params_var)
        elseif tok isa UndirectedTypeEdgeNode
            println("an undirected edge node")
        else
            param_name = get_value(tok)
            params_var[param_name] = (param_count, node_loc)
            param_count += 1

        end
    end
end

function ir_boundary!(ir_builder)
model_boundary = 
    "
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
                node_type node_n = {curr_key, {$type_namespace::Boundary{}, px, py, 0.0}};
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
                node_type node_n = {curr_key, {$type_namespace::Boundary{}, px, py, 0.0}};
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
                node_type node_n = {curr_key, {$type_namespace::Boundary{}, px, py, 0.0}};
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
                node_type node_n = {curr_key, {$type_namespace::Boundary{}, px, py, 0.0}};
                graph.addNode(node_n);
                //connect to previous node
                graph.addEdge(prev_key, curr_key);
                prev_key = curr_key;
            }

            //complete the loop with the first
            graph.addEdge(prev_key, first_key);

        };

    "
    return model_boundary
end

function ir_main(name_space, num_stages, function_table=nothing)

    # Check if any model-loaded functions exist in the function_table
    has_models = false
    func_namespace = ""
    if function_table !== nothing && haskey(function_table, "__func_namespace__")
        has_models = true
        func_namespace = function_table["__func_namespace__"]
    end

    ir_header = [
                "#include <iostream>\n",
                "#include <chrono>\n",
                "#include \"DggFactory.hpp\"\n",
                "#include \"model.h\"\n",
                "#include \"functions.h\"\n",
                "#include \"simdjson.h\"\n",

                "int main(int argc, char **argv) {\n"
                ]

    ir_model_load = []
    if has_models
        push!(ir_model_load, "\t $(func_namespace)::load_models();\n")
    end

    ir_body = [
               "\t if (argc != 2) {\n",
               "\t\t std::cerr << \"Usage: \" << argv[0] << \" <json_file>\" << std::endl;\n",
               "\t\t return 1;\n\t\t}\n",
               "\t std::string filename = argv[1];\n",
            ]

    # Generate sequential stage execution
    for stage_ix in 0:(num_stages-1)
        model_class = "Model_$stage_ix"
        push!(ir_body, "\n\t // ===== Stage $stage_ix =====\n")
        push!(ir_body, "\t {\n")
        push!(ir_body, "\t\t DGGML::SimulatorInterface<$name_space::$model_class> simulator_$stage_ix;\n")
        push!(ir_body, "\t\t $name_space::$model_class model_$stage_ix;\n")
        push!(ir_body, "\t\t simulator_$stage_ix.setModel(model_$stage_ix);\n")
        push!(ir_body, "\t\t simulator_$stage_ix.simulate();\n")
        push!(ir_body, "\t }\n")
    end

    # ir_show_gamma = [
                     # "\t std::cout << Grammar << std::endl;\n",
                     # "\t model_simulator.get_gamma();\n",
                # ]

    ir_footer = [
               "\treturn 0;\n"
                 "}\n"
     ]

    # check_main = join([ir_header; ir_body; ir_show_gamma; ir_footer])
    check_main = join([ir_header; ir_model_load; ir_body; ir_footer])

    check_main
end

function do_simulation()
    """
    Do Simulation — supports multiple sequential RunSimulation stages.
    """

    sim_table = Dict()
    symbol_tables = OrderedDict()
    function_table = OrderedDict()

    # test_folder = "fracture_network"
    # test_folder = "cell_competition"
    # test_folder = "microtubules"
    # test_folder = "lattice"
    # test_folder = "dissolution"
    test_folder = "neural_network"

    test_base = "../tests/"
    base = "../tests/generated_tests/generated_2/"

    # ===== Parse simulation section =====
    tokens_sim = tokenize_file(test_base*"$test_folder/simulation.fflow")
    ast = parse_file!(tokens_sim)

    # ir_simulation_section! returns the namespace name and a list of ir_run_sim dicts (one per RunSimulation call)
    name_space, all_stages = ir_simulation_section!(ast[1], sim_table)
    global type_namespace = name_space
    num_stages = length(all_stages)
    println("Simulation namespace: ", name_space)
    println("Found $num_stages simulation stage(s)")

    # ===== Parse functions (shared across stages) =====
    # TODO: If functions.fflow is missing just make it empty.
    tokens_funct = tokenize_file(test_base*"$test_folder/functions.fflow")
    funct_ast = parse_file!(tokens_funct)
    ir_functions = ir_function_section(funct_ast[1], function_table)

    propensity_table["function_table"] = function_table

    write_file(base*"functions.h", ir_functions)

    # ===== Parse params (shared — use from first stage) =====
    params_file = all_stages[1]["parameters"]
    tokens_params = tokenize_file(test_base*"$test_folder/$params_file")
    ast_params = parse_file!(tokens_params)

    println("Generating Params")
    generated_params = ir_parameter_section(ast_params[1], propensity_table)
    write_file(base*"parameters.h", generated_params)

    # ===== Parse types (shared — use from first stage) =====
    types_file = all_stages[1]["types"]
    tokens_types = tokenize_file(test_base*"$test_folder/$types_file")
    ast_types = parse_file!(tokens_types)

    println("Generating Types")
    generated_types = ir_types_section(ast_types[1], symbol_tables)
    write_file(base*"types.h", generated_types)

    # ===== Per-stage: parse rules, generate rules_N.h, collect rules_tables =====
    # We collect the rules_table for each stage so the model generator
    # can register only that stage's rules in its Model_N class.
    stage_rules_tables = []
    all_rules_includes = []

    for (stage_ix, stage_info) in enumerate(all_stages)
        stage_idx = stage_ix - 1  # 0-based for file naming

        # Reset the global rules_table for this stage
        global rules_table = Dict()
        global propensity_table
        # (propensity_table retains function_table from above)

        rules_file = stage_info["rules"]
        println("Stage $stage_idx: Generating Rules from $rules_file")

        tokens_rules = tokenize_file(test_base*"$test_folder/$rules_file")
        ast_rules = parse_file!(tokens_rules)[1]

        ir_rules_code = ir_rules_section!(ast_rules,
                                          rules_table, symbol_tables,
                                          propensity_table,
                                          type_namespace;
                                          stage_index=stage_idx)

        rules_filename = "rules_$stage_idx.h"
        write_file(base*rules_filename, ir_rules_code)

        # Save a snapshot of this stage's rules_table
        push!(stage_rules_tables, deepcopy(rules_table))
        push!(all_rules_includes, rules_filename)
    end

    # ===== Parse observables (optional — per stage via RunSimulation, or auto-detected) =====
    observable_section = nothing
    obs_file_from_stage = get(all_stages[1], "observables", nothing)
    obs_file = obs_file_from_stage !== nothing ? test_base*"$test_folder/$obs_file_from_stage" :
                   test_base*"$test_folder/observables.fflow"
    if isfile(obs_file)
        tokens_obs = tokenize_file(obs_file)
        ast_obs = parse_file!(tokens_obs)
        if !isempty(ast_obs)
            observable_section = ast_obs[1]
            println("Loaded observables section: $(get_value(observable_section.name))")
        end
    else
        println("No observables file found — tensor attributes will be skipped in VTK output")
    end

    # ===== Generate Model (one Model_N class per stage) =====
    println("Generating Model")
    ir_models = ir_models_section_multistage(all_stages, name_space, symbol_tables,
                                             stage_rules_tables, all_rules_includes;
                                             observable_section=observable_section)
    write_file(base*"model.h", ir_models)

    # ===== Grammar Entry Point =====
    println("Generating Main")
    main_ir = ir_main(name_space, num_stages, function_table)
    write_file(base*"main.cpp", main_ir)

    println("Done with IR Generation")
end

do_simulation()
