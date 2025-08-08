# include("../main_parser_op.jl")
# include("../ast_nodes.jl")
# include("ir_builder.jl")
#
include("ir_rule_generation.jl")
include("ir_parameter_section.jl")
include("ir_types_section.jl")
include("ir_models_section.jl")

import .IRRuleGeneration: ir_rules_section!
import .IRBuildUtils: emit, build, IRBuilder
import ..AstNodes: IntegerNode, FloatNode, IdentifierNode, BinaryOpNode, GroupNode, CallNode
import ..Tokens: IntegerToken, FloatToken, PositionToken, LiteralToken
import .IRUtils: get_value, convert_type_name, write_file
using OrderedCollections

# include("ir_types_section.jl")
# FIXME: Spaces should not be relevant.
# FIXME: need to fix that the rules are updating the wrong nodes.
#    seems like node counts are wrong.

using UUIDs

symbol_tables = OrderedDict()
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
        inner_value = ir_value(ast.expr)
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

        while length(ast.token) > 0
            each_token = popfirst!(ast.token)
            if isa(each_token, IntegerNode)
                int_name = "fflow_"*string(UUIDs.uuid4())[1:6]

                push!(params,"\t\tint $(int_name);\n")
                push!(param_names,"$(int_name)")
                push!(param_types, "int")
                push!(is_list_params, (false, nothing, nothing))

            elseif isa(each_token, FloatNode)

                float_name = "fflow_"*string(UUIDs.uuid4())[1:6]

                push!(params,"\t\tdouble $(float_name);\n")
                push!(param_names,"$(float_name)")
                push!(param_types, "double")
                push!(is_list_params, (false, nothing, nothing))

            elseif isa(each_token, IdentifierNode)

                # TODO: Must have a parameter node afterwards and parse it

                name = get_value(each_token)
                if name == "FixedList"
                    # A fixed list has two parameters
                    # the first one is the size
                    # the second one is the type

                    list_params = popfirst!(ast.token)
                    list_size = list_params.token[1]
                    list_type = list_params.token[2]

                    # (Is a list, size of list, list_type)
                    push!(is_list_params, (true, list_size, list_type))

                    list_name ="fflow_"* string(UUIDs.uuid4())[1:6]

                    dim = list_size.token.position.value

                    push!(params, "\t\tfloat $(list_name)[$(dim)];\n")
                    push!(param_names,"$(list_name)")
                    push!(param_types, list_type)
                else
                    # TODO: Handle if the Identifer is another FoxFlow type
                    # Make sure it also has a parameter node
                end
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

function ir_type_instance(ast, define_type=false)
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

                list_size = is_list[ix][2].token.position.value
                list_type = is_list[ix][3].token.position.value
                list_size_int = parse(Int64, list_size)

                for i in 1:list_size_int
                    symbol_tables[type_name][param_count] = (list_type, name*"[$(i-1)]", is_list[ix])
                    param_count += 1
                end

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

        # It's a literal
        # begin_struct = ["\t"*type_class_name*" $(type_name) "*"= $(ast.value.position.value);"]
        begin_struct = ["\t"*ir_type*" $(type_name) "*"= $(ast.value.position.value);"]
    elseif isa(ast.value, IntegerNode)
        # It's an integer
        # ir_type = "int"

    ir_type = convert_type_name(get_value(ast.type.name))

        begin_struct = ["\t"*ir_type*" $(type_name) "*"= $(ast.value.token.position.value);"]
    elseif isa(ast.value, BinaryOpNode)
        # It's a binary operation
        expression = ir_value(ast.value)

        ir_type = convert_type_name(type_class_name)

        begin_struct = ["\t"*ir_type*" $(type_name) "*"= $(expression);"]

    elseif isa(ast.value, IdentifierNode)
        # It's an identifier
        ir_type = convert_type_name(type_class_name)
        begin_struct = ["\t"*ir_type*" $(type_name) "*"= $(ast.value.token.position.value);"]
    end

    push!(begin_struct,"\n")
    join(begin_struct)
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

        };

    "
    return model_boundary
end

function ir_main(name_space)

    ir_header = [
                "#include <iostream>\n",
                "#include <chrono>\n",
                "#include \"DggFactory.hpp\"\n",
                "#include \"model.h\"\n",
                "#include \"simdjson.h\"\n",

                "int main(int argc, char **argv) {\n"
                ]

    ir_body = [
               "\t if (argc != 2) {\n",
               "\t\t std::cerr << \"Usage: \" << argv[0] << \" <json_file>\" << std::endl;\n",
               "\t\t return 1;\n\t\t}\n",
               "\t std::string filename = argv[1];\n",
               "\t DGGML::SimulatorInterface<$name_space::Model> model_simulator;\n",
               "\t $name_space::Model current_model;\n",
               "\t model_simulator.setModel(current_model);\n",
               "\t model_simulator.simulate();\n",
            ]

    # ir_show_gamma = [
                     # "\t std::cout << Grammar << std::endl;\n",
                     # "\t model_simulator.get_gamma();\n",
                # ]


    ir_footer = [
               "\treturn 0;\n"
                 "}\n"
     ]

    # check_main = join([ir_header; ir_body; ir_show_gamma; ir_footer])
    check_main = join([ir_header; ir_body; ir_footer])

    check_main
end


# NOTE: Solving rule takes in:
# explicit SolvingRule(std::string rname, GraphType& lhs_graph, GraphType& rhs_graph, std::size_t num_eq, initial_condition_t&& ic, solving_t&& ode)  
#

# TODO: Get branching to work.
# TODO: Need to handle settings.json
# # FIXME: Rules are not firing anymore...
function generate_ir()
    """
    Generate IR
    """

    name_space = "Particles"

    # test_folder = "particle_sim"
    test_folder = "particle_sim_branching"

    base = "../tests/generated_tests/generated_2/"

    test_base = "../tests/"

    tokens_params = tokenize_file(test_base*"$test_folder/params.fflow")
    ast_params = parse_file!(tokens_params)

    println("Generating Params")
    generated_params = ir_parameter_section(ast_params[1])

    write_file(base*"parameters.h",
            generated_params)

    tokens_types = tokenize_file(test_base*"$test_folder/types.fflow")
    ast_types = parse_file!(tokens_types)
    println("Generating Types")
    generated_types = ir_types_section(ast_types[1])
    write_file(base*"types.h",
               generated_types)

    tokens_rules = tokenize_file(test_base*"$test_folder/rules.fflow")
    ast_rules = parse_file!(tokens_rules)[1]

    println("Generating Rules")
    ir_rules_section = ir_rules_section!(ast_rules,
                                         rules_table, symbol_tables,
                                         propensity_table, type_namespace)

    write_file(base*"rules.h", ir_rules_section)

    println("Generating Main")
    main_ir = ir_main(name_space)
    write_file(
            base*"main.cpp",
            main_ir
        )

    println("Generating Model")
    ir_models = ir_models_section(nothing)
    write_file(base*"model.h", ir_models)

    println("Done")

end

generate_ir()
