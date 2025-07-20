include("./main_parser_op.jl")

# FIXME: Spaces should not be relevant.
# FIXME: Should undersetand more than 1 rule, but at the moment only understand 1 rule.

using UUIDs

symbol_tables = Dict()
rule_table = Dict()

function get_value(ast)
    ast.token.position.value
end

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

        # Handle if parameter is an array
        # if name == "FixedList"

            # Handle array size 
            # push!(variant_func,"\t\t\tif (index == $(ix-1)) return $(name);\n")
        # else
            # Handle if parameter is a scalar
            # push!(variant_func,"\t\t\tif (index == $(ix-1)) return $(name);\n")
        # end
        # push!(variant_func,"\t\t\tif (index == $(ix-1)) return $(name);\n")
    end

    variant_func = [variant_func;["\t\t\tthrow std::out_of_range(\"Index out of bounds\");\n","\t\t};\n"]]

    return variant_func
end

function convert_type_name(type_name)
    """
    Convert type name to a proper c++ type
    """

    if type_name == "Float"
        return "double"
    elseif type_name == "Integer"
        return "int"
    else
        return type_name
    end
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
        for (ix, name) in enumerate(param_names)

            symbol_tables[type_name][ix] = (param_types[ix], name, is_list[ix])

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
        ir_type = "int"
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

    # ir_value(ast.value)
    # push!(begin_struct,"\n};\n")
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

function ir_types_section(ast)

    section_name = get_value(ast.name)

    section_data = [
        "#ifndef DGGML_$(section_name)_TYPES_HPP\n",
        "#define DGGML_$(section_name)_TYPES_HPP\n",
        "#include \"YAGL_Graph.hpp\" \n",
        "#include \"YAGL_Node.hpp\" \n",
        "#include \"SpatialData3D.hpp\" \n",
        "namespace $(section_name) {\n",
        "\tstruct Type {};\n",
        "\tstruct Boundary {};\n",
        "\tstruct StartType {
            float start_location[2];\n
       };\n"
    ]

    symbol_tables["StartType"] = Dict()
    start_type_info = (
                   true, 
                   IntegerNode(IntegerToken(PositionToken("2", -1, -1))), 
                   FloatNode(FloatToken(PositionToken("Float", -1, -1)))
               )

    symbol_tables["StartType"][1] = ("float", "start_location", start_type_info)

    type_instances_ast = ast.types
    type_instance_names = []
    for type_inst_ast in type_instances_ast

        cur_type_inst = ir_type_instance(type_inst_ast, true)
        cur_type_inst = cur_type_inst*"};"
        type_name = get_value(type_inst_ast.name)
        push!(type_instance_names, type_name)
        push!(section_data, cur_type_inst)
    end

    spatial_node_3d = create_spatial_node_3d(type_instance_names)

    bottom_section_data = ["\tusing graph_type = YAGL::Graph<std::size_t,",
                           spatial_node_3d, ">;\n",
                           "};\n", "#endif"
                        ]

    section_data = [section_data; bottom_section_data]

    final_str = join(section_data)
    final_str
end

function ir_grammar_parameter(ast)
    """
    IR Parameter
    """

    if isa(ast, TypeInstanceNode)
        grammar_parameter = ir_type_instance(ast)
        return grammar_parameter
    else
        println("did not expect this $(ast)")
    end
end

function ir_parameter_section(ast)
    parameter_section_data = [
                            "#ifndef DGGML_PARAMETERS_HPP\n",
                            "#define DGGML_PARAMETERS_HPP\n",
                            "#include <string>\n",
                            "#include \"simdjson.h\"\n"
                        ]

    section_name = get_value(ast.name)

    # TODO: REference Parameter section name
    # push!(parameter_section_data, "struct $(section_name) {\n")
    push!(parameter_section_data, "struct Parameters {\n")
    for (ix, param) in enumerate(ast.parameter_list)
        push!(parameter_section_data, ir_grammar_parameter(param))
    end

    push!(parameter_section_data, "};\n")

    push!(parameter_section_data, "#endif")

    join(parameter_section_data)
end

function write_file(file_name, ans)
    open("$(file_name)", "w+") do file
        write(file, ans)
    end
end

function link_parameter_to_node(param, node_loc, params_var)
    """
    params_var[param_name] = (param_count, node_loc)
    """

    param_count = 1
    for tok in param.token

        if tok isa ParameterNode
            node_loc += 1
            link_parameter_to_node(tok, node_loc, params_var)
        else
            param_name = get_value(tok)
            params_var[param_name] = (param_count, node_loc)
            param_count += 1

        end
    end
end

function ir_rules_section!(ast)
    """
    Generate Intermediate Rules
    for the Section

    # TODO: Handle edges in the graph
    """

    function get_assn_values(
            assn,values,lhs_params_var, lhs_var_track,
            seen, total_assn_ir
        )

        println("assn.value before =========> ", assn.value)
        for val in assn.value

            if (val.position.value in seen) || (!isa(val, IdentifierToken))
                push!(values, val.position.value)
                continue
            end

            # Create a new variable if it does not exist
            # within the current context
            #
            # lhs_node_loc = lhs_params_var[val.position.value][1]
            # lhs_param_loc = lhs_params_var[val.position.value][2]
            lhs_node_loc = lhs_params_var[val.position.value][2]
            lhs_param_loc = lhs_params_var[val.position.value][1]

            println("val.position.value: ", val.position.value)
            println("lhs_params_var: ", lhs_params_var)
            println("lhs_var_track: ", lhs_var_track)
            println("lhs_node_loc: ", lhs_node_loc)
            println("lhs_param_loc: ", lhs_param_loc)
            println("symbol_tables: ", symbol_tables)
            println("check oooooooof: ", lhs_var_track[lhs_node_loc][3]) # -> Goes to symbol tables
            println("oof2: ", symbol_tables[lhs_var_track[lhs_node_loc][3]])

            # Check if param_loc is referencing a parameter that is actually from an array

            # TODO: Handle if the parameter is a list

            # FIXME may 5th 2025: indexing Error 

            # Find where the variable belongs to and type
            param_type = symbol_tables[lhs_var_track[lhs_node_loc][3]][lhs_param_loc][1]
            lhs_node_type = lhs_var_track[lhs_node_loc][3]
            lhs_node_namespace = lhs_var_track[lhs_node_loc][2]

            # Get the attribute name
            param_name = symbol_tables[lhs_var_track[lhs_node_loc][3]][lhs_param_loc][2]

            param_info = symbol_tables[lhs_var_track[lhs_node_loc][3]][lhs_param_loc]

            println("param_info: ", param_info[3])

            total_params = []

            if param_info[3][1] == true
                param_name_arr = param_name*"["*string(lhs_param_loc)*"]"
                var_ir = "\t\t $param_type $(val.position.value) = std::get<$lhs_node_namespace::$lhs_node_type>(lhs[m1[$lhs_node_loc]].data).$param_name_arr;\n"
                # push!(total_params, tmp_var_ir)
            else
                var_ir = "\t\t $param_type $(val.position.value) = std::get<$lhs_node_namespace::$lhs_node_type>(lhs[m1[$lhs_node_loc]].data).$param_name;\n"
            end

            println("var_ir: =========> ", var_ir)

            # println("total_params: ", total_params)

            # exit(0)

            push!(total_assn_ir, var_ir)
            push!(seen, val.position.value)
            push!(values,
                  val.position.value)
        end
    end

    # Associate the user's parameter name with the real parameter name
    # based on location


    type_namespace = "Particles"
    rule_section_name = get_value(ast.name)

    rule_section_header = [
        "#ifndef DGGML_RULES_HPP\n",
        "#define DGGML_RULES_HPP\n",
        "#include \"types.h\"\n",
        "#include \"parameters.h\"\n",
        "namespace $rule_section_name {\n",
        "using GT = $type_namespace::graph_type;\n"
    ]

    # Create lhs and rhs
    # TODO: Support adding edges.
    # TODO: Support ode rules
    # TODO Reference Parameters section name instead of Parameters
    total_rules = []
    for rule in ast.rules_list

        rule_name = get_value(rule.name)
        rule_header = [
           "void $rule_name(DGGML::Grammar<$type_namespace::graph_type> &gamma,
           Particles::graph_type &system_graph,
           Parameters &settings) {\n"
        ]

        rule_table[type_namespace*"::"*rule_name] = true

        lhs = rule.lhs

        lhs_name = "$(rule_name)_lhs"
        lhs_rule_instr = [
          "GT $lhs_name;\n",
        ]
        lhs_node_count = 1

        lhs_var_track = Dict()

        lhs_pos_track = []
        for cur_type in lhs
            node_type = get_value(cur_type.type.name)
            node_name = get_value(cur_type.name)
            add_node = "$(rule_name)_lhs.addNode({$lhs_node_count, {$type_namespace::$node_type{} }});\n"
            lhs_var_track[lhs_node_count] = (node_name, type_namespace, node_type)
            lhs_node_count += 1
            lhs_rule_instr = [lhs_rule_instr; add_node]
        end

        tmp = join(lhs_rule_instr)

        # LHS Parameters
        lhs_param = rule.lhs_parameter
        rhs_param = rule.rhs_parameter

        lhs_params_var = Dict()
        lhs_param_count = 1
        lhs_node_loc = 1
        # Params count and
        # (param_count, node_loc)
        link_parameter_to_node(lhs_param, lhs_node_loc, lhs_params_var)

        rhs_params_var = Dict()
        rhs_param_count = 1
        rhs_node_loc = 1
        link_parameter_to_node(rhs_param, rhs_node_loc, rhs_params_var)

        # Handle parameter of lhs
        rhs = rule.rhs
        rhs_var_track = Dict()
        rhs_name = "$(rule_name)_rhs"

        rhs_rule_instr = [
          "GT $rhs_name;\n",
        ]
        rhs_node_count = 1
        for cur_type in rhs
            node_type = get_value(cur_type.type.name)
            node_name = get_value(cur_type.name)

            println("node_type: $node_type")
            add_node = "$(rule_name)_rhs.addNode({$rhs_node_count, {$type_namespace::$node_type{} }});\n"
            rhs_var_track[rhs_node_count] = (node_name, type_namespace, node_type )

            rhs_node_count += 1
            rhs_rule_instr = [rhs_rule_instr; add_node]
        end

        # Creating the rule
        with_rule_hdr = [
                "DGGML::WithRule<GT> $rule_name(\"$rule_name\", $lhs_name, $rhs_name,
                                                     [&](auto &lhs, auto &m) {\n"
        ]

        # Handling where clause
        with_clause = rule.modify_clause
        where_clause = with_clause.where_clause

        where_header = [
            "[&](auto &lhs, auto &rhs, auto &m1, auto &m2) {\n"
        ]

        # Creating new variables
        total_assn_ir = []
        seen = Set{String}()

        # Goes through the where clause
        while length(where_clause.clause) > 0
            assn = popfirst!(where_clause.clause)

            println("assn===========--->: ", assn)

            println("assn $(length(where_clause.clause)) ===>: ", assn)

            assn_name = get_value(assn.name)
            assn_type = get_value(assn.type.name)

            assn_type_ir = nothing
            if assn_type == "Integer"
                assn_type_ir = "int"
            elseif assn_type == "Float"
                assn_type_ir = "double"
            end

            values = []

            # Gets assignment values
            get_assn_values(assn, values,
                            lhs_params_var, lhs_var_track,
                            seen, total_assn_ir)

            # exit(0)
            
            # Check if val.position.value is on the right hand size
            # Then assign it
            # (param_count, node_loc)
            rhs_assn_param_loc = rhs_params_var[assn_name][1]
            rhs_assn_node_loc = rhs_params_var[assn_name][2]

            # TODO: If the assignment name is not in rhs_params_var
            # then it is a new variable and we are just keeping it to be used
            # later down the line

            # FIXME: NEED TO HANDLE NESTED
            #        Handle array/lists
            # The parameter location needs to reflect z

            println("rhs_params_var: ", rhs_params_var)
            println("assn_name: ", assn_name)
            println("rhs_var_track: ", rhs_var_track)
            println("rhs_assn_node_loc: ", rhs_assn_node_loc)
            println("rhs_assn_param_loc: ", rhs_assn_param_loc)
            println("symbol_tables: ", symbol_tables)

            check = symbol_tables[rhs_var_track[rhs_assn_node_loc][3]]

            println("check ============> ", check)

            oof = symbol_tables[rhs_var_track[rhs_assn_node_loc][3]][rhs_assn_param_loc][3]
            node_assn_ir = ""
            rhs_param_name = symbol_tables[rhs_var_track[rhs_assn_node_loc][3]][rhs_assn_param_loc][2]

            println("oof: ", oof)

            # Check if the parameter is a list
            if oof[1] == true
                println("in list size")

                list_size = parse(Int, oof[2].token.position.value)-1

                ix = 0
                while list_size >= 0

                    if ix != 0
                        assn = popfirst!(where_clause.clause)
                        println("assn.value: ", assn.value)

                        # values = join([x.position.value for x in assn.value])
                        println("next assn: =======> ", assn)
                        values = []
                        get_assn_values(assn, values, lhs_params_var, lhs_var_track,
                                        seen, total_assn_ir
                                       )
                    end

                    println("assn new =====> : ", assn)

                    assn_name = get_value(assn.name)
                    assn_type = get_value(assn.type.name)
                    assn_type_ir = nothing
                    if assn_type == "Integer"
                        assn_type_ir = "int"
                    elseif assn_type == "Float"
                        assn_type_ir = "double"
                    end

                    rhs_assn_node_loc = rhs_params_var[assn_name][2]
                    rhs_param_name = symbol_tables[rhs_var_track[rhs_assn_node_loc][3]][rhs_assn_param_loc][2]

                    # ??? ass_name is going all over the place

                    # This needs to be fixed. At the moment
                    # it does not handle lists.
                    assn_ir = "\t"*assn_type_ir*" "*assn_name*" = "*join(values)*";\n"

                    println("assn_ir =======>: ", assn_ir)
                    println("values  =======>: ", values)

                    # node type
                    rhs_node_type = rhs_var_track[rhs_assn_node_loc][3]
                    rhs_node_type_namespace = rhs_var_track[rhs_assn_node_loc][2]

                    node_assn_ir = "\tstd::get<$rhs_node_type_namespace::$rhs_node_type>(rhs[m2[$rhs_assn_node_loc]].data).$(rhs_param_name)[$ix] = $assn_name; \n"

                    node_pos_ir = "\trhs[m2[$rhs_assn_node_loc]].position[$ix] = $assn_name; \n"

                    # Right now assume that if it is
                    # a list, then it updates position

                    # How to handle position assignment?
                    # How would the user know this?
                    # Maybe default add position to everything

                    push!(total_assn_ir, assn_ir)
                    push!(total_assn_ir, node_assn_ir)
                    push!(total_assn_ir, node_pos_ir)

                    list_size = list_size - 1
                    ix += 1
                end
                # Start popping and adding the values
            else
                assn_ir = "\t"*assn_type_ir*" "*assn_name*" = "*join(values)*";\n"

                rhs_node_type = rhs_var_track[rhs_assn_node_loc][3]
                rhs_node_type_namespace = rhs_var_track[rhs_assn_node_loc][2]

                node_assn_ir = "\tstd::get<$(rhs_node_type_namespace)::$(rhs_node_type)>(rhs[m2[$rhs_assn_node_loc]].data).$rhs_param_name = $assn_name; \n"

                push!(total_assn_ir, assn_ir)
                push!(total_assn_ir, node_assn_ir)
            end

        end

        # Generating propensity
        # Generating assignment data
        propensity_rule = []
        function_node = with_clause.function_node

        function_name_ir = []
        if get_value(function_node.name) == "Heaviside"
            push!(function_name_ir, "DGGML::heaviside")
        else
            continue
        end

        # FIXME: What to do if counter appears multiple times on the lhs? Which one do you
        # grab? Should return an error.
        function_args_ir = []
        for arg in function_node.args

            if arg.position.value == ","
                push!(function_args_ir, ",")
            elseif arg isa IdentifierToken
                attr_name = arg.position.value
                # exit(0)

                (loc_lhs, param_lhs) = get(lhs_params_var, arg.position.value, nothing)

                class_name = lhs_var_track[loc_lhs][3]
                class_name_space = lhs_var_track[loc_lhs][2]

                node_instr = "auto $(attr_name)_node = lhs.findNode(m[$loc_lhs])->second.getData();\n"
                push!(with_rule_hdr, node_instr)

                lhs_param_name = lhs_var_track[loc_lhs][1]

                println("lhs_var_track: ", lhs_var_track)
                println("lhs_node_loc: ", lhs_node_loc)
                println("lhs_params_var: ", lhs_params_var)
                println("symbol_tables: ", symbol_tables)

                param_name = symbol_tables[class_name][param_lhs][2]

                param_type = symbol_tables[class_name][param_lhs][1]

                instr = "$param_type $attr_name = std::get<$class_name_space::$class_name>($(attr_name)_node.data).$param_name;\n"
                # push!(function_args_ir, lhs_param)
                push!(with_rule_hdr, instr)
                push!(function_args_ir, arg.position.value)
            else
                push!(function_args_ir, arg.position.value)
            end

        end

        # TODO: Reference lhs rule
        function_ir = [
            function_name_ir;"(";function_args_ir;")"
        ]

        propensity_rule = [
           "return ";function_ir;";},"
        ]
        
        
        # For each assignment
        # Find where the variable belongs to for each type.
        # And then update them.
        # Need to match
        oof = "});\n"

        rule_creation_ir = [
            with_rule_hdr; propensity_rule; where_header; total_assn_ir; oof
        ]

        # Assign the new variables to the right location.
        rule_end = ["gamma.addRule($rule_name); }\n"]
        total_rule = join(
                          [
                           rule_header; lhs_rule_instr;
                           rhs_rule_instr; rule_creation_ir ;rule_end
                          ]
                         )

        push!(total_rules, total_rule)

    end

    rule_section_footer = ["};\n", "#endif"]
    rule_section_data = [
                             rule_section_header; 
                             total_rules; rule_section_footer
                        ]

    tmp = join(rule_section_data)
    tmp
end

function ir_models_section(ast)
    """
    Generate Intermediate Models
    for the Section
    """

    model_section_name = "Particles"

    # model_section_name = get_value(ast.name)

    # Create boundary
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

    model_create_initial_type = 
"template<typename GraphType, typename CplexType, typename ParamType, typename GenType>
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
            Particles::graph_type tg;
            Particles::StartType tmp = Particles::StartType{};

            // Initial placement of particle node creator
            tmp.start_location[0] = center_x; // X
            tmp.start_location[1] = center_y; // Y
                                                       //
            node_type oof_node = {gen.get_key(),//i*segments,
                                {tmp, 
                                 center_x, center_y, 0.0}
                        };

            graph.addNode(oof_node);

            };"                               
           

    model_section_header = [
        "#ifndef DGGML_MODELS_HPP\n",
        "#define DGGML_MODELS_HPP\n",
        "#include <fstream>\n",
        "#include \"DGGML.h\"\n",
        "#include \"rules.h\"\n",
        "#include \"simdjson.h\"\n",
        "#include \"ExpandedComplex2D.hpp\"\n",
        "#include \"YAGL_Algorithms.hpp\"\n",
        "#include \"parameters.h\"\n",
        "namespace $model_section_name {\n",
        "using graph_grammar_t = DGGML::Grammar<Particles::graph_type>;\n",
        "class Model : public DGGML::Model<graph_grammar_t> {\n",
        "\tpublic:\n",
        "\tParameters settings;\n",
        ""*model_create_initial_type,
    ""*model_boundary,
        "void initialize() override {\n",
        "geoplex2D.init(
                        settings.CELL_NX,
                        settings.CELL_NY,
                        settings.CELL_DX,
                        settings.CELL_DY,
                        settings.MAXIMAL_REACTION_RADIUS
                );
",
    "\t\tparticle_rules::change_to_particle(gamma, this->system_graph, settings);",
        "\n\t\tthis->add_default_type(this->system_graph,
                    geoplex2D,
                    settings,
                    this->gen);",

                    "\n\t\tthis->add_boundary(
                        geoplex2D.reaction_grid,
                        this->system_graph,
                        geoplex2D,
                        settings,
                        this->gen);",
            "}\n"
    ]

    model_section_body = []
    # Add the file writer
    check_point = 
    "void checkpoint(std::size_t step) override {"*
                "std::string results_dir_name = \"my_results\";"*
                "if(step == 0)
                {
                    // Create the local save directory
                    std::filesystem::remove_all(results_dir_name);
                    std::filesystem::create_directory(results_dir_name);

                        DGGML::GridFileWriter grid_writer;
                    grid_writer.save({geoplex2D.reaction_grid,geoplex2D.dim_label},
                                     results_dir_name+\"/expanded_cell_complex\");

                    std::string title = results_dir_name+\"/simulation_step_\";
                    DGGML::VtkFileWriter<graph_type> vtk_writer;
                    vtk_writer.save(system_graph, title+std::to_string(step));
                    collect(step);
                }
                if( step != 0 & step % 5 == 0)
                {
                    DGGML::GridFileWriter grid_writer;
                    grid_writer.save({geoplex2D.reaction_grid,geoplex2D.dim_label},
                                     results_dir_name+\"/expanded_cell_complex\");

                    std::string title = results_dir_name+\"/simulation_step_\";
                    DGGML::VtkFileWriter<graph_type> vtk_writer;
                    vtk_writer.save(system_graph, title+std::to_string(step));
                    collect(step);
                }
                }"
    model_section_footer = [check_point, "};\n};\n", "#endif"]

    model_section_data = [model_section_header; model_section_footer]

    join(model_section_data)
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
                 "}\n"]

    # check_main = join([ir_header; ir_body; ir_show_gamma; ir_footer])
    check_main = join([ir_header; ir_body; ir_footer])

    check_main
end


# TODO: Need to handle settings.json
function generate_ir()
    """
    Generate IR
    """

    name_space = "Particles"

    base = "tests/generated_tests/generated_2/"
    tokens_params = tokenize_file("tests/particle_sim/params.fflow")
    ast_params = parse_file!(tokens_params)

    generated_params = ir_parameter_section(ast_params[1])
    write_file(
            base*"parameters.h",
            generated_params
        )

    tokens_types = tokenize_file("tests/particle_sim/types.fflow")
    ast_types = parse_file!(tokens_types)

    generated_types = ir_types_section(ast_types[1])
    write_file(base*"types.h",
               generated_types)

    println("doing rule stuff")
    tokens_rules = tokenize_file("tests/particle_sim/rules.fflow")

    ast_rules = parse_file!(tokens_rules)[1]
    # println("ast_rules: ", ast_rules)
    write_file(base*"rules.h", ir_rules_section!(ast_rules))

    main_ir = ir_main(name_space)

    write_file(
            base*"main.cpp",
            main_ir
        )

    ir_models = ir_models_section(nothing)
    write_file(base*"model.h", ir_models)
    # NOTE: Sample assignment ; semi-colon is to denote serialized/sequential oeprations
    #
end

# generate_ir()
