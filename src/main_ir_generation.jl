include("./main_parser_op.jl")

# FIXME: Spaces should not be relevant.
#
# FIXME: need to fix that the rules are updating the wrong nodes.
#    seems like node counts are wrong.

using UUIDs

symbol_tables = Dict()
rules_table = Dict()
propensity_table = Dict()

type_namespace = ""
rule_namespace = ""

# Maps ix to the attribute
parameter_table = Dict()

mutable struct IRBuilder
    instructions::Array{Any}
end

function emit(ir_builder::IRBuilder, instruction::String)
    push!(ir_builder.instructions, instruction)
end

function build(ir_builder::IRBuilder)
    return join(ir_builder.instructions, "\n")
end

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
        param_count = 1
        for (ix, name) in enumerate(param_names)

            if is_list[ix][1] == true

                # It's a list
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

    global type_namespace = section_name

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

    symbol_tables["StartType"][1] = ("float", "start_location[0]", start_type_info)
    symbol_tables["StartType"][2] = ("float", "start_location[1]", start_type_info)

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

function build_param_node_link(param, nodes_attr, nodes)

        params_var = Dict()
        param_count = 1
        # node_loc = 1

        # tells us which attribute maps to what parameter number for
        # a single node type
        
        # (param_count, node_loc)
        # node_param_count = Dict()
        # param_stack = copy(param)
        node_count = 1
        global_param_count = 0

        total_param_names = map( (x) -> begin get_value(x) end, param.token)
        total_node_pos = collect(keys(nodes_attr))

        node_pos = 1
        global_param_count_check = 1
        while node_pos <= length(total_node_pos)
            param_count = 1

            node_type = nodes_attr[node_pos]
            node_params = symbol_tables[node_type]
            while param_count <= length(node_params)
                user_param = total_param_names[global_param_count_check]
                params_var[user_param] = (param_count, node_pos)
                param_count += 1
                global_param_count_check += 1
            end

            node_pos += 1
        end

        # Params count and
        # (param_count, node_loc)
    params_var
end

function build_node_pos_to_type(type_nodes)
    """
    From the node position, retrieve
    the type of the node.


    Dict => {
        "node_ix": node_type_name
    }
    """

    # FIXME: For some reason a triple edge node is returning not the correct
    # edge type
    left_name_registry = []
    right_name_registry = []
    solo_name_registry = []

    pos_count = 1

    pos_to_type = Dict()

    # What is the purpose of this function?
    map(
        (info) -> begin
        # node_ix = pos_count
        node = info

        if node isa EdgeNode
            # If it is an edge node, we need to
            # grab the left and right vertex
            # and get their types.
            # node_type = "Edge"

            # This is creating a duplicate vertex if the left
            # vertex is the same.

            # Can we check if the name is the same?

            # Left vertex
            left_vertex = node.left_vertex
            left_name = get_value(left_vertex.name)

            #check if it is in registry
            if !(left_name in left_name_registry) && !(left_name in right_name_registry)
                left_node_type = get_value(left_vertex.type.name)
                pos_to_type[pos_count] = left_node_type
                pos_count += 1
                push!(left_name_registry, left_name)
            else
                # If the name is already registered
                # we should not add it again.
                println("found in registry")
                # exit(0)
                # return
            end

            # Right vertex
            right_vertex = node.right_vertex
            right_name = get_value(right_vertex.name)
            if left_name == right_name
                # If the left and right vertex are the same
                # we should not add it again.

                println("left and right vertex are the same, exiting")
                exit(0)
                return
            end

            # check if it is in registry
            if !(right_name in right_name_registry)
                right_node_type = get_value(right_vertex.type.name)
                pos_to_type[pos_count] = right_node_type
                pos_count += 1
                push!(right_name_registry, right_name)
            else

                println("the node name is already registered: $right_name")
                exit(0)
                # If the name is already registered
                # we should not add it again.
                return
            end
        else

            # It's a solo node
            node_type = get_value(node.type.name)
            node_name = get_value(node.name)

            if !(node_name in solo_name_registry)
                pos_to_type[pos_count] = node_type
                pos_count += 1
                push!(solo_name_registry, node_name)
            else
                # If the name is already registered
                # we should not add it again.
                println("the node name is already registered: $node_name")
                exit(0)
                return
            end

        end
        
        end,
        type_nodes
       )

    pos_to_type
end

function fetch_lhs_params(ir_builder, lhs_param_to_node, lhs_node_ix_to_type)
    map(
        (each_param) -> begin

            param_ix = each_param[1]
            param_val = each_param[2]

            user_param_name = collect(keys(lhs_param_to_node))[param_ix]

            param_count = param_val[1]
            node_count = param_val[2]
            # Now that we have the node_count
            # we should be able to fetch

            # Grabbing node name
            node_type = lhs_node_ix_to_type[node_count]

            node_params = symbol_tables[node_type]

            attr_name = node_params[param_count][2]

            attr_type = convert_type_name(node_params[param_count][1])

            ir = "\t\t $attr_type $user_param_name = std::get<$type_namespace::$node_type>(lhs[m1[$(node_count)]].data).$attr_name;\n"
            emit(ir_builder, ir)

            end,
            enumerate(values(lhs_param_to_node))
    )

end


function ir_propensity!(with_clause, ir_builder)

    # Check to see if the propensity
    # is using a variable. If it is using a propensity.
    
    function_node = with_clause.function_node

    # Check if function node is a function or expression or an identifier

    func_args = nothing
    func_name = nothing

    # TODO: Implement identity function in c++?
    if function_node isa IdentifierNode
        # If it is an identifier, we should
        # fetch the function arguments from the identifier
        func_args = get_value(function_node)
    elseif function_node isa BinaryOpNode
        # TODO finish if it is a binary expression
    else
        func_args = function_node.args
        func_name = get_value(function_node.name)
    end

    # TODO: Check if function arg is from somewhere
    # TODO: Finish creating function

    fetch_ir = []
    arg_str = join(
                   map(

                (arg) -> begin

                if arg.position.value in collect(keys(propensity_table))
                    ir = propensity_table[arg.position.value]
                    emit(ir_builder, ir)
                end

                # push!(fetch_ir, found)

                # TODO: If it is not in the propensity table, it's probably
                # from the where clause. Need to update and handle that.

                return arg.position.value
            end,
        func_args
        )
    )


    # Check the function_node type

    return_ir = nothing
    if function_node isa IdentifierNode
        # If it is an identifier, we should
        # fetch the function arguments from the identifier
        # func_name = get_value(function_node)
        # return_ir = 
        exit(0)
    else
        return_ir = "return DGGML::$func_name($arg_str);},"
    end

    emit(ir_builder, return_ir)
end

# TODO: Handle expressions.
function ir_where_clause!(where_clause, ir_builder,
                        lhs_node_ix_to_type, lhs_param_to_node,
                        rhs_node_ix_to_type, rhs_param_to_node)

    # Build the lhs and rhs nodes
    # Key is the parameter number
    # Value is the parameter name in the
    # class that is generated.
    # If the parameter belongs to a list
    # then the value references the attribute
    # of the list

    # For each lhs and rhs grab the parameters
    # and update their value based on the
    # where clause

    # Do a get on everything on the lhs parameter
    # assign everything on the rhs based on
    # what is inside the where clause or
    # in the rhs parameter definition.

    # Associate rhs parameter location with the
    # appropriate c++ class attribute. For
    # the left hand side we are doing a
    # bunch of gets.

    where_hdr = "[&](auto &lhs, auto &rhs, auto &m1, auto &m2) {"


    println("lhs_parma_to_node: $lhs_param_to_node")

    emit(ir_builder, where_hdr)
    map(
        (each_param) -> begin

            param_ix = each_param[1]
            param_val = each_param[2]

            user_param_name = collect(keys(lhs_param_to_node))[param_ix]

            param_count = param_val[1]
            node_count = param_val[2]
            # Now that we have the node_count
            # we should be able to fetch
            #

            # Grabbing node name
            node_type = lhs_node_ix_to_type[node_count]
            node_params = symbol_tables[node_type]
            attr_name = node_params[param_count][2]

            attr_type = convert_type_name(node_params[param_count][1])

            ir = "$attr_type $user_param_name = std::get<$type_namespace::$node_type>(lhs[m1[$(node_count)]].data).$attr_name;\n"

            ir_prop = "$attr_type $user_param_name = std::get<$type_namespace::$node_type>(lhs[m[$(node_count)]].data).$attr_name;\n"

            propensity_table[user_param_name] = ir_prop

            emit(ir_builder, ir)

            end,
            enumerate(values(lhs_param_to_node))
    )

    # Mapping to assign 
    # Note it can also be mapping to the
    # lhs also...
    node_loc_to_name = []
    map( (each_node) -> begin

            node_name = get_value(each_node.name)

            # NOTE:
            # Check if it is in rhs_param_to_node if it isn't
            # check if it is in the propensity function
            # If it is not in the propensity function, maybe it's just a helper?

            if !(node_name in collect(keys(rhs_param_to_node)))
                if (node_name in collect(keys(propensity_table)))
                    propensity_table[node_name] = each_node
                    # return
                end
                # TODO: Probably just a helper variable, we need to
                # That we have to put in memory.
                return
            end

            assigned_data = rhs_param_to_node[node_name]

            assigned_param = assigned_data[1]
            assigned_node = assigned_data[2]

            println("assigned_node: $assigned_node")


            rhs_node_type = rhs_node_ix_to_type[assigned_node]
            rhs_node_params = symbol_tables[rhs_node_type]

            # Check if the assigned_param is the first two positions
            # If it is, it represents the position of the node
            # on the graph


            rhs_attr = rhs_node_params[assigned_param][2]

            node_name = "rhs_node_"*string(assigned_node)
            if !(node_name in node_loc_to_name)
                push!(node_loc_to_name, node_name)
                assign_ir = "$type_namespace::$rhs_node_type $node_name = std::get<$type_namespace::$rhs_node_type>(rhs[m2[$assigned_node]].data);"
                emit(ir_builder, assign_ir)
            else
            end

            node_type = each_node.type.name
            node_value = each_node.value

            assn_type = convert_type_name(get_value(each_node.type.name))

            local_namespace = type_namespace
            assgn_expr = []
            map( (expr) -> begin

                    ir =  expr.position.value
                    push!(assgn_expr, ir)

                end, node_value)


            joined_expr = join(assgn_expr)
            
            assigned_param_name = assigned_param

            # Update the rhs node data
            ir = "\t\tstd::get<$type_namespace::$rhs_node_type>(rhs[m2[$assigned_node]].data).$rhs_attr = $joined_expr;"

            # Update node position if it is in the first two positions of the
            # new node
            #
        if assigned_param == 1 || assigned_param == 2
            ir_node_pos = "\t\trhs[m2[$assigned_node]].position[$(assigned_param-1)] = $joined_expr;"
            emit(ir_builder, ir_node_pos)
        end

            emit(ir_builder, ir)
        end,
        where_clause
    )

    emit(ir_builder, "}")
end

function ir_rule!(ir_builder, rule_node)

    global propensity_table = Dict()
    rule_name = get_value(rule_node.name)

    # Note we will turn this into a graph eventually
    rules_table[rule_name] = Dict()
    rules_table[rule_name]["lhs"] = rule_node.lhs
    rules_table[rule_name]["rhs"] = rule_node.rhs
    rules_table[rule_name]["lhs_parameter"] = rule_node.lhs_parameter
    rules_table[rule_name]["rhs_parameter"] = rule_node.rhs_parameter
    rules_table[rule_name]["modify_clause"] = rule_node.modify_clause
    rules_table[rule_name]["activated"] = true
    rules_table[rule_name]["namespace"] = rule_namespace

    rule_header = 
       "void $rule_name(DGGML::Grammar<$type_namespace::graph_type> &gamma,
       Particles::graph_type &system_graph,
       Parameters &settings) {\n"

    emit(ir_builder, rule_header)

    # Generating the lhs nodes
    lhs_node_type = rule_node.lhs

    lhs_param = rule_node.lhs_parameter
    lhs_name = "$(rule_name)_lhs"
    emit(ir_builder, "GT $lhs_name;")

    map( (each_node_info) -> begin

        lhs_node_count = each_node_info[1]
        node_type = get_value(each_node_info[2].type.name)

        add_node_ir = "$(rule_name)_lhs.addNode({$lhs_node_count, {$type_namespace::$node_type{} }});\n"

        emit(ir_builder, add_node_ir)
        end,
        enumerate(lhs_node_type)
    )

    # Generating rhs nodes
    rhs_node_type = rule_node.rhs
    rhs_param = rule_node.rhs_parameter

    rhs_name = "$(rule_name)_rhs"
    emit(ir_builder, "GT $rhs_name;")

    rhs_count = 1
    rhs_node_registry = []
    println("rhs_node_type: $rhs_node_type")

    map( (each_node_info) -> begin

        # println("each_node_info: ", each_node_info)
        # Handle if it is an edge node
        # TODO: Handle things like Node 1 --- Node 2 -- Node 3
        # Need to handle node counts with edges

        println("each_node_info: ", each_node_info)
        if each_node_info[2] isa EdgeNode

            vert_0 = each_node_info[2].left_vertex

            node_0_name = get_value(vert_0.name)
            rhs_node_count_0 = rhs_count

            if !(node_0_name in rhs_node_registry)

                # rhs_node_count_0 = each_node_info[1]
                node_type_0 = get_value(vert_0.type.name)

                add_node_ir_0 = "$(rule_name)_rhs.addNode({$rhs_node_count_0, {$type_namespace::$node_type_0{} }});\n"
                rhs_count += 1
                emit(ir_builder, add_node_ir_0)

                push!(rhs_node_registry, node_0_name)
            else
                rhs_node_count_0 = rhs_count - 1
            end

            # rhs_node_count_1 = each_node_info[1]
            vert_1 = each_node_info[2].right_vertex
            rhs_node_count_1 = rhs_count

            node_type_1 = get_value(vert_1.type.name)
            add_node_ir_1 = "$(rule_name)_rhs.addNode({$rhs_node_count_1, {$type_namespace::$node_type_1{} }});\n"
            # rhs_count += 1

            emit(ir_builder, add_node_ir_1)

            # Create an edge
            edge_ir = "$(rule_name)_rhs.addEdge($rhs_node_count_0, $rhs_node_count_1);\n"

            emit(ir_builder, edge_ir)
        else

            # Single Node
            rhs_node_count = rhs_count
            node_type = get_value(each_node_info[2].type.name)

            add_node_ir = "$(rule_name)_rhs.addNode({$rhs_node_count, {$type_namespace::$node_type{} }});\n"
            emit(ir_builder, add_node_ir)

            rhs_count += 1
        end

        end,
        enumerate(rhs_node_type)
    )

    # Generating the rule
    with_rule_hdr = "DGGML::WithRule<GT> $rule_name(\"$rule_name\", $lhs_name, $rhs_name,
                                                     [&](auto &lhs, auto &m) {\n"
    emit(ir_builder, with_rule_hdr)

    lhs_node_ix_to_type = build_node_pos_to_type(lhs_node_type)
    lhs_param_to_node = build_param_node_link(lhs_param, 
                                              lhs_node_ix_to_type,
                                              lhs_node_type)

    rhs_node_ix_to_type = build_node_pos_to_type(rhs_node_type)

    rhs_param_to_node = build_param_node_link(rhs_param,
                                              rhs_node_ix_to_type,
                                              rhs_node_type)

    # Building propensity function
    with_clause = rule_node.modify_clause

    # Parse with clause
    where_ir_builder = IRBuilder([])

    where_clause = with_clause.where_clause.clause

    ir_where_clause!(where_clause, where_ir_builder,
                    lhs_node_ix_to_type,
                    lhs_param_to_node,
                    rhs_node_ix_to_type,
                    rhs_param_to_node)

    prop_ir_builder = IRBuilder([])
    ir_propensity!(with_clause, prop_ir_builder)

    emit(ir_builder, build(prop_ir_builder))
    emit(ir_builder, build(where_ir_builder))

    # Finish propensity
    # Now you need to check that the parameter number
    # matches with the attribute number in the c++ class.
    # fetch the class name
    
    emit(ir_builder, ");")

    emit(ir_builder, "gamma.addRule($rule_name);")
    emit(ir_builder, "};")
end

function ir_rules_section!(ast)
    """
    Generate Intermediate Rules
    for the Section

    # TODO: Handle edges in the graph
    """

    rule_section_name = get_value(ast.name)
    global rule_namespace = rule_section_name

    rule_section_header = [
        "#ifndef DGGML_RULES_HPP",
        "#define DGGML_RULES_HPP",
        "#include \"types.h\"",
        "#include \"parameters.h\"",
        "namespace $rule_section_name {",
        "using GT = $type_namespace::graph_type;"
    ]

    rules_ir = IRBuilder([])

    map(
        (hdr) -> begin
            emit(rules_ir, hdr)
        end,
        rule_section_header
    )


    map(
        (rule) -> begin
            ir_rule!(rules_ir, rule)
        end, ast.rules_list
       )

    emit(rules_ir, "}\n#endif")
    build(rules_ir)
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

function ir_models_section(ast)
    """
    Generate Intermediate Models
    for the Section
    """

    models_ir = IRBuilder([])

    model_boundary = ir_boundary!(models_ir)
    model_section_name = type_namespace

    # model_section_name = get_value(ast.name)
    # Create boundary
    
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

    call_rule_ir = map(
        (rule_name) -> begin
            if rules_table[rule_name]["activated"] == true
                tmp_ir = "\t$(rules_table[rule_name]["namespace"])::$rule_name(gamma, this->system_graph, settings);\n"
                return tmp_ir
            end
        end,
        collect(keys(rules_table))
       )

    initial_rules_ir = join(call_rule_ir)

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
"*initial_rules_ir,
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

    # emit(models_ir, model_section_data)
    # build(models_ir)
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
                 "}\n"
     ]

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

    println("Generating Params")
    generated_params = ir_parameter_section(ast_params[1])

    write_file(base*"parameters.h",
            generated_params)

    tokens_types = tokenize_file("tests/particle_sim/types.fflow")
    ast_types = parse_file!(tokens_types)
    println("Generating Types")
    generated_types = ir_types_section(ast_types[1])
    write_file(base*"types.h",
               generated_types)

    tokens_rules = tokenize_file("tests/particle_sim/rules.fflow")
    ast_rules = parse_file!(tokens_rules)[1]
    println("Generating Rules")
    write_file(base*"rules.h", ir_rules_section!(ast_rules))

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
