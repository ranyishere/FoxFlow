"""
    IRRuleGeneration:
        Module for generating Intermediate Representation (IR) for rules in a grammar.
"""


module IRRuleGeneration
    import ..IRUtils: get_value, convert_type_name

    import ..IRBuildUtils: emit, build, IRBuilder
    import ...AstNodes: RuleNode, IdentifierNode, TypeInstanceNode, UndirectedTypeEdgeNode, BinaryOpNode, ExpressionNode, ParameterNode
    export ir_rule!

    """
    Generates the IR for a rule, including its header, nodes, and edges.
    """
    function ir_rule!(ir_builder, rule_node,
                type_namespace, symbol_tables, propensity_table)

        rule_name = get_value(rule_node.name)

        emit_rule_header(ir_builder, rule_node, type_namespace)
        
        # Generating the lhs nodes
        lhs_node_type = rule_node.lhs
        lhs_param = rule_node.lhs_parameter
        lhs_name = "$(rule_name)_lhs"
        emit(ir_builder, "GT $lhs_name;")
        emit_add_nodes(ir_builder, lhs_name, lhs_node_type, type_namespace)

        # Generating rhs nodes
        rhs_node_type = rule_node.rhs
        rhs_param = rule_node.rhs_parameter
        rhs_name = "$(rule_name)_rhs"
        emit(ir_builder, "GT $rhs_name;")
        emit_add_nodes(ir_builder, rhs_name, rhs_node_type, type_namespace)
        
        # Generating the rule
        with_rule_hdr = "DGGML::WithRule<GT> $rule_name(\"$rule_name\", $lhs_name, $rhs_name,
                                 [&](auto &lhs, auto &m) {\n"
        emit(ir_builder, with_rule_hdr)

        lhs_node_ix_to_type = build_node_pos_to_type(lhs_node_type)
        lhs_param_to_node = build_param_node_link(lhs_param, 
                              lhs_node_ix_to_type,
                              lhs_node_type, symbol_tables)

        rhs_node_ix_to_type = build_node_pos_to_type(rhs_node_type)
        rhs_param_to_node = build_param_node_link(rhs_param,
                              rhs_node_ix_to_type,
                              rhs_node_type, symbol_tables)

        # Building propensity function
        with_clause = rule_node.modify_clause

        # Parse with clause
        where_ir_builder = IRBuilder([])

        where_clause = with_clause.where_clause.clause

        ir_where_clause!(where_clause,
                where_ir_builder,
                lhs_node_ix_to_type,lhs_param_to_node,
                rhs_node_ix_to_type, rhs_param_to_node,
                symbol_tables, type_namespace, propensity_table)

        prop_ir_builder = IRBuilder([])
        ir_propensity!(with_clause, prop_ir_builder, propensity_table)

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


    """
    Creates the header for a rule function in the IR builder.

    Args:
        ir_builder: The IR builder to which the header will be added.
        rule_node: The rule node containing the rule information.
        type_namespace: The namespace for types used in the rule.
    """
    function emit_rule_header(ir_builder, rule_node, type_namespace)

        rule_name = get_value(rule_node.name)
        rule_header = 
           "void $rule_name(DGGML::Grammar<$type_namespace::graph_type> &gamma,
           Particles::graph_type &system_graph,
           Parameters &settings) {\n"

        emit(ir_builder, rule_header)
    end

    function _emit_single_node(ir_builder, graph_name, each_node_info,
                              rhs_count, rhs_node_registry, rhs_node_key,
                              type_namespace)

        node_type = get_value(each_node_info[2].type.name)
        node_name = get_value(each_node_info[2].name)
        rhs_node_count = rhs_count

        if !(node_name in rhs_node_registry)
            add_node_ir = "$graph_name.addNode({$rhs_node_count, {$type_namespace::$node_type{} }});\n"
            emit(ir_builder, add_node_ir)
            rhs_count += 1
            rhs_node_key[node_name] = rhs_node_count
            push!(rhs_node_registry, node_name)
        else
            rhs_node_count = rhs_node_key[node_name]
        end

        return rhs_count
    end

    function _emit_undirected_edge(ir_builder, graph_name, each_node_info, rhs_count,
                                   rhs_node_registry, rhs_node_key, type_namespace)

        vert_0 = each_node_info[2].left_vertex
        node_0_name = get_value(vert_0.name)
        rhs_node_count_0 = rhs_count

        if !(node_0_name in rhs_node_registry)
            node_type_0 = get_value(vert_0.type.name)
            add_node_ir_0 = "$graph_name.addNode({$rhs_node_count_0, {$type_namespace::$node_type_0{} }});\n"
            rhs_count += 1
            emit(ir_builder, add_node_ir_0)
            rhs_node_key[node_0_name] = rhs_node_count_0
            push!(rhs_node_registry, node_0_name)
        else
            rhs_node_count_0 = rhs_node_key[node_0_name]
        end

        rhs_node_count_1 = rhs_count
        vert_1 = each_node_info[2].right_vertex
        node_type_1 = get_value(vert_1.type.name)
        node_1_name = get_value(vert_1.name)

        if !(node_1_name in rhs_node_registry)
            node_type_1 = get_value(vert_1.type.name)
            rhs_node_key[node_1_name] = rhs_node_count_1
            add_node_ir_1 = "$graph_name.addNode({$rhs_node_count_1, {$type_namespace::$node_type_1{} }});\n"
            rhs_count += 1
            emit(ir_builder, add_node_ir_1)
            push!(rhs_node_registry, node_1_name)
        else
            rhs_node_count_1 = rhs_node_key[node_1_name]
        end

        # Create an edge
        edge_ir = "$graph_name.addEdge($rhs_node_count_0, $rhs_node_count_1);\n"
        emit(ir_builder, edge_ir)
        return rhs_count
    end


    """
    Generates the IR for adding nodes and edges to the rule's graph.

    Args:
        ir_builder: The IR builder to which the node addition code will be added.
        graph_name: The name of the graph to which nodes are being added.
        node_type: The type of nodes to be added.
        type_namespace: The namespace for types used in the rule.
    """
    function emit_add_nodes(ir_builder, graph_name, node_types, type_namespace)

        rhs_count = 1
        rhs_node_registry = []
        rhs_node_key = Dict()
        map( (each_node_info) -> begin

            if each_node_info[2] isa UndirectedTypeEdgeNode
                # Undirected Edge
                rhs_count = _emit_undirected_edge(ir_builder, graph_name, each_node_info,
                                      rhs_count, rhs_node_registry, rhs_node_key,
                                      type_namespace)
            else
                # @assert each_node_info[2] isa TypeInstanceNode "Expected TypeInstanceNode, got $(typeof(each_node_info[2]))"

                println("$(typeof(each_node_info[2]))")
                println(typeof(each_node_info[2]) isa Main.UndirectedTypeEdgeNode)

                # Single Node
                rhs_node_count = rhs_count
                node_type = get_value(each_node_info[2].type.name)
                node_name = get_value(each_node_info[2].name)
                add_node_ir = "$graph_name.addNode({$rhs_node_count, {$type_namespace::$node_type{} }});\n"
                emit(ir_builder, add_node_ir)
                rhs_count += 1
            end

        end,
        enumerate(node_types)
        )

    end

    function build_rules_table(rule_node, rules_table, type_namespace)

        rule_name = get_value(rule_node.name)
        rules_table[rule_name] = Dict()
        rules_table[rule_name]["lhs"] = rule_node.lhs
        rules_table[rule_name]["rhs"] = rule_node.rhs
        rules_table[rule_name]["lhs_parameter"] = rule_node.lhs_parameter
        rules_table[rule_name]["rhs_parameter"] = rule_node.rhs_parameter
        rules_table[rule_name]["modify_clause"] = rule_node.modify_clause
        rules_table[rule_name]["activated"] = true
        rules_table[rule_name]["namespace"] = rule_namespace

        return Dict{String, Dict}()
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

            if node isa UndirectedTypeEdgeNode
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

    function build_param_node_link(param, nodes_attr, nodes, symbol_tables)

            params_var = Dict()
            param_count = 1

            # (param_count, node_loc)
            # node_param_count = Dict()
            # param_stack = copy(param)
            node_count = 1
            global_param_count = 0

            total_param_names = map( (x) -> begin get_value(x) end, param.token)
            total_node_pos = collect(keys(nodes_attr))

            println("total_node_pos --->: $total_node_pos")
            println("nodes_attr ---->: $nodes_attr")
            # Note there are actually only 5 nodes here that are unique.
            # but there is a reference to the same node multiple times

            # Map node position
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

    # TODO: Handle expressions.
    """
    Build the lhs and rhs nodes
    Key is the parameter number
    Value is the parameter name in the
    class that is generated.
    If the parameter belongs to a list
    then the value references the attribute
    of the list

    For each lhs and rhs grab the parameters
    and update their value based on the
    where clause

    Do a get on everything on the lhs parameter
    assign everything on the rhs based on
    what is inside the where clause or
    in the rhs parameter definition.

    Associate rhs parameter location with the
    appropriate c++ class attribute. For
    the left hand side we are doing a
    bunch of gets.
    """
    function ir_where_clause!(where_clause, ir_builder,
                            lhs_node_ix_to_type, lhs_param_to_node,
                            rhs_node_ix_to_type, rhs_param_to_node,
                            symbol_tables, type_namespace, propensity_table)

        
        where_hdr = "[&](auto &lhs, auto &rhs, auto &m1, auto &m2) {"

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

function ir_propensity!(with_clause, ir_builder, propensity_table)

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





end
