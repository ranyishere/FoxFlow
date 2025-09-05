"""
    IRRuleGeneration:
        Module for generating Intermediate Representation (IR) for rules in a grammar.
"""

# Strucute Function

module IRRuleGeneration

    import ..IRUtils: get_value, convert_type_name
    import ..IRBuildUtils: emit, build, IRBuilder, build_sameline
    import ...AstNodes: RuleNode, IdentifierNode, TypeInstanceNode,
            UndirectedTypeEdgeNode, BinaryOpNode, ExpressionNode, ParameterNode,
            WithClauseNode, SolveClauseNode, UnaryOpNode, FunctionNode, GroupNode,
            CallNode, LiteralNode, DefinitionNode

    export ir_rules_section!

    using OrderedCollections

    const BUILT_IN_FUNC = [
        "heaviside", "sqrt", "normal_distr",
        "uniform_distr", "cos", "sin", "inverse", "pow", "indicator"
    ]

    function ir_indicator_func(arg_str)
        """
        Handles Indicator function
        """

        # TODO: check if values are in propensity table
        ir = "( $arg_str ? 1.0 : 0.0)"
        return ir
    end


    function traverse_group_node(group_node, variables=Set{String}())
        """
        Traverses a group node and gets the expression inside
        as a string. It should also collect variables used in the expression.
        """

        if group_node.expression isa BinaryOpNode
            bin_node = group_node.expression
            lhs = traverse_group_node(GroupNode(bin_node.lhs), variables)
            rhs = traverse_group_node(GroupNode(bin_node.rhs), variables)
            op = bin_node.expression.position.value
            return "($lhs $op $rhs)"
        elseif group_node.expression isa UnaryOpNode
            un_node = group_node.expression
            operand = traverse_group_node(GroupNode(un_node.operand), variables)
            op = un_node.expression.position.value
            return "($op$operand)"
        elseif group_node.expression isa CallNode

            call_node = group_node.expression
            func_name = get_value(call_node.function_node.name)
            args = call_node.args
            arg_str = join(map( (arg) -> begin
                if arg isa LiteralNode
                    return string(get_value(arg))
                elseif arg isa GroupNode
                    return traverse_group_node(arg, variables)
                elseif arg isa BinaryOpNode
                    return traverse_group_node(GroupNode(arg), variables)
                elseif arg isa UnaryOpNode
                    return traverse_group_node(GroupNode(arg), variables)
                else
                    return arg.token
                end
                end, args), ", ")
            return "$func_name($arg_str)"

        elseif group_node.expression isa LiteralNode
            return string(get_value(group_node.expression))

        elseif group_node.expression isa IdentifierNode
            var_name = get_value(group_node.expression)
            push!(variables, var_name)
            return get_value(group_node.expression)
        elseif group_node isa GroupNode
            return traverse_group_node(group_node.expression, variables)
        else
            throw("Unknown group node expression type: $(typeof(group_node.expression))")
        end
    end

    function ir_builtin_func(func_name, args, ir_builder, propensity_table, prop_body_ir)
        """
        Generates the IR for a built-in function.
        This is a placeholder function that should be
        replaced with actual built-in function handling logic.
        """
        
        variables = Set{String}()
        arg_str = join(map( (arg) -> begin
                if arg isa LiteralNode
                    return string(get_value(arg))
                elseif arg isa GroupNode
                    return traverse_group_node(arg, variables)
                elseif arg isa BinaryOpNode
                    return traverse_group_node(GroupNode(arg), variables)
                elseif arg isa UnaryOpNode
                    return traverse_group_node(GroupNode(arg), variables)
                elseif arg isa IdentifierNode
                    var_name = get_value(arg)
                    push!(variables, var_name)
                    return var_name
                elseif arg isa LiteralNode
                    return string(get_value(arg))
                else
                    return arg.token
                end
        end, args), ", ")

        # Emit out variables
        map( (var) -> begin
            if var in collect(keys(propensity_table))
                ir = propensity_table[var]
                emit(prop_body_ir, ir)
            end
            end,
            collect(variables)
        )

        if func_name == "heaviside"
            # arg_str = parse_func_args(args)
            ir = "DGGML::$func_name($arg_str)"

        elseif func_name == "sqrt"
            ir = "sqrt($arg_str)"
        elseif func_name == "normal_distr"
            ir = "DGGML::normal_distr()"
        elseif func_name == "uniform_distr"
            ir = "DGGML::uniform_distr()"
        elseif func_name == "cos"
            ir = "cos($arg_str)"
        elseif func_name == "sin"
            ir = "sin($arg_str)"
        elseif func_name == "inverse"
            ir = "DGGML::inverse()"
        elseif func_name == "pow"
            ir = "pow()"
        elseif func_name == "abs"
            ir = "abs($arg_str)"
        elseif func_name == "indicator"
            ir = ir_indicator_func(arg_str)
        else
            throw("Unknown built-in function: $func_name")
        end

        emit(ir_builder, ir)
    end

    function ir_rules_section!(
            ast, rules_table, symbol_tables,
            propensity_table, type_namespace
        )
        """
        Generate intermediate rules
        for the section.
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
                build_rules_table(rule,
                                  rules_table,
                                  type_namespace)

        # Reset the random_device each time and the local scope
        var_local_table = Dict{String, Any}()

        var_local_table["random_device"] = false
        var_local_table["rule_lhs"] = Dict()
        var_local_table["rule_lhs"]["declared"] = []

        var_local_table["rule_rhs"] = Dict()

        propensity_table["var_local_table"] = var_local_table

                ir_rule!(rules_ir, rule,
                         type_namespace,
                         symbol_tables,
                         propensity_table)

            end, ast.rules_list
           )

        emit(rules_ir, "}\n#endif")
        build(rules_ir)
    end

    function traverse_ode_expr(expression, ir_builder, var_loc_attr)
        """
        Generates ir for an ODE expression
        by traversing the expression tree.
        """

        
        if !(expression isa BinaryOpNode)

            if (expression isa UnaryOpNode)
                operation = expression.expression
                traverse_ode_expr(expression.operand, ir_builder, var_loc_attr)
                emit(ir_builder, " $(get_value(operation)) ")
            else
                emit(ir_builder, " $(get_value(expression)) ")
            end
            return

        end

        if expression.lhs isa BinaryOpNode
            # If the left hand side is a binary operation
            traverse_ode_expr(expression.lhs, ir_builder, var_loc_attr)
        elseif expression.lhs isa IdentifierNode
            # ix_ir = var_loc_ir[get_value(expression.lhs)]
            ix_ir = "ix_"*get_value(expression.lhs)
            ir = "NV_Ith_S(y, varmap.at(&$ix_ir))"
            emit(ir_builder, ir)
        else
            # If it is a single value, just print it
            emit(ir_builder, get_value(expression.lhs))
        end

        op = expression.expression.position.value
        emit(ir_builder, " $op ")

        if expression.rhs isa BinaryOpNode
            # If the left hand side is a binary operation
            traverse_ode_expr(expression.rhs, ir_builder, var_loc_attr)
        elseif expression.rhs isa IdentifierNode
            ix_ir = "ix_"*get_value(expression.rhs)
            ir = "NV_Ith_S(y, varmap.at(&$ix_ir))"
            emit(ir_builder, ir)
        else

            # What if it is a unary op node
            if expression.rhs isa UnaryOpNode
                operation = expression.rhs.expression
                emit(ir_builder, " $(operation.position.value) ")
                traverse_ode_expr(expression.rhs.operand, ir_builder, var_loc_attr)
                # return
            else
                # If it is a single value, just print it
                emit(ir_builder, get_value(expression.rhs))
            end

        end
    end

    function ir_solve_clause!(solve_clause, lhs_assgn_to_node, 
            rhs_assgn_to_node, ir_builder, type_namespace)
        """
        Solve Clause Node
        """

        var_bind_ir = "[](auto &lhs, auto &m1, auto &varset) {"

        emit(ir_builder, var_bind_ir)


        var_attr_loc = Dict()
        # Takes the binding variable and returns the
        # dependency
        bv_to_dep = Dict()

        # Handle variable binding
        # Bind Variable (bv)
        map((bv_node) -> begin

                # TODO: Generalize for multiple var odes (AKA PDE)
                # Only grabs the first variable.
                dep_vars = get_value(bv_node.value[1])
                bv_name = get_value(bv_node.name)

                bv_to_dep[bv_name] = dep_vars
                # bv_pos = rhs_assgn_to_node[bv_name][1]
                bv_pos = lhs_assgn_to_node[dep_vars][1]
                attr_pos = lhs_assgn_to_node[dep_vars][4][3]
                ir = nothing
                ir_ix_attr_loc = nothing
                if attr_pos > 2

                    bv_attr = lhs_assgn_to_node[dep_vars][4][2]
                    bv_type = lhs_assgn_to_node[dep_vars][3]
                    bv_attr_pos = lhs_assgn_to_node[dep_vars][4][3]

                    ref_name = "node_$(bv_pos)_$bv_attr_pos"
                    ref_fetch = "std::get<$type_namespace::$bv_type>(
                        lhs[m1[$bv_pos]].data
                    ).$bv_attr"

                    # Fetching attr
                    ir_fetch = "auto &$ref_name = $ref_fetch;"
                    emit(ir_builder, ir_fetch)

                    # attr_ir = "&lhs[m1[$bv_pos]].$bv_attr"
                    # ir = "varset.insert(&lhs[m1[$bv_pos]].$bv_attr);"
                    ir = "varset.insert(&$ref_name);"

                    ir_ix_attr_loc = "$ref_fetch"
                else
                    ir = "varset.insert(&lhs[m1[$bv_pos]].position[$(attr_pos-1)]);"

                    attr_ir = "lhs[m1[$bv_pos]].position[$(attr_pos-1)]"
                    ir_ix_attr_loc = "$attr_ir"
                end

                var_attr_loc[dep_vars] = ir_ix_attr_loc

                emit(ir_builder, ir)
            end,
            solve_clause.variables
           )

        emit(ir_builder, "},")

        emit(ir_builder, "[&](auto &lhs, auto &m1, auto y, auto ydot, auto &varmap) {")

        # Fetch all the dependent variables first
        map(
            dep_var -> begin
                var_loc_ir = var_attr_loc[dep_var]
                ir = "auto &ix_$dep_var = $var_loc_ir;"
                emit(ir_builder, ir)
            end, collect(keys(var_attr_loc))
       )

        map( assgn_node -> begin
            ode_name = get_value(assgn_node.name)
            ode_value = assgn_node.value

            bv_name = ode_name
            bv_pos = rhs_assgn_to_node[bv_name][1]
            attr_pos = rhs_assgn_to_node[bv_name][4][3]

            # FIXME: It should be grabbing from the data attribute.
            ir_attr = nothing
            ref = nothing
            if attr_pos > 2
                bv_attr = rhs_assgn_to_node[bv_name][4][2]
                ir_attr = bv_attr
                ref = "ix_$(bv_to_dep[bv_name])"
            else
                ir_attr = "position[$(attr_pos-1)]"
                ref = "ix_$(bv_to_dep[bv_name])"
            end

            # ir = "NV_Ith_S(ydot, varmap[&lhs[m1[$bv_pos]].$ir_attr]) += "
            ir = "NV_Ith_S(ydot, varmap[&$ref]) += "

            expr_ir = IRBuilder([])
            traverse_ode_expr(ode_value, expr_ir, var_attr_loc)
            build_expr_ir = join(expr_ir.instructions)
            ir = ir * build_expr_ir * ";"
            emit(ir_builder, ir)

            end, solve_clause.clause
        )

        emit(ir_builder, "}")
        # check = build(ir_builder)
    end

    function ir_rule!(ir_builder, rule_node,
        type_namespace, symbol_tables, propensity_table)
        """
        Generates the IR for a rule, including its header, nodes, and edges.
        """

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

        # Building propensity function
        modify_clause = rule_node.modify_clause
        
        rule_hdr = nothing
        # Generating the rule
        if modify_clause isa WithClauseNode
            rule_hdr = "DGGML::WithRule<GT> $rule_name(\"$rule_name\", $lhs_name, $rhs_name,"
        elseif modify_clause isa SolveClauseNode
            rule_hdr = "DGGML::SolvingRule<GT> $rule_name(\"$rule_name\", $lhs_name, $lhs_name,"
        else
            throw("Unknown modify clause type: $(typeof(modify_clause))")
        end
        emit(ir_builder, rule_hdr)

        lhs_names, lhs_types, lhs_order_of_nodes = build_node_pos_to_type(lhs_node_type)
        lhs_param_to_node, lhs_assgn_to_node = link_param_to_nodes(lhs_param,
                                                lhs_types,
                                                symbol_tables,
                                                lhs_names)

        # rhs_node_ix_to_type, rhs_names, rhs_types = build_node_pos_to_type(rhs_node_type)
        rhs_names, rhs_types, rhs_order_of_nodes = build_node_pos_to_type(rhs_node_type)

        rhs_param_to_node, rhs_assgn_to_node = link_param_to_nodes(rhs_param,
                                                rhs_types,
                                                symbol_tables,
                                                rhs_names)

        if modify_clause isa WithClauseNode
            with_clause = modify_clause
            # Parse with clause
            where_ir_builder = IRBuilder([])

            where_clause = with_clause.where_clause.clause

            ir_where_clause!(where_clause,
                    where_ir_builder,
                    lhs_param_to_node, lhs_assgn_to_node,
                    rhs_param_to_node, rhs_assgn_to_node,
                    symbol_tables, type_namespace, propensity_table)

            prop_ir_builder = IRBuilder([])
            prop_body_ir = IRBuilder([])

            prop_hdr = "[&](auto &lhs, auto &m) {\n"
            emit(prop_body_ir, prop_hdr)
            emit(prop_ir_builder, "return ")

            ir_propensity!(with_clause, prop_ir_builder, propensity_table, prop_body_ir)
            emit(prop_ir_builder, ";},")

            # Build the body
            emit(ir_builder, build(prop_body_ir))
            # Build the return statement
            emit(ir_builder, build(prop_ir_builder))
            emit(ir_builder, build(where_ir_builder))

            # Finish propensity
            # Now you need to check that the parameter number
            # matches with the attribute number in the c++ class.
            # fetch the class name
        elseif modify_clause isa SolveClauseNode
            solve_clause = modify_clause
            num_vars = length(solve_clause.variables)
            emit(ir_builder, "$num_vars,")
            solve_ir_builder = IRBuilder([])
            ir_solve_clause!(solve_clause,
                             lhs_assgn_to_node,
                             rhs_assgn_to_node,
                             solve_ir_builder, type_namespace)
            check = build(solve_ir_builder)
            emit(ir_builder, build(solve_ir_builder))
        end

        emit(ir_builder, ");")
        emit(ir_builder, "gamma.addRule($rule_name);")
        emit(ir_builder, "};")

    end
    
    function emit_rule_header(ir_builder, rule_node, type_namespace)
        """
        Creates the header for a rule function in the IR builder.

        Args:
            ir_builder: The IR builder to which the header will be added.
            rule_node: The rule node containing the rule information.
            type_namespace: The namespace for types used in the rule.
        """

        rule_name = get_value(rule_node.name)
        rule_header = 
           "void $rule_name(DGGML::Grammar<$type_namespace::graph_type> &gamma,
           $type_namespace::graph_type &system_graph,
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

    function emit_add_nodes(ir_builder, graph_name, node_types, type_namespace)
        """
        Generates the IR for adding nodes and edges to the rule's graph.

        Args:
            ir_builder: The IR builder to which the node addition code will be added.
            graph_name: The name of the graph to which nodes are being added.
            node_type: The type of nodes to be added.
            type_namespace: The namespace for types used in the rule.
        """

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

        function extract_edge_vertex_type(node)
            if node isa UndirectedTypeEdgeNode
                left_vertex = node.left_vertex
                right_vertex = node.right_vertex
                left_type_name = get_value(left_vertex.type.name)
                right_type_name = get_value(right_vertex.type.name)
                return (left_type_name,
                        right_type_name)
            else
                return (get_value(node.type.name),)
            end
        end

        function extract_edge_vertex_name(node)
            if node isa UndirectedTypeEdgeNode
                left_vertex = node.left_vertex
                right_vertex = node.right_vertex
                return (get_value(left_vertex.name), get_value(right_vertex.name))
            else
                return (get_value(node.name),)
            end
        end

        function get_node_order(node_names, node_types)
            
            node_type_map = OrderedDict{String, Any}()

            order_of_nodes = []
            count = 1
            for (name_tuple, type_tuple) in zip(node_names, node_types)
                for (name, ntype) in zip(name_tuple, type_tuple)
                    if !haskey(node_type_map, name)
                        push!(order_of_nodes, name)
                        # node_type_map[name] = [ntype]
                        node_type_map[name] = count
                        count += 1
                    else
                        # push!(node_type_map[name], ntype)
                    end
                end
            end
            # Display result
            # for (name, ntype) in node_type_map
                # println("$name => $ntype")
            # end
            # node_type_map
            # println("node_type_map: $node_type_map\n")
            # order_of_nodes
            node_type_map
        end


        node_states = Dict()

        zipped_together = []

        # If the node is an undirected edge and the next node is an undirected edge.
        # Check if the left vertex is the same as the right vertex.

        node_types = []
        node_names = []

        count = 0
        types = map(
                    (node) -> begin

                    extracted_types = extract_edge_vertex_type(node)
                    extracted_names = extract_edge_vertex_name(node)

                    push!(node_types, extracted_types)
                    push!(node_names, extracted_names)
            end,
            type_nodes
        )

        unique_vert = []
        unique_types = []
        for (ix, edge) in enumerate(node_names)
            if length(edge) < 2
                # If it is a single node, we can just push it
                push!(unique_vert, edge[1])
                push!(unique_types, node_types[ix][1])
            else

                if isempty(unique_vert)
                    push!(unique_vert, edge[1])
                    push!(unique_types, node_types[ix][1])
                elseif edge[1] != unique_vert[end]
                    # If the left vertex is not the same as the last unique vertex
                    push!(unique_vert, edge[1])
                    push!(unique_types, node_types[ix][1])
                end

                if edge[2] != edge[1]
                    push!(unique_vert, edge[2])
                    push!(unique_types, node_types[ix][2])
                end

            end
        end

        order_of_nodes = get_node_order(node_names, node_types)

        # verts, node_names, node_types
        unique_vert, unique_types, order_of_nodes
    end

    function build_node_pos_to_type_old(type_nodes)
        """
        From the node position, retrieve
        the type of the node.


        Dict => {
            "node_ix": node_type_name
        }
        """

        # TODO: Count all nodes even duplicates? Right now 

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
                pos_count += 1
                pos_to_type[pos_count] = node_type
                push!(solo_name_registry, node_name)

            end
            
            end,
            type_nodes
           )

        pos_to_type
    end

    function link_param_to_nodes(param, node_types, symbol_tables, node_names)
        """
        Links each parameter to its corresponding node.
        Returns a tuple of two dictionaries:
        - params_var_count: Maps parameter names to their counts in nodes.
        - params_var_loc: Maps parameter names to their node locations.

        Args:
            param: The parameter node containing the parameters.
            nodes_attr: The attributes of the nodes.
            nodes: The nodes in the graph.
            symbol_tables: The symbol tables for the nodes.

        Returns:
            A tuple containing two dictionaries.
        """

        # println("param: $param\n")
        user_param_names = map(
            (x) -> begin get_value(x) end,
            param.token
        )

        println("user_param_names ========================>: ", user_param_names, "\n")

        ix = 0
        node_attrs = []

        for node_type in node_types

            if !(node_type in collect(keys(symbol_tables)))
                println("Node type $node_type not found in symbol tables.")
                continue
            end

            total_node_attr = symbol_tables[node_type]

            # Map this to user_param_names by position
            tmp_attr = []
            for ix_node_attr in 1:length(total_node_attr)
                node_attr = total_node_attr[ix_node_attr]
                attr_type = convert_type_name(node_attr[1])
                attr_name = node_attr[2]

                push!(tmp_attr, (attr_type, attr_name, ix_node_attr))
            end

            push!(node_attrs, tmp_attr)

        end

        flatten_node_attrs = map(
                                 (x) -> begin
                                    map(
                                        (y) -> begin
                                            y
                                        end,
                                        x
                                    )
                                 end, node_attrs
                            )

        flatten_node_attrs = collect(Iterators.flatten(flatten_node_attrs))

        println("flatten_node_attrs ========================>: ", flatten_node_attrs, "\n")

        node_to_params = []
        node_to_types = []
        node_to_loc = []
        unique_count = 0
        seen = Dict()
        # Map parameters to their node locations
        for (ix, node_name) in enumerate(node_names)

            cur_count = nothing
            if !(node_name in collect(keys(seen)))
                unique_count += 1
                seen[node_name] = unique_count
                cur_count = unique_count
            else
                cur_count = seen[node_name]
            end
            cur_attr = node_attrs[ix]
            push!(node_to_params, fill(node_name, length(cur_attr)))
            push!(node_to_types, fill(node_types[ix], length(cur_attr)))
            push!(node_to_loc, fill(cur_count, length(cur_attr)))
        end

        node_to_params = collect(Iterators.flatten(node_to_params))
        node_to_types = collect(Iterators.flatten(node_to_types))
        node_to_loc = collect(Iterators.flatten(node_to_loc))

        zipped = collect(zip(node_to_loc, node_to_params, node_to_types,
                             flatten_node_attrs, user_param_names))

        param_name_to_node = Dict()
        for (node_loc, node_param, node_type, attr, user_param_name) in zipped
            param_name_to_node[user_param_name] = (node_loc, node_param, node_type, attr)
        end

        keys(param_name_to_node) == Set(user_param_names) || throw("Parameter names do not match node parameters.")

        zipped, param_name_to_node
    end

    function ir_where_left_clause!(lhs, type_namespace, symbol_tables,
                                   propensity_table, ir_builder)

        for (ix, param) in enumerate(lhs)

            node_loc = param[1]
            node_type = param[3]
            attr_type = param[4][1]
            attr_name = param[4][2]
            user_param_name = param[end]

            ir = "$attr_type $user_param_name = std::get<$type_namespace::$node_type>(lhs[m1[$node_loc]].data).$attr_name;\n"

            # emit(ir_builder, ir)

            ir_prop = "$attr_type $user_param_name = std::get<$type_namespace::$node_type>(lhs[m[$node_loc]].data).$attr_name;\n"
            propensity_table[user_param_name] = ir_prop
            propensity_table["var_local_table"]["rule_lhs"][user_param_name] = ir
        end

    end

    function ir_where_right_clause!(
                rhs, type_namespace, symbol_tables,
                propensity_table, ir_builder
            )

        rhs_table = Dict()

        for (ix, param) in enumerate(rhs)

            node_loc = param[1]
            node_type = param[3]
            attr_type = param[4][1]
            attr_name = param[4][2]
            user_param_name = param[end]

            node_name = "rhs_node_"*string(node_loc)

            # NOTE: Do I need this??
            if !(node_name in collect(keys(rhs_table)))
                # assign_ir = "$type_namespace::$node_type $node_name = std::get<$type_namespace::$node_type>(rhs[m2[$node_loc]].data);"
                assign_ir = "$type_namespace::$node_type $node_name = std::get<$type_namespace::$node_type>(rhs[m2[$node_loc]].data);"
                rhs_table[node_name] = assign_ir
                # emit(ir_builder, assign_ir)
            end

            ir_prop = "$attr_type $user_param_name = std::get<$type_namespace::$node_type>(lhs[m[$node_loc]].data).$attr_name;\n"

            propensity_table["var_local_table"]["rule_rhs"][user_param_name] = ir_prop
        end

    rhs_table
    end

    function ir_distribution!(node, ir_builder, propensity_table, propensity)
        """
        Generates the IR for a distribution.
        This is a placeholder function that should be
        replaced with actual distribution handling logic.
        """


        # NOTE: Let's just for now assume that we are only dealing with
        # arguments that have identifiers already loaded.
        # Probably should add a scope dictionary...
        distribution = node.function_node.name
        args = node.args

        vars = Set{String}()
        parsed_args = []
        for arg in args
            arg_ir = IRBuilder([])
            ir_prop_expr!(arg, arg_ir, propensity_table, arg_ir, propensity)
            push!(parsed_args, join(arg_ir.instructions))
        end

        args = join(parsed_args, ", ")

        # FIXME: Handle local variable scoping?
        # Check to see if a random device is already been defined
        if propensity_table["var_local_table"]["random_device"] == false
            propensity_table["random_device"] = true
            ir0 = "std::random_device random_device;\n"
            emit(ir_builder, ir0)

            ir1 = "std::mt19937 random_engine(random_device());\n"
            emit(ir_builder, ir1)
            propensity_table["var_local_table"]["random_device"] = true

        end
        
        # TODO: Check dist_type and return the correct one here instead.
        ir2 = "std::uniform_real_distribution<double>("*args *")(random_engine)"
        return ir2
    end

    function ir_where_clause!(where_clause, ir_builder,
        lhs_param_to_node, lhs_assgn_to_node, rhs_param_to_node,
        rhs_assgn_to_node, symbol_tables, type_namespace, propensity_table
    )
        """
        Handles where clause
        """

        where_hdr = "[&](auto &lhs, auto &rhs, auto &m1, auto &m2) {"
        emit(ir_builder, where_hdr)


        # These functions will populate the propensity table
        # with the local variables.
        ir_where_left_clause!(lhs_param_to_node, type_namespace,
                              symbol_tables, propensity_table, ir_builder)

        rhs_table = ir_where_right_clause!(rhs_param_to_node, type_namespace,
                                   symbol_tables, propensity_table, ir_builder)

        map(
            (assign_node) -> begin

                ir_where_assignment!(assign_node, type_namespace, rhs_assgn_to_node,
                                     propensity_table, ir_builder)

                
                    # exit(0)

                    # println(

                    # Need to check for intermediate expressions.

                    # This handles the where clause expressions.
                    # Lets go ahead and handle special cases here.
                    # FIXME: This should be done already at the parser level.
                    # We should not have to do this here. but I am going to because
                    # I am lazy at this moment.
                    # assgn_expr = []
                    # expr_count = 1
                    # while expr_count <= length(node_value)
                        # expr = node_value[expr_count]
                        # ir =  expr.position.value
                        # if ir == "~"
                            # expr_count += 1
                            # dist_type = node_value[expr_count].position.value
                            # expr_count += 1
                            # Remove (
                            # expr_count += 1
                            # if propensity_table["var_local_table"]["random_device"] == false
                                # # We need to define the random device and engine
                                # # only once.
                                # ir0 = "std::random_device random_device;\n"
                                # emit(ir_builder, ir0)
                                # ir1 = "std::mt19937 random_engine(random_device());\n"
                                # emit(ir_builder, ir1)
                                # propensity_table["var_local_table"]["random_device"] = true
                            # end
                            # args = []
                            # while length(node_value) > expr_count
                                # arg = node_value[expr_count].position.value
                                # if arg == ")"
                                    # break
                                # end
                                # push!(args, arg)
                                # expr_count += 1
                            # end
                            # args = join(args)

                            # TODO: Check dist_type and return the correct one here instead.
                            # ir2 = "std::uniform_real_distribution<double>("*args *")(random_engine)"
                            # ir = ir2
                            # push!(assgn_expr, ir)
                            # expr_count += 1
                            # ir = ir_distribution!(node_value, ir_builder)
                        # else
                            # push!(assgn_expr, ir)
                            # expr_count += 1
                        # end
                    # end

                    # joined_expr = join(assgn_expr)

                    # Add to propensity function if it is not used
                    # in node assignment.
                    # assigned_param = rhs_assgn_node[4][3]
                    # assigned_node = rhs_assgn_node[1]
                    # rhs_attr = rhs_assgn_node[4][2]
                    # node_type = rhs_assgn_node[3]

                    # assigned_param_name = assigned_param
                    # ir = "\t\tstd::get<$type_namespace::$node_type>(rhs[m2[$assigned_node]].data).$rhs_attr = $joined_expr;"
                    # if assigned_param == 1 || assigned_param == 2
                        # ir_node_pos = "\t\trhs[m2[$assigned_node]].position[$(assigned_param-1)] = $joined_expr;"
                        # emit(ir_builder, ir_node_pos)
                    # end

                    # emit(ir_builder, ir)
            end, 
            where_clause
        )

        emit(ir_builder, "}")
    end


    function ir_prop_expr!(with_clause, ir_builder, propensity_table, prop_body_ir, propensity)

        function_node = with_clause

        # Check if function node is a function or expression or an identifier
        func_args = nothing
        func_name = nothing

        # TODO: Implement identity function in c++?
        if function_node isa IdentifierNode

            # If it is an identifier, we should
            # fetch the function arguments from the identifier
            arg = function_node.token

            # println("propensity_table: ", propensity_table)

            # FIXME: Check to see if it is in the rule_lhs also..
            if arg.position.value in collect(keys(propensity_table)) && propensity == true
                ir = propensity_table[arg.position.value]
                # Loads it into the body
                emit(prop_body_ir, ir)
                emit(ir_builder, arg.position.value)
            elseif arg.position.value in collect(keys(propensity_table["parameter_table"]))
                # ir = propensity_table[arg.position.value]
                ir = "settings." * arg.position.value
                emit(ir_builder, ir)

            elseif arg.position.value in collect(keys(propensity_table["var_local_table"]["rule_rhs"]))
                ir = propensity_table["var_local_table"]["rule_rhs"][arg.position.value]
                emit(ir_builder, ir)

            elseif arg.position.value in collect(keys(propensity_table["var_local_table"]["rule_lhs"]))

                # Check if it is already emitted
                # If it is already emitted, we should not emit it again.

                if !(arg.position.value in propensity_table["var_local_table"]["rule_lhs"]["declared"])
                    ir = propensity_table["var_local_table"]["rule_lhs"][arg.position.value]
                    emit(prop_body_ir, ir)
                    push!(propensity_table["var_local_table"]["rule_lhs"]["declared"], arg.position.value)
                end

                emit(ir_builder, arg.position.value)

            else

                println("propensity_table keys: ", collect(keys(propensity_table)))
                println("rhs keys: ", collect(keys(propensity_table["var_local_table"]["rule_rhs"])))

                # TODO: Check if it is inside the settings also (global variable)
                # TODO: Check where clause
                throw("Error: Propensity variable $(arg.position.value) not found in propensity table.")
            end

        elseif function_node isa BinaryOpNode
            lhs = function_node.lhs
            ir_prop_expr!(lhs, ir_builder, propensity_table, prop_body_ir, propensity)
            emit(ir_builder, " $(function_node.expression.position.value) ")
            rhs = function_node.rhs
            ir_prop_expr!(rhs, ir_builder, propensity_table, prop_body_ir, propensity)

        elseif function_node isa CallNode
            call_args = function_node.args
            function_node = function_node.function_node
            func_name = get_value(function_node.name)
            func_args = function_node.args
            ir_builtin_func(func_name, func_args,
                            ir_builder, propensity_table, prop_body_ir)

        # TODO: check if it is just a regular digit or something
        elseif function_node isa GroupNode
            # If it is a group node, we need to
            # traverse the expression inside the group
            inner_expr = function_node.expression

            # Should emit parenthesis too but I don't know where to, if it is in the prop body
            # or in the ir builder. 

            emit(ir_builder, "(")
            ir_prop_expr!(inner_expr, ir_builder, propensity_table, prop_body_ir, propensity)
            emit(ir_builder, ")")

        elseif function_node isa LiteralNode
            # If it is a literal, we can just emit it
            emit(ir_builder, get_value(function_node))
        elseif function_node isa UnaryOpNode
            operation = function_node.expression

            # println("funfction_node unar op ===============> ", function_node)
            if operation.position.value == "~"
                ir = ir_distribution!(function_node.operand, prop_body_ir, propensity_table, propensity)
            elseif operation.position.value == "-"
                println("operation: ", operation)
                ir = "-" * get_value(function_node.operand)
            else
                throw("Error: Unary operation $(operation.position.value) not recognized.")
            end
            emit(ir_builder, ir)
        else
            # func_args = function_node.args
            # func_name = get_value(function_node.name)
            throw("Error: $(function_node) not recognized.")
        end

    end

    function ir_propensity!(with_clause, ir_builder, propensity_table, prop_body_ir)
        """
        Handles the propensity function.
        """
        
        # Check to see if the propensity
        # is using a variable. If it is using a propensity.
        function_node = with_clause.function_node

        ir_prop_expr!(function_node, ir_builder, propensity_table, prop_body_ir, true)
    end

    function ir_where_assignment!(assign_node, type_namespace,
                                 rhs_assgn_to_node,
                                 propensity_table, ir_builder)
        """
        Handles the assignment in the where clause.
        """

        if assign_node isa DefinitionNode
            name = get_value(assign_node.name)
            type = convert_type_name(get_value(assign_node.type.name))
            value = assign_node.value

            if !(value isa GroupNode)
                value = GroupNode(value)
            end

            # Generate ir
            ir_value = IRBuilder([])
            ir_prop_expr!(value, ir_value, propensity_table, ir_builder, false)
            ir_value = build_sameline(ir_value)

            # Adds the definition to propensity table local scope.
            println("assigning name: $name with type: $type and value: $ir_value")
            propensity_table["var_local_table"]["rule_rhs"][name] = "$name"

            ir = "$type $name = $ir_value;"
            emit(ir_builder, ir)
        else

            # Check if it is an assign node or just an intermediate expression
            name = get_value(assign_node.name)
            type = convert_type_name(get_value(assign_node.type.name))
            value = assign_node.value

            println("rhs_assgn_to_node: ", collect(keys(rhs_assgn_to_node)))

            rhs_assgn_node = rhs_assgn_to_node[name]
            node_value = assign_node.value
            value = assign_node.value

            local_namespace = type_namespace

            println("assign_node: $assign_node")
            println("Processing assignment for node: $name with value: $node_value")
            if !(value isa GroupNode)
                value = GroupNode(value)
            end

            # Generate ir
            ir_value = IRBuilder([])
            ir_prop_expr!(value, ir_value, propensity_table, ir_builder, false)
            ir_value = build_sameline(ir_value)

            # Adds the definition to propensity table local scope.
            println("assigning name: $name with type: $type and value: $ir_value")
            propensity_table["var_local_table"]["rule_rhs"][name] = "$name"
            new_ir = "$type $name = $ir_value;"

            emit(ir_builder, new_ir)

            assigned_param = rhs_assgn_node[4][3]
            assigned_node = rhs_assgn_node[1]
            rhs_attr = rhs_assgn_node[4][2]
            node_type = rhs_assgn_node[3]

            # assigned_param_name = assigned_param
            ir = "\t\tstd::get<$type_namespace::$node_type>(rhs[m2[$assigned_node]].data).$rhs_attr = $name;"

            if assigned_param == 1 || assigned_param == 2
                ir_node_pos = "\t\trhs[m2[$assigned_node]].position[$(assigned_param-1)] = $name;"
                emit(ir_builder, ir_node_pos)
            end

            emit(ir_builder, ir)
        end
    end

end
