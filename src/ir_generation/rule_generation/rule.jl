# Top-level per-rule IR generation (ir_rule!) and graph construction.
# Included into module IRRuleGeneration (see ../ir_rule_generation.jl).

    """
    Generates the IR for a rule, including its header, nodes, and edges.
    """
    function ir_rule!(ir_builder, rule_node,
        type_namespace, symbol_tables, propensity_table)

        """
        Visits each node in the rule and builds the graph
        for the lhs and rhs of the rule
        """
        function visit_rule_node!(each_node, graph_table, ir_builder,
            graph_name, type_namespace=type_namespace)
            
            # NOTE: Shouldn't this include adding
            # the symbol table?

            if each_node isa UndirectedTypeEdgeNode
                left_node = each_node.left_vertex
                right_node = each_node.right_vertex

                visit_rule_node!(left_node, graph_table,
                                ir_builder, graph_name, type_namespace)

                visit_rule_node!(right_node, graph_table,
                                ir_builder, graph_name, type_namespace)

                rhs_node_count_0 = graph_table[(get_value(left_node.name), get_value(left_node.type.name))]
                rhs_node_count_1 = graph_table[(get_value(right_node.name), get_value(right_node.type.name))]

                edge_ir = "$graph_name.addEdge($rhs_node_count_0, $rhs_node_count_1);\n"
                # edge_ir = "$graph_name.addEdge($lhs_ix, $ix);\n"

                emit(ir_builder, edge_ir)

                # emit edge

            elseif each_node isa TypeInstanceNode
                node_type = get_value(each_node.type.name)
                node_name = get_value(each_node.name)

                # if haskey(graph_table, (node_name, node_type))
                    # throw("Duplicate node name found in rule $graph_name: $node_name, type: $node_type")
                # end
                # Just get the max value in graph_table
                vals  = values(graph_table)

                res = nothing
                if length(vals) == 0
                    res = 0
                else
                    res = maximum(collect(vals))
                end

                if !((node_name, node_type) in keys(graph_table))

                    node_names = [x[1] for x in keys(graph_table)]

                    if node_name in node_names
                        throw("Duplicate node name found in rule $graph_name: $node_name, type: $node_type. Are you sure you want to do this?")
                    end

                    res = res + 1
                    graph_table[(node_name, node_type)] = res
                    # Grab from symbol table
                    # cur_count = graph_table[(node_name, node_type)]
                    # add_node_ir = "$graph_name.addNode({$rhs_node_count, {$type_namespace::$node_type{} }});\n"
                    # Emit it since it hasn't been seen yet.
                    add_node_ir = "$graph_name.addNode({$res, {$type_namespace::$node_type{} }});\n"
                    emit(ir_builder, add_node_ir)

                else
                    res = graph_table[(node_name, node_type)]
                end

                else
                    println("Visiting other node type: ", typeof(each_node))
            end
        end

        rule_name = get_value(rule_node.name)

        println("Rule name ====================> ", rule_name)

        emit_rule_header(ir_builder, rule_node, type_namespace)
        
        # Generating the lhs nodes
        lhs_node_type = rule_node.lhs

        lhs_param = rule_node.lhs_parameter
        lhs_name = "$(rule_name)_lhs"
        emit(ir_builder, "GT $lhs_name;")

        # Generating lhs graph
        lhs_symbol = Dict()
        # lhs_count = 1
        for each_node in lhs_node_type
            visit_rule_node!(each_node, lhs_symbol, ir_builder, lhs_name, type_namespace)
        end

        # emit_add_nodes(ir_builder, lhs_name, lhs_node_type, type_namespace, lhs_symbol)

        # Generating rhs graph
        rhs_node_type = rule_node.rhs
        rhs_param = rule_node.rhs_parameter

        rhs_node_type = rule_node.rhs
        rhs_symbol = Dict()
        # rhs_count = 1
        rhs_name = "$(rule_name)_rhs"
        emit(ir_builder, "GT $rhs_name;")
        for each_node in rhs_node_type
            visit_rule_node!(each_node, rhs_symbol,
                            ir_builder, rhs_name, type_namespace)
        end
        # rhs_symbol = visit_rule_node!(rhs_node_type, rhs_symbol)
        # emit_add_nodes(ir_builder, rhs_name, rhs_node_type, type_namespace, rhs_symbol)

        # Building propensity function
        modify_clause = rule_node.modify_clause
        
        rule_hdr = nothing
        # Generating the rule
        if modify_clause isa WithClauseNode
            rule_hdr = "DGGML::WithRule<GT> $rule_name(\"$rule_name\", $lhs_name, $rhs_name,"
        elseif modify_clause isa EveryClauseNode
            validate_every_rule(rule_name, rule_node, lhs_symbol, rhs_symbol)
            rule_hdr = "DGGML::EveryRule<GT> $rule_name(\"$rule_name\", $lhs_name, $rhs_name,"
        elseif modify_clause isa SolveClauseNode
            rule_hdr = "DGGML::SolvingRule<GT> $rule_name(\"$rule_name\", $lhs_name, $lhs_name,"
        else
            throw("Unknown modify clause type: $(typeof(modify_clause))")
        end

        emit(ir_builder, rule_hdr)

        # println("lhs_types: ", lhs_types)

        lhs_names, lhs_types, lhs_order_of_nodes = build_node_pos_to_type(lhs_node_type)


        println("rule_name: ", rule_name)
        lhs_param_to_node, lhs_assgn_to_node = link_param_to_nodes(lhs_param,
                                                    lhs_types,
                                                    symbol_tables,
                                                    lhs_names)
        # exit(0)

        # rhs_node_ix_to_type, rhs_names, rhs_types = build_node_pos_to_type(rhs_node_type)
        rhs_names, rhs_types, rhs_order_of_nodes = build_node_pos_to_type(rhs_node_type)

        # println("rhs_name: ", rhs_name)
        rhs_param_to_node, rhs_assgn_to_node = link_param_to_nodes(rhs_param,
                                                rhs_types,
                                                symbol_tables,
                                                rhs_names)

        # Map every RHS pattern node name to its (m2 loc, type) using the SAME
        # first-appearance ordering link_param_to_nodes uses, so implicit carry-over
        # can resolve a persistent node even when ALL of its RHS attributes are
        # omitted (and therefore absent from rhs_param_to_node).
        rhs_node_locs = Dict()
        let c = 0
            for (ix, nm) in enumerate(rhs_names)
                if !haskey(rhs_node_locs, nm)
                    c += 1
                    rhs_node_locs[nm] = (c, rhs_types[ix])
                end
            end
        end

        if modify_clause isa WithClauseNode
            with_clause = modify_clause
            # Parse with clause
            where_ir_builder = IRBuilder([])

            where_clause = with_clause.where_clause.clause

            ir_where_clause!(where_clause,
                    where_ir_builder,
                    lhs_param_to_node, lhs_assgn_to_node,
                    rhs_param_to_node, rhs_assgn_to_node,
                    symbol_tables, type_namespace, propensity_table,
                    rhs_node_locs)

            prop_ir_builder = IRBuilder([])
            prop_body_ir = IRBuilder([])

            prop_hdr = "[&](auto &lhs, auto &m1) {\n"
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
        elseif modify_clause isa EveryClauseNode
            every_clause = modify_clause
            where_ir_builder = IRBuilder([])

            ir_where_clause!(every_clause.where_clause.clause,
                    where_ir_builder,
                    lhs_param_to_node, lhs_assgn_to_node,
                    rhs_param_to_node, rhs_assgn_to_node,
                    symbol_tables, type_namespace, propensity_table,
                    rhs_node_locs)

            emit(ir_builder, build(where_ir_builder))

            n_r = _every_schedule_expr(every_clause.n_r, propensity_table, rule_name, "n_r")
            p_r = _every_schedule_expr(every_clause.p_r, propensity_table, rule_name, "p_r")
            t_clock = _every_schedule_expr(every_clause.t_clock, propensity_table, rule_name, "T_clock")
            n_r = n_r === nothing ? "1" : "static_cast<int>($n_r)"
            p_r = p_r === nothing ? "0.0" : "static_cast<double>($p_r)"
            t_clock = t_clock === nothing ? "0.0" : "static_cast<double>($t_clock)"
            emit(ir_builder, ",\n$n_r, $p_r, $t_clock")

            if every_clause.additive
                principal_type = first(k[2] for (k, v) in lhs_symbol if v == 1)
                emit(ir_builder, ",\n" * ir_every_combine(type_namespace, principal_type, rhs_assgn_to_node))
            end
        elseif modify_clause isa SolveClauseNode
            solve_clause = modify_clause

            # ── Collect read-only ODE vars for num_vars computation ──
            # First, gather all vars referenced in ODE RHS expressions
            # AND DefinitionNode values (temporary definitions within the
            # solve clause that may reference dynamic vertex attributes).
            all_ode_ref = Set{String}()
            for assgn_node in solve_clause.clause
                if assgn_node isa ODENode
                    union!(all_ode_ref,
                           collect_ode_referenced_vars(assgn_node.value, lhs_assgn_to_node))
                elseif assgn_node isa DefinitionNode
                    union!(all_ode_ref,
                           collect_ode_referenced_vars(assgn_node.value, lhs_assgn_to_node))
                end
            end

            # FixedList attributes are stack arrays and cannot be bound into the
            # ODE solver; only scalars and TFTensor attributes can.
            let ode_vars = Set{String}(all_ode_ref)
                for bv_node in solve_clause.variables
                    ode_var = bv_node.value[1]
                    push!(ode_vars, ode_var isa IndexAccessNode ? get_value(ode_var.object) : get_value(ode_var))
                end
                for v in ode_vars
                    if haskey(lhs_assgn_to_node, v) && is_tf_array_type(lhs_assgn_to_node[v].cpp_var.type)
                        throw("ODE variable '$v' is a FixedList (stack array); ODE solving supports " *
                              "scalars and TFTensor attributes only. Declare the attribute as TFTensor<<...>>.")
                    end
                end
            end

            # Compute num_vars: solving vars + read-only vars.
            # For tensor bindings, count tensor_size + position sync elements.
            # For scalar bindings, count 1 as before.
            # For indexed bindings (dpos[0] := D(im_pos[0], t)), avoid double-counting.
            num_vars = 0
            seen_dep_vars = Set()
            # Count solving variables
            for bv_node in solve_clause.variables
                ode_var = bv_node.value[1]
                dep_vars_name = ode_var isa IndexAccessNode ? get_value(ode_var.object) : get_value(ode_var)
                if dep_vars_name in seen_dep_vars
                    continue
                end
                push!(seen_dep_vars, dep_vars_name)
                if haskey(lhs_assgn_to_node, dep_vars_name)
                    rp = lhs_assgn_to_node[dep_vars_name]
                    tsize = rp.cpp_var.tensor_size
                    if rp.cpp_var.type == "torch::Tensor" && tsize > 0
                        if rp.cpp_var.index == 1
                            # Position solving tensor now binds ONLY the deep-copied
                            # `.position[]` array (min(3,tsize) vars) — the aliasable
                            # torch::Tensor pointer is no longer registered.
                            num_vars += min(3, tsize)
                        else
                            # tensor elements + position sync
                            num_vars += tsize + min(3, tsize)
                        end
                    else
                        num_vars += 1
                    end
                else
                    num_vars += 1
                end
            end
            # Count read-only variables (referenced in ODE RHS but not solving)
            for ro_var in setdiff(all_ode_ref, seen_dep_vars)
                if haskey(lhs_assgn_to_node, ro_var)
                    rp = lhs_assgn_to_node[ro_var]
                    tsize = rp.cpp_var.tensor_size
                    if rp.cpp_var.type == "torch::Tensor" && tsize > 0
                        if rp.cpp_var.index == 1
                            # Position read-only tensor still binds its elements.
                            num_vars += tsize
                        else
                            # Non-position read-only tensor (e.g. Direction) is read
                            # as a live constant and is NOT bound into the solver, so
                            # it contributes no SUNDIALS variables. Must mirror the
                            # binding logic in ir_solve_variable_binding or the
                            # declared num_vars will exceed the varset.insert count
                            # and mislocate the tau/propensity equations -> heap
                            # corruption.
                        end
                    else
                        num_vars += 1
                    end
                end
            end

            emit(ir_builder, "$num_vars,")
            solve_ir_builder = IRBuilder([])
            ir_solve_clause!(solve_clause,
                             lhs_assgn_to_node,
                             rhs_assgn_to_node,
                             solve_ir_builder, type_namespace,
                             propensity_table, symbol_tables, lhs_param_to_node)

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

    function _emit_undirected_edge_old(ir_builder, graph_name, each_node_info, rhs_count,
            rhs_node_registry, rhs_node_key, type_namespace, graph_symbol=Dict())

        vert_0 = each_node_info[2].left_vertex
        node_0_name = get_value(vert_0.name)
        # rhs_node_count_0 = rhs_count
        rhs_node_count_0 = graph_symbol[(node_0_name, node_type_0)]

        if !(node_0_name in rhs_node_registry)
            node_type_0 = get_value(vert_0.type.name)
            add_node_ir_0 = "$graph_name.addNode({$rhs_node_count_0, {$type_namespace::$node_type_0{} }});\n"
            rhs_count += 1
            emit(ir_builder, add_node_ir_0)
            rhs_node_key[node_0_name] = rhs_node_count_0
            push!(rhs_node_registry, node_0_name)
        else
            # rhs_node_count_0 = rhs_node_key[node_0_name]
            # rhs_node_count_0 = rhs_node_key[node_0_name]
            rhs_node_count_0 = graph_symbol[(node_0_name, node_type_0)]
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


    function _emit_undirected_edge(ir_builder, graph_name, each_node_info, rhs_count,
            rhs_node_registry, rhs_node_key, type_namespace, graph_symbol=Dict())

        vert_0 = each_node_info[2].left_vertex
        node_0_name = get_value(vert_0.name)
        node_type_0 = get_value(vert_0.type.name)
        rhs_node_count_0 = rhs_count
        if !(node_0_name in rhs_node_registry)
            add_node_ir_0 = "$graph_name.addNode({$rhs_node_count_0, {$type_namespace::$node_type_0{} }});\n"
            rhs_count += 1
            emit(ir_builder, add_node_ir_0)
            rhs_node_key[node_0_name] = rhs_node_count_0
            push!(rhs_node_registry, node_0_name)
        end

        rhs_node_count_1 = rhs_count
        vert_1 = each_node_info[2].right_vertex
        node_type_1 = get_value(vert_1.type.name)
        node_1_name = get_value(vert_1.name)
        rhs_node_count_1 = graph_symbol[(node_1_name, node_type_1)]
        if !(node_1_name in rhs_node_registry)
            node_type_1 = get_value(vert_1.type.name)
            # rhs_node_key[node_1_name] = rhs_node_count_1
            add_node_ir_1 = "$graph_name.addNode({$rhs_node_count_1, {$type_namespace::$node_type_1{} }});\n"
            rhs_count += 1
            emit(ir_builder, add_node_ir_1)
            push!(rhs_node_registry, node_1_name)
        end

        # Create an edge
        edge_ir = "$graph_name.addEdge($rhs_node_count_0, $rhs_node_count_1);\n"
        emit(ir_builder, edge_ir)
        return rhs_count
    end

    function emit_add_nodes(ir_builder, graph_name, node_types, type_namespace, graph_symbol=Dict())
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

            # node_type = get_value(each_node_info[2].type.name)
            # node_name = get_value(each_node_info[2].name)

            if !((node_name, node_type) in rhs_node_registry)
                # Single Node
                rhs_node_count = rhs_count
                node_type = get_value(each_node_info[2].type.name)
                node_name = get_value(each_node_info[2].name)

                # Grab from symbol table
                cur_count = graph_symbol[(node_name, node_type)]
                # add_node_ir = "$graph_name.addNode({$rhs_node_count, {$type_namespace::$node_type{} }});\n"
                add_node_ir = "$graph_name.addNode({$cur_count, {$type_namespace::$node_type{} }});\n"

                # Add to registry
                push!(rhs_node_registry, node_name)

                emit(ir_builder, add_node_ir)
                rhs_count += 1
            end

            if each_node_info[2] isa UndirectedTypeEdgeNode
                # Undirected Edge
                rhs_count = _emit_undirected_edge(ir_builder, graph_name, each_node_info,
                                      rhs_count, rhs_node_registry, rhs_node_key,
                                      type_namespace, graph_symbol)
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
