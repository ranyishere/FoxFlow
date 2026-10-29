# Rule parameter bindings: CPPVariable/RuleParam and node/attribute linking.
# Included into module IRRuleGeneration (see ../ir_rule_generation.jl).

    struct CPPVariable
        type::String
        name::String
        index::Int # Index is the attr pos as it appears in c++
        tensor_size::Int # Number of elements if tensor, 0 otherwise
    end

    CPPVariable(type, name, index) = CPPVariable(type, name, index, 0)

    struct RuleParam
        index::Int # the index as it appear sin the lhs
        name::String
        type::String
        cpp_var::CPPVariable
    end

    function create_cpp_var(rule_param::RuleParam, side="rhs", match_arr="m2")
        cpp_var = rule_param.cpp_var
        ns = _type_namespace[]
        cpp_var_ir = "std::get<$(ns)::$(rule_param.type)>($(side)[$(match_arr)[ $(rule_param.index) ]].data).$(cpp_var.name)"
        return cpp_var_ir
    end

    function create_pos_cpp_var(rule_param::RuleParam, side="rhs", match_arr="m2")
        # pos = cpp_var.index
        cpp_var = rule_param.cpp_var
        pos = cpp_var.index

        param_ix = rule_param.index

        cpp_var_ir = "$(side)[$(match_arr)[ $(param_ix) ]].position[$(pos-1)]"
        return cpp_var_ir
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
                else
                    # Same vertex as last — skip to avoid duplicate
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

    function link_param_to_nodes_positional(param, node_types, symbol_tables, node_names)
        """
        POSITIONAL binding (legacy / fallback).

        Links each parameter to its corresponding node by POSITION: it flattens
        every field of every pattern node in declaration order and zips that list
        1:1 against the declared parameters. This requires each node to list ALL of
        its fields, in order, and ignores the `::Attribute` annotation for slot
        resolution. Kept as a fallback for the name-based `link_param_to_nodes`.

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

        println("param: ", param)

        tmp_count = 0
        user_param_names = map(
            (x) -> begin

                if (x isa NamedParameterNode)
                    tmp_count += 1
                    (get_value(x.name), tmp_count)
                else
                    tmp_count += 1
                    (get_value(x), tmp_count)
                end

            end,
            param.token
        )

        # println("user_param_names: ", user_param_names)
        # exit(0)

        ix = 0
        node_attrs = []

        for node_type in node_types

            if !(node_type in collect(keys(symbol_tables)))
                println("Node type $node_type not found in symbol tables.")
                throw("Node type $node_type not found in symbol tables.")
                continue
            end

            # Fetches all the attributes for this node type from the symbol table
            total_node_attr = symbol_tables[node_type]

            # Map this to user_param_names by position
            tmp_attr = []
            for ix_node_attr in 1:length(total_node_attr)

                node_attr = total_node_attr[ix_node_attr]

                # println("node_attr: ", node_attr)
                # println("ix_node_attr: ", ix_node_attr)
                # exit(0)

                is_tensor = is_list_type_tensor(node_attr[1])
                attr_type = convert_type_name(node_attr[1])

                attr_name = node_attr[2]

                # Extract tensor size from symbol table is_list info
                # is_list_info[2] can be:
                #   - A single node (e.g., IntegerNode("3")) for 1D FixedList<<3, Float>>
                #   - A list of strings (e.g., ["3", "4"]) for multi-dim FixedList<<(3,4), Float>>
                tensor_size = 0
                if is_tensor && length(node_attr) >= 3
                    is_list_info = node_attr[3]
                    if is_list_info[1] == true
                        dim_info = is_list_info[2]
                        if dim_info isa AbstractVector || dim_info isa AbstractArray
                            # Multi-dimensional: product of all dimensions
                            tensor_size = prod([parse(Int64, d) for d in dim_info])
                        else
                            # Single dimension
                            tensor_size = parse(Int64, get_value(dim_info))
                        end
                    end
                end

                push!(tmp_attr, (attr_type, attr_name, ix_node_attr, tensor_size))
            end

            push!(node_attrs, tmp_attr)

        end

        # Validate: user must declare exactly as many parameters per node
        #    as the type has fields.  A mismatch causes a silent misalignment
        #    in the positional zip below (fields from node N bleed into node N+1).
        total_type_fields = sum(length(a) for a in node_attrs)
        total_user_params = length(user_param_names)
        if total_type_fields != total_user_params
            # Build a per-node breakdown for a helpful error message
            err_details = []
            for (ix, ntype) in enumerate(node_types)
                nfields = length(node_attrs[ix])
                nname = ix <= length(node_names) ? node_names[ix] : "node_$ix"
                field_names = [a[2] for a in node_attrs[ix]]
                push!(err_details,
                    "  node '$(nname)' (type $(ntype)): type has $(nfields) field(s) [$(join(field_names, ", "))]")
            end
            throw(
                "Parameter count mismatch in rule pattern: " *
                "type definitions expect $(total_type_fields) parameter(s) total " *
                "but the rule declares $(total_user_params).\n" *
                "Each node in the pattern must list ALL fields of its type.\n" *
                join(err_details, "\n") * "\n" *
                "Declared parameters: [$(join([p[1] for p in user_param_names], ", "))]"
            )
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
        count = 0
        for (node_loc, node_param, node_type, attr, user_param_name) in zipped

            param_name_to_node[user_param_name] = (node_loc, node_param, node_type, attr)
            count += 1
        end

        setdiff(Set(user_param_names), Set(keys(param_name_to_node))) != Set() &&
            println("Missing parameters: ",
                    setdiff(Set(user_param_names), Set(keys(param_name_to_node))))

        oof0 = keys(param_name_to_node)
        oof1 = Set(user_param_names)

        # Where is this error from?
        keys(param_name_to_node) == Set(user_param_names) || throw("Parameter names do not match node parameters.")

        linked_params = Dict()
        linked_params_new = Dict()
        for param in param_name_to_node
            key = param[1][1]
            # cpp_var_name = param[2][end][1]
            #
            # println("key: ", key)
            # println("param[2]: ", param[2])
            #
            attr_tuple = param[2][end]
            # attr_tuple is now (attr_type, attr_name, ix_node_attr, tensor_size)
            tensor_size = length(attr_tuple) >= 4 ? attr_tuple[4] : 0
            cpp_var = CPPVariable(attr_tuple[1], attr_tuple[2], attr_tuple[3], tensor_size)

            rule_param = RuleParam(param[2][1], param[2][2], param[2][3], cpp_var)

            println("rule_param ========> ", rule_param)

            linked_params[key] = param[2]
            linked_params_new[key] = rule_param
        end

        # zipped, param_name_to_node
        zipped, linked_params_new

    end

    function link_param_to_nodes(param, node_types, symbol_tables, node_names)
        """
        NAME-BASED binding.

        Resolves each `alias : node::Attribute` parameter to its attribute slot by
        NAME (node name + attribute name) rather than by position. This lets a rule
        omit fields it does not touch: only the attributes actually referenced need
        be listed, in any order. Produces the exact same output structures as
        link_param_to_nodes_positional (zipped, linked_params_new) so all downstream
        codegen is unchanged; when a rule lists every field in declaration order the
        emitted output is byte-identical to the positional path.

        Args:
            param: The ParameterNode containing the `<< ... >>` bindings.
            node_types: Types of the pattern nodes (parallel to node_names).
            symbol_tables: Type -> field list symbol tables.
            node_names: Names of the pattern nodes, in pattern order.

        Returns:
            (zipped, linked_params_new)
        """

        println("param: ", param)

        # Assign each unique pattern node a 1-based location (its m1/m2 index)
        #    in first-appearance order, and remember its type.
        node_loc_of = Dict{Any, Int}()
        type_of_node = Dict{Any, Any}()
        unique_count = 0
        for (ix, node_name) in enumerate(node_names)
            if !haskey(node_loc_of, node_name)
                unique_count += 1
                node_loc_of[node_name] = unique_count
                type_of_node[node_name] = node_types[ix]
            end
        end

        # Build, per node type, a lookup: attr_name -> (attr_type, attr_name,
        #    ix_node_attr, tensor_size).  Mirrors the tensor-size extraction in the
        #    positional path so downstream cpp_var info is identical. ──
        attr_table_cache = Dict{Any, Any}()
        function attr_table_for(node_type)
            if haskey(attr_table_cache, node_type)
                return attr_table_cache[node_type]
            end
            if !(node_type in collect(keys(symbol_tables)))
                throw("Node type $node_type not found in symbol tables.")
            end
            total_node_attr = symbol_tables[node_type]
            tbl = Dict{Any, Any}()
            order = String[]
            for ix_node_attr in 1:length(total_node_attr)
                node_attr = total_node_attr[ix_node_attr]
                is_tensor = is_list_type_tensor(node_attr[1])
                attr_type = convert_type_name(node_attr[1])
                attr_name = node_attr[2]

                tensor_size = 0
                if is_tensor && length(node_attr) >= 3
                    is_list_info = node_attr[3]
                    if is_list_info[1] == true
                        dim_info = is_list_info[2]
                        if dim_info isa AbstractVector || dim_info isa AbstractArray
                            tensor_size = prod([parse(Int64, d) for d in dim_info])
                        else
                            tensor_size = parse(Int64, get_value(dim_info))
                        end
                    end
                end

                tbl[attr_name] = (attr_type, attr_name, ix_node_attr, tensor_size)
                push!(order, attr_name)
            end
            tbl["__order__"] = order
            attr_table_cache[node_type] = tbl
            return tbl
        end

        # ── Resolve each parameter by name. ──
        zipped = []
        param_name_to_node = Dict()
        seen_binding = Dict()  # (node_loc, attr_name) -> alias, to catch double binds
        tmp_count = 0
        for x in param.token
            tmp_count += 1

            if !(x isa NamedParameterNode) || x.parameter === nothing ||
               !(x.parameter isa IdentifierNode) || x.parameter.namespace === nothing
                throw("Name-based binding requires the `alias : node::Attribute` form; " *
                      "encountered a parameter without a node::attribute reference. " *
                      "(Swap the caller to link_param_to_nodes_positional to use the " *
                      "old positional syntax.)")
            end

            alias = get_value(x.name)
            node_name = x.parameter.namespace.position.value
            attr_name = get_value(x.parameter)

            if !haskey(node_loc_of, node_name)
                throw("Rule binding '$alias : $node_name::$attr_name' references node " *
                      "'$node_name', which is not one of the pattern nodes " *
                      "[$(join(unique(node_names), ", "))].")
            end
            node_loc = node_loc_of[node_name]
            node_type = type_of_node[node_name]
            tbl = attr_table_for(node_type)
            if !haskey(tbl, attr_name)
                throw("Rule binding '$alias : $node_name::$attr_name' references " *
                      "attribute '$attr_name', which is not a field of type " *
                      "'$node_type' [$(join(tbl["__order__"], ", "))].")
            end
            attr_tuple = tbl[attr_name]

            if haskey(seen_binding, (node_loc, attr_name))
                throw("Attribute '$node_name::$attr_name' is bound more than once in " *
                      "the same rule pattern (aliases " *
                      "'$(seen_binding[(node_loc, attr_name)])' and '$alias').")
            end
            seen_binding[(node_loc, attr_name)] = alias

            user_param_name = (alias, tmp_count)
            push!(zipped, (node_loc, node_name, node_type, attr_tuple, user_param_name))
            param_name_to_node[user_param_name] = (node_loc, node_name, node_type, attr_tuple)
        end

        # ── Build linked_params_new (alias -> RuleParam), identical shape to the
        #    positional path. ──
        linked_params_new = Dict()
        for param_pair in param_name_to_node
            key = param_pair[1][1]  # alias string
            attr_tuple = param_pair[2][end]
            tensor_size = length(attr_tuple) >= 4 ? attr_tuple[4] : 0
            cpp_var = CPPVariable(attr_tuple[1], attr_tuple[2], attr_tuple[3], tensor_size)
            rule_param = RuleParam(param_pair[2][1], param_pair[2][2], param_pair[2][3], cpp_var)
            linked_params_new[key] = rule_param
        end

        zipped, linked_params_new

    end

    function link_param_to_nodes_old(param, node_types, symbol_tables, node_names)
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

        tmp_count = 0
        user_param_names = map(
            (x) -> begin
                tmp_count += 1
                (get_value(x), tmp_count)
            end,
            param.token
        )

        ix = 0
        node_attrs = []

        for node_type in node_types

            if !(node_type in collect(keys(symbol_tables)))
                println("Node type $node_type not found in symbol tables.")
                throw("Node type $node_type not found in symbol tables.")
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
        count = 0
        for (node_loc, node_param, node_type, attr, user_param_name) in zipped

            param_name_to_node[user_param_name] = (node_loc, node_param, node_type, attr)
            count += 1
        end

        setdiff(Set(user_param_names), Set(keys(param_name_to_node))) != Set() &&
            println("Missing parameters: ",
                    setdiff(Set(user_param_names), Set(keys(param_name_to_node))))

        oof0 = keys(param_name_to_node)
        oof1 = Set(user_param_names)

        # Where is this error from?
        keys(param_name_to_node) == Set(user_param_names) || throw("Parameter names do not match node parameters.")

        linked_params = Dict()
        linked_params_new = Dict()
        for param in param_name_to_node
            key = param[1][1]
            # cpp_var_name = param[2][end][1]
            #
            println("key: ", key)
            println("param[2]: ", param[2])

            cpp_var = CPPVariable(param[2][end][1], param[2][end][2], param[2][end][3])
            rule_param = RuleParam(param[2][1], param[2][2], param[2][3], cpp_var)

            linked_params[key] = param[2]
            linked_params_new[key] = rule_param
        end

        # zipped, param_name_to_node
        zipped, linked_params_new
    end
