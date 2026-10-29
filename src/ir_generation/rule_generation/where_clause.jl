# Where-clause IR generation (variable declarations, autograd pre-scan).
# Included into module IRRuleGeneration (see ../ir_rule_generation.jl).

    function ir_where_left_clause!(lhs, type_namespace, symbol_tables,
                                   propensity_table, ir_builder)
        """
        IR Where left clause adds nodes attributes to propensity table.
        From the left hand side of the rule, we extract the node attributes and add them to the propensity table.
        """

        for (ix, param) in enumerate(lhs)

            node_loc = param[1]
            node_type = param[3]
            attr_type = param[4][1]
            attr_name = param[4][2]
            attr_loc = param[4][3]

            user_param_name = param[end][1]

            ir = nothing
            if is_list_type_tensor(attr_type) && attr_name == "Position"
                # The torch::Tensor Position goes STALE after ODE integration:
                # CVODE writes only the mirrored .position[] double array, never
                # the tensor (which stays frozen at the node's birth/seed value).
                # Bind a SHALLOW alias to the node's tensor (NOT a clone) so that
                # intra-block read-after-write still works (a where block may write
                # <var>[i] and later read it back — e.g. creation computing a unit
                # vector from a position it just set), then sync the aliased tensor
                # from the live .position[] array so reads of a matched/preserved
                # node use its current integrated position instead of the frozen
                # birth/seed value. Fixes boundary-catastrophe teleport-to-seed
                # without breaking creation's direction vectors.
                ir = "$attr_type $user_param_name = std::get<$type_namespace::$node_type>(lhs[m1[$node_loc]].data).$attr_name;\n" *
                     "for (int _i = 0; _i < 3 && _i < $user_param_name.numel(); _i++) { $user_param_name[_i] = lhs[m1[$node_loc]].position[_i]; }\n"
            elseif is_list_type_tensor(attr_type)
                # Non-Position tensor attribute (e.g. Direction) — read from .data
                ir = "$attr_type $user_param_name = std::get<$type_namespace::$node_type>(lhs[m1[$node_loc]].data).$attr_name;\n"
            elseif attr_loc <= 3 && attr_name == "Position"
                # Only the spatial Position attribute is mirrored in .position[]
                ir = "$attr_type $user_param_name = lhs[m1[$node_loc]].position[$(attr_loc-1)];"
            else
                # All other scalar attributes live in .data
                ir = "$attr_type $user_param_name = std::get<$type_namespace::$node_type>(lhs[m1[$node_loc]].data).$attr_name;\n"
            end

            # ir_prop = "$attr_type $user_param_name = std::get<$type_namespace::$node_type>(lhs[m1[$node_loc]].data).$attr_name;\n"
            propensity_table[user_param_name] = ir
            propensity_table["var_local_table"]["rule_lhs"][user_param_name] = ir

        end

    end

    function ir_where_right_clause!(
                rhs, type_namespace, symbol_tables,
                propensity_table, ir_builder
            )
        """
        IR Where right clause adds nodes attributes to propensity table.
        to be used in the where clause of the rule.
        """

        rhs_table = Dict()

        for (ix, param) in enumerate(rhs)

            node_loc = param[1]
            node_type = param[3]
            attr_type = param[4][1]
            attr_name = param[4][2]
            attr_loc = param[4][3]

            user_param_name = param[end][1]

            node_name = "rhs_node_"*string(node_loc)

            if !(node_name in collect(keys(rhs_table)))
                assign_ir = "$type_namespace::$node_type $node_name = std::get<$type_namespace::$node_type>(rhs[m2[$node_loc]].data);"
                rhs_table[node_name] = assign_ir
            end

            ir = nothing
            if is_list_type_tensor(attr_type) && attr_name == "Position"
                # See ir_where_left_clause!: bind a SHALLOW alias to the node's
                # tensor (NOT a clone) so intra-block read-after-write works (e.g.
                # creation writes pos_pos[i] then reads it back to build a unit
                # vector), then sync it from the live .position[] array so reads of
                # a preserved node use its current integrated position.
                ir = "$attr_type $user_param_name = std::get<$type_namespace::$node_type>(rhs[m2[$node_loc]].data).$attr_name;\n" *
                     "for (int _i = 0; _i < 3 && _i < $user_param_name.numel(); _i++) { $user_param_name[_i] = rhs[m2[$node_loc]].position[_i]; }\n"
            elseif is_list_type_tensor(attr_type)
                # Non-Position tensor attribute (e.g. Direction) — read from .data
                ir = "$attr_type $user_param_name = std::get<$type_namespace::$node_type>(rhs[m2[$node_loc]].data).$attr_name;\n"
            elseif attr_loc <= 3 && attr_name == "Position"
                # Only the spatial Position attribute is mirrored in .position[]
                ir = "$attr_type $user_param_name = rhs[m2[$node_loc]].position[$(attr_loc-1)];"
            else
                # All other scalar attributes live in .data
                ir = "$attr_type $user_param_name = std::get<$type_namespace::$node_type>(rhs[m2[$node_loc]].data).$attr_name;\n"
            end

            # ir_prop = "$attr_type $user_param_name = std::get<$type_namespace::$node_type>(lhs[m[$node_loc]].data).$attr_name;\n"
            # propensity_table["var_local_table"]["rule_rhs"][user_param_name] = ir_prop
            propensity_table["var_local_table"]["rule_rhs"][user_param_name] = ir

        end

        rhs_table
    end

    """
        collect_grad_vars(where_clause) -> Set{String}

    Pre-scan pass: walk all nodes in a where clause body and collect every
    variable name passed as an argument to a `grad(x)` call.  These are the
    tensors that need `requires_grad_(true)` injected at their declaration site
    so that the computation graph is built correctly before `backward()` is
    called.
    """
    function collect_grad_vars(where_clause)
        result = Set{String}()
        for node in where_clause
            _scan_grad_calls!(node, result)
        end
        return result
    end

    function _scan_grad_calls!(node, result)
        if node isa DefinitionNode
            _scan_grad_calls!(node.value, result)
        elseif node isa CallNode
            println("Scanning CallNode for grad calls: ", node.function_node.token)
            fname = get_value(node.function_node)
            if fname == "grad" && !isempty(node.args)
                arg = node.args[1]
                if arg isa IdentifierNode
                    push!(result, get_value(arg))
                end
            else
                for a in node.args
                    _scan_grad_calls!(a, result)
                end
            end
        elseif node isa BinaryOpNode
            _scan_grad_calls!(node.lhs, result)
            _scan_grad_calls!(node.rhs, result)
        elseif node isa UnaryOpNode
            _scan_grad_calls!(node.operand, result)
        elseif node isa GroupNode
            _scan_grad_calls!(node.expression, result)
        elseif node isa IndexAccessNode
            _scan_grad_calls!(node.object, result)
        elseif node isa ArrayLiteralNode
            for e in node.elements
                _scan_grad_calls!(e, result)
            end
        end
        # Leaf nodes (IdentifierNode, LiteralNode, etc.) — nothing to do
    end

    function ir_where_clause!(where_clause, ir_builder,
            lhs_param_to_node, lhs_assgn_to_node, rhs_param_to_node,
            rhs_assgn_to_node, symbol_tables, type_namespace, propensity_table,
            rhs_node_locs=nothing
        )
        """
        Handles where clause
        """

        where_hdr = "[&](auto &lhs, auto &rhs, auto &m1, auto &m2) {"
        emit(ir_builder, where_hdr)

        # ── Pre-scan: find all variables passed to grad(...) ──
        grad_vars = collect_grad_vars(where_clause)
        if !isempty(grad_vars)
            println("  Autograd pre-scan — grad vars detected: ", grad_vars)
        end

        # These functions will populate the propensity table
        # with the local variables.
        ir_where_left_clause!(lhs_param_to_node, type_namespace,
                              symbol_tables, propensity_table, ir_builder)

        # For LHS-bound attributes that need grad, emit requires_grad_ immediately
        # after they are declared by ir_where_left_clause!.
        for (_, _, _, _, user_param_name) in lhs_param_to_node
            vname = user_param_name[1]
            if vname in grad_vars
                emit(ir_builder, "$vname.requires_grad_(true);")
                println("  Injected requires_grad_(true) for LHS attr: ", vname)
            end
        end


        rhs_table = ir_where_right_clause!(rhs_param_to_node, type_namespace,
                                           symbol_tables, propensity_table, ir_builder)

        # NOTE: This is where i might add if the user doesn't define an attribute
        # We can just add it to the propensity table with the default value from the symbol table.
        lhs_user_param_names = map(x -> x[end][1], lhs_param_to_node)
        rhs_user_param_names = map(x -> x[end][1], rhs_param_to_node)

        # if the rhs name is equal to the lhs, we should just assign it to that automatically
        # This is because in the where clause, the user might want to use the same variable name for both the lhs and rhs, and we can just assume that they are the same variable. We can check if there are any variable names that are the same in both lhs and rhs, and if so, we can just assign them to be the same in the propensity table.

        # Build a set of shared names between LHS and RHS
        shared_names = intersect(Set(lhs_user_param_names), Set(rhs_user_param_names))

        if length(shared_names) > 0
            # Iterate over every RHS parameter entry (from the zipped list) so that
            # duplicate names across multiple RHS nodes are all handled, not just
            # the last one stored in the rhs_assgn_to_node dict.
            for (node_loc, node_param, node_type, attr, user_param_name) in rhs_param_to_node
                var_name = user_param_name[1]
                if !(var_name in shared_names)
                    continue
                end

                # Build a RuleParam for this specific RHS entry
                attr_tuple = attr
                tensor_size = length(attr_tuple) >= 4 ? attr_tuple[4] : 0
                rhs_cpp_var = CPPVariable(attr_tuple[1], attr_tuple[2], attr_tuple[3], tensor_size)
                rhs_rule_param = RuleParam(node_loc, node_param, node_type, rhs_cpp_var)

                rhs_ir = create_cpp_var(rhs_rule_param)
                lhs_ir = create_cpp_var(lhs_assgn_to_node[var_name], "lhs", "m1")

                # setting them equal. For non-position torch::Tensor attributes
                # (Direction, Count) deep-copy with .clone(): a plain `rhs = lhs`
                # shallow-copies the tensor so both nodes share storage. When DGGML
                # later mutates one node's tensor (e.g. a growing MT tip's Direction)
                # every node sharing that storage changes too -> tips end up sharing
                # one direction that is overwritten each step (vibration instead of
                # straight growth). Position (index 1) is intentionally left shallow:
                # the deep-copied spatial `.position[]` array is the source of truth,
                # and later `pos_pos`-style snapshots rely on aliasing to observe the
                # in-place indexed writes when computing unit vectors.
                if rhs_rule_param.cpp_var.type == "torch::Tensor" && rhs_rule_param.cpp_var.index != 1
                    ir = "$rhs_ir = $lhs_ir.clone();\n"
                else
                    ir = "$rhs_ir = $lhs_ir;\n"
                end
                emit(ir_builder, ir)

                if is_list_type_tensor(rhs_rule_param.cpp_var.type) && rhs_rule_param.cpp_var.index == 1
                    # For Position tensor attribute, also copy spatial .position from LHS to RHS
                    lhs_node_idx = lhs_assgn_to_node[var_name].index
                    rhs_node_idx = rhs_rule_param.index
                    pos_ir = "std::copy(std::begin(lhs[m1[$lhs_node_idx]].position), std::end(lhs[m1[$lhs_node_idx]].position), std::begin(rhs[m2[$rhs_node_idx]].position));\n"
                    emit(ir_builder, pos_ir)
                elseif !is_list_type_tensor(rhs_rule_param.cpp_var.type) && rhs_rule_param.cpp_var.index <= 3
                    rhs_pos = create_pos_cpp_var(rhs_rule_param, "rhs", "m2")
                    lhs_pos = create_pos_cpp_var(lhs_assgn_to_node[var_name], "lhs", "m1")
                    pos_ir = "$rhs_pos = $lhs_pos;\n"
                    emit(ir_builder, pos_ir)
                end
            end
        end

        # Implicit LHS to RHS carry-over
        # A node that persists across the rewrite (same node name on both the LHS and
        # RHS patterns; a retyped node keeps only attributes both types declare) keeps
        # any attribute the user listed on the LHS
        # but did NOT rebind on the RHS. RHS nodes are freshly constructed, so an
        # omitted field would otherwise be default-initialised and the persistent
        # node would silently lose that state. Auto-copying it from the matched LHS
        # node lets a rule omit unchanged RHS bindings instead of re-listing every
        # preserved field. Attributes rebound on the RHS pattern (handled by the
        # shared-name block above and/or the where-clause assignments) are skipped,
        # so a rule that lists every field emits exactly the same code as before.

        # node_name -> (rhs_loc, rhs_type). Prefer the full pattern-derived map
        # (covers nodes with no RHS bindings); fall back to bindings if absent.
        rhs_node_of = Dict()
        if rhs_node_locs !== nothing
            rhs_node_of = rhs_node_locs
        else
            for (node_loc, node_name, node_type, _, _) in rhs_param_to_node
                if !haskey(rhs_node_of, node_name)
                    rhs_node_of[node_name] = (node_loc, node_type)
                end
            end
        end

        # (node_name, attr_name) pairs already bound on the RHS pattern
        rhs_bound_attrs = Set{Tuple{Any,Any}}()
        for (_, node_name, _, attr_tuple, _) in rhs_param_to_node
            push!(rhs_bound_attrs, (node_name, attr_tuple[2]))
        end

        carried_over = Set{Tuple{Any,Any}}()  # guard against duplicate emits
        for (lhs_loc, lhs_name, lhs_type, lhs_attr, _) in lhs_param_to_node
            attr_name = lhs_attr[2]
            haskey(rhs_node_of, lhs_name) || continue
            (rhs_loc, rhs_type) = rhs_node_of[lhs_name]
            # Skip attributes the RHS already binds (shared-name / where handled them).
            ((lhs_name, attr_name) in rhs_bound_attrs) && continue
            ((lhs_name, attr_name) in carried_over) && continue

            tensor_size = length(lhs_attr) >= 4 ? lhs_attr[4] : 0
            cpp_var = CPPVariable(lhs_attr[1], lhs_attr[2], lhs_attr[3], tensor_size)
            rhs_cpp_var = cpp_var
            if rhs_type != lhs_type
                # Retyped node (same name, different type): carry over only the
                # attributes the new type also declares, with the same C++ type.
                # The slot index is resolved in the RHS type's own field list.
                haskey(symbol_tables, rhs_type) || continue
                rhs_ix = 0
                rhs_fields = symbol_tables[rhs_type]
                for ix in 1:length(rhs_fields)
                    if rhs_fields[ix][2] == attr_name &&
                       convert_type_name(rhs_fields[ix][1]) == lhs_attr[1]
                        rhs_ix = ix
                        break
                    end
                end
                rhs_ix == 0 && continue
                rhs_cpp_var = CPPVariable(lhs_attr[1], lhs_attr[2], rhs_ix, tensor_size)
            end
            push!(carried_over, (lhs_name, attr_name))

            lhs_rp = RuleParam(lhs_loc, lhs_name, lhs_type, cpp_var)
            rhs_rp = RuleParam(rhs_loc, lhs_name, rhs_type, rhs_cpp_var)

            rhs_ir = create_cpp_var(rhs_rp)               # rhs[m2[rhs_loc]].data.Attr
            lhs_ir = create_cpp_var(lhs_rp, "lhs", "m1")  # lhs[m1[lhs_loc]].data.Attr

            # Deep-copy non-position tensors so persistent nodes don't share storage
            # (see the shared-name block above for why aliasing corrupts growth).
            if cpp_var.type == "torch::Tensor" && cpp_var.index != 1
                emit(ir_builder, "$rhs_ir = $lhs_ir.clone();\n")
            else
                emit(ir_builder, "$rhs_ir = $lhs_ir;\n")
            end

            if is_list_type_tensor(cpp_var.type) && cpp_var.index == 1
                # Position tensor: also carry the authoritative spatial .position[].
                pos_ir = "std::copy(std::begin(lhs[m1[$lhs_loc]].position), std::end(lhs[m1[$lhs_loc]].position), std::begin(rhs[m2[$rhs_loc]].position));\n"
                emit(ir_builder, pos_ir)
            elseif !is_list_type_tensor(cpp_var.type) && cpp_var.index <= 3
                rhs_pos = create_pos_cpp_var(rhs_rp, "rhs", "m2")
                lhs_pos = create_pos_cpp_var(lhs_rp, "lhs", "m1")
                emit(ir_builder, "$rhs_pos = $lhs_pos;\n")
            end
        end

        # Tracks torch::Tensor node attributes that have already been detached
        # (given independent storage via .clone()) in this rule body, so we only
        # emit the detach once before the first in-place indexed write.
        detached_tensors = Set{String}()

        map(
            (assign_node) -> begin

                ir_where_assignment!(assign_node, type_namespace, rhs_assgn_to_node,
                                     propensity_table, ir_builder, grad_vars, detached_tensors)
                            end, 
            where_clause
        )

        emit(ir_builder, "}")
    end
