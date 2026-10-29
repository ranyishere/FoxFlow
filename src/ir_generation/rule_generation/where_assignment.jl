# Where-clause assignment IR generation and local function helpers.
# Included into module IRRuleGeneration (see ../ir_rule_generation.jl).

    function is_self_referential(code_line::String)
        # ^(\w+)      : Capture the variable name at the start
        # \s*=        : Match the equals sign
        # .*?         : Non-greedy match for any characters in between
        # \b\1\b      : Match the exact same word from the capture group
        pattern = r"^(\w+)\s*=.*?\b\1\b"

        return occursin(pattern, strip(code_line))
    end


    function _ir_local_fn_scalar_type(type_str)
        if type_str == "Float"       return "double"
        elseif type_str == "Integer" return "int64_t"
        elseif type_str == "FixedList" return "torch::Tensor"
        else return type_str
        end
    end

    function _ir_local_fn_return_type(ret_type_node)
        # ret_type_node is a TypeClassNode; its .name is an IdentifierNode
        return _ir_local_fn_scalar_type(get_value(ret_type_node.name))
    end

    function _ir_local_fn_param(arg_node)
        # FunctionArgNode has .name (IdentifierNode) and .type (TypeClassNode)
        name     = get_value(arg_node.name)
        type_str = get_value(arg_node.type.name)
        return "$(_ir_local_fn_scalar_type(type_str)) $name"
    end

    function ir_where_assignment!(assign_node, type_namespace,
                                 rhs_assgn_to_node,
                                 propensity_table, ir_builder,
                                 grad_vars=Set{String}(),
                                 detached_tensors=Set{String}())
        """
        Handles the assignment in the where clause.
        """

        # FIXME: if you assign a node to another node who is just assigned it will not work.

        # Bare call statement (e.g. backward(loss))
        # Parsed as a CallNode with no LHS.  Emit as a side-effecting statement.
        if assign_node isa CallNode
            ir_value = IRBuilder([])
            ir_definition!(GroupNode(assign_node), ir_value, propensity_table, ir_builder)
            ir_str = build_sameline(ir_value)
            emit(ir_builder, "$ir_str;")
            return
        end

        # Local function definition (emitted as a C++ lambda)
        if assign_node isa FunctionDefinitionNode
            fn_name    = get_value(assign_node.name)
            sig        = assign_node.signature
            ret_cpp    = _ir_local_fn_return_type(sig.output)
            params_cpp = join([_ir_local_fn_param(a) for a in sig.args], ", ")

            emit(ir_builder, "auto $fn_name = [&]($params_cpp) -> $ret_cpp {")

            # body statements
            if assign_node.body !== nothing
                for stmt in assign_node.body.expressions
                    stmt_name = get_value(stmt.name)
                    stmt_type = _ir_local_fn_scalar_type(get_value(stmt.type))
                    stmt_ir   = IRBuilder([])
                    ir_definition!(GroupNode(stmt.value), stmt_ir, propensity_table, ir_builder)
                    stmt_val  = build_sameline(stmt_ir)
                    emit(ir_builder, "    $stmt_type $stmt_name = $stmt_val;")
                end
            end

            # return statement
            if assign_node.fun_return !== nothing
                ret_ir  = IRBuilder([])
                ir_definition!(GroupNode(assign_node.fun_return.value), ret_ir, propensity_table, ir_builder)
                ret_str = build_sameline(ret_ir)
                emit(ir_builder, "    return $ret_str;")
            end

            emit(ir_builder, "};")
            return
        end

        if assign_node isa DefinitionNode
            name = get_value(assign_node.name)
            # FixedList<<N, Float>> -> TFArray<N>; TFTensor -> torch::Tensor;
            # FixedList without a size -> auto (deduced from the value).
            type = convert_type_node(assign_node.type; dimless="auto")
            def_is_array = get_value(assign_node.type.name) == "FixedList"
            value = assign_node.value

            # autodiff(expr, W) means "differentiate expr w.r.t. W".
            # We emit an immediately-invoked lambda so that requires_grad_(true)
            # is set on W *before* expr is evaluated, ensuring W is in the graph:
            #
            #   auto dW = ([&]() {
            #       W.requires_grad_(true);
            #       auto _fflow_ad_tmp = <expr>;
            #       _fflow_ad_tmp.backward();
            #       return W.grad();
            #   })();
            raw_value = value isa GroupNode ? value.expression : value
            if raw_value isa CallNode && get_value(raw_value.function_node) == "autodiff"
                if length(raw_value.args) != 2
                    throw("autodiff requires exactly 2 arguments: autodiff(expr, var)")
                end
                ad_expr_node = raw_value.args[1]
                ad_var_node  = raw_value.args[2]
                if !(ad_var_node isa IdentifierNode)
                    throw("Second argument to autodiff must be a plain variable name")
                end
                ad_var = get_value(ad_var_node)

                # Emit the forward expression into a temp string
                ad_expr_ir = IRBuilder([])
                ir_definition!(GroupNode(ad_expr_node), ad_expr_ir, propensity_table, ir_builder)
                ad_expr_str = build_sameline(ad_expr_ir)

                push!(propensity_table["var_local_table"]["rule_rhs"]["declared"], "$name")
                emit(ir_builder, "$type $name = ([&]() {")
                emit(ir_builder, "    $ad_var.requires_grad_(true);")
                emit(ir_builder, "    auto _fflow_ad_tmp = $ad_expr_str;")
                emit(ir_builder, "    _fflow_ad_tmp.backward();")
                emit(ir_builder, "    return $ad_var.grad();")
                emit(ir_builder, "})();")
                if name in grad_vars
                    emit(ir_builder, "$name.requires_grad_(true);")
                end
                return
            end

            if !(value isa GroupNode)
                value = GroupNode(value)
            end

            # Cant allow definition using same name found in propensity table
            if name in collect(keys(propensity_table["var_local_table"]["rule_rhs"]))
                throw("Error: Variable '$name' already defined in rhs propensity table. Please choose a different name.")
                exit(0)
            end

            if name in collect(keys(propensity_table["var_local_table"]["rule_lhs"]))
                throw("Error: Variable '$name' already defined in lhs propensity table. Please choose a different name.")
                exit(0)
            end

            println("processing definition node: ", name)

            # Generate ir
            ir_value = IRBuilder([])

            ARRAY_CTX[] = def_is_array
            try
                ir_definition!(value, ir_value, propensity_table, ir_builder)
            finally
                ARRAY_CTX[] = false
            end

            ir_value = build_sameline(ir_value)

            # Adds the definition to propensity table local scope.
            # propensity_table["var_local_table"]["rule_rhs"]["declared"][name] = "$name"
            push!(propensity_table["var_local_table"]["rule_rhs"]["declared"], "$name")
            println("added $name to propensity table declared variables: ", propensity_table["var_local_table"]["rule_rhs"]["declared"])

            # This seems like it broke something. Its no longer being output to the ir builder.
            # I think it might be because of the context that is being passed in.
            # I need to check if the context is being modified correctly inside ir_definition.

            ir = "$type $name = $ir_value;"
            emit(ir_builder, ir)
            # ── Autograd: inject requires_grad_(true) if this var is differentiated over ──
            if name in grad_vars
                emit(ir_builder, "$name.requires_grad_(true);")
                println("  Injected requires_grad_(true) for local def: ", name)
            end
        elseif assign_node.name isa IndexAccessNode
            # Indexed assignment: count[0] = 1
            # Extract the base identifier and index expressions
            idx_node = assign_node.name
            base_name = get_value(idx_node.object)

            if !(base_name in keys(rhs_assgn_to_node))
                throw("Error: Variable '$base_name' in indexed where clause assignment not found in rhs parameters.")
            end

            rhs_assgn_node = rhs_assgn_to_node[base_name]

            # Build the index expression string
            # Build the index expression strings, handling slices
            idx_parts = String[]
            for idx in idx_node.indices
                if idx isa SliceNode
                    push!(idx_parts, ir_slice_expr(idx, ir_builder, propensity_table, ir_builder))
                else
                    idx_ir = IRBuilder([])
                    ir_prop_expr!(idx isa GroupNode ? idx : GroupNode(idx), idx_ir, propensity_table, ir_builder, false)
                    push!(idx_parts, build_sameline(idx_ir))
                end
            end

            # Build the RHS value expression
            value = assign_node.value
            if !(value isa GroupNode)
                value = GroupNode(value)
            end
            ir_value = IRBuilder([])
            ir_prop_expr!(value, ir_value, propensity_table, ir_builder, false)
            ir_value_str = build_sameline(ir_value)

            # Emit assignment. Use .index_put_ for slices or multi-index, plain [] for single scalar index.
            cpp_var = create_cpp_var(rhs_assgn_node)

            # Detach tensor storage before an in-place indexed write. torch::Tensor
            # is shallow-copied when DGGML duplicates a node, so several nodes can
            # share one underlying buffer; writing X[i] in place would mutate all
            # of them. Cloning once gives this node independent storage (harmless
            # on a freshly created node whose tensor is still default zeros).
            # Position (index 1) is EXCLUDED: its indexed writes must remain visible
            # through existing shallow snapshots (e.g. `pos_pos`) so unit-vector
            # computations in the same body read the freshly written coordinates.
            if rhs_assgn_node.cpp_var.type == "torch::Tensor" && rhs_assgn_node.cpp_var.index != 1 && !(cpp_var in detached_tensors)
                emit(ir_builder, "$cpp_var = $cpp_var.clone();")
                push!(detached_tensors, cpp_var)
            end

            if has_slice(idx_node.indices)
                ir = "$cpp_var.index_put_({$(join(idx_parts, ", "))}, $ir_value_str);"
            elseif length(idx_parts) == 1
                ir = "$cpp_var[$(idx_parts[1])] = $ir_value_str;"
            else
                ir = "$cpp_var.index_put_({$(join(idx_parts, ", "))}, $ir_value_str);"
            end
            emit(ir_builder, ir)

            # Also sync .position for single element if it's the Position attribute (index 1)
            if is_list_type_tensor(rhs_assgn_node.cpp_var.type) && length(idx_parts) == 1 && rhs_assgn_node.cpp_var.index == 1
                rhs_node_idx = rhs_assgn_node.index
                pos_sync = "if ($(idx_parts[1]) < 3) { rhs[m2[$rhs_node_idx]].position[$(idx_parts[1])] = static_cast<double>($ir_value_str); }"
                emit(ir_builder, pos_sync)
            end
        else

            # Check if it is an assign node or just an intermediate expression
            name = get_value(assign_node.name)
            
            type = rhs_assgn_to_node[name].cpp_var.type

            value = assign_node.value

            if !(name in keys(rhs_assgn_to_node))
                throw("Error: Variable '$name' in where clause assignment not found in rhs parameters.
                      Please make sure to reference a variable defined in the rhs parameters.")
                exit(0)
            end

            # TODO: When doing this, maybe we add a guard rail to say did you mean to reference the node that you just defined in the where clause?
            # Because that is a common mistake that I can see happening.
            rhs_assgn_node = rhs_assgn_to_node[name]

            node_value = assign_node.value
            value = assign_node.value

            local_namespace = type_namespace

            if !(value isa GroupNode)
                value = GroupNode(value)
            end

            # Generate ir
            ir_value = IRBuilder([])
            ARRAY_CTX[] = is_tf_array_type(type)
            try
                ir_prop_expr!(value, ir_value, propensity_table, ir_builder, false)
            finally
                ARRAY_CTX[] = false
            end
            ir_value = build_sameline(ir_value)

            # Adds the definition to propensity table local scope.
            propensity_table["var_local_table"]["rule_rhs"][name] = "$name"

            # if occursin(name, ir_value)
            if is_self_referential(ir_value)
                throw(
                      "Error: Self-referential assignment detected for variable '$name'. This is not supported."
                      )
            end

            new_ir = "$type $name = $ir_value;"

            emit(ir_builder, new_ir)
            # exit(0)

            # assigned_param = rhs_assgn_node[4][3]
            # assigned_param = rhs_assgn_node.cpp_var.index

            # assigned_node = rhs_assgn_node[1]
            # assigned_node = rhs_assgn_node.index

            # rhs_attr = rhs_assgn_node[4][2]
            # rhs_attr = rhs_assgn_node.cpp_var[2]

            # node_type = rhs_assgn_node[3]

            # assigned_param_name = assigned_param
            # ir = "\t\tstd::get<$type_namespace::$node_type>(rhs[m2[$assigned_node]].data).$rhs_attr = $name;"

            # For non-position torch::Tensor attributes, deep-copy the assigned
            # value with .clone() so this node owns independent storage instead of
            # aliasing the local (or a value shared with another node). See the
            # indexed / carry-copy cases above: shared tensor storage makes
            # duplicated nodes mutate together (the MT-tip "vibration" bug).
            # Position (index 1) stays a plain assignment to preserve snapshot
            # aliasing used by unit-vector computations.
            if rhs_assgn_node.cpp_var.type == "torch::Tensor" && rhs_assgn_node.cpp_var.index != 1
                ir = "$(create_cpp_var(rhs_assgn_node)) = $name.clone();"
            else
                ir = "$(create_cpp_var(rhs_assgn_node)) = $name;"
            end

            assigned_param = rhs_assgn_node.cpp_var.index

            if is_list_type_tensor(rhs_assgn_node.cpp_var.type) && rhs_assgn_node.cpp_var.index == 1
                # For the Position attribute (index 1), sync spatial .position from tensor value
                rhs_node_idx = rhs_assgn_node.index
                pos_sync = "for (int _i = 0; _i < 3 && _i < $name.numel(); _i++) { rhs[m2[$rhs_node_idx]].position[_i] = $name[_i].template item<double>(); }"
                emit(ir_builder, pos_sync)
            elseif !is_list_type_tensor(rhs_assgn_node.cpp_var.type) && (assigned_param == 1 || assigned_param == 2 || assigned_param == 3)
                # For scalar attributes at positions 1-3, sync to .position
                ir_node_pos = create_pos_cpp_var(rhs_assgn_node, "rhs", "m2") * " = $name;"
                emit(ir_builder, ir_node_pos)
            end

            emit(ir_builder, ir)
        end
    end
