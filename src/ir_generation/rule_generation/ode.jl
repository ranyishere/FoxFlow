# ODE (solve clause) IR generation for SolvingRule.
# Included into module IRRuleGeneration (see ../ir_rule_generation.jl).

    """
        collect_ode_referenced_vars(expression, lhs_assgn_to_node) -> Set{String}

    Walk an ODE RHS expression tree and collect all IdentifierNode names that
    reference LHS-matched node attributes (keys in lhs_assgn_to_node).
    These are the variables the ODE RHS reads — both solving and read-only.
    """
    function collect_ode_referenced_vars(expression, lhs_assgn_to_node)
        result = Set{String}()
        _collect_ode_vars!(expression, lhs_assgn_to_node, result)
        return result
    end

    function _collect_ode_vars!(node, lhs_assgn_to_node, result)
        if node isa BinaryOpNode
            _collect_ode_vars!(node.lhs, lhs_assgn_to_node, result)
            _collect_ode_vars!(node.rhs, lhs_assgn_to_node, result)
        elseif node isa UnaryOpNode
            _collect_ode_vars!(node.operand, lhs_assgn_to_node, result)
        elseif node isa GroupNode
            _collect_ode_vars!(node.expression, lhs_assgn_to_node, result)
        elseif node isa CallNode
            for arg in node.args
                _collect_ode_vars!(arg, lhs_assgn_to_node, result)
            end
        elseif node isa IndexAccessNode
            if node.object isa IdentifierNode
                name = get_value(node.object)
                if haskey(lhs_assgn_to_node, name)
                    push!(result, name)
                end
            end
            for idx in node.indices
                _collect_ode_vars!(idx, lhs_assgn_to_node, result)
            end
        elseif node isa IdentifierNode
            name = get_value(node)
            if haskey(lhs_assgn_to_node, name)
                push!(result, name)
            end
        elseif node isa ArrayLiteralNode
            for el in node.elements
                _collect_ode_vars!(el, lhs_assgn_to_node, result)
            end
        end
        # LiteralNode, IntegerNode, FloatNode etc. — nothing to collect
    end

    # Symbolic ODE pretty-printer 
    # Walks the same ODE AST that traverse_ode_expr uses, but instead of
    # emitting C++ IR it returns a human-readable algebraic string like
    #   "11.11 * (P0 - P1)"
    # This is used to emit comments and runtime debug prints showing the
    # symbolic form of each ODE equation.
    function symbolic_ode_expr(expression, assgn_info, dep_vars, bv_to_dep)::String

        if expression isa BinaryOpNode
            l = symbolic_ode_expr(expression.lhs, assgn_info, dep_vars, bv_to_dep)
            op = expression.expression.position.value
            r = symbolic_ode_expr(expression.rhs, assgn_info, dep_vars, bv_to_dep)
            return "$l $op $r"

        elseif expression isa GroupNode
            inner = symbolic_ode_expr(expression.expression, assgn_info, dep_vars, bv_to_dep)
            return "($inner)"

        elseif expression isa UnaryOpNode
            op = get_value(expression.expression)
            operand = symbolic_ode_expr(expression.operand, assgn_info, dep_vars, bv_to_dep)
            return "$op$operand"

        elseif expression isa CallNode
            func_name = get_value(expression.function_node)
            ns = expression.function_node.namespace
            prefix = ""
            if !(ns isa Nothing)
                prefix = ns.position.value * "::"
            end
            arg_strs = [symbolic_ode_expr(a, assgn_info, dep_vars, bv_to_dep)
                        for a in expression.args]
            return "$prefix$func_name(" * join(arg_strs, ", ") * ")"

        elseif expression isa IndexAccessNode
            obj = symbolic_ode_expr(expression.object, assgn_info, dep_vars, bv_to_dep)
            idx_strs = [symbolic_ode_expr(i, assgn_info, dep_vars, bv_to_dep)
                        for i in expression.indices]
            return "$obj[" * join(idx_strs, ", ") * "]"

        elseif expression isa IdentifierNode
            name = get_value(expression)
            # If it's a binding variable that maps to a dep var, show the binding name
            # (which is the user-facing name from the .fflow file)
            return name

        elseif expression isa IntegerNode || expression isa FloatNode
            return string(get_value(expression))

        elseif expression isa ArrayLiteralNode
            elts = [symbolic_ode_expr(e, assgn_info, dep_vars, bv_to_dep)
                    for e in expression.elements]
            return "[" * join(elts, ", ") * "]"

        else
            # Fallback: try get_value
            try
                return string(get_value(expression))
            catch
                return "<expr>"
            end
        end
    end

    function traverse_ode_expr(expression, ir_builder, var_loc_attr,
            assgn_info, dep_vars, propensity_table, context)
        """
        Generates ir for an ODE expression
        by traversing the expression tree.
        """
        
        if !(expression isa BinaryOpNode)

            if (expression isa UnaryOpNode)
                operation = expression.expression
                traverse_ode_expr(expression.operand, ir_builder, var_loc_attr,
                  assgn_info, dep_vars, propensity_table, context)
                emit(ir_builder, " $(get_value(operation)) ")
            elseif (expression isa GroupNode)
                emit(ir_builder, "(")
                traverse_ode_expr(expression.expression, ir_builder, var_loc_attr,
                  assgn_info, dep_vars, propensity_table, context)
                emit(ir_builder, ")")
            elseif (expression isa UnaryOpNode)
                operation = expression.expression
                emit(ir_builder, " $(get_value(operation)) ")
                traverse_ode_expr(expression.operand, ir_builder, var_loc_attr,
                  assgn_info, dep_vars, propensity_table, context)
            elseif (expression isa CallNode)
                # Render a function call inside an ODE RHS. Arguments MUST be
                # rendered ODE-aware (via traverse_ode_expr) so that references
                # to the solved dependent variables read the LIVE integrator
                # state through NV_Ith_S(y, varmap...) instead of a stale tensor
                # snapshot. Routing args through ir_builtin_func emitted
                # `pos_pos[0].template item<double>()` — the value at solve
                # start — which is why e.g. boundary_check never saw the tip
                # move and could not stop growth at the wall.
                func_node = expression.function_node
                func_name = get_value(func_node)
                namespace = func_node.namespace
                if !(namespace isa Nothing)
                    namespace = namespace.position.value
                end

                arg_strs = String[]
                for a in expression.args
                    a_ir = IRBuilder([])
                    traverse_ode_expr(a, a_ir, var_loc_attr,
                        assgn_info, dep_vars, propensity_table, context)
                    push!(arg_strs, build_sameline(a_ir))
                end
                arg_str = join(arg_strs, ", ")

                # Map the function name to its C++ form, mirroring the tail of
                # ir_builtin_func but reusing the ODE-aware arg_str above.
                call_ir = _builtin_to_cpp(func_name, arg_str)
                if call_ir === nothing
                    if func_name == "indicator"
                        call_ir = ir_indicator_func(arg_str)
                    elseif func_name == "normal_distr"
                        call_ir = "DGGML::normal_distr()"
                    elseif func_name == "uniform_distr"
                        call_ir = "DGGML::uniform_distr()"
                    elseif func_name == "inverse"
                        call_ir = "DGGML::inverse()"
                    elseif func_name == "gamma_distr"
                        call_ir = "std::gamma_distribution<double>($arg_str)"
                    else
                        # User-defined function (e.g. HELP::boundary_check) or a
                        # DGGML built-in referenced by namespace.
                        call_ir = add_namespace_identifier(func_name, namespace, arg_str)
                    end
                end
                emit(ir_builder, call_ir)

            elseif (expression isa IndexAccessNode)
                # Check if this is an indexed access on a tensor dep var (e.g., position[1])
                # If so, read the current state from SUNDIALS via NV_Ith_S
                if expression.object isa IdentifierNode
                    obj_name = get_value(expression.object)
                    if haskey(var_loc_attr, obj_name)
                        ptr_name = var_loc_attr[obj_name]
                        # Check if this is a tensor dep var (ptr_name[idx] appears in dep_vars)
                        test_ref = "$(ptr_name)[0]"
                        if test_ref in dep_vars
                            # Tensor dep var: emit NV_Ith_S(y, varmap.at(&ptr[idx]))
                            if length(expression.indices) == 1
                                idx_ir = IRBuilder([])
                                traverse_ode_expr(expression.indices[1], idx_ir, var_loc_attr,
                                    assgn_info, dep_vars, propensity_table, context)
                                idx_str = build_sameline(idx_ir)
                                ref = "$(ptr_name)[$idx_str]"
                                emit(ir_builder, "NV_Ith_S(y, varmap.at(&$ref))")
                            else
                                throw("Multi-index tensor dep var access not supported in ODE expressions.")
                            end
                            return
                        end
                    end

                    # Check if this is a non-dep LHS tensor param (e.g., p_unit[0])
                    # These are declared via ir_where_left_clause! into propensity_table
                    # but are NOT SUNDIALS dep vars, so we just read them directly.
                    if haskey(assgn_info, obj_name) && assgn_info[obj_name].cpp_var.type == "torch::Tensor"
                        # Ensure the tensor variable is declared in the ODE lambda
                        lhs_table = propensity_table["var_local_table"]["rule_lhs"]
                        if haskey(lhs_table, obj_name) && !(obj_name in lhs_table["declared"])
                            emit(context, lhs_table[obj_name])
                            push!(lhs_table["declared"], obj_name)
                        end
                        # Emit the variable name; the indexing ([N].template item<double>())
                        # will be appended by the generic index-access code below.
                        emit(ir_builder, " $obj_name ")
                        # Fall through to generic index-access handling below
                        # (skip the recursive traverse on the object)
                        if has_slice(expression.indices)
                            idx_parts = String[]
                            for idx in expression.indices
                                if idx isa SliceNode
                                    push!(idx_parts, ir_slice_expr(idx, ir_builder, propensity_table, context))
                                else
                                    idx_ir = IRBuilder([])
                                    traverse_ode_expr(idx, idx_ir, var_loc_attr,
                                        assgn_info, dep_vars, propensity_table, context)
                                    push!(idx_parts, build_sameline(idx_ir))
                                end
                            end
                            emit(ir_builder, ".index({" * join(idx_parts, ", ") * "})")
                        elseif length(expression.indices) == 1
                            idx_ir = IRBuilder([])
                            traverse_ode_expr(expression.indices[1], idx_ir, var_loc_attr,
                                assgn_info, dep_vars, propensity_table, context)
                            emit(ir_builder, "[" * build_sameline(idx_ir) * "].template item<double>()")
                        else
                            idx_parts = String[]
                            for idx in expression.indices
                                idx_ir = IRBuilder([])
                                traverse_ode_expr(idx, idx_ir, var_loc_attr,
                                    assgn_info, dep_vars, propensity_table, context)
                                push!(idx_parts, build_sameline(idx_ir))
                            end
                            emit(ir_builder, ".index({" * join(idx_parts, ", ") * "}).template item<double>()")
                        end
                        return
                    end
                end

                # Non-tensor IndexAccessNode: original behavior, with slice support
                traverse_ode_expr(expression.object, ir_builder, var_loc_attr,
                    assgn_info, dep_vars, propensity_table, context)
                if has_slice(expression.indices)
                    # Slice present: use .index({...}) with torch::indexing types
                    idx_parts = String[]
                    for idx in expression.indices
                        if idx isa SliceNode
                            push!(idx_parts, ir_slice_expr(idx, ir_builder, propensity_table, context))
                        else
                            idx_ir = IRBuilder([])
                            traverse_ode_expr(idx, idx_ir, var_loc_attr,
                                assgn_info, dep_vars, propensity_table, context)
                            push!(idx_parts, build_sameline(idx_ir))
                        end
                    end
                    emit(ir_builder, ".index({" * join(idx_parts, ", ") * "})")
                elseif length(expression.indices) == 1
                    idx_ir = IRBuilder([])
                    traverse_ode_expr(expression.indices[1], idx_ir, var_loc_attr,
                        assgn_info, dep_vars, propensity_table, context)
                    emit(ir_builder, "[" * build_sameline(idx_ir) * "].template item<double>()")
                else
                    idx_parts = String[]
                    for idx in expression.indices
                        idx_ir = IRBuilder([])
                        traverse_ode_expr(idx, idx_ir, var_loc_attr,
                            assgn_info, dep_vars, propensity_table, context)
                        push!(idx_parts, build_sameline(idx_ir))
                    end
                    emit(ir_builder, ".index({" * join(idx_parts, ", ") * "}).template item<double>()")
                end

            elseif expression isa ArrayLiteralNode
                emit(ir_builder, ir_array_literal(expression, propensity_table, context))

            elseif expression isa IdentifierNode
                # Bare identifier (e.g. a function argument in an ODE RHS).
                # Resolve it the same way binary-operand identifiers are: a
                # solved scalar reads the live SUNDIALS state, a rule-bound
                # attribute reads from the match, a parameter reads from
                # settings.*, and anything else (e.g. a solve-local `auto`
                # variable) is emitted by name.
                id_name = get_value(expression)
                if id_name in dep_vars
                    emit(ir_builder, "NV_Ith_S(y, varmap.at(&ix_$id_name))")
                elseif haskey(assgn_info, id_name)
                    cur_info = assgn_info[id_name]
                    if cur_info.cpp_var.index <= 3 && cur_info.cpp_var.type != "torch::Tensor"
                        emit(ir_builder, create_pos_cpp_var(cur_info, "lhs", "m1"))
                    else
                        emit(ir_builder, create_cpp_var(cur_info, "lhs", "m1"))
                    end
                elseif id_name in collect(keys(propensity_table["parameter_table"]))
                    emit(ir_builder, "settings.$id_name")
                else
                    emit(ir_builder, " $id_name ")
                end

            else
                emit(ir_builder, " $(get_value(expression)) ")
            end
            return

        end

        if expression.lhs isa BinaryOpNode
            # If the left hand side is a binary operation
            traverse_ode_expr(expression.lhs, ir_builder, var_loc_attr,
                  assgn_info, dep_vars, propensity_table, context)

        elseif expression.lhs isa IdentifierNode

            if get_value(expression.lhs) in dep_vars
                ix_ir = "ix_"*get_value(expression.lhs)
                ir = "NV_Ith_S(y, varmap.at(&$ix_ir))"
                emit(ir_builder, ir)
            elseif get_value(expression.lhs) in collect(keys(assgn_info))

                cur_info = assgn_info[get_value(expression.lhs)]

                println("assgn_info for $(get_value(expression.lhs)): ", assgn_info[get_value(expression.lhs)])

                pos = cur_info.cpp_var.index
                # pos = cur_info[4][3]
                # node_pos = cur_info[1]

                # NOTE: This is where we need to differentiate between position attributes and regular attributes.
                # If it's a position attribute, we grab it from the position vector.
                # If it's a regular attribute, we grab it from the data variant.
                # Tensor attributes are never spatial position attributes.
                if pos <= 3 && cur_info.cpp_var.type != "torch::Tensor"
                    ir_node_pos = create_pos_cpp_var(cur_info, "lhs", "m1")
                    # ir_node_pos = "\t\tlhs[m1[$node_pos]].position[$(pos-1)]"
                    emit(ir_builder, ir_node_pos)
                else

                    # Non-positional, non-tensor scalar attribute — grab from data variant
                    ir = create_cpp_var(cur_info, "lhs", "m1")

                    emit(ir_builder, ir)
                    println("Variable $(get_value(expression.lhs)) is a regular attribute.")
                end

            elseif get_value(expression.lhs) in collect(keys(propensity_table["parameter_table"]))
                # emit(ir_builder, get_value(expression.lhs))
                emit(ir_builder, "settings."*get_value(expression.lhs))

            else
                throw("Variable $(get_value(expression.lhs)) not found in dependent variables or assignment info.")
            end
            
        elseif expression.lhs isa UnaryOpNode
            operation = expression.lhs.expression
            emit(ir_builder, " $(operation.position.value) ")
            traverse_ode_expr(expression.lhs.operand, ir_builder, var_loc_attr,
                  assgn_info, dep_vars, propensity_table, context)
        elseif expression.lhs isa GroupNode
            emit(ir_builder, "(")
            traverse_ode_expr(expression.lhs.expression, ir_builder, var_loc_attr,
                  assgn_info, dep_vars, propensity_table, context)
            emit(ir_builder, ")")
            # TODO implement call node
        elseif expression.lhs isa CallNode

            # Delegate to the non-binary dispatch so the call's arguments are
            # rendered ODE-aware: references to solved (dependent) variables
            # read the live SUNDIALS state via NV_Ith_S(y, varmap...) instead
            # of a stale tensor snapshot.
            traverse_ode_expr(expression.lhs, ir_builder, var_loc_attr,
                assgn_info, dep_vars, propensity_table, context)

        elseif expression.lhs isa IndexAccessNode
            traverse_ode_expr(expression.lhs, ir_builder, var_loc_attr,
                  assgn_info, dep_vars, propensity_table, context)

        elseif expression.lhs isa ArrayLiteralNode
            emit(ir_builder, ir_array_literal(expression.lhs, propensity_table, context))

        else
            # If it is a single value, just print it
            emit(ir_builder, get_value(expression.lhs))
        end

        op = expression.expression.position.value
        emit(ir_builder, " $op ")

        if expression.rhs isa BinaryOpNode
            # If the left hand side is a binary operation
            traverse_ode_expr(expression.rhs, ir_builder, var_loc_attr,
                  assgn_info, dep_vars, propensity_table, context)
        elseif expression.rhs isa IdentifierNode
            if get_value(expression.rhs) in dep_vars
                ix_ir = "ix_"*get_value(expression.rhs)
                ir = "NV_Ith_S(y, varmap.at(&$ix_ir))"
                emit(ir_builder, ir)

            elseif get_value(expression.rhs) in collect(keys(assgn_info))

                # FIXME: Depending on if its a position or regular attribute
                # ir_attr = assgn_info[get_value(expression.rhs)][4][2]
                # ir_type = assgn_info[get_value(expression.rhs)][3]
                # pos = assgn_info[get_value(expression.rhs)][4][3]
                # node_pos = assgn_info[get_value(expression.rhs)][1]

                rule_param = assgn_info[get_value(expression.rhs)]
                pos = rule_param.cpp_var.index

                if pos <= 3 && rule_param.cpp_var.type != "torch::Tensor"
                    ir_node_pos = create_pos_cpp_var(assgn_info[get_value(expression.rhs)], "lhs", "m1")
                    # ir_node_pos = "\t\tlhs[m1[$node_pos]].position[$(pos-1)]"
                    emit(ir_builder, ir_node_pos)
                else
                    ir = create_cpp_var(rule_param, "lhs", "m1")
                    # namespace = "Microtubule"
                    # ir = "std::get<$namespace::$ir_type>(lhs[m1[ $(assgn_info[get_value(expression.rhs)][1]) ]].data).$ir_attr"
                    emit(ir_builder, ir)
                end

            elseif get_value(expression.rhs) in collect(keys(propensity_table["parameter_table"]))
                # emit(ir_builder, get_value(expression.rhs))
                emit(ir_builder, "settings."*get_value(expression.rhs))

            else
                throw("Variable $(get_value(expression.rhs)) not found in dependent variables or assignment info.")
            end

        else

            # What if it is a unary op node
            if expression.rhs isa UnaryOpNode
                operation = expression.rhs.expression
                emit(ir_builder, " $(operation.position.value) ")
                traverse_ode_expr(expression.rhs.operand, ir_builder,
                  var_loc_attr, assgn_info, dep_vars, propensity_table, context)
                # return
            elseif expression.rhs isa GroupNode
                emit(ir_builder, "(")
                traverse_ode_expr(expression.rhs.expression, ir_builder, var_loc_attr,
                  assgn_info, dep_vars, propensity_table, context)
                emit(ir_builder, ")")
            elseif expression.rhs isa IndexAccessNode
                traverse_ode_expr(expression.rhs, ir_builder, var_loc_attr,
                  assgn_info, dep_vars, propensity_table, context)
            elseif expression.rhs isa ArrayLiteralNode
                emit(ir_builder, ir_array_literal(expression.rhs, propensity_table, context))
            elseif expression.rhs isa CallNode

                # Delegate to the non-binary dispatch so the call's arguments
                # are rendered ODE-aware: references to solved (dependent)
                # variables read the live SUNDIALS state via
                # NV_Ith_S(y, varmap...) instead of a stale tensor snapshot.
                traverse_ode_expr(expression.rhs, ir_builder, var_loc_attr,
                    assgn_info, dep_vars, propensity_table, context)
            else
                # If it is a single value, just print it
                emit(ir_builder, get_value(expression.rhs))
            end

        end
    end


    function ir_solve_variable_binding(ir_builder, solve_clause,
            lhs_assgn_to_node, type_namespace, readonly_vars=Set{String}())
        """
        Does variable binding for solve clause.
        Registers scalar variables directly with varset.insert(&var).
        For tensor variables (torch::Tensor), gets a double* to the underlying data
        and registers each element individually: varset.insert(&ptr[i]).
        Also registers read-only variables (readonly_vars) that appear in ODE RHS
        expressions but are not solving variables — these must be in varset so the
        ODE lambda can read them from the SUNDIALS y vector via varmap.
        """

        var_bind_ir = "[](auto &lhs, auto &m1, auto &varset) {"
        emit(ir_builder, var_bind_ir)
        var_attr_loc = Dict()

        # Takes the binding variable and returns the
        # dependency. For tensors, also tracks per-element info.

        bv_to_dep = Dict()
        # For tensor bindings, store the tensor size and per-element pointer name
        # key: dep_var_name => (tensor_size, pointer_name)
        tensor_binding_info = Dict()

        # Handle variable binding
        # Bind Variable (bv)
        # Track which base tensor variables have already been set up
        tensor_setup_done = Set()

        map((bv_node) -> begin

                # TODO: Generalize for multiple var odes
                # Only grabs the first variable.
                # Handle indexed ODE variables: D(im_pos[0], t) -> base name "im_pos"
                ode_var = bv_node.value[1]
                if ode_var isa IndexAccessNode
                    dep_vars = get_value(ode_var.object)
                else
                    dep_vars = get_value(ode_var)
                end

                # Handle indexed binding name: dpos[0] -> base name "dpos"
                bv_name_node = bv_node.name
                if bv_name_node isa IndexAccessNode
                    bv_name = get_value(bv_name_node.object)
                else
                    bv_name = get_value(bv_name_node)
                end

                bv_to_dep[bv_name] = dep_vars

                bv_pos = lhs_assgn_to_node[dep_vars].index
                attr_pos = lhs_assgn_to_node[dep_vars].cpp_var.index
                bv_type_str = lhs_assgn_to_node[dep_vars].cpp_var.type
                tensor_size = lhs_assgn_to_node[dep_vars].cpp_var.tensor_size

                if bv_type_str == "torch::Tensor" && tensor_size > 0
                    # === TENSOR VARIABLE BINDING ===
                    # For indexed bindings (dpos[0] := D(im_pos[0], t)) we still
                    # need to set up the full tensor pointer once, since SUNDIALS
                    # needs all elements registered.
                    if !(dep_vars in tensor_setup_done)
                        push!(tensor_setup_done, dep_vars)

                        bv_attr = lhs_assgn_to_node[dep_vars].cpp_var.name
                        bv_node_type = lhs_assgn_to_node[dep_vars].type
                        bv_attr_pos = lhs_assgn_to_node[dep_vars].cpp_var.index

                        if bv_attr_pos == 1
                            # === POSITION SOLVING TENSOR ===
                            # Bind ONLY the deep-copied per-node `.position[]` array as the
                            # solver variable, do NOT bind the torch::Tensor Position data
                            # pointer. The Position tensor is a shallow (ref-counted) copy
                            # that is SHARED across duplicated graph nodes, so two distinct
                            # MT tips can point at one storage buffer. Because varset is a
                            # std::set<double*>, that shared pointer is de-duplicated and
                            # both nodes collapse onto a SINGLE solver variable — their
                            # positions then snap together and lock ("vertices combine").
                            # `.position[]` is deep-copied per node (SpatialNode3D) and is
                            # never aliased, so each node gets its own solver variable.
                            pos_count = min(3, tensor_size)
                            for i in 0:(pos_count - 1)
                                emit(ir_builder, "varset.insert(&lhs[m1[$bv_pos]].position[$i]);")
                            end
                            pos_accessor = "lhs[m1[$bv_pos]].position"
                            # Store the .position accessor as the "pointer"; downstream RHS
                            # emits &pos_accessor[i] = &lhs[m1[i]].position[i].
                            tensor_binding_info[dep_vars] = (pos_count, pos_accessor, bv_pos)
                            var_attr_loc[dep_vars] = pos_accessor
                        else
                            # === NON-POSITION SOLVING TENSOR (e.g. neural-ODE weights) ===
                            # Fetch the tensor reference
                            tensor_ref = "std::get<$type_namespace::$bv_node_type>(lhs[m1[$bv_pos]].data).$bv_attr"
                            ptr_name = "tensor_ptr_$(bv_pos)_$(bv_attr_pos)"

                            # Get raw double* pointer from tensor
                            emit(ir_builder, "auto &tensor_ref_$(bv_pos)_$(bv_attr_pos) = $tensor_ref;")
                            emit(ir_builder, "double* $ptr_name = tensor_ref_$(bv_pos)_$(bv_attr_pos).template data_ptr<double>();")

                            # Register each element with varset
                            for i in 0:(tensor_size - 1)
                                emit(ir_builder, "varset.insert(&$(ptr_name)[$i]);")
                            end

                            # Store tensor info for use in ODE lambda
                            tensor_binding_info[dep_vars] = (tensor_size, ptr_name, bv_pos)

                            # var_attr_loc maps dep_var to a per-element accessor pattern
                            var_attr_loc[dep_vars] = ptr_name
                        end
                    end

                elseif attr_pos > 3
                    # === NON-TENSOR, NON-POSITION SCALAR ATTRIBUTE ===
                    bv_attr = lhs_assgn_to_node[dep_vars].cpp_var.name
                    bv_node_type = lhs_assgn_to_node[dep_vars].type
                    bv_attr_pos = lhs_assgn_to_node[dep_vars].cpp_var.index

                    ref_name = "node_$(bv_pos)_$bv_attr_pos"
                    ref_fetch = "std::get<$type_namespace::$bv_node_type>(lhs[m1[$bv_pos]].data).$bv_attr"

                    # Fetching attr
                    ir_fetch = "auto &$ref_name = $ref_fetch;"
                    emit(ir_builder, ir_fetch)

                    ir = "varset.insert(&$ref_name);"
                    ir_ix_attr_loc = "$ref_fetch"

                    var_attr_loc[dep_vars] = ir_ix_attr_loc
                    emit(ir_builder, ir)
                else
                    # === POSITION SCALAR ATTRIBUTE ===
                    ir = "varset.insert(&lhs[m1[$bv_pos]].position[$(attr_pos-1)]);"

                    attr_ir = "lhs[m1[$bv_pos]].position[$(attr_pos-1)]"
                    ir_ix_attr_loc = "$attr_ir"

                    println("ir_ix_attr_loc for $bv_name: ", ir_ix_attr_loc)

                    var_attr_loc[dep_vars] = ir_ix_attr_loc
                    emit(ir_builder, ir)
                end

            end,
            solve_clause.variables
           )

        # === REGISTER READ-ONLY ODE VARIABLES ===
        # These are LHS-matched variables referenced in ODE RHS expressions
        # but NOT solving variables. They must be in varset so the ODE lambda
        # reads their current values from the SUNDIALS y vector (via varmap)
        # instead of stale graph memory during RK intermediate stages.
        readonly_registered = Set{String}()
        for ro_var in readonly_vars
            if haskey(lhs_assgn_to_node, ro_var) && !haskey(var_attr_loc, ro_var)
                ro_info = lhs_assgn_to_node[ro_var]
                ro_pos = ro_info.index
                ro_attr_pos = ro_info.cpp_var.index
                ro_type_str = ro_info.cpp_var.type
                ro_tensor_size = ro_info.cpp_var.tensor_size

                if ro_type_str == "torch::Tensor" && ro_tensor_size > 0 && ro_attr_pos != 1
                    # === NON-POSITION READ-ONLY TENSOR ===
                    # Do NOT bind this tensor into the SUNDIALS solver. Binding it
                    # would capture a raw double* into a shallow, ref-counted,
                    # movable torch::Tensor at reinit and then write back through
                    # that pointer in copy_back() every solve a use-after-free
                    # hazard if a graph rewrite reallocates/frees the tensor storage
                    # between reinit and copy_back, and the same std::set<double*>
                    # dedup aliasing that collapses shared-storage nodes. A
                    # read-only tensor is a CONSTANT during integration (its ydot is
                    # never assigned), so it is read directly as a live value via
                    # `pos_dir[i].template item<double>()` in the RHS instead. This
                    # matches CMA's ode_mt_growth, which reads unit_vec as a plain
                    # constant (num_eq=3). Leaving it out of var_attr_loc lets
                    # traverse_ode_expr fall through to its tensor-param read branch,
                    # which declares the tensor local on demand from the where-clause
                    # binding table. Nothing to register here.

                elseif ro_type_str == "torch::Tensor" && ro_tensor_size > 0
                    # === POSITION READ-ONLY TENSOR (attr_pos == 1) ===
                    # A reference position (e.g. the Intermediate's position used
                    # inside a distance() gate) that is NOT solved but must be read
                    # at its live value during integration. Bind the deep-copied
                    # per-node `.position[]` array — NEVER the shared, ref-counted
                    # torch Position tensor pointer — for the same aliasing-safety
                    # reasons as the solving-position branch, and crucially so the
                    # binding is addressable in BOTH the varset lambda AND the solve
                    # lambda. A `tensor_ptr_*` declared only in the varset lambda is
                    # out of scope in the solve lambda (`tensor_ptr_2_1 was not
                    # declared in this scope`). Registered with a zero derivative,
                    # it stays constant at its solve-start value.
                    pos_count = min(3, ro_tensor_size)
                    for i in 0:(pos_count - 1)
                        emit(ir_builder, "varset.insert(&lhs[m1[$ro_pos]].position[$i]);")
                    end
                    pos_accessor = "lhs[m1[$ro_pos]].position"
                    tensor_binding_info[ro_var] = (pos_count, pos_accessor, ro_pos)
                    var_attr_loc[ro_var] = pos_accessor
                    push!(readonly_registered, ro_var)

                elseif ro_attr_pos > 3
                    # Non-positional scalar read-only variable
                    ro_attr = ro_info.cpp_var.name
                    ro_node_type = ro_info.type
                    ref_name = "node_$(ro_pos)_$ro_attr_pos"
                    ref_fetch = "std::get<$type_namespace::$ro_node_type>(lhs[m1[$ro_pos]].data).$ro_attr"
                    emit(ir_builder, "auto &$ref_name = $ref_fetch;")
                    emit(ir_builder, "varset.insert(&$ref_name);")
                    var_attr_loc[ro_var] = ref_fetch
                    push!(readonly_registered, ro_var)

                else
                    # Position scalar read-only variable
                    emit(ir_builder, "varset.insert(&lhs[m1[$ro_pos]].position[$(ro_attr_pos-1)]);")
                    var_attr_loc[ro_var] = "lhs[m1[$ro_pos]].position[$(ro_attr_pos-1)]"
                    push!(readonly_registered, ro_var)
                end
            end
        end

        emit(ir_builder, "},")


        var_attr_loc, bv_to_dep, tensor_binding_info, readonly_registered
    end


    function ir_solve_clause!(solve_clause, lhs_assgn_to_node, rhs_assgn_to_node,
        ir_builder, type_namespace, propensity_table, symbol_tables, lhs_param_to_node)
        """
        Solve Clause Node.
        Handles both scalar and tensor (FixedList) dependent variables.
        For tensors, uses double* pointers into torch::Tensor data for SUNDIALS integration.
        Supports:
          - dx[i] : ODE = expr   (indexed, per-element ODE for tensor binding)
          - dx : ODE = expr      (unindexed, broadcast to all tensor elements)
          - dx : ODE = expr      (scalar, original behavior)
        """

        # Collect all variables referenced in ODE RHS expressions
        # Walk ALL nodes in the solve clause body (ODENode values AND
        # DefinitionNode values) so that tensor vertex attributes referenced
        # in temporary definitions are also registered as read-only vars
        # and read from the SUNDIALS y vector instead of stale graph memory.
        all_ode_referenced = Set{String}()
        for assgn_node in solve_clause.clause
            if assgn_node isa ODENode
                union!(all_ode_referenced,
                       collect_ode_referenced_vars(assgn_node.value, lhs_assgn_to_node))
            elseif assgn_node isa DefinitionNode
                union!(all_ode_referenced,
                       collect_ode_referenced_vars(assgn_node.value, lhs_assgn_to_node))
            end
        end

        # Compute the set of solving variable names (from binding vars)
        solving_var_names = Set{String}()
        for bv_node in solve_clause.variables
            ode_var = bv_node.value[1]
            dep_name = ode_var isa IndexAccessNode ? get_value(ode_var.object) : get_value(ode_var)
            push!(solving_var_names, dep_name)
        end

        # Read-only vars = referenced in ODE RHS but not solving variables
        readonly_vars = setdiff(all_ode_referenced, solving_var_names)
        println("  ODE referenced vars: ", all_ode_referenced)
        println("  Solving vars: ", solving_var_names)
        println("  Read-only ODE vars: ", readonly_vars)

        var_attr_loc, bv_to_dep, tensor_binding_info, readonly_registered = ir_solve_variable_binding(
                ir_builder, solve_clause, lhs_assgn_to_node,
                type_namespace, readonly_vars
            )


        println("solving var names: ", solving_var_names)
        println("readonly vars: ", readonly_vars)
        println("all ode referenced vars: ", all_ode_referenced)

        # Defining ODE and temporary definitions
        emit(ir_builder, "[&](auto &lhs, auto &m1, auto y, auto ydot, auto &varmap) {")

        ir_where_left_clause!(lhs_param_to_node, type_namespace,
                              symbol_tables, propensity_table, ir_builder)

        dep_vars = []
        # Collect symbolic equations for debug printing
        symbolic_equations = String[]
        # Fetch all the dependent variables first (solving + read-only).
        # For tensor dep vars, create a double* pointer and per-element refs.
        # For scalar dep vars, create a single auto& ref as before.
        # dep_vars stores RAW variable names so traverse_ode_expr can match
        # get_value(expression) against them and emit NV_Ith_S(y, varmap[...]).
        for dep_var in collect(keys(var_attr_loc))
            if haskey(tensor_binding_info, dep_var)
                # Tensor: fetch the double* pointer in the ODE lambda
                tsize, ptr_name, bv_pos = tensor_binding_info[dep_var]
                if lhs_assgn_to_node[dep_var].cpp_var.index == 1
                    # Position solving tensor is bound directly to the deep-copied
                    # `.position[]` array (ptr_name is the `.position` accessor, not a
                    # torch pointer). Do NOT emit a `double* ...` declaration; just
                    # register the per-element refs so any RHS read resolves to
                    # NV_Ith_S(y, varmap.at(&lhs[m1[i]].position[k])).
                    for i in 0:(tsize - 1)
                        push!(dep_vars, "$(ptr_name)[$i]")
                    end
                else
                    bv_attr = lhs_assgn_to_node[dep_var].cpp_var.name
                    bv_node_type = lhs_assgn_to_node[dep_var].type
                    tensor_ref = "std::get<$type_namespace::$bv_node_type>(lhs[m1[$bv_pos]].data).$bv_attr"
                    emit(ir_builder, "double* $ptr_name = $tensor_ref.template data_ptr<double>();")
                    # Register per-element refs so traverse_ode_expr can read them via NV_Ith_S
                    for i in 0:(tsize - 1)
                        push!(dep_vars, "$(ptr_name)[$i]")
                    end
                end
            else

                # Scalar: create auto& reference for varmap lookup
                var_loc_ir = var_attr_loc[dep_var]

                println("Binding scalar dep var '$dep_var' with varmap ref: ", var_loc_ir)

                ir = "auto &ix_$dep_var = $var_loc_ir;"
                # Store the RAW name so traverse_ode_expr's
                # `get_value(node) in dep_vars` check matches
                push!(dep_vars, dep_var)
                emit(ir_builder, ir)
            end
        end

        # Process each node in the solve clause body
        map( assgn_node -> begin

            if assgn_node isa DefinitionNode

                def_name = get_value(assgn_node.name)
                def_type = assgn_node.type

                ir_value = IRBuilder([])

                ir_definition!(assgn_node.value, ir_value,
                               propensity_table, ir_builder,
                               true, false)

                ir = "auto $def_name = " * build(ir_value) * ";"

                println("ir: ", ir)

                emit(ir_builder, ir)

                # Register this solve-local variable so later expressions
                # (e.g. built-in / user function calls like boundary_check(...))
                # recognize it as already declared and don't re-emit or error.
                if !(def_name in propensity_table["var_local_table"]["rule_rhs"]["declared"])
                    push!(propensity_table["var_local_table"]["rule_rhs"]["declared"], def_name)
                end

            elseif assgn_node isa ODENode

                ode_name_node = assgn_node.name
                ode_value = assgn_node.value

                if ode_name_node isa IndexAccessNode
                    # === INDEXED TENSOR ODE: dx[i] : ODE = expr ===
                    bv_name = get_value(ode_name_node.object)
                    idx_val = get_value(ode_name_node.indices[1])  # integer index

                    dep_var_name = bv_to_dep[bv_name]

                    if !haskey(tensor_binding_info, dep_var_name)
                        throw("Indexed ODE $bv_name[$idx_val] but $dep_var_name is not a tensor binding.")
                    end

                    tsize, ptr_name, bv_pos = tensor_binding_info[dep_var_name]
                    ref = "$(ptr_name)[$idx_val]"

                    ir = "NV_Ith_S(ydot, varmap[&$ref]) += "

                    expr_ir = IRBuilder([])
                    traverse_ode_expr(ode_value, expr_ir, var_attr_loc,
                                      lhs_assgn_to_node, dep_vars, propensity_table, ir_builder)

                    build_expr_ir = join(expr_ir.instructions)
                    ir = ir * build_expr_ir * ";"
                    emit(ir_builder, ir)

                    # Collect symbolic form
                    sym_rhs = symbolic_ode_expr(ode_value, lhs_assgn_to_node, dep_vars, bv_to_dep)
                    push!(symbolic_equations, "d($bv_name[$idx_val])/dt += $sym_rhs")

                elseif ode_name_node isa IdentifierNode
                    bv_name = get_value(ode_name_node)

                    if haskey(bv_to_dep, bv_name) && haskey(tensor_binding_info, bv_to_dep[bv_name])
                        # === UNINDEXED TENSOR ODE (Neural ODE): dw : ODE = f(weights) ===
                        # Pattern:
                        #   1. Reconstruct the tensor from SUNDIALS y (so it reflects ODE state)
                        #   2. Evaluate the RHS once into a temp torch::Tensor
                        #   3. Scatter elements into ydot element-by-element
                        dep_var_name = bv_to_dep[bv_name]
                        tsize, ptr_name, bv_pos = tensor_binding_info[dep_var_name]
                        bv_attr     = lhs_assgn_to_node[dep_var_name].cpp_var.name
                        bv_node_type = lhs_assgn_to_node[dep_var_name].type

                        # === from_blob zero-copy pattern ===
                        # SUNDIALS stores this tensor's elements contiguously starting at
                        # varmap[&ptr[0]].  We wrap y/ydot memory directly into LibTorch
                        # tensors — no copy in either direction.
                        emit(ir_builder, "{  // Neural ODE: $bv_name")
                        emit(ir_builder, "const int _base_$bv_name = varmap.at(&$(ptr_name)[0]);")
                        emit(ir_builder, "// Zero-copy views into SUNDIALS y / ydot memory")
                        emit(ir_builder, "torch::Tensor _y_$bv_name = torch::from_blob(")
                        emit(ir_builder, "    N_VGetArrayPointer(y) + _base_$bv_name, {$tsize}, torch::kFloat64);")
                        emit(ir_builder, "torch::Tensor _ydot_$bv_name = torch::from_blob(")
                        emit(ir_builder, "    N_VGetArrayPointer(ydot) + _base_$bv_name, {$tsize}, torch::kFloat64);")

                        # Rebind the user-visible attribute name to the y view
                        emit(ir_builder, "auto& $bv_attr = _y_$bv_name;")

                        # Evaluate RHS once; flush any preamble declarations first
                        rhs_preamble = IRBuilder([])
                        rhs_expr    = IRBuilder([])
                        traverse_ode_expr(ode_value, rhs_expr, var_attr_loc,
                                          lhs_assgn_to_node, dep_vars, propensity_table, rhs_preamble)
                        for preamble_line in rhs_preamble.instructions
                            emit(ir_builder, preamble_line)
                        end
                        rhs_str = join(rhs_expr.instructions)
                        emit(ir_builder, "torch::Tensor _rhs_$bv_name = $rhs_str;")

                        # Accumulate into ydot via the zero-copy view (no scatter loop)
                        emit(ir_builder, "_ydot_$(bv_name) += _rhs_$(bv_name);")
                        emit(ir_builder, "}  // end Neural ODE: $bv_name")

                        # Collect symbolic form
                        sym_rhs = symbolic_ode_expr(ode_value, lhs_assgn_to_node, dep_vars, bv_to_dep)
                        push!(symbolic_equations, "d($bv_name[0..$( tsize-1 )])/dt += $sym_rhs")

                    elseif haskey(bv_to_dep, bv_name)
                        # === SCALAR ODE: dx : ODE = expr (original behavior) ===
                        dep_var_name = bv_to_dep[bv_name]

                        # Look up from lhs_assgn_to_node since this is a binding var
                        bv_pos = lhs_assgn_to_node[dep_var_name].index
                        attr_pos = lhs_assgn_to_node[dep_var_name].cpp_var.index

                        ref = nothing
                        if attr_pos > 3
                            ref = "ix_$dep_var_name"
                        else
                            ref = "ix_$dep_var_name"
                        end

                        ir = "NV_Ith_S(ydot, varmap[&$ref]) += "

                        expr_ir = IRBuilder([])
                        traverse_ode_expr(ode_value, expr_ir, var_attr_loc,
                                          lhs_assgn_to_node, dep_vars, propensity_table, ir_builder)

                        build_expr_ir = join(expr_ir.instructions)
                        ir = ir * build_expr_ir * ";"
                        emit(ir_builder, ir)

                        # Collect symbolic form
                        sym_rhs = symbolic_ode_expr(ode_value, lhs_assgn_to_node, dep_vars, bv_to_dep)
                        push!(symbolic_equations, "d($bv_name)/dt += $sym_rhs")

                    else
                        throw("ODE name '$bv_name' not found in binding variables.")
                    end
                else
                    throw("Unknown ODE name node type: $(typeof(ode_name_node))")
                end

            else
                throw("Unknown node type in solve clause: $(typeof(assgn_node))")
            end

            end, solve_clause.clause
        )

        # Emit position sync ODEs for tensor bindings.
        # SUNDIALS tracks .position[i] as independent variables, so their ydot
        # must mirror the corresponding tensor element's ydot to stay in lockstep.
        # Only sync for the Position attribute (first FixedList, attr_pos == 1).
        #
        # IMPORTANT: sync ONLY solving (D-target) variables. Read-only tensor
        # variables are registered in varset for READING their current value, but
        # their node.position[] array is deliberately NOT registered in the ic
        # varset. Emitting a sync line for a read-only var would reference an
        # unregistered pointer; `varmap[...]` (operator[]) then default-inserts it
        # with index 0, and copy_back() later writes y[0] (a different node's
        # integrated value) into that node's position[] — corrupting a read-only
        # node's coordinates (e.g. an Intermediate's z jumping off the z=0 plane).
        for (dep_var, info) in tensor_binding_info
            tsize, ptr_name, bv_pos = info
            attr_pos = lhs_assgn_to_node[dep_var].cpp_var.index
            if attr_pos == 1 && (dep_var in solving_var_names) && startswith(ptr_name, "tensor_ptr_")
                pos_count = min(3, tsize)
                for i in 0:(pos_count - 1)
                    tensor_ref = "$(ptr_name)[$i]"
                    pos_ref = "lhs[m1[$bv_pos]].position[$i]"
                    emit(ir_builder, "NV_Ith_S(ydot, varmap[&$pos_ref]) += NV_Ith_S(ydot, varmap[&$tensor_ref]);")
                end
            end
        end

        #Emit symbolic ODE equations as C++ comments
        if !isempty(symbolic_equations)
            # Print to Julia stdout during codegen
            println("=== Symbolic ODE system ===")
            for eq in symbolic_equations
                println("    $eq")
            end

            emit(ir_builder, "// == Symbolic ODE system ==")
            for eq in symbolic_equations
                emit(ir_builder, "// $eq")
            end
            # Runtime print gated by ODE_DUMP env var, fires once per rule instance
            emit(ir_builder, "{ static bool _sym_dumped = false;")
            emit(ir_builder, "  if (!_sym_dumped && std::getenv(\"ODE_DUMP\")) { _sym_dumped = true;")
            emit(ir_builder, "    std::cout << \"  ── Symbolic ODE ──\" << std::endl;")
            for eq in symbolic_equations
                # Escape any backslashes/quotes for C++ string literal
                escaped = replace(eq, "\\" => "\\\\")
                escaped = replace(escaped, "\"" => "\\\"")
                emit(ir_builder, "    std::cout << \"    $escaped\" << std::endl;")
            end
            emit(ir_builder, "  }")
            emit(ir_builder, "}")
        end

        emit(ir_builder, "}")
    end
