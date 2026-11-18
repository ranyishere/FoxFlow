# Propensity expression IR generation for `with` rules.
# Included into module IRRuleGeneration (see ../ir_rule_generation.jl).

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
            #

            # FIXME: Check to see if it is in the rule_lhs also..
            if arg.position.value in collect(keys(propensity_table)) && propensity == true

                ir = propensity_table[arg.position.value]


                # FIXME: Loads it into the body depending if its a position or an attribute
                if !(arg.position.value in propensity_table["var_local_table"]["propensity"]["declared"])
                    emit(prop_body_ir, ir)
                    push!(propensity_table["var_local_table"]["propensity"]["declared"], arg.position.value)
                end

                emit(ir_builder, arg.position.value)

            # It's in the parameters file
            elseif arg.position.value in collect(keys(propensity_table["parameter_table"]))

                # if propensity == true
                    # emit(prop_body_ir, "settings."*arg.position.value)
                # end

                # ir = propensity_table[arg.position.value]
                ir = "settings." * arg.position.value
                emit(ir_builder, ir)

            # It's in the local variable table for the rhs of the rule
            elseif arg.position.value in collect(
                                     keys(propensity_table["var_local_table"]["rule_rhs"])
                                )

                # Dedup against ALL THREE declared sets (see the rule_lhs branch
                # below): the same RHS-bound tensor may already have been declared
                # into this same body by a HELP:: call that recorded it under
                # propensity["declared"] instead of rule_rhs["declared"].
                if !(arg.position.value in propensity_table["var_local_table"]["rule_rhs"]["declared"]) &&
                   !(arg.position.value in propensity_table["var_local_table"]["rule_lhs"]["declared"]) &&
                   !(arg.position.value in propensity_table["var_local_table"]["propensity"]["declared"])
                    ir = propensity_table["var_local_table"]["rule_rhs"][arg.position.value]
                    emit(prop_body_ir, ir)
                    push!(propensity_table["var_local_table"]["rule_rhs"]["declared"], arg.position.value)
                end

                emit(ir_builder, arg.position.value)

            elseif arg.position.value in propensity_table["var_local_table"]["rule_rhs"]["declared"]

                emit(ir_builder, arg.position.value)

            # It is in the lhs of the rule
            elseif arg.position.value in collect(keys(propensity_table["var_local_table"]["rule_lhs"]))

                # Dedup against ALL THREE declared sets. The same LHS tensor may
                # already have been declared into this same body by another code
                # path that tracks dedup in a different set: a HELP:: call inside a
                # where-block definition goes through ir_builtin_func with
                # propensity=true and records the tensor under
                # propensity["declared"], whereas this arithmetic path only pushed
                # to rule_lhs["declared"]. Checking only rule_lhs["declared"] here
                # re-declared the same torch::Tensor and produced a C++
                # "redeclaration of at::Tensor" error (e.g. crossover's
                # `s1_raw := HELP::intersection_s1(..., pos_dir)` followed by
                # `new_pos_pos[0] = ... pos_dir[0]`).
                if !(arg.position.value in propensity_table["var_local_table"]["rule_lhs"]["declared"]) &&
                   !(arg.position.value in propensity_table["var_local_table"]["rule_rhs"]["declared"]) &&
                   !(arg.position.value in propensity_table["var_local_table"]["propensity"]["declared"])
                    ir = propensity_table["var_local_table"]["rule_lhs"][arg.position.value]
                    emit(prop_body_ir, ir)
                    push!(propensity_table["var_local_table"]["rule_lhs"]["declared"], arg.position.value)
                end

                emit(ir_builder, arg.position.value)

            else
                # TODO: Check if it is inside the settings also (global variable)
                # TODO: Check where clause
                # print("Declared: ", 
                      # propensity_table["var_local_table"]["rule_rhs"]["declared"])
		# println("arg: ", arg)
		throw("Error: Propensity variable $(arg.position.value) not found in propensity table. $(arg)")
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

            func_name = get_value(function_node)

            # Function namespace
            namespace = function_node.namespace

            if !(namespace isa Nothing)
                namespace = namespace.position.value
            end

            # func_args = function_node.args
            func_args = call_args

            # Loop through func args and emit them t
            ir_builtin_func(func_name, func_args, namespace,
                            ir_builder, propensity_table, prop_body_ir, propensity)

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
                ir = ir_distribution!(function_node.operand, ir_builder, prop_body_ir, propensity_table, propensity)
            elseif operation.position.value == "-"
                value = function_node.operand

                # If the operand is a literal (IntegerNode, FloatNode, LiteralNode),
                # just emit the negated value directly without propensity table lookup.
                if value isa LiteralNode
                    ir = "-" * get_value(value)
                elseif value isa IndexAccessNode || value isa CallNode || value isa BinaryOpNode || value isa GroupNode
                    # Complex operand: use ir_prop_expr! to generate its IR,
                    # then prepend the negation sign.
                    inner_ir = IRBuilder([])
                    ir_prop_expr!(value, inner_ir, propensity_table, prop_body_ir, propensity)
                    ir = "-" * build_sameline(inner_ir)
                else
                    # Simple identifier: fetch from propensity table
                    fetched = find_and_fetch_propensity_var(value, propensity_table)
                    if fetched != nothing
                        emit(prop_body_ir, fetched)
                    end

                    arg_val = get_value(function_node.operand)

                    if arg_val in collect(keys(propensity_table["parameter_table"]))
                        arg_val = "settings." * arg_val
                    end

                    ir = "-" * arg_val
                end

            else
                throw("Error: Unary operation $(operation.position.value) not recognized.")
            end

            emit(ir_builder, ir)

        elseif function_node isa IndexAccessNode
            # Tensor index access in propensity expression, with slice support
            ir_prop_expr!(function_node.object, ir_builder, propensity_table, prop_body_ir, propensity)
            if has_slice(function_node.indices)
                idx_parts = String[]
                for idx in function_node.indices
                    if idx isa SliceNode
                        push!(idx_parts, ir_slice_expr(idx, ir_builder, propensity_table, ir_builder))
                    else
                        idx_ir = IRBuilder([])
                        ir_prop_expr!(idx, idx_ir, propensity_table, prop_body_ir, propensity)
                        push!(idx_parts, build_sameline(idx_ir))
                    end
                end
                emit(ir_builder, ".index({" * join(idx_parts, ", ") * "})")
            else
                indices_ir = IRBuilder([])
                for (i, idx) in enumerate(function_node.indices)
                    if i > 1
                        emit(indices_ir, ", ")
                    end
                    ir_prop_expr!(idx, indices_ir, propensity_table, prop_body_ir, propensity)
                end
                if length(function_node.indices) == 1
                    emit(ir_builder, "[" * build_sameline(indices_ir) * "].template item<double>()")
                else
                    emit(ir_builder, ".index({" * build_sameline(indices_ir) * "}).template item<double>()")
                end
            end

        elseif function_node isa ArrayLiteralNode
            # For array literals that may contain complex expressions (index access,
            # function calls, etc.), we pre-evaluate each element and store in temp
            # variables, then construct the tensor from those temps.
            # This avoids issues with brace-initializer lists and template parsing,
            # and properly triggers variable declarations for referenced LHS/RHS vars.
            n_elems = length(function_node.elements)
            has_complex = any(el -> !(el isa LiteralNode), function_node.elements)

            if has_complex
                # Generate unique temp names using global counter
                elem_names = String[]
                for (i, el) in enumerate(function_node.elements)
                    el_ir = IRBuilder([])
                    ir_prop_expr!(el, el_ir, propensity_table, prop_body_ir, propensity)
                    el_str = build_sameline(el_ir)
                    tmp_name = "_arr_tmp_$(_arr_tmp_counter[])"
                    _arr_tmp_counter[] += 1
                    emit(prop_body_ir, "double $tmp_name = $el_str;")
                    push!(elem_names, tmp_name)
                end
                if ARRAY_CTX[]
                    emit(ir_builder, "TFArray<$(n_elems)>{{" * join(elem_names, ", ") * "}}")
                else
                    emit(ir_builder, "torch::tensor({" * join(elem_names, ", ") * "}, torch::kFloat64)")
                end
            else
                # All elements are simple literals — safe to inline directly
                elem_strs = String[]
                for el in function_node.elements
                    el_ir = IRBuilder([])
                    ir_prop_expr!(el, el_ir, propensity_table, prop_body_ir, propensity)
                    push!(elem_strs, build_sameline(el_ir))
                end
                if ARRAY_CTX[]
                    emit(ir_builder, "TFArray<$(n_elems)>{{" * join(elem_strs, ", ") * "}}")
                else
                    emit(ir_builder, "torch::tensor({" * join(elem_strs, ", ") * "}, torch::kFloat64)")
                end
            end

        else
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

        # reset the declared variables for the propensity function
        propensity_table["var_local_table"]["propensity"]["declared"] = []

        ir_prop_expr!(function_node, ir_builder, propensity_table, prop_body_ir, true)
    end
