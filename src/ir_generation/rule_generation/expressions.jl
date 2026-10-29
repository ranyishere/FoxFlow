# Expression IR generation: definitions, built-in calls, distributions.
# Included into module IRRuleGeneration (see ../ir_rule_generation.jl).

    function find_and_fetch_propensity_var(arg, propensity_table)
        """
        Finds the propensity variable in the propensity table and returns its IR.
        """

        val = get_value(arg)

        if val in collect(keys(propensity_table))
            ir = propensity_table[val]
            return ir
        elseif val in keys(propensity_table["var_local_table"]["rule_rhs"])
            ir = propensity_table["var_local_table"]["rule_rhs"][val]
            return ir
        elseif val in propensity_table["var_local_table"]["rule_rhs"]["declared"]
            return nothing
        elseif val in collect(keys(propensity_table["var_local_table"]["rule_lhs"]))
            ir = propensity_table["var_local_table"]["rule_lhs"][arg.position.value]
            return ir
        elseif val in propensity_table["var_local_table"]["rule_lhs"]["declared"]
            return nothing
        else
            throw("Error: Propensity variable $(val) not found in propensity table.")
        end
    end


    function add_namespace_identifier(func_name, namespace, arg_str)
        ir = nothing
        if namespace != nothing
            ir = "$(namespace)::$func_name($arg_str)"
        else
            ir = "$func_name($arg_str)"
        end
        ir
    end

    function ir_definition!(expression, ir_builder,
            propensity_table, context, propensity=false, where_clause=true)
        """
        Handles definition node, which is used for defining variables
        in the where clause of a with statement.
        """

        # Check if function node is a function or expression or an identifier
        func_args = nothing
        func_name = nothing

        # TODO: Implement identity function in c++?
        if expression isa IdentifierNode

            # If it is an identifier, we should
            # fetch the function arguments from the identifier
            arg = expression.token

            # It's in the parameters file
            if arg.position.value in collect(keys(propensity_table["parameter_table"]))

                # ir = propensity_table[arg.position.value]
                ir = "settings." * arg.position.value
                emit(ir_builder, ir)

            # It is in the lhs of the rule
            elseif arg.position.value in collect(keys(propensity_table["var_local_table"]["rule_lhs"]))

                # Dedup against ALL THREE declared sets. The same LHS tensor can
                # be first read from three different code paths that each track
                # dedup in a different set but emit into the same `context`
                # buffer:
                #   * this arithmetic path            -> rule_lhs["declared"]
                #   * ir_builtin_func where_clause=true -> rule_rhs["declared"]
                #   * ir_builtin_func propensity=true   -> propensity["declared"]
                # A HELP:: call inside a where-block definition reaches
                # ir_builtin_func with propensity=true (see the CallNode branch
                # below), so it records the tensor under propensity["declared"].
                # Checking only rule_lhs/rule_rhs here re-declared the same
                # torch::Tensor and produced a C++ "redeclaration of at::Tensor"
                # error (e.g. `s1 := HELP::intersection_s1(im4_pos, ...)` followed
                # by `jpx := im4_pos[0] + ...`).


                # why do we need to check the propensity?
                # well that's because the propensity table is used to store the variables that are used in the propensity function
                if !(arg.position.value in propensity_table["var_local_table"]["rule_lhs"]["declared"]) &&
                   !(arg.position.value in propensity_table["var_local_table"]["rule_rhs"]["declared"]) &&
                   !(arg.position.value in propensity_table["var_local_table"]["propensity"]["declared"])
                    ir = propensity_table["var_local_table"]["rule_lhs"][arg.position.value]
                    emit(context, ir)
                    push!(propensity_table["var_local_table"]["rule_lhs"]["declared"], arg.position.value)
                end

                emit(ir_builder, arg.position.value)

            elseif arg.position.value in propensity_table["var_local_table"]["rule_rhs"]["declared"]
                emit(ir_builder, arg.position.value)

            else
                # TODO: Check if it is inside the settings also (global variable)
                # TODO: Check where clause
                println("declared: ", propensity_table["var_local_table"]["rule_rhs"]["declared"])
                throw("Error: variable $(arg.position.value) not found in variable table. $(arg)")
            end

        elseif expression isa BinaryOpNode
            lhs = expression.lhs
            ir_definition!(lhs, ir_builder, propensity_table, context)
            emit(ir_builder, " $(expression.expression.position.value) ")
            rhs = expression.rhs
            ir_definition!(rhs, ir_builder, propensity_table, context)

        elseif expression isa CallNode
            call_args = expression.args
            function_node = expression.function_node

            func_name = get_value(function_node)

            # Function namespace
            namespace = function_node.namespace

            if !(namespace isa Nothing)
                namespace = namespace.position.value
            end

            func_args = call_args

            ir_value = IRBuilder([])

            # Loop through func args and emit them t
            ir_builtin_func(func_name, func_args, namespace,
                            # ir_builder, propensity_table, context, false, true)
                            ir_value, propensity_table, context, true, false)

            println("ir_builder after ir_builtin_func for definition: ", ir_builder.instructions)
            println("context after ir_builtin_func for definition: ", context.instructions)
            println("ir_value after ir_builtin_func for definition: ", build(ir_value))

            emit(ir_builder, build(ir_value))

            # ir_builtin_func(func_name, func_args, namespace,
                            # context, propensity_table, ir_builder)

        # TODO: check if it is just a regular digit or something
        elseif expression isa GroupNode
            # If it is a group node, we need to
            # traverse the expression inside the group
            inner_expr = expression.expression

            # Should emit parenthesis too but I don't know where to, if it is in the prop body
            # or in the ir builder. 

            emit(ir_builder, "(")
            ir_definition!(inner_expr, ir_builder, propensity_table, context)
            emit(ir_builder, ")")

        elseif expression isa LiteralNode
            # If it is a literal, we can just emit it
            emit(ir_builder, get_value(expression))

        elseif expression isa UnaryOpNode
            operation = expression.expression

            if operation.position.value == "~"
                ir = ir_distribution!(expression.operand, ir_builder, context, propensity_table, false, true)
                println("debug: ir after distribution: ", ir)
                emit(ir_builder, ir)

            elseif operation.position.value == "-"
                # Recurse so parameters/locals in the operand are resolved
                # (e.g. `-WOBBLE_ANGLE` -> `-settings.WOBBLE_ANGLE`).
                emit(ir_builder, "-")
                ir_definition!(expression.operand, ir_builder, propensity_table, context)

            else
                throw("Error: Unary operation $(operation.position.value) not recognized.")
            end

        elseif expression isa IndexAccessNode
            # Tensor index access: supports plain indexing and slicing
            ir_index_access(expression, ir_builder, propensity_table, context)

        elseif expression isa ArrayLiteralNode
            # Inline tensor literal: [1,2,3] -> torch::tensor({1,2,3}, torch::kFloat64)
            emit(ir_builder, ir_array_literal(expression, propensity_table, context))

        elseif expression isa StringNode
            # String literal \u2014 emit as quoted C++ string (used e.g. in einsum equations)
            emit(ir_builder, "\"$(expression.token.position.value)\"")

        else
            throw("Error: $(expression) not recognized.")
        end

    end

    function ir_indicator_func(arg_str)
        """
        Handles Indicator function
        """

        # TODO: check if values are in propensity table
        ir = "( $arg_str ? 1.0 : 0.0)"
        return ir
    end

    function ir_distribution_func(distribution_node, args)
        distribution = get_value(distribution_node)
        if distribution == "UniformDistribution"
            return "std::uniform_real_distribution<double>"
        elseif distribution == "NormalDistribution"
            return "std::normal_distribution<double>"
        elseif distribution == "GammaDistribution"
            return "std::gamma_distribution<double>"
        elseif distribution == "LogNormalDistribution"
            return "std::lognormal_distribution<double>"
        else
            throw("Unknown distribution: $distribution")
        end
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
            func_name = get_value(call_node.function_node)

            namespace = call_node.function_node.namespace

            args = call_node.args

            # tensordot special case
            # torch::tensordot(self, other, IntArrayRef dims_self, IntArrayRef dims_other)
            # Dim args must be rendered as {0,1} (IntArrayRef), not torch::tensor(...).
            # FoxFlow syntax: tensordot(A, B, [dims_self...], [dims_other...])
            if func_name == "tensordot"
                if length(args) != 4
                    throw("tensordot requires 4 arguments: tensordot(A, B, [dims_self], [dims_other])")
                end
                t1 = traverse_group_node(GroupNode(args[1]), variables)
                t2 = traverse_group_node(GroupNode(args[2]), variables)
                d1 = _arr_to_intarrayref(args[3])
                d2 = _arr_to_intarrayref(args[4])
                return "torch::tensordot($t1, $t2, $d1, $d2)"
            end

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
                elseif arg isa IndexAccessNode
                    return traverse_group_node(GroupNode(arg), variables)
                elseif arg isa CallNode
                    return traverse_group_node(GroupNode(arg), variables)
                elseif arg isa ArrayLiteralNode
                    return ir_array_literal(arg)
                elseif arg isa StringNode
                    return "\"$(arg.token.position.value)\""
                else
                    return arg.token
                end
                end, args), ", ")

            # println("doing call node in group node traversal: ", func_name, " with args: ", arg_str)
            # println("doing call: ", call_node)

            if !(namespace isa Nothing)
                namespace = namespace.position.value
            end

            # Route built-in functions through the same mapping used by ir_builtin_func
            # so that nested calls like rand_matrix(4,4) or transpose(W2) inside
            # mat_add / mat_mul args get their C++ equivalents emitted correctly.
            if func_name in BUILT_IN_FUNC
                ir = _builtin_to_cpp(func_name, arg_str)
            else
                ir = add_namespace_identifier(func_name, namespace, arg_str)
            end
            return ir

        elseif group_node.expression isa LiteralNode
            return string(get_value(group_node.expression))

        elseif group_node.expression isa IdentifierNode
            var_name = get_value(group_node.expression)
            push!(variables, var_name)
            return get_value(group_node.expression)
        elseif group_node.expression isa IndexAccessNode
            idx_node = group_node.expression
            obj_str = traverse_group_node(GroupNode(idx_node.object), variables)
            idx_strs = map(idx -> begin
                if idx isa SliceNode
                    # Generate slice expression inline
                    s = "torch::indexing::None"
                    e = "torch::indexing::None"
                    st = nothing
                    if !isnothing(idx.start)
                        s = traverse_group_node(GroupNode(idx.start), variables)
                    end
                    if !isnothing(idx.stop)
                        e = traverse_group_node(GroupNode(idx.stop), variables)
                    end
                    if !isnothing(idx.step)
                        st = traverse_group_node(GroupNode(idx.step), variables)
                    end
                    if isnothing(idx.start) && isnothing(idx.stop) && isnothing(idx.step)
                        "torch::indexing::Slice()"
                    elseif isnothing(st)
                        "torch::indexing::Slice($s, $e)"
                    else
                        "torch::indexing::Slice($s, $e, $st)"
                    end
                else
                    traverse_group_node(GroupNode(idx), variables)
                end
            end, idx_node.indices)
            if has_slice(idx_node.indices)
                return "$(obj_str).index({$(join(idx_strs, ", "))})"
            elseif length(idx_strs) == 1
                return "$(obj_str)[$(idx_strs[1])].template item<double>()"
            else
                return "$(obj_str).index({$(join(idx_strs, ", "))}).template item<double>()"
            end
        elseif group_node.expression isa ArrayLiteralNode
            return ir_array_literal(group_node.expression)
        elseif group_node isa GroupNode
            return traverse_group_node(group_node.expression, variables)
        else
            throw("Unknown group node expression type: $(typeof(group_node.expression))")
        end
    end

    """
        _arr_to_intarrayref(node) -> String

    Renders an ArrayLiteralNode (or integer IdentifierNode/IntegerNode) as a
    C++ initializer-list suitable for IntArrayRef, e.g. `{0, 1, 2}`.
    Used by tensordot and permute to avoid torch::tensor(...).
    """
    function _arr_to_intarrayref(node)
        if node isa ArrayLiteralNode
            inner = _ir_array_literal_inner(node)
            return inner  # already produces {0, 1, ...}
        elseif node isa IntegerNode
            return "{$(get_value(node))}"
        elseif node isa IdentifierNode
            return get_value(node)
        else
            return "{$(get_value(node))}"
        end
    end

    """
        _builtin_to_cpp(func_name, arg_str) -> String

    Pure string-level mapping from a FoxFlow built-in function name and its
    already-rendered argument string to the corresponding C++ expression.
    Used by both ir_builtin_func and traverse_group_node so that nested
    built-in calls (e.g. rand_matrix inside mat_add) are always expanded.
    Returns `nothing` for functions that are NOT pure string mappings
    (e.g. indicator, distribution functions that need extra context).
    """
    function _builtin_to_cpp(func_name, arg_str)
        if func_name == "heaviside"
            return "DGGML::$func_name($arg_str)"
        elseif func_name == "sqrt"
            return "sqrt($arg_str)"
        elseif func_name == "cos"
            return "cos($arg_str)"
        elseif func_name == "arccos"
            return "acos($arg_str)"
        elseif func_name == "sin"
            return "sin($arg_str)"
        elseif func_name == "abs"
            return "abs($arg_str)"
        elseif func_name == "pow"
            return "pow($arg_str)"
        elseif func_name == "zeros_matrix"
            return _array_ctx_fill("zeros", arg_str, "torch::zeros({$arg_str}, torch::kFloat64)")
        elseif func_name == "ones_matrix"
            return _array_ctx_fill("ones", arg_str, "torch::ones({$arg_str}, torch::kFloat64)")
        elseif func_name == "dot"
            return "dot($arg_str).template item<double>()"
        elseif func_name == "to_tensor"
            return "tf::to_tensor($arg_str)"
        elseif func_name == "to_array"
            parts = split(arg_str, ", ", limit=2)
            length(parts) == 2 || throw("to_array(tensor, N) expects two arguments.")
            return "tf::to_array<$(parts[2])>($(parts[1]))"
        elseif func_name == "rand_matrix"
            return "torch::rand({$arg_str}, torch::kFloat64)"
        elseif func_name == "eye_matrix"
            return "torch::eye($arg_str, torch::kFloat64)"
        elseif func_name == "mat_add"
            parts = split(arg_str, ", ", limit=2)
            return "($(parts[1]) + $(parts[2]))"
        elseif func_name == "mat_mul"
            parts = split(arg_str, ", ", limit=2)
            return "torch::mm($(parts[1]), $(parts[2]))"
        elseif func_name == "mat_dot"
            parts = split(arg_str, ", ", limit=2)
            return "torch::mv($(parts[1]), $(parts[2]))"
        elseif func_name == "transpose"
            return "($arg_str).t()"
        elseif func_name == "tensordot"
            # tensordot(A, B, {dims_self}, {dims_other})
            # arg_str contains 4 parts: tensor1, tensor2, intarrayref1, intarrayref2
            parts = split(arg_str, ", ", limit=4)
            return "torch::tensordot($(parts[1]), $(parts[2]), $(parts[3]), $(parts[4]))"
        elseif func_name == "einsum"
            parts = split(arg_str, ", ", limit=2)
            tensor_list = join(split(parts[2], ", "), ", ")
            return "torch::einsum($(parts[1]), {$tensor_list})"
        elseif func_name == "permute"
            parts = split(arg_str, ", ", limit=2)
            dims = join(split(parts[2], ", "), ", ")
            return "$(parts[1]).permute({$dims})"
        elseif func_name == "backward"
            return "$arg_str.backward()"
        elseif func_name == "grad"
            return "$arg_str.grad()"
        else
            return nothing  # needs propensity_table context — handled by ir_builtin_func
        end
    end

    function ir_builtin_func(func_name, args, namespace, ir_builder,
            propensity_table, prop_body_ir, propensity=false, where_clause=false)
        """
        Generates the IR for a built-in function.

        TODO: Add gamma function and lognormal distribution
        """

        # FIXME: Why is this namespace here?
        # name_space = "FractureNetwork"
        
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

                    tmp_namespace = arg.namespace

                    var_name = nothing
                    if !(tmp_namespace isa Nothing)
                        tmp_namespace = tmp_namespace.position.value
                        var_name = "$(tmp_namespace)::$(get_value(arg))"
                    else
                        var_name = get_value(arg)
                    end

                    push!(variables, var_name)

                    return var_name
                elseif arg isa LiteralNode
                    return string(get_value(arg))
                elseif arg isa CallNode
                    return traverse_group_node(GroupNode(arg), variables)
                elseif arg isa IndexAccessNode
                    return traverse_group_node(GroupNode(arg), variables)
                elseif arg isa ArrayLiteralNode
                    return ir_array_literal(arg)
                elseif arg isa StringNode
                    # Emit as a quoted C++ string literal for e.g. einsum equations
                    return "\"$(arg.token.position.value)\""
                else
                    return arg.token
                end
        end, args), ", ")

        # Emit out variables
        # I should probably emit something if it's from the where clause
        # if we detect a variable that needs to be used, we declare it out
        # innto the context
        map( (var) -> begin
            println("variable in builtin func: ", var)
            if (((var in propensity_table["var_local_table"]["rule_lhs"]["declared"])
                || (var in propensity_table["var_local_table"]["rule_rhs"]["declared"]))
                && (where_clause == true || propensity == false))
                    println("rule_lhs declared vars: ", propensity_table["var_local_table"]["rule_lhs"]["declared"])
                    nothing
            else
                # And its not already inside prop_body_ir
                println("collect(keys(propensity_table)): ", collect(keys(propensity_table)))
                if var in collect(keys(propensity_table))

                    ir = propensity_table[var]

                    if propensity == true
                        println("prop_body_ir.instructions for $(var): ", prop_body_ir.instructions)

                        # if (ir in prop_body_ir.instructions) || (var in propensity_table["var_local_table"]["rule_rhs"]["declared"])
                        if (ir in prop_body_ir.instructions) || (var in propensity_table["var_local_table"]["propensity"]["declared"])
                            nothing
                        else
                            # has been declared, assume it's in the context of the where clause, so we need to declare it in the prop body ir.
                            push!(propensity_table["var_local_table"]["propensity"]["declared"], var)
                            emit(prop_body_ir, ir)
                        end

                    else
                        push!(propensity_table["var_local_table"]["rule_rhs"]["declared"], var)

                        # Emit the declaration as a preamble side-effect so it never
                        # lands inside the expression string being built in ir_builder
                        # (e.g. `x = a + indicator(b > a)` with a, b first read here).
                        emit(prop_body_ir, ir)
                    end

                elseif var in collect(keys(propensity_table["var_local_table"]["rule_rhs"]))
                    # RHS-only binding read for the first time (e.g. a freshly
                    # created node attribute referenced inside a where-clause
                    # expression). Emit its declaration from the RHS match
                    # (rhs[m2[...]]) into the preamble so the variable is in
                    # scope; previously this was a no-op and only compiled by
                    # accident when a same-named LHS binding leaked in from a
                    # prior rule. Mirrors the propensity_table branch above.
                    ir = propensity_table["var_local_table"]["rule_rhs"][var]

                    if strip(string(ir)) == var
                        # The stored IR is just the variable's own name: the
                        # variable was assigned earlier in the where-body (e.g.
                        # `F = D2Q9::feq(...)`) and is already a declared local.
                        # Emitting the name would produce a bare, semicolon-less
                        # line in the generated C++.
                        nothing
                    elseif propensity == true
                        if (ir in prop_body_ir.instructions) || (var in propensity_table["var_local_table"]["propensity"]["declared"])
                            nothing
                        else
                            push!(propensity_table["var_local_table"]["propensity"]["declared"], var)
                            emit(prop_body_ir, ir)
                        end
                    else
                        push!(propensity_table["var_local_table"]["rule_rhs"]["declared"], var)
                        emit(prop_body_ir, ir)
                    end

                # NOTE: It's already been declared do nothing.
                elseif var in propensity_table["var_local_table"]["rule_rhs"]["declared"]
                    nothing

                elseif var in collect(keys(propensity_table["parameter_table"]))

                    ir = "auto $var =  settings.$var;\n"

                    if ir in prop_body_ir.instructions
                        println("skipping")
                        nothing
                    else
                        println("Emitting variable from parameter table: ", var)
                        emit(prop_body_ir, ir)
                    end

                else
                    println("Variable $var not found in variable table, propensity table, or parameter table for built-in function $func_name.")
                    println("var local table lhs: ", collect(keys(propensity_table["var_local_table"]["rule_lhs"]["declared"])))
                    println("var local table lhs: ", collect(keys(propensity_table["var_local_table"])))
                    throw("Variable $var not found in var table for built-in function $func_name.")
                end
                end
            end,
            collect(variables)
        )

        # Try the pure string mapping first
        ir = _builtin_to_cpp(func_name, arg_str)
        if ir === nothing
            # Functions that need propensity_table / extra context
            if func_name == "normal_distr"
                ir = "DGGML::normal_distr()"
            elseif func_name == "uniform_distr"
                ir = "DGGML::uniform_distr()"
            elseif func_name == "inverse"
                ir = "DGGML::inverse()"
            elseif func_name == "indicator"
                ir = ir_indicator_func(arg_str)
            elseif func_name == "gamma_distr"
                ir = "std::gamma_distribution<double>($arg_str)"
            elseif func_name in collect(keys(propensity_table["function_table"]))
                ir = add_namespace_identifier(func_name, namespace, arg_str)
            else
                throw("Unknown built-in function: $func_name")
            end
        end

        emit(ir_builder, ir)
    end

    function ir_distribution!(node, ir_builder, ir_context, propensity_table, propensity=false, where_clause=false)
        """
        Generates the IR for a distribution.
        """

        # NOTE: Let's just for now assume that we are only dealing with
        # arguments that have identifiers already loaded.
        # Probably should add a scope dictionary...

        distribution = node.function_node
        args = node.args

        vars = Set{String}()
        parsed_args = []
        for arg in args
            arg_ir = IRBuilder([])

            if where_clause == true
                ir_definition!(arg, arg_ir, propensity_table, ir_context)
                # ir_prop_expr!(arg, arg_ir, propensity_table, arg_ir, propensity)
                # ir_prop_expr!(arg, arg_ir, propensity_table, arg_ir, propensity)
            else
                ir_prop_expr!(arg, arg_ir, propensity_table, arg_ir, propensity)
            end

            push!(parsed_args, join(arg_ir.instructions))
        end

        args = join(parsed_args, ", ")

        # FIXME: Handle local variable scoping?
        # Check to see if a random device is already been defined
        if propensity_table["var_local_table"]["random_device"] == false
            propensity_table["random_device"] = true
            ir0 = "std::random_device random_device;\n"
            # emit(ir_builder, ir0)
            #

            println("IR_DISTR node: ", node)
            println("emitting in distribution: ", ir0)
            emit(ir_context, ir0)

            ir1 = "std::mt19937 random_engine(random_device());\n"
            # emit(ir_builder, ir1)
            println("emitting in distribution: ", ir1)
            emit(ir_context, ir1)
            propensity_table["var_local_table"]["random_device"] = true

        end

        func_name = ir_distribution_func(distribution, args)
        
        # TODO: Check dist_type and return the correct one here instead.
        # ir2 = "std::uniform_real_distribution<double>("*args *")(random_engine)"
        ir2 = "$func_name("*args *")(random_engine)"
        return ir2
    end
