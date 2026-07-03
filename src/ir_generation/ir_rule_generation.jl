"""
    IRRuleGeneration:
        Module for generating Intermediate Representation (IR) for rules in a grammar.
"""

module IRRuleGeneration

    import ..IRUtils: get_value, convert_type_name, is_list_type_tensor
    import ..IRBuildUtils: emit, build, IRBuilder, build_sameline
    import ...AstNodes: RuleNode, IdentifierNode, TypeInstanceNode,
            UndirectedTypeEdgeNode, BinaryOpNode, ExpressionNode, ParameterNode,
            WithClauseNode, SolveClauseNode, UnaryOpNode, FunctionNode, GroupNode,
            CallNode, LiteralNode, DefinitionNode, ODENode, BindingVariableNode, WhereClauseNode, NamedParameterNode,
            IndexAccessNode, ArrayLiteralNode, SliceNode,
            IntegerNode, FloatNode, StringNode,
            FunctionDefinitionNode, FunctionSignatureNode, FunctionArgNode,
            FunctionBodyExpressionNode, FunctionDefinitionExpressionNode, ReturnNode

    export ir_rules_section!

    using OrderedCollections

    const BUILT_IN_FUNC = [
        "heaviside", "sqrt", "normal_distr",
        "uniform_distr", "cos", "sin", "inverse", "pow", "indicator",
        # Matrix / tensor operations (FixedList<<M,N,Float>> treated as 2-D tensor)
        "zeros_matrix", "ones_matrix", "rand_matrix", "eye_matrix",
        "mat_add", "mat_mul", "mat_dot", "transpose",
        # N-D tensor operations
        "tensordot", "einsum", "permute",
        # Autograd
        "backward", "grad",
        # Single-call automatic differentiation: autodiff(expr, var)
        "autodiff"
    ]

    # Module-level type namespace, set by ir_rules_section! at entry
    _type_namespace = Ref("")

    # Global counter for generating unique array literal temp variable names
    _arr_tmp_counter = Ref(0)

    """
        ir_array_literal(node::ArrayLiteralNode) -> String

    Recursively converts an ArrayLiteralNode into a C++ torch::tensor literal.
    - 1D: [1, 2, 3]         -> torch::tensor({1, 2, 3}, torch::kFloat64)
    - 2D: [[1,2], [3,4]]    -> torch::tensor({{1, 2}, {3, 4}}, torch::kFloat64)
    """
    function ir_array_literal(node, propensity_table=nothing, context=nothing)
        inner = _ir_array_literal_inner(node, propensity_table, context)
        return "torch::tensor($inner, torch::kFloat64)"
    end

    function _ir_array_literal_inner(node, propensity_table=nothing, context=nothing)
        if node isa ArrayLiteralNode
            parts = [_ir_array_literal_inner(el, propensity_table, context) for el in node.elements]
            return "{" * join(parts, ", ") * "}"
        elseif node isa LiteralNode
            return string(get_value(node))
        elseif node isa UnaryOpNode
            op = node.expression.position.value
            return op * string(get_value(node.operand))
        elseif node isa IdentifierNode
            return get_value(node)
        elseif node isa BinaryOpNode
            lhs = _ir_array_literal_inner(node.lhs, propensity_table, context)
            op = node.expression.position.value
            rhs = _ir_array_literal_inner(node.rhs, propensity_table, context)
            return "($lhs $op $rhs)"
        elseif node isa GroupNode
            inner_expr = node.expression
            lhs = _ir_array_literal_inner(inner_expr, propensity_table, context)
            return lhs
        elseif node isa IndexAccessNode
            # Generate C++ for indexed expressions like cb_count[0]
            obj_str = _ir_array_literal_inner(node.object, propensity_table, context)
            idx_strs = map(idx -> begin
                if idx isa SliceNode
                    s = isnothing(idx.start) ? "torch::indexing::None" : _ir_array_literal_inner(idx.start, propensity_table, context)
                    e = isnothing(idx.stop) ? "torch::indexing::None" : _ir_array_literal_inner(idx.stop, propensity_table, context)
                    st = isnothing(idx.step) ? nothing : _ir_array_literal_inner(idx.step, propensity_table, context)
                    if isnothing(idx.start) && isnothing(idx.stop) && isnothing(idx.step)
                        "torch::indexing::Slice()"
                    elseif isnothing(st)
                        "torch::indexing::Slice($s, $e)"
                    else
                        "torch::indexing::Slice($s, $e, $st)"
                    end
                else
                    _ir_array_literal_inner(idx, propensity_table, context)
                end
            end, node.indices)
            if has_slice(node.indices)
                return "$(obj_str).index({$(join(idx_strs, ", "))})"
            elseif length(idx_strs) == 1
                return "$(obj_str)[$(idx_strs[1])].template item<double>()"
            else
                return "$(obj_str).index({$(join(idx_strs, ", "))}).template item<double>()"
            end
        elseif node isa CallNode
            # Handle function calls like sqrt(x), sin(y), etc.
            func_name = get_value(node.function_node)
            arg_strs = [_ir_array_literal_inner(arg, propensity_table, context) for arg in node.args]
            return "$(func_name)($(join(arg_strs, ", ")))"
        else
            return string(get_value(node))
        end
    end


    """
        ir_index_element(idx, ir_builder, propensity_table, context)

    Generate C++ for a single index element. Returns a string.
    - For SliceNode: emits torch::indexing::Slice(start, stop, step)
    - For other nodes: emits the expression value directly
    """
    function ir_index_element(idx, ir_builder, propensity_table, context)
        if idx isa SliceNode
            return ir_slice_expr(idx, ir_builder, propensity_table, context)
        else
            idx_ir = IRBuilder([])
            ir_definition!(idx, idx_ir, propensity_table, context)
            return build_sameline(idx_ir)
        end
    end

    """
        ir_slice_expr(slice::SliceNode, ...) -> String

    Generates C++ for a SliceNode:
      SliceNode(nothing, nothing, nothing) -> torch::indexing::Slice()
      SliceNode(start, stop, nothing)      -> torch::indexing::Slice(start, stop)
      SliceNode(start, stop, step)         -> torch::indexing::Slice(start, stop, step)
      SliceNode(nothing, stop, nothing)    -> torch::indexing::Slice(torch::indexing::None, stop)
      SliceNode(start, nothing, nothing)   -> torch::indexing::Slice(start, torch::indexing::None)
      etc.
    """
    function ir_slice_expr(slice, ir_builder, propensity_table, context)
        start_str = "torch::indexing::None"
        stop_str = "torch::indexing::None"
        step_str = nothing

        if !isnothing(slice.start)
            s_ir = IRBuilder([])
            ir_definition!(slice.start, s_ir, propensity_table, context)
            start_str = build_sameline(s_ir)
        end
        if !isnothing(slice.stop)
            s_ir = IRBuilder([])
            ir_definition!(slice.stop, s_ir, propensity_table, context)
            stop_str = build_sameline(s_ir)
        end
        if !isnothing(slice.step)
            s_ir = IRBuilder([])
            ir_definition!(slice.step, s_ir, propensity_table, context)
            step_str = build_sameline(s_ir)
        end

        # Bare colon: Slice() means select all
        if isnothing(slice.start) && isnothing(slice.stop) && isnothing(slice.step)
            return "torch::indexing::Slice()"
        elseif isnothing(step_str)
            return "torch::indexing::Slice($start_str, $stop_str)"
        else
            return "torch::indexing::Slice($start_str, $stop_str, $step_str)"
        end
    end

    """
        has_slice(indices) -> Bool

    Returns true if any index in the list is a SliceNode.
    """
    function has_slice(indices)
        return any(idx -> idx isa SliceNode, indices)
    end

    """
        ir_index_access(expression, ir_builder, propensity_table, context)

    Generate C++ for an IndexAccessNode. Uses torch::indexing API when any
    index is a slice, otherwise uses the simpler [...] or .index({...}) syntax.
    """
    function ir_index_access(expression, ir_builder, propensity_table, context)
        ir_definition!(expression.object, ir_builder, propensity_table, context)
        indices = expression.indices

        if has_slice(indices)
            # Any slice present: must use .index({...}) with torch::indexing types
            idx_parts = String[]
            for idx in indices
                push!(idx_parts, ir_index_element(idx, ir_builder, propensity_table, context))
            end
            emit(ir_builder, ".index({" * join(idx_parts, ", ") * "})")
        elseif length(indices) == 1
            idx_str = ir_index_element(indices[1], ir_builder, propensity_table, context)
            emit(ir_builder, "[$idx_str].template item<double>()")
        else
            idx_parts = String[]
            for idx in indices
                push!(idx_parts, ir_index_element(idx, ir_builder, propensity_table, context))
            end
            emit(ir_builder, ".index({" * join(idx_parts, ", ") * "}).template item<double>()")
        end
    end

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

                if !(arg.position.value in propensity_table["var_local_table"]["rule_lhs"]["declared"])
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

            elseif operation.position.value == "-"
                ir = "-" * get_value(expression.operand)

            else
                throw("Error: Unary operation $(operation.position.value) not recognized.")
            end
            emit(ir_builder, ir)

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

            # ── tensordot special case ──────────────────────────────────────────────
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
            # ────────────────────────────────────────────────────────────────────────

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
            return "torch::zeros({$arg_str}, torch::kFloat64)"
        elseif func_name == "ones_matrix"
            return "torch::ones({$arg_str}, torch::kFloat64)"
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
                || (var in propensity_table["var_local_table"]["rule_rhs"]["declared"])) && where_clause == true)
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

                        if where_clause == true
                            # Emit declaration as a preamble side-effect so it doesn't
                            # pollute the expression string (e.g. in Neural ODE RHS).
                            emit(prop_body_ir, ir)
                        else
                            emit(ir_builder, ir)
                        end
                    end

                elseif var in collect(keys(propensity_table["var_local_table"]["rule_rhs"]))
                    nothing

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

    function ir_rules_section!(
            ast, rules_table, symbol_tables,
            propensity_table, type_namespace;
            stage_index::Union{Int, Nothing}=nothing
        )
        """
        Generate intermediate rules
        for the section.
        """

        rule_section_name = get_value(ast.name)
        global rule_namespace = rule_section_name
        _type_namespace[] = type_namespace

        # Generate a unique header guard per stage
        guard_suffix = stage_index !== nothing ? "_STAGE_$(stage_index)" : ""
        guard_name = "DGGML_RULES$(guard_suffix)_HPP"

        rule_section_header = [
            "#ifndef $guard_name",
            "#define $guard_name",
            "#include \"types.h\"",
            "#include \"parameters.h\"",
            "#include \"functions.h\"",
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
        var_local_table["propensity"] = Dict()
        var_local_table["propensity"]["declared"] = []

        var_local_table["rule_rhs"] = Dict()
        var_local_table["rule_rhs"]["declared"] = []

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

    # ── Symbolic ODE pretty-printer ──────────────────────────────────────
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
                ir_value = IRBuilder([])

                func_node = expression.function_node

                func_name = get_value(func_node)
                namespace = func_node.namespace
                if !(namespace isa Nothing)
                    namespace = namespace.position.value
                end
                args = expression.args

                ir_builtin_func(func_name, args, namespace,
                            ir_value, propensity_table,
                            context, false, true)

                emit(ir_builder, build(ir_value))

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

            ir_value = IRBuilder([])

            func_node = expression.lhs.function_node
            func_name = get_value(func_node)
            namespace = func_node.namespace
            if !(namespace isa Nothing)
                namespace = namespace.position.value
            end
            args = expression.lhs.args

            ir_builtin_func(func_name, args, namespace, ir_value,
                propensity_table, context, true, false)

            # FIXME
            emit(ir_builder, build(ir_value))

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

                        # Only sync position[0..2] if this is the Position attribute (first FixedList, index 1)
                        if bv_attr_pos == 1
                            pos_count = min(3, tensor_size)
                            for i in 0:(pos_count - 1)
                                emit(ir_builder, "varset.insert(&lhs[m1[$bv_pos]].position[$i]);")
                            end
                        end

                        # Store tensor info for use in ODE lambda
                        tensor_binding_info[dep_vars] = (tensor_size, ptr_name, bv_pos)

                        # var_attr_loc maps dep_var to a per-element accessor pattern
                        var_attr_loc[dep_vars] = ptr_name
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

                if ro_type_str == "torch::Tensor" && ro_tensor_size > 0
                    # Tensor read-only variable
                    ro_attr = ro_info.cpp_var.name
                    ro_node_type = ro_info.type
                    tensor_ref = "std::get<$type_namespace::$ro_node_type>(lhs[m1[$ro_pos]].data).$ro_attr"
                    ptr_name = "tensor_ptr_$(ro_pos)_$(ro_attr_pos)"
                    emit(ir_builder, "auto &tensor_ref_$(ro_pos)_$(ro_attr_pos) = $tensor_ref;")
                    emit(ir_builder, "double* $ptr_name = tensor_ref_$(ro_pos)_$(ro_attr_pos).template data_ptr<double>();")
                    for i in 0:(ro_tensor_size - 1)
                        emit(ir_builder, "varset.insert(&$(ptr_name)[$i]);")
                    end
                    tensor_binding_info[ro_var] = (ro_tensor_size, ptr_name, ro_pos)
                    var_attr_loc[ro_var] = ptr_name
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

#     function ir_solve_ode_expr(expression, ir_builder,
#           var_loc_attr, assgn_info, dep_vars, propensity_table)
#         """
#         Generates ir for an ODE expression
#         by traversing the expression tree.
#         """
# 
#         traverse_ode_expr(expression, ir_builder, var_loc_attr,
#             assgn_info, dep_vars, propensity_table)
#     end

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

        # ── Collect all variables referenced in ODE RHS expressions ──
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
                bv_attr = lhs_assgn_to_node[dep_var].cpp_var.name
                bv_node_type = lhs_assgn_to_node[dep_var].type
                tensor_ref = "std::get<$type_namespace::$bv_node_type>(lhs[m1[$bv_pos]].data).$bv_attr"
                emit(ir_builder, "double* $ptr_name = $tensor_ref.template data_ptr<double>();")
                # Register per-element refs so traverse_ode_expr can read them via NV_Ith_S
                for i in 0:(tsize - 1)
                    push!(dep_vars, "$(ptr_name)[$i]")
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
        for (dep_var, info) in tensor_binding_info
            tsize, ptr_name, bv_pos = info
            attr_pos = lhs_assgn_to_node[dep_var].cpp_var.index
            if attr_pos == 1
                pos_count = min(3, tsize)
                for i in 0:(pos_count - 1)
                    tensor_ref = "$(ptr_name)[$i]"
                    pos_ref = "lhs[m1[$bv_pos]].position[$i]"
                    emit(ir_builder, "NV_Ith_S(ydot, varmap[&$pos_ref]) += NV_Ith_S(ydot, varmap[&$tensor_ref]);")
                end
            end
        end

        # ── Emit symbolic ODE equations as C++ comments ──
        if !isempty(symbolic_equations)
            # Print to Julia stdout during codegen
            println("  ── Symbolic ODE system ──")
            for eq in symbolic_equations
                println("    $eq")
            end

            emit(ir_builder, "// ── Symbolic ODE system ──")
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

    function ir_rule!(ir_builder, rule_node,
        type_namespace, symbol_tables, propensity_table)
        """
        Generates the IR for a rule, including its header, nodes, and edges.
        """

        function visit_rule_node!(each_node, graph_table, ir_builder,
            graph_name, type_namespace=type_namespace)
            """
            Visits each node in the rule and builds the graph
            for the lhs and rhs of the rule
            """

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
        elseif modify_clause isa SolveClauseNode
            rule_hdr = "DGGML::SolvingRule<GT> $rule_name(\"$rule_name\", $lhs_name, $lhs_name,"
        else
            throw("Unknown modify clause type: $(typeof(modify_clause))")
        end

        emit(ir_builder, rule_hdr)

        # println("lhs_types: ", lhs_types)

        lhs_names, lhs_types, lhs_order_of_nodes = build_node_pos_to_type(lhs_node_type)


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
                        # tensor elements + position sync
                        num_vars += tsize + min(3, tsize)
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
                        num_vars += tsize
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

        # ── Validate: user must declare exactly as many parameters per node
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

    function ir_where_left_clause!(lhs, type_namespace, symbol_tables,
                                   propensity_table, ir_builder)
        """
        IR Where left clause adds nodes attributes to propensity table.
        """

        for (ix, param) in enumerate(lhs)

            node_loc = param[1]
            node_type = param[3]
            attr_type = param[4][1]
            attr_name = param[4][2]
            attr_loc = param[4][3]

            user_param_name = param[end][1]

            ir = nothing
            if attr_type == "torch::Tensor"
                # Tensor attribute — always read from .data
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
            if attr_type == "torch::Tensor"
                # Tensor attribute — always read from .data
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
            rhs_assgn_to_node, symbol_tables, type_namespace, propensity_table
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

                # setting them equal
                ir = "$rhs_ir = $lhs_ir;\n"
                emit(ir_builder, ir)

                if rhs_rule_param.cpp_var.type == "torch::Tensor" && rhs_rule_param.cpp_var.index == 1
                    # For Position tensor attribute, also copy spatial .position from LHS to RHS
                    lhs_node_idx = lhs_assgn_to_node[var_name].index
                    rhs_node_idx = rhs_rule_param.index
                    pos_ir = "std::copy(std::begin(lhs[m1[$lhs_node_idx]].position), std::end(lhs[m1[$lhs_node_idx]].position), std::begin(rhs[m2[$rhs_node_idx]].position));\n"
                    emit(ir_builder, pos_ir)
                elseif rhs_rule_param.cpp_var.type != "torch::Tensor" && rhs_rule_param.cpp_var.index <= 3
                    rhs_pos = create_pos_cpp_var(rhs_rule_param, "rhs", "m2")
                    lhs_pos = create_pos_cpp_var(lhs_assgn_to_node[var_name], "lhs", "m1")
                    pos_ir = "$rhs_pos = $lhs_pos;\n"
                    emit(ir_builder, pos_ir)
                end
            end
        end

        map(
            (assign_node) -> begin

                ir_where_assignment!(assign_node, type_namespace, rhs_assgn_to_node,
                                     propensity_table, ir_builder, grad_vars)
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

                if !(arg.position.value in propensity_table["var_local_table"]["rule_rhs"]["declared"])
                    ir = propensity_table["var_local_table"]["rule_rhs"][arg.position.value]
                    emit(prop_body_ir, ir)
                    push!(propensity_table["var_local_table"]["rule_rhs"]["declared"], arg.position.value)
                end

                emit(ir_builder, arg.position.value)

            elseif arg.position.value in propensity_table["var_local_table"]["rule_rhs"]["declared"]

                emit(ir_builder, arg.position.value)

            # It is in the lhs of the rule
            elseif arg.position.value in collect(keys(propensity_table["var_local_table"]["rule_lhs"]))

                # Check if it is already emitted
                # If it is already emitted, we should not emit it again.
                #

                if !(arg.position.value in propensity_table["var_local_table"]["rule_lhs"]["declared"])
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
                emit(ir_builder, "torch::tensor({" * join(elem_names, ", ") * "}, torch::kFloat64)")
            else
                # All elements are simple literals — safe to inline directly
                elem_strs = String[]
                for el in function_node.elements
                    el_ir = IRBuilder([])
                    ir_prop_expr!(el, el_ir, propensity_table, prop_body_ir, propensity)
                    push!(elem_strs, build_sameline(el_ir))
                end
                emit(ir_builder, "torch::tensor({" * join(elem_strs, ", ") * "}, torch::kFloat64)")
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

    function is_self_referential(code_line::String)
        # ^(\w+)      : Capture the variable name at the start
        # \s*=        : Match the equals sign
        # .*?         : Non-greedy match for any characters in between
        # \b\1\b      : Match the exact same word from the capture group
        pattern = r"^(\w+)\s*=.*?\b\1\b"

        return occursin(pattern, strip(code_line))
    end

    # ── helpers for local (where-body) function definitions ──────────────────

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

    # ─────────────────────────────────────────────────────────────────────────

    function ir_where_assignment!(assign_node, type_namespace,
                                 rhs_assgn_to_node,
                                 propensity_table, ir_builder,
                                 grad_vars=Set{String}())
        """
        Handles the assignment in the where clause.
        """

        # FIXME: if you assign a node to another node who is just assigned it will not work.

        # ── Bare call statement (e.g. backward(loss)) ──
        # Parsed as a CallNode with no LHS.  Emit as a side-effecting statement.
        if assign_node isa CallNode
            ir_value = IRBuilder([])
            ir_definition!(GroupNode(assign_node), ir_value, propensity_table, ir_builder)
            ir_str = build_sameline(ir_value)
            emit(ir_builder, "$ir_str;")
            return
        end

        # ── Local function definition (emitted as a C++ lambda) ──
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
            type = convert_type_name(get_value(assign_node.type.name))
            value = assign_node.value

            # ── Special case: autodiff(expr, var) ──────────────────────────────
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
            # ── end autodiff special case ───────────────────────────────────────

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

            ir_definition!(value, ir_value, propensity_table, ir_builder)

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
            if has_slice(idx_node.indices)
                ir = "$cpp_var.index_put_({$(join(idx_parts, ", "))}, $ir_value_str);"
            elseif length(idx_parts) == 1
                ir = "$cpp_var[$(idx_parts[1])] = $ir_value_str;"
            else
                ir = "$cpp_var.index_put_({$(join(idx_parts, ", "))}, $ir_value_str);"
            end
            emit(ir_builder, ir)

            # Also sync .position for single element if it's the Position attribute (index 1)
            if rhs_assgn_node.cpp_var.type == "torch::Tensor" && length(idx_parts) == 1 && rhs_assgn_node.cpp_var.index == 1
                rhs_node_idx = rhs_assgn_node.index
                pos_sync = "if ($(idx_parts[1]) < 3) { rhs[m2[$rhs_node_idx]].position[$(idx_parts[1])] = static_cast<double>($ir_value_str); }"
                emit(ir_builder, pos_sync)
            end
        else

            # Check if it is an assign node or just an intermediate expression
            name = get_value(assign_node.name)

            # println("assign_node: ", assign_node)
            # println("assign_node: ", assign_node.name)
            # println("assign_node: ", assign_node.value)
            # println("rhs_assgn_to_node: ", rhs_assgn_to_node)

            type = rhs_assgn_to_node[name].cpp_var.type

            # Need to fetch type from types table.

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
            ir_prop_expr!(value, ir_value, propensity_table, ir_builder, false)
            ir_value = build_sameline(ir_value)

            # Adds the definition to propensity table local scope.
            propensity_table["var_local_table"]["rule_rhs"][name] = "$name"

            # if occursin(name, ir_value)
            if is_self_referential(ir_value)
                println("propensity_table: ", propensity_table["var_local_table"]["rule_rhs"])
                println("ir_value: ", ir_value)

                println(assign_node)

                throw(
                      "Error: Self-referential assignment detected for variable '$name'. This is not supported.")
            end

            new_ir = "$type $name = $ir_value;"

            emit(ir_builder, new_ir)

            println("rhs_assgn_node: ", rhs_assgn_node)
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

            ir = "$(create_cpp_var(rhs_assgn_node)) = $name;"

            assigned_param = rhs_assgn_node.cpp_var.index

            if rhs_assgn_node.cpp_var.type == "torch::Tensor" && rhs_assgn_node.cpp_var.index == 1
                # For the Position attribute (index 1), sync spatial .position from tensor value
                rhs_node_idx = rhs_assgn_node.index
                pos_sync = "for (int _i = 0; _i < 3 && _i < $name.numel(); _i++) { rhs[m2[$rhs_node_idx]].position[_i] = $name[_i].template item<double>(); }"
                emit(ir_builder, pos_sync)
            elseif rhs_assgn_node.cpp_var.type != "torch::Tensor" && (assigned_param == 1 || assigned_param == 2 || assigned_param == 3)
                # For scalar attributes at positions 1-3, sync to .position
                ir_node_pos = create_pos_cpp_var(rhs_assgn_node, "rhs", "m2") * " = $name;"
                emit(ir_builder, ir_node_pos)
            end

            emit(ir_builder, ir)
        end
    end

end
