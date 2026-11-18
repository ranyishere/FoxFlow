# Array-literal and index/slice IR generation for rule expressions.
# Included into module IRRuleGeneration (see ../ir_rule_generation.jl).

    """
        ir_array_literal(node::ArrayLiteralNode) -> String

    Recursively converts an ArrayLiteralNode into a C++ torch::tensor literal.
    - 1D: [1, 2, 3]         -> torch::tensor({1, 2, 3}, torch::kFloat64)
    - 2D: [[1,2], [3,4]]    -> torch::tensor({{1, 2}, {3, 4}}, torch::kFloat64)
    """
    # True while emitting a value whose declared type is a FixedList (TFArray<N>):
    # 1-D array literals and ones_matrix/zeros_matrix become stack arrays.
    const ARRAY_CTX = Ref(false)

    # zeros_matrix(N)/ones_matrix(N) as a stack array in array context, else `fallback`.
    function _array_ctx_fill(kind, arg_str, fallback)
        ARRAY_CTX[] || return fallback
        if occursin(r"^\s*\d+\s*$", arg_str)
            return "tf::$(kind)<$(strip(arg_str))>()"
        end
        throw("$(kind)_matrix in a FixedList context needs a literal integer size, " *
              "got '$(arg_str)'. Use TFTensor for tensors with computed shapes.")
    end

    function ir_array_literal(node, propensity_table=nothing, context=nothing)
        if ARRAY_CTX[]
            if any(el -> el isa ArrayLiteralNode, node.elements)
                throw("Nested array literals are not supported for FixedList; " *
                      "use TFTensor for multi-dimensional data.")
            end
            parts = ["static_cast<double>($(_ir_array_literal_inner(el, propensity_table, context)))" for el in node.elements]
            return "TFArray<$(length(parts))>{{" * join(parts, ", ") * "}}"
        end
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

    function ir_slice_expr(slice, ir_builder, propensity_table, context)
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

    function has_slice(indices)
    """
        has_slice(indices) -> Bool
    Returns true if any index in the list is a SliceNode.
    """
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
