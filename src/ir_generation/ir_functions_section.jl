import ...AstNodes: RuleNode, IdentifierNode, TypeInstanceNode,
            UndirectedTypeEdgeNode, BinaryOpNode, ExpressionNode, ParameterNode,
            WithClauseNode, SolveClauseNode, UnaryOpNode, FunctionNode, GroupNode,
            CallNode, LiteralNode, DefinitionNode, IndexAccessNode, ArrayLiteralNode,
            ModelLoadNode, StringNode


WHITE_SPACE = " "
COMMA = ","
LEFT_PAREN = "("
RIGHT_PAREN = ")"
LEFT_CURLY = "{"
RIGHT_CURLY = "}"
SEMI_COLON = ";"
RETURN = "return"
EQUAL = "="
NEW_LINE = "\n"

# Registry of model-loaded functions: (name, filepath) pairs
MODEL_REGISTRY = []

"""
    ir_array_literal(node::ArrayLiteralNode) -> String

Recursively converts an ArrayLiteralNode into a C++ torch::tensor literal.
"""
function ir_array_literal(node)
    inner = _ir_array_literal_inner(node)
    return "torch::tensor($inner, torch::kFloat64)"
end

function _ir_array_literal_inner(node)
    if node isa ArrayLiteralNode
        parts = [_ir_array_literal_inner(el) for el in node.elements]
        return "{" * join(parts, ", ") * "}"
    elseif node isa LiteralNode
        return string(get_value(node))
    elseif node isa UnaryOpNode
        op = node.expression.position.value
        return op * string(get_value(node.operand))
    elseif node isa IdentifierNode
        return get_value(node)
    elseif node isa BinaryOpNode
        op = node.expression.position.value
        lhs = _ir_array_literal_inner(node.lhs)
        rhs = _ir_array_literal_inner(node.rhs)
        return "$lhs $op $rhs"
    elseif node isa GroupNode
        inner_expr = node.expression
        return "(" * _ir_array_literal_inner(inner_expr) * ")"
    else
        return string(get_value(node))
    end
end

"""
    builtin_to_cpp(func_name, arg_str) -> String

Maps a FoxFlow built-in function name and its already-rendered C++ argument
string to the corresponding C++ expression.  Returns `nothing` when the name
is not a known built-in (caller should fall back to a raw call).
"""
function builtin_to_cpp(func_name, arg_str)
    if func_name == "heaviside"
        return "DGGML::heaviside($arg_str)"
    elseif func_name == "sqrt"
        return "sqrt($arg_str)"
    elseif func_name == "normal_distr"
        return "DGGML::normal_distr()"
    elseif func_name == "uniform_distr"
        return "DGGML::uniform_distr()"
    elseif func_name == "indicator"
        parts = split(arg_str, ", ", limit=2)
        return length(parts) == 2 ? "(($(parts[1])) > 0 ? ($(parts[2])) : 0.0)" : "($arg_str)"
    elseif func_name == "zeros_matrix"
        return "torch::zeros({$arg_str}, torch::kFloat64)"
    elseif func_name == "ones_matrix"
        return "torch::ones({$arg_str}, torch::kFloat64)"
    elseif func_name == "rand_matrix"
        return "torch::rand({$arg_str}, torch::kFloat64)"
    elseif func_name == "eye_matrix"
        return "torch::eye($arg_str, torch::kFloat64)"
    elseif func_name == "mat_dot"
        return "torch::mm($arg_str)"
    elseif func_name == "mat_add"
        parts = split(arg_str, ", ", limit=2)
        return length(parts) == 2 ? "($(parts[1]) + $(parts[2]))" : "($arg_str)"
    elseif func_name == "mat_mul"
        parts = split(arg_str, ", ", limit=2)
        return length(parts) == 2 ? "($(parts[1]) * $(parts[2]))" : "($arg_str)"
    elseif func_name == "transpose"
        return "$arg_str.t()"
    elseif func_name == "einsum"
        return "torch::einsum($arg_str)"
    elseif func_name == "permute"
        parts = split(arg_str, ", ", limit=2)
        return length(parts) == 2 ? "$(parts[1]).permute({$(parts[2])})" : "($arg_str)"
    elseif func_name == "sum"
        return "$arg_str.sum()"
    elseif func_name == "argmax"
        return "$arg_str.argmax()"
    elseif func_name == "softmax"
        parts = split(arg_str, ", ", limit=2)
        return length(parts) == 2 ? "torch::softmax($(parts[1]), $(parts[2]))" : "torch::softmax($arg_str, 0)"
    elseif func_name == "relu"
        return "torch::relu($arg_str)"
    elseif func_name == "forward"
        parts = split(arg_str, ", ", limit=2)
        return length(parts) == 2 ? "$(parts[1]).forward($(parts[2]))" : "$arg_str.forward()"
    else
        return nothing  # unknown — caller emits raw call
    end
end

function traverse_group_node_expr(group_node, variables=Set{String}())
    """
    Traverses a group node and gets the expression inside
    as a string. It should also collect variables used in the expression.
    """

    if group_node.expression isa BinaryOpNode
        bin_node = group_node.expression
        lhs = traverse_group_node_expr(GroupNode(bin_node.lhs), variables)
        rhs = traverse_group_node_expr(GroupNode(bin_node.rhs), variables)
        op = bin_node.expression.position.value

        # Handle powered expression
        if op == "^"
            return "std::pow($lhs, $rhs)"
        else
            return "($lhs $op $rhs)"
        end

    elseif group_node.expression isa UnaryOpNode
        un_node = group_node.expression
        operand = traverse_group_node_expr(GroupNode(un_node.operand), variables)
        op = un_node.expression.position.value
        return "($op$operand)"
    elseif group_node.expression isa CallNode

        call_node = group_node.expression
        func_name = get_value(call_node.function_node)
        args = call_node.args
        arg_str = join(map( (arg) -> begin
            if arg isa LiteralNode
                return string(get_value(arg))
            elseif arg isa GroupNode
                return traverse_group_node_expr(arg, variables)
            elseif arg isa BinaryOpNode
                return traverse_group_node_expr(GroupNode(arg), variables)
            elseif arg isa UnaryOpNode
                return traverse_group_node_expr(GroupNode(arg), variables)
            elseif arg isa IdentifierNode
                var_name = get_value(arg)
                push!(variables, var_name)
                return var_name
            elseif arg isa IndexAccessNode
                return traverse_group_node_expr(GroupNode(arg), variables)
            elseif arg isa ArrayLiteralNode
                return ir_array_literal(arg)
            elseif arg isa CallNode || arg isa BinaryOpNode
                return traverse_group_node_expr(GroupNode(arg), variables)
            elseif arg isa FloatNode || arg isa IntegerNode
                return string(get_value(arg))
            elseif arg isa StringNode
                return "\"$(arg.token.position.value)\""
            else
                return traverse_group_node_expr(GroupNode(arg), variables)
            end
            end, args), ", ")
        cpp = builtin_to_cpp(func_name, arg_str)
        return cpp !== nothing ? cpp : "$func_name($arg_str)"

    elseif group_node.expression isa LiteralNode
        return string(get_value(group_node.expression))

    elseif group_node.expression isa IdentifierNode

        var_name = get_value(group_node.expression)
        push!(variables, var_name)
        return get_value(group_node.expression)

    elseif group_node.expression isa IndexAccessNode
        idx_node = group_node.expression
        obj_str = traverse_group_node_expr(GroupNode(idx_node.object), variables)
        idx_strs = map(idx -> traverse_group_node_expr(GroupNode(idx), variables), idx_node.indices)
        if length(idx_strs) == 1
            return "$(obj_str)[$(idx_strs[1])].template item<double>()"
        else
            return "$(obj_str).index({$(join(idx_strs, ", "))}).template item<double>()"
        end

    elseif group_node.expression isa ArrayLiteralNode
        return ir_array_literal(group_node.expression)

    elseif group_node isa GroupNode
        return traverse_group_node_expr(group_node.expression, variables)

    else
        throw("Unknown group node expression type: $(typeof(group_node.expression))")
    end

end


function ir_expr_builtin_func!(func_name, args, ir_builder, scope)
    """
    Generates the IR for a built-in function.
    """
    
    variables = Set{String}()
    arg_str = join(map( (arg) -> begin
        if arg isa LiteralNode
            return string(get_value(arg))
        elseif arg isa GroupNode
            return traverse_group_node_expr(arg, variables)
        elseif arg isa BinaryOpNode
            return traverse_group_node_expr(GroupNode(arg), variables)
        elseif arg isa UnaryOpNode
            return traverse_group_node_expr(GroupNode(arg), variables)
        elseif arg isa IdentifierNode
            var_name = get_value(arg)
            push!(variables, var_name)
            return var_name
        elseif arg isa LiteralNode
            return string(get_value(arg))
        elseif arg isa IndexAccessNode
            return traverse_group_node_expr(GroupNode(arg), variables)
        elseif arg isa ArrayLiteralNode
            return ir_array_literal(arg)
        elseif arg isa CallNode
            return traverse_group_node_expr(GroupNode(arg), variables)
        elseif arg isa FloatNode || arg isa IntegerNode
            return string(get_value(arg))
        elseif arg isa StringNode
            return "\"$(arg.token.position.value)\""
        else
            return traverse_group_node_expr(GroupNode(arg), variables)
        end
    end, args), ", ")

    # check to see variables are correctly defined
    map( (var) -> begin
            if var in collect(keys(scope))
            else
                throw("Variable $var not found in propensity table for built-in function $func_name.")
            end
        end,
        collect(variables)
    )

    if func_name == "cos"
        ir = "cos($arg_str)"
    elseif func_name == "arccos"
        ir = "acos($arg_str)"
    elseif func_name == "sin"
        ir = "sin($arg_str)"
    elseif func_name == "inverse"
        ir = "DGGML::inverse()"
    elseif func_name == "pow"
        ir = "pow()"
    elseif func_name == "abs"
        ir = "abs($arg_str)"
    elseif func_name == "argmax"
        # returns an integer index
        ir = "$arg_str.argmax().template item<int64_t>()"
    else
        cpp = builtin_to_cpp(func_name, arg_str)
        if cpp !== nothing
            ir = cpp
        else
            throw("Unknown built-in function: $func_name")
        end
    end

    emit(ir_builder, ir)
end

function ir_expr_distribution!(node, ir_builder, scope)
    """
    Generates the IR for a distribution.
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
        println("arg: ", arg)
        exit(0)
        # ir_prop_expr!(arg, arg_ir, propensity_table, arg_ir, propensity)
        push!(parsed_args, join(arg_ir.instructions))
    end

    args = join(parsed_args, ", ")
    if scope["random_device"] == false
        scope["random_device"] = true
        ir0 = "std::random_device random_device;\n"
        emit(ir_builder, ir0)

        ir1 = "std::mt19937 random_engine(random_device());\n"
        emit(ir_builder, ir1)
        scope["random_device"] = true
    end

    func_name = ir_distribution_func(distribution, args)
    ir2 = "$func_name("*args *")(random_engine)"
    return ir2
end

function ir_func_expr!(expr_value, scope, ir_builder)
    """
    Function Expression
    """

    function_node = expr_value
    func_args = nothing
    func_name = nothing

    if function_node isa IdentifierNode

        arg = function_node.token
        if arg.position.value in collect(keys(scope))
            # Loads it into the body
            emit(ir_builder, arg.position.value)
        else
            throw("Error: Propensity variable $(arg.position.value) not found in function scope.")
        end

    elseif function_node isa BinaryOpNode
        lhs = function_node.lhs
        ir_func_expr!(lhs, scope, ir_builder)
        emit(ir_builder, " $(function_node.expression.position.value) ")
        rhs = function_node.rhs
        ir_func_expr!(rhs, scope, ir_builder)

    elseif function_node isa CallNode
        call_args = function_node.args
        function_node = function_node.function_node
        func_name = get_value(function_node)
        println("func_name ---> ", func_name)
        # func_args = function_node.args
        ir_expr_builtin_func!(func_name, call_args,
                        ir_builder, scope)
    elseif function_node isa GroupNode
        # If it is a group node, we need to
        # traverse the expression inside the group
        inner_expr = function_node.expression

        emit(ir_builder, "(")
        ir_func_expr!(inner_expr, scope, ir_builder)
        emit(ir_builder, ")")

    elseif function_node isa LiteralNode
        emit(ir_builder, get_value(function_node))
    elseif function_node isa IndexAccessNode
        # Tensor index access in function expression
        ir_func_expr!(function_node.object, scope, ir_builder)
        indices_ir = []
        for idx in function_node.indices
            idx_builder = IRBuilder([])
            ir_func_expr!(idx, scope, idx_builder)
            push!(indices_ir, build_sameline(idx_builder))
        end
        if length(indices_ir) == 1
            emit(ir_builder, "[" * indices_ir[1] * "].template item<double>()")
        else
            emit(ir_builder, ".index({" * join(indices_ir, ", ") * "}).template item<double>()")
        end

    elseif function_node isa ArrayLiteralNode
        emit(ir_builder, ir_array_literal(function_node))

    elseif function_node isa UnaryOpNode

        operation = function_node.expression

        if operation.position.value == "~"
            ir = ir_expr_distribution!(function_node.operand, ir_builder, scope)
        elseif operation.position.value == "-"
            ir = "-" * get_value(function_node.operand)
        else
            throw("Error: Unary operation $(operation.position.value) not recognized.")
        end

        emit(ir_builder, ir)
    else
        throw("Error: $(function_node) not recognized.")
    end

end


function ir_function!(ast, functions_table, ir_builder_func)
    """
    Handles Single Function Definition
    """

    name = ast.name
    name = get_value(name)

    signature = ast.signature

    functions_table[name] = Dict()

    output_type = convert_type_name(get_value(signature.output.name))

    # Check if this is a model-loaded function
    if ast.model_load !== nothing
        model_path = ast.model_load.filepath.token.position.value

        # Register the model for load_models() generation
        push!(MODEL_REGISTRY, (name, model_path))

        # Mark this function as a model in the functions_table
        functions_table[name]["__is_model__"] = true
        functions_table[name]["__model_path__"] = model_path

        # Emit static torch::jit::Module declaration
        emit(ir_builder_func, "static torch::jit::Module $(name)_module;")
        emit(ir_builder_func, NEW_LINE)

        # Emit wrapper function
        emit(ir_builder_func, output_type)
        emit(ir_builder_func, WHITE_SPACE)
        emit(ir_builder_func, name)

        # Function arguments
        emit(ir_builder_func, LEFT_PAREN)
        arg_names = []
        sig_copy = deepcopy(signature)
        while length(sig_copy.args) > 0
            current_arg = popfirst!(sig_copy.args)
            arg_name = get_value(current_arg.name)
            arg_type = convert_type_name(get_value(current_arg.type.name))

            functions_table[name][arg_name] = arg_type
            push!(arg_names, (arg_name, arg_type))

            emit(ir_builder_func, arg_type)
            emit(ir_builder_func, WHITE_SPACE)
            emit(ir_builder_func, arg_name)

            if length(sig_copy.args) > 0
                emit(ir_builder_func, COMMA)
                emit(ir_builder_func, WHITE_SPACE)
            end
        end

        emit(ir_builder_func, RIGHT_PAREN)
        emit(ir_builder_func, LEFT_CURLY)

        # Build the forward call body
        # Pack inputs into a vector of IValues
        emit(ir_builder_func, "std::vector<torch::jit::IValue> inputs;")

        for (arg_name, arg_type) in arg_names
            if arg_type == "torch::Tensor"
                emit(ir_builder_func, "inputs.push_back($(arg_name));")
            else
                # Wrap scalar as a tensor
                emit(ir_builder_func, "inputs.push_back(torch::tensor($(arg_name)));")
            end
        end

        # Call forward and return
        if output_type == "torch::Tensor"
            emit(ir_builder_func, "return $(name)_module.forward(inputs).toTensor();")
        elseif output_type == "double"
            emit(ir_builder_func, "return $(name)_module.forward(inputs).toTensor().item<double>();")
        elseif output_type == "int"
            emit(ir_builder_func, "return $(name)_module.forward(inputs).toTensor().item<int>();")
        else
            # Fallback: assume tensor output
            emit(ir_builder_func, "return $(name)_module.forward(inputs).toTensor();")
        end

        emit(ir_builder_func, RIGHT_CURLY)

        check = build_sameline(ir_builder_func)
        return check
    end

    # Regular function definition (not model-loaded)
    body = ast.body
    ret = ast.fun_return

    function ir_function_body!(ast, func_scope, ir_builder)
        """
        Handles the function definition body.
        """

        while length(ast.expressions) > 0
            ir_expr_build = IRBuilder([])
            expr = popfirst!(ast.expressions)

            expr_name = get_value(expr.name)
            expr_type = convert_type_name(get_value(expr.type))

            emit(ir_expr_build, expr_type)
            emit(ir_expr_build, WHITE_SPACE)
            emit(ir_expr_build, expr_name)

            expr_value = expr.value
            if !(expr_value isa GroupNode)
                expr_value = GroupNode(expr_value)
            end

            emit(ir_expr_build, WHITE_SPACE)
            emit(ir_expr_build, EQUAL)
            emit(ir_expr_build, WHITE_SPACE)

            # Add function expression to scope
            ir_func_expr!(expr_value, func_scope, ir_expr_build)

            emit(ir_expr_build, SEMI_COLON)
            final_expr = build_sameline(ir_expr_build)
            emit(ir_builder, final_expr)

            func_scope[expr_name] = expr_type
        end

    end

    output_type_reg = convert_type_name(get_value(ast.signature.output.name))

    # Output Type
    emit(ir_builder_func, output_type_reg)
    emit(ir_builder_func, WHITE_SPACE)

    # Function Name
    emit(ir_builder_func, name)

    # Adding function arguments
    emit(ir_builder_func, LEFT_PAREN)
    while length(signature.args) > 0
        current_arg = popfirst!(signature.args)
        arg_name = get_value(current_arg.name)
        arg_type = convert_type_name(get_value(current_arg.type.name))

        functions_table[name][arg_name] = arg_type

        emit(ir_builder_func, arg_type)
        emit(ir_builder_func, WHITE_SPACE)
        emit(ir_builder_func, arg_name)

        if length(signature.args) > 0
            emit(ir_builder_func, COMMA)
            emit(ir_builder_func, WHITE_SPACE)
        end

    end

    emit(ir_builder_func, RIGHT_PAREN)
    emit(ir_builder_func, LEFT_CURLY)

    # Body Here
    # Get the current variables in scope
    func_scope = functions_table[name]
    ir_function_body!(body, func_scope, ir_builder_func)

    # Emit Return Statement
    ret_ir_builder = IRBuilder([])
    emit(ret_ir_builder, RETURN)
    emit(ret_ir_builder, WHITE_SPACE)
    ret_value = ret.value
    ir_func_expr!(ret_value, func_scope, ret_ir_builder)
    emit(ir_builder_func, build_sameline(ret_ir_builder))
    emit(ir_builder_func, SEMI_COLON)

    emit(ir_builder_func, RIGHT_CURLY)

    check = build_sameline(ir_builder_func)
end

function ir_function_list!(ast, functions_table, ir_builder)
    """
    Handles multiple function definitions in a single
    name space.
    """

    ir_total_functions = IRBuilder([])
    while length(ast) > 0
        ir_builder_func = IRBuilder([])
        func_ast = popfirst!(ast)
        ir_builder_func = ir_function!(func_ast, functions_table, ir_builder_func)
        emit(ir_total_functions, ir_builder_func)
        emit(ir_total_functions, SEMI_COLON)
        emit(ir_total_functions, NEW_LINE)
    end

    emit(ir_builder, build(ir_total_functions))
end

function ir_func_sec_hdr!(ir_builder, section_name)
    """
    Adds ir function section header information
    """

    hdr = "#ifndef DGGML_FUNCTIONS_$(section_name)_HPP\n#define DGGML_FUNCTIONS_$(section_name)_HPP\n#include<cmath>\n#include <torch/script.h>\n namespace $section_name {\n"
    emit(ir_builder, hdr)
    return ir_builder
end

function ir_func_sec_foot!(ir_builder)
    footer = "}\n#endif"
    emit(ir_builder, footer)
    return ir_builder
end

function ir_load_models!(ir_builder)
    """
    Generates a load_models() function that initializes all
    torch::jit::Module static variables from their file paths.
    """

    if isempty(MODEL_REGISTRY)
        return
    end

    emit(ir_builder, "void load_models() {")
    for (model_name, model_path) in MODEL_REGISTRY
        emit(ir_builder, "    $(model_name)_module = torch::jit::load(\"$(model_path)\");")
        emit(ir_builder, "    $(model_name)_module.eval();")
    end
    emit(ir_builder, "}")
    emit(ir_builder, NEW_LINE)
end

function ir_function_section(ast, functions_table)
    """
    Handles Function Section
    """

    # Clear model registry for this section
    empty!(MODEL_REGISTRY)

    ir_builder = IRBuilder([])
    section_name = get_value(ast.name)

    ir_func_sec_hdr!(ir_builder, section_name)

    functions_table[section_name] = Dict()

    ir_function_list!(ast.functions, functions_table, ir_builder)

    # Generate load_models() if any model-loaded functions were found
    ir_load_models!(ir_builder)

    # Store namespace metadata if models were registered
    if !isempty(MODEL_REGISTRY)
        functions_table["__func_namespace__"] = section_name
    end

    ir_func_sec_foot!(ir_builder)
    build(ir_builder)
end
