import ...AstNodes: RuleNode, IdentifierNode, TypeInstanceNode,
            UndirectedTypeEdgeNode, BinaryOpNode, ExpressionNode, ParameterNode,
            WithClauseNode, SolveClauseNode, UnaryOpNode, FunctionNode, GroupNode,
            CallNode, LiteralNode, DefinitionNode


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
        func_name = get_value(call_node.function_node.name)
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
            else
                return arg.token
            end
            end, args), ", ")
        return "$func_name($arg_str)"

    elseif group_node.expression isa LiteralNode
        return string(get_value(group_node.expression))

    elseif group_node.expression isa IdentifierNode

        var_name = get_value(group_node.expression)
        push!(variables, var_name)
        return get_value(group_node.expression)

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
        else
            return arg.token
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

    if func_name == "heaviside"
        # arg_str = parse_func_args(args)
        ir = "DGGML::$func_name($arg_str)"

    elseif func_name == "sqrt"
        ir = "sqrt($arg_str)"
    elseif func_name == "normal_distr"
        ir = "DGGML::normal_distr()"
    elseif func_name == "uniform_distr"
        ir = "DGGML::uniform_distr()"
    elseif func_name == "cos"
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
    elseif func_name == "indicator"
        ir = ir_indicator_func(arg_str)
    else
        throw("Unknown built-in function: $func_name")
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

    name = ast.name
    signature = ast.signature
    body = ast.body
    ret = ast.fun_return

    functions_table[name] = Dict()

    output_type = convert_type_name(get_value(signature.output.name))

    # Output Type
    emit(ir_builder_func, output_type)
    emit(ir_builder_func, WHITE_SPACE)

    # Function Name
    emit(ir_builder_func, get_value(name))

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

function ir_func_sec_hdr!(ir_builder)
    """
    Adds ir function section header informatoin
    """

    hdr = "#ifndef DGGML_FUNCTIONS_HPP\n#define DGGML_FUNCTIONS_HPP\n#include<cmath>\n namespace FractureNetwork {\n"
    emit(ir_builder, hdr)
    return ir_builder
end

function ir_func_sec_foot!(ir_builder)
    footer = "}\n#endif"
    emit(ir_builder, footer)
    return ir_builder
end

function ir_function_section(ast, functions_table)
    """
    Handles Function Section
    """

    ir_builder = IRBuilder([])
    ir_func_sec_hdr!(ir_builder)

    section_name = get_value(ast.name)

    functions_table[section_name] = Dict()

    ir_function_list!(ast.functions, functions_table, ir_builder)
    # functions_table[]
    ir_func_sec_foot!(ir_builder)
    build(ir_builder)
end
