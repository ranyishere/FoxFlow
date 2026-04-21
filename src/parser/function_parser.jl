
function convert_token_to_node(token)
    if token isa FloatToken
        return FloatNode(token)
    elseif token isa IdentifierToken
        return IdentifierNode(token)
    elseif token isa IntegerToken
        return IntegerNode(token)
    else
        throw("Unsupported token type for conversion to node: $token")
    end
end

function get_value_safe(node)
    """Helper to safely get the string value from a node."""
    if node isa IdentifierNode
        return node.token.position.value
    elseif node isa FloatNode
        return node.token.position.value
    elseif node isa IntegerNode
        return node.token.position.value
    else
        return ""
    end
end

function process_arg_token(arg_token)

    if lookahead(arg_token) isa IdentifierToken == false
        throw("Expected argument name identifier, found: $(lookahead(arg_token))")
    end

    arg_name = IdentifierNode(popfirst!(arg_token))

    if (lookahead(arg_token) isa SingleColonToken) == false
        throw("Expected ':' after argument name in function definition got : $(lookahead(arg_token))")
    end

    popfirst!(arg_token)  # Consume ':'

    if ((lookahead(arg_token) isa IdentifierToken) || (lookahead(arg_token) isa FloatToken) ) == false
        throw("Expected argument type identifier after ':' in function definition got : $(lookahead(arg_token))")
    end

    arg_node = nothing
    arg_node = convert_token_to_node(popfirst!(arg_token))

    # Check if this is a FixedList<<...>> type
    param_nodes = ParameterNode(nothing)
    if get_value_safe(arg_node) == "FixedList" && !isempty(arg_token) && (lookahead(arg_token) isa LeftAngleBracketToken)
        popfirst!(arg_token)  # Consume '<<'
        # Collect the inner parameters until we hit '>>'
        inner_params = []
        while !isempty(arg_token) && !(lookahead(arg_token) isa RightAngleBracketToken)
            if (lookahead(arg_token) isa PunctuationToken) || (lookahead(arg_token) isa EndLineToken)
                popfirst!(arg_token)
                continue
            end
            push!(inner_params, convert_token_to_node(popfirst!(arg_token)))
        end
        if !isempty(arg_token) && (lookahead(arg_token) isa RightAngleBracketToken)
            popfirst!(arg_token)  # Consume '>>'
        end
        param_nodes = ParameterNode(inner_params)
    end

    arg_type = TypeClassNode(arg_node, param_nodes)

    return FunctionArgNode(arg_name, arg_type)
end

function parse_function_args!(tokens)

    function_args = []
    while isempty(tokens) == false && (lookahead(tokens) isa RightAngleBracketToken) == false

        if ((lookahead(tokens) isa PunctuationToken)
            || (lookahead(tokens) isa EndLineToken)
            || (lookahead(tokens) isa RightParenthesisToken)
            || (lookahead(tokens) isa LeftParenthesisToken))
            popfirst!(tokens)
            continue
        end
        push!(function_args, process_arg_token(tokens))
       end

    if isempty(tokens)
        throw("Unexpected end of tokens while parsing function arguments, expected '>>'")
    end

    # Placeholder for function arguments parsing logic
    return function_args
end

function parse_function_type_signature!(tokens)

    if (lookahead(tokens) isa  FunctionToken) == false
        throw("Expected keyword 'Function' after ':' in function definition")
    end

    popfirst!(tokens)  # Consume 'Function'

    if (lookahead(tokens) isa LeftAngleBracketToken) == false
        throw("Expected '<<' after 'Function' keyword")
    end
    popfirst!(tokens)  # Consume '<<'

    args = parse_function_args!(tokens)

    if (lookahead(tokens) isa RightAngleBracketToken) == false
        throw("Expected '>>' after 'Function' keyword")
    end
    popfirst!(tokens)  # Consume '>>'

    if (lookahead(tokens) isa RightArrowToken) == false
        throw("Expected '->' after function arguments")
    end
    popfirst!(tokens)  # Consume '->'

    # Should be another type here
    if (lookahead(tokens) isa Union{IdentifierToken, FloatToken}) == false
        throw("Expected return type identifier after '->' got $(lookahead(tokens))")
    end
    ret_type = popfirst!(tokens)  # Consume return type identifier

    ret_type_node = convert_token_to_node(ret_type)

    # Check if return type is FixedList<<...>>
    ret_param_nodes = ParameterNode(nothing)
    if get_value_safe(ret_type_node) == "FixedList" && !isempty(tokens) && (lookahead(tokens) isa LeftAngleBracketToken)
        popfirst!(tokens)  # Consume '<<'
        inner_params = []
        while !isempty(tokens) && !(lookahead(tokens) isa RightAngleBracketToken)
            if (lookahead(tokens) isa PunctuationToken) || (lookahead(tokens) isa EndLineToken)
                popfirst!(tokens)
                continue
            end
            push!(inner_params, convert_token_to_node(popfirst!(tokens)))
        end
        if !isempty(tokens) && (lookahead(tokens) isa RightAngleBracketToken)
            popfirst!(tokens)  # Consume '>>'
        end
        ret_param_nodes = ParameterNode(inner_params)
    end

    ret_type_class = TypeClassNode(ret_type_node, ret_param_nodes)
    return FunctionSignatureNode(args, ret_type_class)
end

function parse_function_body!(tokens)


    if lookahead(tokens) isa LeftBracketToken
        popfirst!(tokens)  # Consume '{'
    else
        throw("Expected '{' at the beginning of function body got : $(lookahead(tokens))")
    end

    body_expressions = []
    while isempty(tokens) == false && (lookahead(tokens) isa ReturnToken) == false
        if lookahead(tokens) isa ReturnToken
            break
        end

        if lookahead(tokens) isa EndLineToken
            popfirst!(tokens)  # Consume endline tokens
            continue
        end

        if lookahead(tokens) isa IdentifierToken == false
            throw("Expected expression name identifier in function body got : $(lookahead(tokens))")
        end

        # Parse the line name and type
        expr_name = IdentifierNode(popfirst!(tokens))

        # Pop :
        if (lookahead(tokens) isa SingleColonToken) == false
            throw("Expected ':' after expression name in function body got : $(lookahead(tokens))")
        end

        popfirst!(tokens)  # Consume ':'

        # Get type
        if ((lookahead(tokens) isa IdentifierToken) || (lookahead(tokens) isa FloatToken) ) == false
            throw("Expected expression type identifier after ':' in function body got : $(lookahead(tokens))")
        end

        expr_type = popfirst!(tokens)
        expr_type = convert_token_to_node(expr_type)

        # Pop :=
        if (lookahead(tokens) isa DefineToken) == false
            throw("Expected ':=' after expression type in function body got : $(lookahead(tokens))")
        end
        popfirst!(tokens)  # Consume ':='

        expression = parse_expression!(tokens)

        func_expr = FunctionDefinitionExpressionNode(expr_name, expr_type, expression)
        push!(body_expressions, func_expr)
    end

    if lookahead(tokens) isa ReturnToken == false
        throw("Expected 'return' token at the end of function body got : $(lookahead(tokens))")
    end

    return FunctionBodyExpressionNode(body_expressions)
end

function parse_function_definition!(tokens)

    if (lookahead(tokens) isa IdentifierToken) == false
        throw("Expected function name identifier, found: $(lookahead(tokens))")
    end

    function_name_token = IdentifierNode(popfirst!(tokens))

    if (lookahead(tokens) isa SingleColonToken) == false
        throw("Expected ':' after function name")
    end

    popfirst!(tokens)  # Consume ':'

    function_signature = parse_function_type_signature!(tokens)
    if (lookahead(tokens) isa DefineToken) == false
        throw("Expected ':=' after function type signature got : $(lookahead(tokens))")
    end
    popfirst!(tokens)  # Consume ':='

    # Check if this is a model load (load("file.pt")) or a regular function body
    if lookahead(tokens) isa LoadFileToken
        # This is a neural network model loaded from a file
        popfirst!(tokens)  # Consume 'load'

        if (lookahead(tokens) isa LeftParenthesisToken) == false
            throw("Expected '(' after 'load' in model function definition")
        end
        popfirst!(tokens)  # Consume '('

        if (lookahead(tokens) isa StringToken) == false
            throw("Expected string file path in load() for model function definition, got: $(lookahead(tokens))")
        end
        model_path = StringNode(popfirst!(tokens))

        if (lookahead(tokens) isa RightParenthesisToken) == false
            throw("Expected ')' after file path in load() for model function definition")
        end
        popfirst!(tokens)  # Consume ')'

        # Consume any trailing endline tokens
        while !isempty(tokens) && lookahead(tokens) isa EndLineToken
            popfirst!(tokens)
        end

        return FunctionDefinitionNode(function_name_token, function_signature, nothing, nothing, ModelLoadNode(model_path))
    end

    function_body = parse_function_body!(tokens)

    # Should expect to see a return token at the end
    if (lookahead(tokens) isa ReturnToken) == false || isempty(tokens)
        throw("Expected 'return' token followed by '}' at the end of functions list")
    end

    popfirst!(tokens)  # Consume 'return'
    ret_val = ReturnNode(parse_expression!(tokens))

    # remove any trailing endline tokens
    while lookahead(tokens) isa EndLineToken
        popfirst!(tokens) # Consume endline tokens
    end

    if lookahead(tokens) isa RightBracketToken
        popfirst!(tokens)  # Consume '}'
    end

    return FunctionDefinitionNode(function_name_token, function_signature, function_body, ret_val, nothing)
end

function parse_functions_list!(tokens)

    functions = []
    while isempty(tokens) == false && (lookahead(tokens) isa RightBracketToken) == false
        while lookahead(tokens) isa EndLineToken
            popfirst!(tokens) # Consume endline tokens
        end
        func = parse_function_definition!(tokens)
        push!(functions, func)
        while lookahead(tokens) isa EndLineToken
            popfirst!(tokens) # Consume endline tokens
        end
    end

    if isempty(tokens)
        throw("Unexpected end of tokens while parsing functions list expected '}'")
    end

    return functions
end


function parse_functions_section!(tokens)

    functions_token = lookahead(tokens)

    if (functions_token isa FunctionSectionToken) == false
        throw("Expected 'functions' token, found: $functions_token")
    end

    popfirst!(tokens)  # Consume 'functions' token

    if (lookahead(tokens) isa IdentifierToken) == false
        throw("Expected identifier after 'functions' token got : $(lookahead(tokens))")
    end

    function_section_name = IdentifierNode(popfirst!(tokens))

    if (lookahead(tokens) isa LeftBracketToken) == false
        throw("Expected '{' after function section name")
    end
    popfirst!(tokens)  # Consume '{'

    user_defined_functions = parse_functions_list!(tokens)

    if (lookahead(tokens) isa RightBracketToken) == false
        throw("Expected '}' at the end of function section")
    end
    popfirst!(tokens)  # Consume '}'

    return FunctionSectionNode(function_section_name, user_defined_functions)
end
