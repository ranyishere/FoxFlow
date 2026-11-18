# Core expression parsing: index/slice elements, factors, terms, primaries, binary operators.

"""
    parse_index_element!(tokens)

Parse a single index element inside [...]. This can be:
  - A plain expression:        tensor[i]
  - A slice with start:stop:    tensor[0:3]
  - A slice with start:stop:step: tensor[1:5:2]
  - A bare colon (select all):  tensor[:]
  - A colon with step (::step):  tensor[::2]

Returns either a normal expression Node or a SliceNode.
"""
function parse_index_element!(tokens)
    # Check for bare colon first: [:] or [::step] or [:stop] or [:stop:step]
    if !isempty(tokens) && lookahead(tokens) isa SingleColonToken
        popfirst!(tokens)  # consume ':'
        start_expr = nothing
        # Check for second colon (::step pattern) or stop expression
        if !isempty(tokens) && lookahead(tokens) isa SingleColonToken
            # ::step
            popfirst!(tokens)  # consume second ':'
            if !isempty(tokens) && !(lookahead(tokens) isa RightSquareBracketToken) && !(lookahead(tokens) isa PunctuationToken)
                step_expr = parse_expression!(tokens)
                return SliceNode(nothing, nothing, step_expr)
            else
                return SliceNode(nothing, nothing, nothing)
            end
        elseif !isempty(tokens) && !(lookahead(tokens) isa RightSquareBracketToken) && !(lookahead(tokens) isa PunctuationToken)
            # :stop or :stop:step
            stop_expr = parse_expression!(tokens)
            if !isempty(tokens) && lookahead(tokens) isa SingleColonToken
                popfirst!(tokens)  # consume ':'
                step_expr = parse_expression!(tokens)
                return SliceNode(nothing, stop_expr, step_expr)
            end
            return SliceNode(nothing, stop_expr, nothing)
        else
            # bare : (select all)
            return SliceNode(nothing, nothing, nothing)
        end
    end

    # Parse the first expression (could be start of a slice or a plain index)
    expr = parse_expression!(tokens)

    # Check if followed by colon -> slice
    if !isempty(tokens) && lookahead(tokens) isa SingleColonToken
        popfirst!(tokens)  # consume ':'
        start_expr = expr
        # Check for second colon right away (start::step)
        if !isempty(tokens) && lookahead(tokens) isa SingleColonToken
            popfirst!(tokens)  # consume second ':'
            if !isempty(tokens) && !(lookahead(tokens) isa RightSquareBracketToken) && !(lookahead(tokens) isa PunctuationToken)
                step_expr = parse_expression!(tokens)
                return SliceNode(start_expr, nothing, step_expr)
            else
                return SliceNode(start_expr, nothing, nothing)
            end
        elseif !isempty(tokens) && !(lookahead(tokens) isa RightSquareBracketToken) && !(lookahead(tokens) isa PunctuationToken)
            # start:stop possibly followed by :step
            stop_expr = parse_expression!(tokens)
            if !isempty(tokens) && lookahead(tokens) isa SingleColonToken
                popfirst!(tokens)  # consume ':'
                step_expr = parse_expression!(tokens)
                return SliceNode(start_expr, stop_expr, step_expr)
            end
            return SliceNode(start_expr, stop_expr, nothing)
        else
            # start: (from start to end)
            return SliceNode(start_expr, nothing, nothing)
        end
    end

    return expr
end

function parse_factor!(tokens)
    """
    Parse Factor
    """

    # Remove EndLineTokens
    while isa(lookahead(tokens), EndLineToken)
        popfirst!(tokens)
    end

    if lookahead(tokens) isa LeftParenthesisToken

        popfirst!(tokens)
        expr = parse_expression!(tokens)

        lookahead_token = lookahead(tokens)
        if lookahead_token isa RightParenthesisToken
            popfirst!(tokens)
            return GroupNode(expr)
        else
            return expr
        end
        # expect_token!(tokens, RightParenthesisToken)

    elseif lookahead(tokens) isa IdentifierToken

        id_token = popfirst!(tokens)
        namespace = nothing

        if (length(tokens) > 0)
            if (lookahead(tokens) isa DoubleColonToken)
                # popping Double Colon
                popfirst!(tokens)
                namespace = id_token
                if lookahead(tokens) isa IdentifierToken
                    id_token = popfirst!(tokens)
                else
                    throw("Error expected an IdentifierToken after '::'")
                end
            end
        end

        if !isempty(tokens) && lookahead(tokens) isa LeftParenthesisToken

            
            popfirst!(tokens)  # consume '('
            args = Node[]

            # Tolerate newlines after '(' in multi-line calls
            while !isempty(tokens) && lookahead(tokens) isa EndLineToken
                popfirst!(tokens)
            end

            while !(lookahead(tokens) isa RightParenthesisToken)
                push!(args, parse_expression!(tokens))

                # Skip newlines before a separator or the closing ')'
                while !isempty(tokens) && lookahead(tokens) isa EndLineToken
                    popfirst!(tokens)
                end

                if lookahead(tokens) isa PunctuationToken
                    popfirst!(tokens)

                    # Skip newlines after a separator
                    while !isempty(tokens) && lookahead(tokens) isa EndLineToken
                        popfirst!(tokens)
                    end
                end

            end
            popfirst!(tokens)  # consume ')'

            # FIXME: Namespace not appearing with variable.
            result = CallNode(IdentifierNode(id_token, namespace), args)
        else
            result = IdentifierNode(id_token)
        end

        # Postfix index access: expr[i] or expr[i, j] or expr[0:3] etc.
        while !isempty(tokens) && lookahead(tokens) isa LeftSquareBracketToken
            popfirst!(tokens)  # consume '['
            indices = Node[]
            push!(indices, parse_index_element!(tokens))
            while !isempty(tokens) && lookahead(tokens) isa PunctuationToken
                popfirst!(tokens)  # consume ','
                push!(indices, parse_index_element!(tokens))
            end
            if isempty(tokens) || !(lookahead(tokens) isa RightSquareBracketToken)
                error("Expected ']' after index expression")
            end
            popfirst!(tokens)  # consume ']'
            result = IndexAccessNode(result, indices)
        end

        return result

    elseif lookahead(tokens) isa StringToken
        return StringNode(popfirst!(tokens))

    elseif lookahead(tokens) isa IntegerToken
        return IntegerNode(popfirst!(tokens))

    elseif lookahead(tokens) isa FloatToken
        return FloatNode(popfirst!(tokens))

    elseif lookahead(tokens) isa RightParenthesisToken
        return popfirst!(tokens)

    elseif lookahead(tokens) isa MinusToken
        minus_token = popfirst!(tokens)  # consume '-'
        right = parse_factor!(tokens)
        return UnaryOpNode(minus_token, right)
    elseif lookahead(tokens) isa PlusToken
        plus_token = popfirst!(tokens)  # consume '+'
        right = parse_factor!(tokens)
        return UnaryOpNode(plus_token, right)
    elseif lookahead(tokens) isa NotToken
        not_token = popfirst!(tokens)  # consume '!'
        right = parse_factor!(tokens)
        return UnaryOpNode(not_token, right)
    elseif lookahead(tokens) isa SampleToken
        sample_token = popfirst!(tokens)  # consume '~'
        right = parse_factor!(tokens)
        return UnaryOpNode(sample_token, right)

    elseif lookahead(tokens) isa LeftSquareBracketToken
        # Array/tensor literal: [expr, expr, ...] or [[expr,...], [expr,...], ...]
        popfirst!(tokens)  # consume '['
        elements = Node[]
        if !(lookahead(tokens) isa RightSquareBracketToken)
            push!(elements, parse_expression!(tokens))
            while !isempty(tokens) && lookahead(tokens) isa PunctuationToken
                popfirst!(tokens)  # consume ','
                push!(elements, parse_expression!(tokens))
            end
        end
        if isempty(tokens) || !(lookahead(tokens) isa RightSquareBracketToken)
            error("Expected ']' after array literal")
        end
        popfirst!(tokens)  # consume ']'
        return ArrayLiteralNode(elements)

    # elseif lookahead(tokens) isa RightAngleBracketToken
        # right_angle_token = popfirst!(tokens)  # consume '>'
        # right = parse_factor!(tokens)
        # return UnaryOpNode(right_angle_token, right)
    else
        println("tokens: $tokens")
        error("Unexpected token in factor: $(lookahead(tokens))")
    end

end


function parse_term!(tokens)
    left = parse_factor!(tokens)

    while !isempty(tokens) && (lookahead(tokens) isa AsteriskToken || lookahead(tokens) isa SlashToken)
        op = popfirst!(tokens)
        right = parse_factor!(tokens)
        left = BinaryOpNode(op, left, right)
    end

    return left
end

# <expression> ::= <term>
                # | <expression> "+" <term> 
                # | <expression> "-" <term>
function parse_primary!(tokens)
    token = popfirst!(tokens)

    if token isa IntegerToken
        return IntegerNode(token)
    elseif token isa FloatToken
        return FloatNode(token)
    elseif token isa IdentifierToken
        if !isempty(tokens) && tokens[1] isa LeftParenthesisToken
            popfirst!(tokens)  # consume '('
            args = Node[]
            while !(tokens[1] isa RightParenthesisToken)
                push!(args, parse_expression!(tokens))
                if tokens[1] isa CommaToken
                    popfirst!(tokens)
                end
            end
            popfirst!(tokens)  # consume ')'
            return CallNode(IdentifierNode(token), args)
        else
            return IdentifierNode(token)
        end
    elseif token isa LeftParenthesisToken
        expr = parse_expression!(tokens)
        expect_token!(tokens, RightParenthesisToken)
        return GroupNode(expr)
    else
        error("Unexpected token in expression: $token")
    end
end

function parse_binary_op!(tokens, min_prec)
    left = parse_primary!(tokens)

    while !isempty(tokens)
        op_token = tokens[1]
        op_str = string(typeof(op_token))
        op_str = replace(op_str, r".*\\.(\\w+)Token" => s"\\1")

        if !haskey(PRECEDENCE, op_str)
            break
        end

        prec = PRECEDENCE[op_str]
        if prec < min_prec
            break
        end

        popfirst!(tokens)  # consume operator
        right = parse_binary_op!(tokens, prec + 1)
        left = BinaryOpNode(op_token, left, right)
    end

    return left
end
