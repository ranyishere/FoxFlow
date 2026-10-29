# Type-section parsing: type signatures, assignments, updates, declarations, function types.

# Parse symbol parameters `<<param>>` or empty
function parse_symbol_parameters!(tokens)
    params = []

    lookahead_token = lookahead(tokens)
    
    if !isempty(tokens) && isa(lookahead(tokens), LeftAngleBracketToken)

        popfirst!(tokens) # Consume `<<`

        while !isempty(tokens) && !isa(lookahead(tokens), RightAngleBracketToken)

            # Handling nested parameters
            if isa(lookahead(tokens), LeftAngleBracketToken)
                nested_params = parse_symbol_parameters!(tokens)
                push!(params, nested_params)  # Add nested params to params array

            elseif isa(lookahead(tokens), LeftParenthesisToken) || isa(lookahead(tokens), RightParenthesisToken)
                popfirst!(tokens)
            elseif isa(lookahead(tokens), EndLineToken)
                # Pop EndLineToken and ignore
                popfirst!(tokens)

            elseif isa(lookahead(tokens), IdentifierToken) && length(tokens) > 1 && isa(tokens[2], SingleColonToken)
                # Handle cases like `param: Type` or `alias: symbol::attribute`
                param_name = parse_symbol_name!(tokens)
                popfirst!(tokens) # Consume `:`

                if isa(lookahead(tokens), IdentifierToken) && length(tokens) > 1 && isa(tokens[2], DoubleColonToken)
                    # Vertex attribute binding: `alias : symbol::attribute`
                    namespace_token = popfirst!(tokens) # symbol
                    popfirst!(tokens)                   # Consume `::`
                    if !isa(lookahead(tokens), IdentifierToken)
                        throw("Expected an identifier after `::` in attribute binding")
                    end
                    attr_token = popfirst!(tokens)      # attribute
                    param_value = IdentifierNode(attr_token, namespace_token)
                    push!(params, NamedParameterNode(param_name, param_value))
                else
                    param_type = parse_type_signature_list!(tokens)

                    # println("param_name: $param_name, param_type: $param_type")

                    push!(params, NamedParameterNode(param_name, param_type))
                    # push!(params, (param_name, param_type))  # Add parameter name and type as a tuple to params array
                end

            else

                # It could jsut be a regular identifer token. Is this valid?
                cur_expr = parse_expression!(tokens)
                push!(params, cur_expr)  # Add parsed expression to params array

            end

            if !isempty(tokens) && isa(lookahead(tokens), PunctuationToken)
                popfirst!(tokens) # Consume `,`
            end

        end

        if isempty(tokens) || !isa(popfirst!(tokens), RightAngleBracketToken)
            throw("Expected `>>` to close parameter list")
        end

    end
    return ParameterNode(params)
end

function parse_type_symbol_parameters!(tokens)
    params = []

    lookahead_token = lookahead(tokens)
    
    if !isempty(tokens) && isa(lookahead(tokens), LeftAngleBracketToken)

        popfirst!(tokens) # Consume `<<`

        while !isempty(tokens) && !isa(lookahead(tokens), RightAngleBracketToken)

            if !(lookahead(tokens) isa IdentifierToken)
                throw("Expected IdentifierToken in type parameter list, got $(lookahead(tokens))")
            end

            param_name = parse_symbol_name!(tokens)
            if !(lookahead(tokens) isa SingleColonToken)
                throw("Expected SingleColonToken after parameter name in type parameter list, got $(lookahead(tokens))")
            end
            popfirst!(tokens) # Consume `:`

            # Handling nested parameters
            if isa(lookahead(tokens), LeftAngleBracketToken)
                nested_params = parse_symbol_parameters!(tokens)
                push!(params, nested_params)  # Add nested params to params array

            elseif isa(lookahead(tokens), LeftParenthesisToken) || isa(lookahead(tokens), RightParenthesisToken)
                popfirst!(tokens)
            elseif isa(lookahead(tokens), EndLineToken)
                # Pop EndLineToken and ignore
                popfirst!(tokens)
            else
                cur_expr = parse_expression!(tokens)
                push!(params, cur_expr)  # Add parsed expression to params array
            end

            if !isempty(tokens) && isa(lookahead(tokens), PunctuationToken)
                popfirst!(tokens) # Consume `,`
            end

        end

        if isempty(tokens) || !isa(popfirst!(tokens), RightAngleBracketToken)
            throw("Expected `>>` to close parameter list")
        end

    end
    return ParameterNode(params)
end


function parse_type_signature_list!(tokens)

    typename = nothing
    if isa(lookahead(tokens), IdentifierToken)
        type_name = parse_symbol_name!(tokens)
    elseif typeof(lookahead(tokens)) in [IntegerToken, FloatToken]
        type_name = IdentifierNode(popfirst!(tokens))
    elseif isa(lookahead(tokens), TypeToken)
        type_name = IdentifierNode(popfirst!(tokens))
    else
        cur_token = popfirst!(tokens)
        throw("Error couldnt determine type got: $cur_token")
    end

    if !isempty(tokens) && isa(lookahead(tokens), RightArrowToken)
        popfirst!(tokens)  # Consume `->`
        return TypeClassNode(type_name, parse_type_signature_list!(tokens))
    elseif !isempty(tokens) && isa(lookahead(tokens), LeftAngleBracketToken)
        # popfirst!(tokens)  # Consume `->`
        param_node = parse_symbol_parameters!(tokens)
        return TypeClassNode(type_name, param_node)
    else
        return TypeClassNode(type_name, ParameterNode([]))
    end
end

function parse_type_assignment!(tokens)
    """
    Parse Type Assignment
    """

    symbol_name = parse_symbol_name!(tokens)
    symbol_parameters = parse_symbol_parameters!(tokens)

    if !isa(popfirst!(tokens), SingleColonToken)
        throw("Expected `:` after symbol name got $(lookahead(tokens))")
    end

    # Parse Type Signature
    type_signature_list = parse_type_signature_list!(tokens)

    # println("type_signature_list: $type_signature_list")

    if !isempty(tokens) && isa(lookahead(tokens), DefineToken)

        popfirst!(tokens)  # Consume `=`

        if !isempty(tokens) && isa(lookahead(tokens), LeftBracketToken)
            popfirst!(tokens) # Consume `{`
            type_declarations = parse_type_declaration_list!(tokens)

            if isempty(tokens) || !isa(popfirst!(tokens), RightBracketToken)
                throw("Expected `}` to close type declaration list")
            end

            return TypeInstanceNode(symbol_name, symbol_parameters, type_signature_list, type_declarations)
        else

            # Assume it's a literal/expression
            tmp = parse_expression!(tokens)
            # literal = popfirst!(tokens)  
            #
            return TypeInstanceNode(symbol_name, symbol_parameters, type_signature_list, tmp)
        end
    end


    return TypeInstanceNode(symbol_name,
                            symbol_parameters,
                            type_signature_list,
                            nothing)

end

function parse_type!(tokens)
    """
    Parse Type
    """

    symbol_name = parse_symbol_name!(tokens)

    if !isa(popfirst!(tokens), SingleColonToken)
        throw("Expected `:` in type declaration")
    end

    type_signature_list = parse_type_signature_list!(tokens)

    return TypeInstanceNode(symbol_name, nothing, type_signature_list, nothing)
end

function collect_statement_tokens!(tokens)
    # Collect the tokens that make up a single statement's right-hand-side
    # expression.  A statement normally ends at a newline, but an expression
    # may span several physical lines when it is wrapped in unbalanced
    # parentheses or square brackets (e.g. a multi-line function call).  We
    # therefore only treat a newline (or a `}` closing the where body) as the
    # end of the statement when we are at bracket depth 0.  Interior newlines
    # that appear inside brackets are dropped so the expression parser sees a
    # clean token stream.
    #
    # In addition, a depth-0 newline does NOT end the statement when the
    # expression is syntactically incomplete across the line break, i.e. either
    # the last collected token is a binary operator (trailing-operator
    # continuation, `... a +` ⏎ `b ...`) or the next non-newline token is a
    # binary operator (leading-operator continuation, `... a` ⏎ `+ b ...`).
    # A new where-body statement always begins with an identifier, never with
    # an operator, so this is unambiguous.
    expression = Any[]
    depth = 0
    while !isempty(tokens)
        tok = lookahead(tokens)
        if depth == 0 && isa(tok, EndLineToken)
            # Trailing-operator continuation: expression is unfinished.
            if !isempty(expression) && isa(expression[end], OperatorToken)
                popfirst!(tokens)  # consume newline, keep collecting
                continue
            end
            # Leading-operator continuation: look past blank lines for the next
            # meaningful token; if it is a binary operator, keep collecting.
            j = 1
            while j <= length(tokens) && isa(tokens[j], EndLineToken)
                j += 1
            end
            if j <= length(tokens) && isa(tokens[j], OperatorToken)
                popfirst!(tokens)  # consume this newline, keep collecting
                continue
            end
            break
        end
        if depth == 0 && isa(tok, RightBracketToken)
            break
        end
        tok = popfirst!(tokens)
        if isa(tok, LeftParenthesisToken) || isa(tok, LeftSquareBracketToken)
            depth += 1
        elseif isa(tok, RightParenthesisToken) || isa(tok, RightSquareBracketToken)
            depth -= 1
        end
        if isa(tok, EndLineToken)
            # Interior newline inside brackets — skip it.
            continue
        end
        push!(expression, tok)
    end
    return expression
end

function parse_type_update!(tokens)
    # Updates to existing types done by a rule or sees a definition of a tmp variable

    # ── Bare call statement: e.g. backward(loss) ──
    # If the next two tokens are Identifier then '(', parse it as a bare
    # expression statement (no assignment, no type annotation).  This lets
    # users write side-effecting calls like `backward(loss)` directly in a
    # where body without needing `_ = ...`.
    if !isempty(tokens) && lookahead(tokens) isa IdentifierToken
        saved = copy(tokens)
        peek_name = popfirst!(tokens)  # consume identifier
        if !isempty(tokens) && lookahead(tokens) isa LeftParenthesisToken
            # It's a bare call — put the name back and parse as expression
            prepend!(tokens, [peek_name])
            expr = parse_expression!(tokens)
            return expr   # CallNode — ir_where_assignment! will handle it
        else
            # Not a bare call — restore and fall through to normal parsing
            prepend!(tokens, [peek_name])
        end
    end

    symbol_name = parse_symbol_name!(tokens)

    # Check for indexed LHS: name[idx] = expr
    lhs_node = symbol_name
    if !isempty(tokens) && lookahead(tokens) isa LeftSquareBracketToken
        while !isempty(tokens) && lookahead(tokens) isa LeftSquareBracketToken
            popfirst!(tokens)  # consume '['
            indices = Node[]
            push!(indices, parse_index_element!(tokens))
            while !isempty(tokens) && lookahead(tokens) isa PunctuationToken
                popfirst!(tokens)  # consume ','
                push!(indices, parse_index_element!(tokens))
            end
            if isempty(tokens) || !(lookahead(tokens) isa RightSquareBracketToken)
                error("Expected ']' after index expression in LHS")
            end
            popfirst!(tokens)  # consume ']'
            lhs_node = IndexAccessNode(lhs_node, indices)
        end
    end

    symbol_parameters = parse_symbol_parameters!(tokens)
    type_signature_list = nothing

    # Its a definition node
    if lookahead(tokens) isa SingleColonToken
        popfirst!(tokens)  # Consume `:`

        # Local function definition: name : Function <<params>> -> RetType := { ... }
        if !isempty(tokens) && lookahead(tokens) isa FunctionToken
            return parse_function_definition_named!(tokens, symbol_name)
        end

        type_signature_list = parse_type_signature_list!(tokens)

        # This is an intermediate value and not a type update
        if isa(lookahead(tokens), DefineToken)
            popfirst!(tokens)  # Consume `:=`
            expression = collect_statement_tokens!(tokens)

            expression_nodes = parse_expression!(expression)
            return DefinitionNode(symbol_name, type_signature_list, expression_nodes)
        end


    end
    
    if !isempty(tokens) && isa(lookahead(tokens), EqualToken)

        popfirst!(tokens)  # Consume `=`

        if !isempty(tokens) && isa(lookahead(tokens), LeftBracketToken)

            popfirst!(tokens) # Consume `{`
            type_declarations = parse_type_declaration_list!(tokens)

            if isempty(tokens) || !isa(popfirst!(tokens), RightBracketToken)
                throw("Expected `}` to close type declaration list")
            end


            # return TypeInstanceNode(symbol_name, symbol_parameters, type_signature_list, type_declarations)
            return TypeInstanceUpdateNode(lhs_node, type_declarations)
        else

            expression = collect_statement_tokens!(tokens)

            expression_nodes = parse_expression!(expression)
            # return TypeInstanceNode(symbol_name, symbol_parameters, type_signature_list, expression_nodes)
            return TypeInstanceUpdateNode(lhs_node, expression_nodes)

        end
    end

    return TypeInstanceNode(symbol_name, symbol_parameters, type_signature_list, nothing)
end


# Parse a type declaration (non-assignment case)
function parse_type_declaration!(tokens)

    # Nucleator
    type_name = parse_symbol_name!(tokens)

    # << param >>
    symbol_parameters = parse_symbol_parameters!(tokens)

    if !isa(lookahead(tokens), SingleColonToken)
        throw("Expected `:` in type declaration got $(lookahead(tokens))")
    else
        popfirst!(tokens)  # Consume `:`
    end

    type_signature_list = parse_type_signature_list!(tokens)

    return TypeInstanceNode(type_name, symbol_parameters, type_signature_list, nothing)
end

# Parse a list of type declarations recursively
function parse_type_declaration_list!(tokens)
    declarations = []

    while !isempty(tokens) && isa(lookahead(tokens), IdentifierToken)
        push!(declarations, parse_type_declaration!(tokens))

        if !isempty(tokens) && isa(lookahead(tokens), EndLineToken)
            popfirst!(tokens)  # Consume newline
        end
    end

    return declarations
end

function parse_types_section!(tokens)

    if !isa(lookahead(tokens), TypeSectionToken)
        throw("Expected `types` section")
    end

    popfirst!(tokens)  # Consume `types`

    section_name = parse_symbol_name!(tokens)

    if !isa(popfirst!(tokens), LeftBracketToken)
        throw("Expected `{` to open type section")
    end

    if isa(lookahead(tokens), EndLineToken)
        popfirst!(tokens)
    end
    type_declarations = parse_type_declaration_list!(tokens)

    if !isa(popfirst!(tokens), RightBracketToken)
        throw("Expected `}` to close type section")
    end

    return TypeSectionNode(section_name, type_declarations)
end

function parse_function_type!(tokens)

    function_name = parse_symbol_name!(tokens)

    cur_token = popfirst!(tokens)
    if !isa(cur_token, LeftParenthesisToken)
        throw("Expected LeftParenthesisToken got: $cur_token")
    end

    function_args = []
    while !isa(lookahead(tokens), RightParenthesisToken)

        # TODO: Have this handle expression
        cur_arg = popfirst!(tokens)
        # cur_arg = parse_symbol_name!(tokens)
        push!(function_args, cur_arg)
        if isa(lookahead(tokens), CommaToken)
            popfirst!(tokens)
        end
    end

    cur_token = popfirst!(tokens)
    if !isa(cur_token, RightParenthesisToken)
        throw("Expected RightParenthesisToken got: $cur_token")
    end

    FunctionNode(function_name, function_args)
end

# function parse_type_assignment_list!(tokens)
#     """
#     Parse Type Assignment List
#     """
# 
#     type_assignments = []
# 
#     while isa(lookahead(tokens), EndLineToken)
#         popfirst!(tokens)
#     end
# 
#     cur_type_assignment = parse_type_assignment!(tokens)
# 
#     push!(type_assignments, cur_type_assignment)
# 
#     if isa(lookahead(tokens), EndLineToken)
#         tmp = parse_type_assignment_list!(tokens)
#         type_assignments = [type_assignments;tmp]
#     end
# 
#     return type_assignments
# end

function parse_type_update_list!(tokens)

    type_updates = []

    while isa(lookahead(tokens), EndLineToken)
        popfirst!(tokens)
    end

    if isa(lookahead(tokens), RightBracketToken)
        type_updates
    else

        while isa(lookahead(tokens), EndLineToken)
            popfirst!(tokens)
        end

        cur_type_update = parse_type_update!(tokens)

        push!(type_updates, cur_type_update)

        if isa(lookahead(tokens), EndLineToken)
            popfirst!(tokens)
            tmp = parse_type_update_list!(tokens)
            type_updates = [type_updates;tmp]
        end

        return type_updates
    end
end
