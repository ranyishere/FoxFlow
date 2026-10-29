# Rule-section parsing: where/with/solve clauses, binding variables, rules and the rules section.

function parse_where_clause!(tokens)
    """
    Parse Where Clause
    """


    WhereClauseNode(
        parse_type_update_list!(tokens)
   )
end

function parse_with_clause!(tokens)
    """
    Parse With Clause
    """

    # cur_token = popfirst!(tokens)
    if !isa(lookahead(tokens), LeftParenthesisToken)
        println("Expected LeftParenthesis got $cur_token")
    end

    propensity = parse_expression!(tokens)

    # remove endline tokens
    while isa(lookahead(tokens), EndLineToken)
        popfirst!(tokens)
    end


    cur_token = popfirst!(tokens)

    if !isa(cur_token, WhereToken)
        throw("Expected WhereToken got $cur_token")
    end

    cur_token = popfirst!(tokens)

    if !isa(cur_token, LeftBracketToken)
        throw("Expected LeftBracketToken got $cur_token")
    end

    where_clause = parse_where_clause!(tokens)

    while isa(lookahead(tokens), EndLineToken)
        popfirst!(tokens)
    end

    cur_token = popfirst!(tokens)
    if !isa(cur_token, RightBracketToken)
        throw("Expected RightBracketToken got $cur_token")
    end

    # if !isa(cur_token, RightBracketToken)
        # throw("Expected RightBracketToken got $cur_token")
    # end

    WithClauseNode(propensity, where_clause)
end

function parse_every_clause!(tokens)
    """
    Parse Every Clause (the `every` identifier is already consumed):

        every (n_r, p_r, T_clock) [additive] where { ... }
    """

    if !isa(lookahead(tokens), LeftParenthesisToken)
        throw("Expected LeftParenthesis after `every` got $(lookahead(tokens))")
    end
    popfirst!(tokens)

    args = Node[]
    while !isa(lookahead(tokens), RightParenthesisToken)
        if isa(lookahead(tokens), PunctuationToken) || isa(lookahead(tokens), CommaToken)
            popfirst!(tokens)
        else
            push!(args, parse_expression!(tokens))
        end
    end
    popfirst!(tokens)  # `)`

    if isempty(args) || length(args) > 3
        throw("`every` takes 1 to 3 arguments (n_r, p_r, T_clock), got $(length(args))")
    end

    n_r = args[1]
    p_r = length(args) >= 2 ? args[2] : nothing
    t_clock = length(args) >= 3 ? args[3] : nothing

    additive = false
    if isa(lookahead(tokens), IdentifierToken) && lookahead(tokens).position.value == "additive"
        popfirst!(tokens)
        additive = true
    end

    while isa(lookahead(tokens), EndLineToken)
        popfirst!(tokens)
    end

    cur_token = popfirst!(tokens)
    if !isa(cur_token, WhereToken)
        throw("Expected WhereToken got $cur_token")
    end

    cur_token = popfirst!(tokens)
    if !isa(cur_token, LeftBracketToken)
        throw("Expected LeftBracketToken got $cur_token")
    end

    where_clause = parse_where_clause!(tokens)

    while isa(lookahead(tokens), EndLineToken)
        popfirst!(tokens)
    end

    cur_token = popfirst!(tokens)
    if !isa(cur_token, RightBracketToken)
        throw("Expected RightBracketToken got $cur_token")
    end

    EveryClauseNode(n_r, p_r, t_clock, additive, where_clause)
end

function parse_solving_content(tokens)
    """
    Parse Solving Content
    """

    content = []

    while !isempty(tokens) && !isa(lookahead(tokens), EndLineToken)
        if isa(lookahead(tokens), EndLineToken)
            popfirst!(tokens)
        elseif isa(lookahead(tokens), LeftParenthesisToken) || isa(lookahead(tokens), RightParenthesisToken)
            popfirst!(tokens)
        else
            push!(content, parse_expression!(tokens))
        end
    end

    return content
end

function parse_binding_variable!(tokens)
    """
    Parse Binding Variable
    """

    # Should be identifier name, possibly with index access like dpos[0]
    var_name = popfirst!(tokens)
    if !isa(var_name, IdentifierToken)
        throw("Expected IdentifierToken got $var_name")
    end

    # Check for index access: dpos[0]
    var_node = IdentifierNode(var_name)
    if !isempty(tokens) && lookahead(tokens) isa LeftSquareBracketToken
        popfirst!(tokens)  # consume '['
        indices = Node[]
        push!(indices, parse_expression!(tokens))
        while !isempty(tokens) && lookahead(tokens) isa PunctuationToken
            popfirst!(tokens)  # consume ','
            push!(indices, parse_expression!(tokens))
        end
        if isempty(tokens) || !(lookahead(tokens) isa RightSquareBracketToken)
            error("Expected ']' after index in binding variable")
        end
        popfirst!(tokens)  # consume ']'
        var_node = IndexAccessNode(var_node, indices)
    end

    # Should pop := define symbol
    cur_token = popfirst!(tokens)
    if !(cur_token isa DefineToken)
        throw("Expected DefineToken got $cur_token")
    end

    D_token = popfirst!(tokens)
    if D_token.position.value != "D"
        throw("Expected Derivative Operator D(..) got $D_token")
    end

    left_paren = popfirst!(tokens)
    if !isa(left_paren, LeftParenthesisToken)
        throw("Expected LeftParenthesisToken got $left_paren")
    end

    ode_variables = []
    while !(lookahead(tokens) isa RightParenthesisToken)

        if lookahead(tokens) isa PunctuationToken
            popfirst!(tokens)
        else
            ode_variable = parse_symbol_name!(tokens)
            # Check for index access on ODE variable: D(im_pos[0], t)
            if !isempty(tokens) && lookahead(tokens) isa LeftSquareBracketToken
                popfirst!(tokens)  # consume '['
                indices = Node[]
                push!(indices, parse_expression!(tokens))
                while !isempty(tokens) && lookahead(tokens) isa PunctuationToken
                    popfirst!(tokens)  # consume ','
                    push!(indices, parse_expression!(tokens))
                end
                if isempty(tokens) || !(lookahead(tokens) isa RightSquareBracketToken)
                    error("Expected ']' after index in ODE variable")
                end
                popfirst!(tokens)  # consume ']'
                ode_variable = IndexAccessNode(ode_variable, indices)
            end
            push!(ode_variables, ode_variable)
        end

    end

    right_paren = popfirst!(tokens)
    if !isa(right_paren, RightParenthesisToken)
        throw("Expected RightParenthesisToken got $right_paren")
    end

    BindingVariableNode(var_node, ode_variables)
end

# TODO: Add support for DefinitionNode
function parse_solve_clause!(tokens)
    """
    Parses Solve Clause and returns an ODENode.
    Handles both simple names (dx : ODE = expr) and
    indexed names (dx[0] : ODE = expr) for tensor component ODEs.
    """

    var_name = popfirst!(tokens)

    if !(var_name isa IdentifierToken)
        throw("Expected Identifier token got $var_name")
    end

    var_name = IdentifierNode(var_name)

    # Check for index access: dx[i] : ODE = expr or dx[0:3] : ODE = expr
    if !isempty(tokens) && lookahead(tokens) isa LeftSquareBracketToken
        popfirst!(tokens)  # consume '['
        indices = Node[]
        push!(indices, parse_index_element!(tokens))
        while !isempty(tokens) && lookahead(tokens) isa PunctuationToken
            popfirst!(tokens)  # consume ','
            push!(indices, parse_index_element!(tokens))
        end
        if isempty(tokens) || !(lookahead(tokens) isa RightSquareBracketToken)
            error("Expected ']' after index expression in ODE name")
        end
        popfirst!(tokens)  # consume ']'
        var_name = IndexAccessNode(var_name, indices)
    end

    colon_token = popfirst!(tokens)
    if !(colon_token isa SingleColonToken)
        throw("Error expected single colon token")
    end

    # Should be an ode token (Type)
    ode_token = popfirst!(tokens)
    if !((ode_token isa ODEToken) || (ode_token isa FloatToken) || (ode_token isa IntegerToken))
        throw("Error expected an ODE type or Float or Integer. Got: $ode_token.")
    end

    eq_token = popfirst!(tokens)
    if !(eq_token isa EqualToken || eq_token isa DefineToken)
        throw("Error expected '=' or ':=' got $eq_token")
    end

    expression = parse_expression!(tokens)

    if eq_token isa EqualToken
        return ODENode(var_name, expression)
    elseif eq_token isa DefineToken
        return DefinitionNode(var_name, TypeClassNode(IdentifierNode(ode_token), ParameterNode([])), expression)
    else
        throw("Unexpected token in solve clause: $eq_token")
    end
end

function parse_rule_solve!(tokens)
    """
    Parses the Binding Variables
    """

    # Pops left parenthesis
    cur_token = popfirst!(tokens)
    if !isa(cur_token, LeftParenthesisToken)
        throw("Expected LeftParenthesisToken got $cur_token")
    end

    total_binding_variables = []
    while !(lookahead(tokens) isa RightParenthesisToken)

        if (isa(lookahead(tokens), EndLineToken) 
            || isa(lookahead(tokens), PunctuationToken))
            popfirst!(tokens)
        else
            push!(total_binding_variables, parse_binding_variable!(tokens))
        end

    end

    # Pops right parenthesis
    cur_token = popfirst!(tokens)
    if !isa(cur_token, RightParenthesisToken)
        throw("Expected RightParenthesisToken got $cur_token")
    end

    # Need to parse the  { ... } at the end

    # { token
    cur_token = popfirst!(tokens)
    if !isa(cur_token, LeftBracketToken)
        throw("Expected LeftBracketToken got $cur_token")
    end

    solve_clause = []
    while !isempty(tokens) && !isa(lookahead(tokens), RightBracketToken)

        if isa(lookahead(tokens), EndLineToken)
            popfirst!(tokens)
        elseif isa(lookahead(tokens), LeftParenthesisToken) || isa(lookahead(tokens), RightParenthesisToken)
            popfirst!(tokens)
            # Check its a definition token
        else
            ode_node = parse_solve_clause!(tokens)
            # push!(solve_clause, parse_expression!(tokens))
            push!(solve_clause, ode_node)
        end
    end

    cur_token = popfirst!(tokens)
    if !isa(cur_token, RightBracketToken)
        throw("Expected RightBracketToken got $cur_token")
    end

    SolveClauseNode(
        total_binding_variables,
        solve_clause
    )
end

function parse_rule!(tokens)
    """
    Parse Rule
    """

    rule_name = parse_symbol_name!(tokens)
    cur_token = popfirst!(tokens)

    if !isa(cur_token, DefineToken)
        throw("Expected assignment Token ':=' got $cur_token")
    end

    lhs = []
    while !isempty(tokens) && !isa(lookahead(tokens), LeftAngleBracketToken) && !isa(lookahead(tokens), RightArrowToken)

        if isa(lookahead(tokens), EndLineToken)
            popfirst!(tokens)

        elseif isa(lookahead(tokens), LeftParenthesisToken) || isa(lookahead(tokens), RightParenthesisToken)
            popfirst!(tokens)

        # TODO: Handle directed edges
        elseif isa(lookahead(tokens), EdgeToken)

            vert_0 = pop!(lhs)

            if vert_0 isa UndirectedTypeEdgeNode
                # If the last node was an edge node, we need to pop it
                # tmp = vert_0.left_vertex
                tmp = vert_0.right_vertex
                push!(lhs, vert_0)
                vert_0 = tmp
            end

            # Popping edge token
            popfirst!(tokens)

            # Popping left parenthesis token
            popfirst!(tokens)

            vert_1 = parse_type_assignment!(tokens)

            # Popping Right parenthesis token
            popfirst!(tokens)

            edge_node = UndirectedTypeEdgeNode(vert_0, vert_1)
            push!(lhs, edge_node)

        else
            tmp = parse_type_assignment!(tokens)
            push!(lhs, tmp)
            # push!(lhs, popfirst!(tokens))
        end
    end

    left_rule_parameter_node = ParameterNode([])
    if isa(lookahead(tokens), LeftAngleBracketToken)
        left_rule_parameter_node = parse_symbol_parameters!(tokens)
    end

    # println("left rule parameter node: $left_rule_parameter_node")

    # Gets rid of new lines
    while isa(lookahead(tokens), EndLineToken)
        popfirst!(tokens)
    end

    cur_token = popfirst!(tokens)

    if isempty(tokens) || !isa(cur_token, RightArrowToken)
        throw("Expected `->` in rule definition got $cur_token")
    end

    # Parsing RHS
    rhs = []
    while !isempty(tokens) && !isa(lookahead(tokens), LeftAngleBracketToken)
        if isa(lookahead(tokens), EndLineToken)
            popfirst!(tokens)
        elseif isa(lookahead(tokens), LeftParenthesisToken) || isa(lookahead(tokens), RightParenthesisToken)
            popfirst!(tokens)

        # TODO: Handle directed edges
        elseif isa(lookahead(tokens), EdgeToken)

            # TODO: Handle things like Node -- Node -- Node
            vert_0 = pop!(rhs)

            if vert_0 isa UndirectedTypeEdgeNode
                # If the last node was an edge node, we need to pop it
                tmp = vert_0.right_vertex
                push!(rhs, vert_0)
                vert_0 = tmp
            end

            # Popping edge token
            popfirst!(tokens)

            # Popping left parenthesis token
            popfirst!(tokens)

            vert_1 = parse_type_assignment!(tokens)

            # Popping Right parenthesis token
            popfirst!(tokens)

            edge_node = UndirectedTypeEdgeNode(vert_0, vert_1)
            push!(rhs, edge_node)

        else
            tmp = parse_type_assignment!(tokens)
            push!(rhs, tmp)
            # push!(lhs, popfirst!(tokens))
        end
        # push!(rhs, popfirst!(tokens))
    end

    right_rule_parameter_node = ParameterNode([])
    if isa(lookahead(tokens), LeftAngleBracketToken)
        right_rule_parameter_node = parse_symbol_parameters!(tokens)
    end

    # TODO: check parameters
    modify_clause = nothing

    # Removes Endlines Before Solving
    while(isa(lookahead(tokens), EndLineToken))
        popfirst!(tokens)  # Consume EndLineToken
    end

    if (!isempty(tokens) 
        && (isa(lookahead(tokens), WithToken) 
            || isa(lookahead(tokens), SolvingToken))
       )

        # Consume `with` or `solving`
        clause_token = popfirst!(tokens) 
        clause_content = []
        
        if isa(clause_token, WithToken)

            with_clause = parse_with_clause!(tokens)

            modify_clause = with_clause
            # modify_clause = WithClauseNode(clause_token, clause_content)
        else

            modify_clause = parse_rule_solve!(tokens)
        end
    elseif (length(tokens) > 1
            && isa(tokens[1], IdentifierToken) && tokens[1].position.value == "every"
            && isa(tokens[2], LeftParenthesisToken))

        # `every` is a contextual keyword, so it is not reserved by the lexer.
        popfirst!(tokens)
        modify_clause = parse_every_clause!(tokens)
    end


    # Find where token
    return RuleNode(
                    rule_name, 
                    lhs, left_rule_parameter_node,
                    rhs, right_rule_parameter_node,
                    modify_clause
                   )
end

function parse_rules_list!(tokens)
    """
    parse rules list
    """

    rules = []

    # Remove extra lines
    while isa(lookahead(tokens),EndLineToken)
        popfirst!(tokens)
    end

    # Parse rule
    rule = parse_rule!(tokens)
    push!(rules, rule)

    # Remove nedline tokens and return the rules
    # if there are no more tokens
    while isa(lookahead(tokens),EndLineToken)
        popfirst!(tokens)
        if isempty(tokens)
            return rules
        end
    end

    # If there are no more endline tokens
    # There are still rules left.
    if length(tokens) != 0
        # Still another rule
        if isa(lookahead(tokens), IdentifierToken)
        rules = [rules;parse_rules_list!(tokens)]
        end
    end

    return rules

end

function parse_rules_section!(tokens)
    """
    Rule Section
    """

    cur_token = popfirst!(tokens)
    section_name = nothing
    rules_list = []
    if isa(cur_token, RuleSectionToken)
        section_name = parse_symbol_name!(tokens)

        cur_token = popfirst!(tokens)
        if !isa(cur_token, LeftBracketToken)
            throw("Expected Left Bracket Token got $cur_token ")
        end

        rules_list = parse_rules_list!(tokens)
        cur_token = popfirst!(tokens)

        if !isa(cur_token, RightBracketToken)
            throw("Expected Right Bracket Token got: $cur_token")
        end

    else
        throw("Error expected Rule Section Token")
    end

    # Check if tokens are empty
    if isempty(tokens)
        return RuleSectionNode(section_name, rules_list)
    else

        while !isempty(tokens)
            cur_token = popfirst!(tokens)
            if isa(cur_token, EndLineToken)
                continue
            end
            throw("Expected empty token list, got: $(tokens)")
        end
        return RuleSectionNode(section_name, rules_list)
    end
end
