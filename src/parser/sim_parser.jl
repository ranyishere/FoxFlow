
function remove_endlines_func_args!(tokens)
    """
    Removes endliens in function arguments
    """
    while lookahead(tokens) isa EndLineToken
        popfirst!(tokens)  # pop EndlineToken
    end
end

function parse_number!(tokens)
    if lookahead(tokens) isa IntegerToken
        # Pop type
        popfirst!(tokens)  # pop IntegerToken
        popfirst!(tokens)  # pop DefineToken
        int_token = popfirst!(tokens) # pop IntegerToken
        return IntegerNode(int_token)
    elseif lookahead(tokens) isa FloatToken
        popfirst!(tokens)  # pop FloatToken
        popfirst!(tokens)  # pop DefineToken
        float_token = popfirst!(tokens) # pop FloatToken
        return FloatNode(float_token)
    else
        throw(ErrorException("Expected IntegerToken or FloatToken got $(lookahead(tokens))"))
    end
end

function parse_string!(tokens)

    if !(lookahead(tokens) isa QuoteToken)
        throw(ErrorException("Expected QuoteToken got $(lookahead(tokens))"))
    end
    cur_token = popfirst!(tokens) # pop QuoteToken
    position = cur_token.position

    file_token = ""
    while !(lookahead(tokens) isa QuoteToken)

        if !(
             (lookahead(tokens) isa IdentifierToken) || (lookahead(tokens) isa DotToken)
             || (lookahead(tokens) isa SlashToken) ||(lookahead(tokens) isa BackslashToken)
            )
            throw(ErrorException("Expected IdentifierToken got $(lookahead(tokens))"))
        end

        file_token *= popfirst!(tokens).position.value
    end

    if !(lookahead(tokens) isa QuoteToken)
        throw(ErrorException("Expected QuoteToken got $(lookahead(tokens))"))
    end
    popfirst!(tokens) # pop QuoteToken
    return StringNode( StringToken(PositionToken(file_token, position.line_no, position.col_no)) )
end

function parse_sim_types!(tokens)
    if lookahead(tokens) isa DefineToken
        popfirst!(tokens)  # pop StateToken
    else
        throw(ErrorException("Expected := got $(lookahead(tokens))"))
    end

    if lookahead(tokens) isa LoadFileToken
        type_node = parse_load_file!(tokens)
        return SimulationTypesNode(type_node)
    else
        throw(ErrorException("Expected LoadFileToken got $(lookahead(tokens))"))
    end
end

function parse_sim_rules!(tokens)
    if !(lookahead(tokens) isa DefineToken)
        throw(ErrorException("Expected DefineToken got $(lookahead(tokens))"))
    end
    popfirst!(tokens)  # pop DefineToken
    if !(lookahead(tokens) isa LoadFileToken)
        throw(ErrorException("Expected LoadFileToken got $(lookahead(tokens))"))
    end
    load_node = parse_load_file!(tokens)

    return SimulationRulesNode(load_node)
end

function parse_sim_state!(tokens)
    if lookahead(tokens) isa DefineToken
        popfirst!(tokens)  # pop StateToken
    else
        throw(ErrorException("Expected := got $(lookahead(tokens))"))
    end

    if lookahead(tokens) isa LoadFileToken
        state_node = parse_load_file!(tokens)
        return SimulationStateNode(state_node)
    else
        throw(ErrorException("Expected LoadFileToken got $(lookahead(tokens))"))
    end
end

function parse_load_file!(tokens)
    """
    Parses a load file statement from the token stream.
    """

    if !(lookahead(tokens) isa LoadFileToken)
        throw(ErrorException("Expected LoadFileToken got $(lookahead(tokens))"))
    end

    popfirst!(tokens) # pop LoadFileToken

    # Should be (
    if lookahead(tokens) isa LeftParenthesisToken
        popfirst!(tokens) # pop LeftParenToken
    else
        throw(ErrorException("Expected LeftParenToken got $(lookahead(tokens))"))
    end

    file_token = nothing
    if lookahead(tokens) isa IdentifierToken
        file_token = popfirst!(tokens) # pop StringToken
    elseif lookahead(tokens) isa QuoteToken
        file_token = parse_string!(tokens)
    elseif lookahead(tokens) isa StringToken
        file_token = StringNode(popfirst!(tokens)) # pop StringToken
    else
        throw(ErrorException("Expected StringToken got $(lookahead(tokens))"))
    end

    if lookahead(tokens) isa RightParenthesisToken
        popfirst!(tokens) # pop EndLineToken
    else
        throw(ErrorException("Expected RightParenToken got $(lookahead(tokens))"))
    end

    # exit(0)
    LoadNode(file_token)
end

function parse_sim_param!(tokens)

    if !(lookahead(tokens) isa DefineToken)
        throw(ErrorException("Expected Identifier token got $(lookahead(tokens))"))
    end

    popfirst!(tokens) # pop DefineToken

    if lookahead(tokens) isa LoadFileToken
        load_node = parse_load_file!(tokens)
        return SimulationParametersNode(load_node)
    else
        throw(ErrorException("Expected LoadFileToken got $(lookahead(tokens))"))
    end
end

function parse_run_simulation!(tokens)

    if !(lookahead(tokens) isa DefineToken)
        throw(ErrorException("Expected Define Token got $(lookahead(tokens))"))
    end

    popfirst!(tokens)  # pop DefineToken

    if !(lookahead(tokens) isa RunSimulationToken)
        throw(ErrorException("Expected RunToken got $(lookahead(tokens))"))
    end
    popfirst!(tokens)  # pop RunSimulationToken

    if !(lookahead(tokens) isa LeftParenthesisToken)
        throw(ErrorException("Expected LeftParenthesisToken got $(lookahead(tokens))"))
    end
    popfirst!(tokens)  # pop LeftParenthesisToken

    while lookahead(tokens) isa EndLineToken
        popfirst!(tokens)  # pop EndlineToken
    end

    # Grab initial state
    initial_state = nothing
    if lookahead(tokens) isa IdentifierToken
        initial_state = IdentifierNode(popfirst!(tokens)) # pop IdentifierToken
    else
        throw(ErrorException("Expected IdentifierToken got $(lookahead(tokens))"))
    end

    if lookahead(tokens) isa PunctuationToken
        popfirst!(tokens)  # pop CommaToken
    else
        throw(ErrorException("Expected CommaToken got $(lookahead(tokens))"))
    end

    remove_endlines_func_args!(tokens)

    # Grab parameter
    sim_paramters = nothing
    if lookahead(tokens) isa IdentifierToken
        sim_paramters = IdentifierNode(popfirst!(tokens)) # pop IdentifierToken
    else
        throw(ErrorException("Expected IdentifierToken got $(lookahead(tokens))"))
    end

    if lookahead(tokens) isa PunctuationToken
        popfirst!(tokens)  # pop CommaToken
    else
        throw(ErrorException("Expected CommaToken got $(lookahead(tokens))"))
    end
    remove_endlines_func_args!(tokens)

    # Grab Rules
    sim_rules = nothing
    if lookahead(tokens) isa IdentifierToken
        sim_rules = IdentifierNode(popfirst!(tokens)) # pop IdentifierToken
    else
        throw(ErrorException("Expected IdentifierToken got $(lookahead(tokens))"))
    end

    if lookahead(tokens) isa PunctuationToken
        popfirst!(tokens)  # pop CommaToken
    else
        throw(ErrorException("Expected CommaToken got $(lookahead(tokens))"))
    end
    remove_endlines_func_args!(tokens)

    # Grab Types
    sim_types = nothing
    if lookahead(tokens) isa IdentifierToken
        sim_types = IdentifierNode(popfirst!(tokens)) # pop IdentifierToken
    else
        throw(ErrorException("Expected IdentifierToken got $(lookahead(tokens))"))
    end
    if lookahead(tokens) isa PunctuationToken
        popfirst!(tokens)  # pop CommaToken
    else
        throw(ErrorException("Expected CommaToken got $(lookahead(tokens))"))
    end
    remove_endlines_func_args!(tokens)

    # Grab Iterations
    sim_iterations = nothing
    if lookahead(tokens) isa IntegerToken
        sim_iterations = parse_number!(tokens)
    elseif lookahead(tokens) isa FloatToken
        sim_iterations = parse_number!(tokens)
    elseif lookahead(tokens) isa IdentifierToken
        sim_iterations = IdentifierNode(popfirst!(tokens)) # pop IdentifierToken
    else
        throw(ErrorException("Expected IntegerToken got $(lookahead(tokens))"))
    end

    if !(lookahead(tokens) isa RightParenthesisToken)
        throw(ErrorException("Expected RightParenthesisToken got $(lookahead(tokens))"))
    end
    popfirst!(tokens)  # pop RightParenthesisToken

    return RunSimulationNode(initial_state, sim_paramters, sim_rules, sim_types, sim_iterations)
end

function parse_simulation_declaration!(tokens)
    """
    Parses a single simulation declaration from the token stream.
    """

    if !(lookahead(tokens) isa IdentifierToken)
        throw(ErrorException("Expected Identifier token got $(lookahead(tokens))"))
    end

    node_name = IdentifierNode(popfirst!(tokens)) # pop IdentifierToken

    # : Token
    if !(lookahead(tokens) isa SingleColonToken)
        throw(ErrorException("Expected ColonToken got $(lookahead(tokens))"))
    end

    popfirst!(tokens)  # pop EqualToken

    node_value = nothing

    # Type
    if (lookahead(tokens) isa SimulationParametersToken)
        popfirst!(tokens)  # pop ParameterToken
        node_value = parse_sim_param!(tokens)
    elseif (lookahead(tokens) isa StateToken)
        popfirst!(tokens)  # pop StateToken
        node_value = parse_sim_state!(tokens)
    elseif (lookahead(tokens) isa IntegerToken || lookahead(tokens) isa FloatToken)
        node_value = parse_number!(tokens)
    elseif (lookahead(tokens) isa RulesToken)
        popfirst!(tokens)  # pop RulesToken
        node_value = parse_sim_rules!(tokens)
    elseif (lookahead(tokens) isa SimulationTypesToken)
        popfirst!(tokens)
        node_value = parse_sim_types!(tokens)
    elseif lookahead(tokens) isa SimulationToken
        popfirst!(tokens)  # pop SimulationToken
        node_value = parse_run_simulation!(tokens)
    else
        throw(ErrorException("Expected SimulationParameterToken, StateToken, RulesToken, IntegerToken or FloatToken got $(lookahead(tokens))"))
    end

    while lookahead(tokens) isa EndLineToken
        popfirst!(tokens)  # pop EndlineToken
    end

    decl_node = SimDeclarationNode(node_name, node_value)
end

function parse_sim_declarations!(tokens)
    """
    Declarations of simulations inside a simulations section.
    """

    declarations = []
    while lookahead(tokens) != RightBracketToken
        if length(tokens) == 0
            throw(ErrorException("Unexpected end of file while parsing simulation declarations"))
        end

        if lookahead(tokens) isa EndLineToken
            popfirst!(tokens)  # pop EndlineToken
            continue
        end

        if lookahead(tokens) isa RightBracketToken
            break
        end

        decl = parse_simulation_declaration!(tokens)
        push!(declarations, decl)
    end

    return declarations
end

function parse_simulations_section!(tokens)
    """
    Parses a simulations section from the token stream and updates the simulation table.
    """

    if !(lookahead(tokens) isa SimulationSectionToken)
        throw(ErrorException("Expected SimulationSectionToken got $(lookahead(tokens))"))
    end
    popfirst!(tokens)  # pop SimulationSectionToken

    if !(lookahead(tokens) isa IdentifierToken)
        throw(ErrorException("Expected Identifier token got $(lookahead(tokens))"))
    end

    sim_name = IdentifierNode(popfirst!(tokens)) # pop IdentifierToken

    if !(lookahead(tokens) isa LeftBracketToken)
        throw(ErrorException("Expected LeftBracketToken got $(lookahead(tokens))"))
    end
    popfirst!(tokens)  # pop LeftBracketToken

    declarations = parse_sim_declarations!(tokens)

    if !(lookahead(tokens) isa RightBracketToken)
        throw(ErrorException("Expected RightBracketToken got $(lookahead(tokens))"))
    end
    popfirst!(tokens)  # pop RightBracketToken

    SimulationSectionNode(sim_name, declarations)
end
