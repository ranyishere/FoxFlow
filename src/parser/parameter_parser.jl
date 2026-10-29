# Parameters-section parsing.

function parse_parameter!(tokens)

    # Ignoring new lines
    while isa(lookahead(tokens), EndLineToken)
        popfirst!(tokens)
    end

    # Here it is parsing type assignment.
    # I should be parsing either
    # type declaration
    # type assignment.

    cur_parameter = parse_type_assignment!(tokens)

    # Ignoring new lines
    while isa(lookahead(tokens), EndLineToken)
        popfirst!(tokens)
    end
    return cur_parameter
end

function parse_parameter_list!(tokens)
    """
    Parse a list of grammar parameters
    """

    parameter_list = []

    parameter = parse_parameter!(tokens)

    push!(parameter_list, parameter)

    while isa(lookahead(tokens), EndLineToken)
        popfirst!(tokens)
    end

    if isa(lookahead(tokens), RightBracketToken)
        return parameter_list
    else
        inner_parameter = parse_parameter_list!(tokens)
        parameter_list = [parameter_list;inner_parameter]
    end

    return parameter_list
end

function parse_parameters_section!(tokens)

    cur_token = popfirst!(tokens)

    if !isa(cur_token, ParameterSectionToken)
        throw(
              "Error expected ParameterSectionToken keyword got $(cur_token)"
             )
    end

    symbol_name = parse_symbol_name!(tokens)

    cur_token = popfirst!(tokens)

    if isa(cur_token, LeftBracketToken)
        parameters = parse_parameter_list!(tokens)
    else
        throw("Error expected LeftBracketToken keyword got $(cur_token)")
    end

    cur_token = popfirst!(tokens)
    if !isa(cur_token, RightBracketToken)
        throw("Error expected RightBracketToken keyword got $(cur_token)")
    end
    ParameterSectionNode(symbol_name, parameters)
end
