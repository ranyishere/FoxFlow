# Observables-section parsing.

"""
Parse a destructured parameter:  name : TypeName << alias1 : T1, alias2 : T2, ... >>
The token stream starts just after the opening '('.
Consumes up to (but not including) the closing ')'.
"""
function parse_observable_destructured_param!(tokens)
    while isa(lookahead(tokens), EndLineToken); popfirst!(tokens); end
    name = IdentifierNode(expect_token!(tokens, IdentifierToken))
    expect_token!(tokens, SingleColonToken)
    while isa(lookahead(tokens), EndLineToken); popfirst!(tokens); end
    type_tok = expect_token!(tokens, IdentifierToken)
    type_name = IdentifierNode(type_tok)

    # Optional field destructuring: << alias1 : T1, alias2 : T2, ... >>
    fields = FunctionArgNode[]
    while isa(lookahead(tokens), EndLineToken); popfirst!(tokens); end
    if isa(lookahead(tokens), LeftAngleBracketToken)
        popfirst!(tokens)  # consume <<
        # parse_function_args! reads until RightAngleBracketToken, skipping commas/parens/newlines
        fields = parse_function_args!(tokens)
        expect_token!(tokens, RightAngleBracketToken)  # >>
    end

    return ObservableDestructuredParamNode(name, type_name, fields)
end

"""
Parse a local Function definition inside an Observable body:
    name : Function << (p1 : TypeName << alias1 : T1, ... >>) >> -> RetType := { body; return expr }
"""
function parse_observable_local_fn!(tokens)
    while isa(lookahead(tokens), EndLineToken); popfirst!(tokens); end
    name = IdentifierNode(expect_token!(tokens, IdentifierToken))
    expect_token!(tokens, SingleColonToken)
    expect_token!(tokens, FunctionToken)

    # << (p1 : TypeName) << fields... >> >>
    while isa(lookahead(tokens), EndLineToken); popfirst!(tokens); end
    expect_token!(tokens, LeftAngleBracketToken)  # outer <<
    while isa(lookahead(tokens), Union{EndLineToken, LeftParenthesisToken}); popfirst!(tokens); end
    # Parse name : TypeName (fields come after the closing ')')
    fn_param_name = IdentifierNode(expect_token!(tokens, IdentifierToken))
    expect_token!(tokens, SingleColonToken)
    while isa(lookahead(tokens), EndLineToken); popfirst!(tokens); end
    fn_param_type = IdentifierNode(expect_token!(tokens, IdentifierToken))
    # consume closing ')'
    while isa(lookahead(tokens), Union{EndLineToken, RightParenthesisToken}); popfirst!(tokens); end
    # Optional fields: << alias1 : T1, alias2 : T2, ... >>
    fn_fields = FunctionArgNode[]
    while isa(lookahead(tokens), EndLineToken); popfirst!(tokens); end
    if isa(lookahead(tokens), LeftAngleBracketToken)
        popfirst!(tokens)  # consume <<
        fn_fields = parse_function_args!(tokens)
        expect_token!(tokens, RightAngleBracketToken)  # >>
    end
    while isa(lookahead(tokens), EndLineToken); popfirst!(tokens); end
    expect_token!(tokens, RightAngleBracketToken)  # outer >>
    param = ObservableDestructuredParamNode(fn_param_name, fn_param_type, fn_fields)

    # -> RetType
    while isa(lookahead(tokens), EndLineToken); popfirst!(tokens); end
    expect_token!(tokens, RightArrowToken)
    while isa(lookahead(tokens), EndLineToken); popfirst!(tokens); end
    ret_tok = popfirst!(tokens)
    ret_type = IdentifierNode(IdentifierToken(ret_tok.position))

    # :=
    while isa(lookahead(tokens), EndLineToken); popfirst!(tokens); end
    expect_token!(tokens, DefineToken)

    # { body; return expr; }
    fn_body = parse_function_body!(tokens)
    expect_token!(tokens, ReturnToken)
    body_ret = ReturnNode(parse_expression!(tokens))
    while isa(lookahead(tokens), EndLineToken); popfirst!(tokens); end
    expect_token!(tokens, RightBracketToken)

    return ObservableLocalFunctionNode(name, param, ret_type, fn_body, body_ret)
end

"""
Parse an observables section:

    observables <Name> {
        obs_name : Observable << (p1 : TypeName << alias1: T1, ... >>) >> -> RetType := {
            fn_name : Function << (p1 : TypeName << alias1: T1, ... >>) >> -> RetType := {
                return expr
            }
            return expr
        }
    }
"""
function parse_observables_section!(tokens)
    expect_token!(tokens, ObservableSectionToken)
    while isa(lookahead(tokens), EndLineToken); popfirst!(tokens); end
    section_name = IdentifierNode(expect_token!(tokens, IdentifierToken))
    while isa(lookahead(tokens), EndLineToken); popfirst!(tokens); end
    expect_token!(tokens, LeftBracketToken)

    definitions = ObservableDefinitionNode[]

    while true
        while isa(lookahead(tokens), EndLineToken); popfirst!(tokens); end
        if isa(lookahead(tokens), RightBracketToken)
            popfirst!(tokens); break
        end

        # obs_name : Observable << (p1 : TypeName << fields... >>) >> -> RetType := { ... }
        obs_name = IdentifierNode(expect_token!(tokens, IdentifierToken))
        expect_token!(tokens, SingleColonToken)
        expect_token!(tokens, ObservableToken)

        # << (p1 : TypeName) << fields... >> >>
        while isa(lookahead(tokens), EndLineToken); popfirst!(tokens); end
        expect_token!(tokens, LeftAngleBracketToken)  # outer <<
        while isa(lookahead(tokens), Union{EndLineToken, LeftParenthesisToken}); popfirst!(tokens); end
        # Parse name : TypeName (fields come after the closing ')')
        param_name = IdentifierNode(expect_token!(tokens, IdentifierToken))
        expect_token!(tokens, SingleColonToken)
        while isa(lookahead(tokens), EndLineToken); popfirst!(tokens); end
        param_type = IdentifierNode(expect_token!(tokens, IdentifierToken))
        # consume closing ')'
        while isa(lookahead(tokens), Union{EndLineToken, RightParenthesisToken}); popfirst!(tokens); end
        # Optional fields: << alias1 : T1, alias2 : T2, ... >>
        obs_fields = FunctionArgNode[]
        while isa(lookahead(tokens), EndLineToken); popfirst!(tokens); end
        if isa(lookahead(tokens), LeftAngleBracketToken)
            popfirst!(tokens)  # consume <<
            obs_fields = parse_function_args!(tokens)
            expect_token!(tokens, RightAngleBracketToken)  # >>
        end
        while isa(lookahead(tokens), EndLineToken); popfirst!(tokens); end
        expect_token!(tokens, RightAngleBracketToken)  # outer >>
        param = ObservableDestructuredParamNode(param_name, param_type, obs_fields)
        println("================================")

        # -> RetType
        while isa(lookahead(tokens), EndLineToken); popfirst!(tokens); end
        expect_token!(tokens, RightArrowToken)
        while isa(lookahead(tokens), EndLineToken); popfirst!(tokens); end
        ret_tok = popfirst!(tokens)
        ret_type = IdentifierNode(IdentifierToken(ret_tok.position))

        # :=
        while isa(lookahead(tokens), EndLineToken); popfirst!(tokens); end
        expect_token!(tokens, DefineToken)
        while isa(lookahead(tokens), EndLineToken); popfirst!(tokens); end
        expect_token!(tokens, LeftBracketToken)

        # Body: optional local Function defs, then return expr
        local_fns = ObservableLocalFunctionNode[]
        while true
            while isa(lookahead(tokens), EndLineToken); popfirst!(tokens); end
            if isa(lookahead(tokens), ReturnToken) || isa(lookahead(tokens), RightBracketToken)
                break
            end
            push!(local_fns, parse_observable_local_fn!(tokens))
        end

        expect_token!(tokens, ReturnToken)
        body_ret = ReturnNode(parse_expression!(tokens))
        while isa(lookahead(tokens), EndLineToken); popfirst!(tokens); end
        expect_token!(tokens, RightBracketToken)

        push!(definitions, ObservableDefinitionNode(obs_name, param, ret_type, local_fns, body_ret))
    end

    return ObservableSectionNode(section_name, definitions)
end
