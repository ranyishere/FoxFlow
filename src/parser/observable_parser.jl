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
Parse a global observable (the `name : Observable` prefix is already consumed):

    := <reduction> <pattern> [within <r>] [<< alias : node::Attribute, ... >>] [return ( expr )]

The pattern uses the rule-pattern syntax, e.g. `(a : Intermediate) -- (b : Positive)`. It must be
connected, except for a pair observable: `(a : T) (b : U) within r` has exactly two nodes and no
edges, and matches every ordered pair of distinct nodes closer than `r` (found with a grid, so
it stays cheap on large graphs). `expr` is the deterministic value a matched configuration returns; it is
evaluated once per match and `reduction` (sum, mean, count, min, max) combines the
values. `count` needs no expression (it counts matches; with an
expression it counts the matches where the expression is nonzero).
"""
function parse_global_observable!(tokens, name::IdentifierNode)
    skip_eols!() = while !isempty(tokens) && isa(lookahead(tokens), EndLineToken); popfirst!(tokens); end
    # Index of the next non-newline token (without consuming anything).
    function peek_ix()
        j = 1
        while j <= length(tokens) && tokens[j] isa EndLineToken; j += 1; end
        return j
    end

    name_str = String(name.token.position.value)
    pos = name.token.position

    skip_eols!()
    expect_token!(tokens, DefineToken)
    skip_eols!()
    red_tok = expect_token!(tokens, IdentifierToken)
    reduction = String(red_tok.position.value)
    if !(reduction in OBSERVABLE_REDUCTIONS)
        throw("Observable '$name_str': unknown reduction '$reduction' (expected one of " *
              join(OBSERVABLE_REDUCTIONS, ", ") * ")")
    end

    # ---- pattern ----
    lhs = Any[]
    while true
        j = peek_ix()
        j > length(tokens) && break
        t = tokens[j]
        if t isa LeftParenthesisToken || t isa EdgeToken
            for _ in 1:(j - 1); popfirst!(tokens); end
            if t isa EdgeToken
                isempty(lhs) && throw("Observable '$name_str': pattern cannot start with an edge")
                vert_0 = pop!(lhs)
                if vert_0 isa UndirectedTypeEdgeNode
                    tmp = vert_0.right_vertex
                    push!(lhs, vert_0)
                    vert_0 = tmp
                end
                popfirst!(tokens)                         # --
                skip_eols!()
                expect_token!(tokens, LeftParenthesisToken)
                vert_1 = parse_type_assignment!(tokens)
                expect_token!(tokens, RightParenthesisToken)
                push!(lhs, UndirectedTypeEdgeNode(vert_0, vert_1))
            else
                popfirst!(tokens)                         # (
                push!(lhs, parse_type_assignment!(tokens))
                expect_token!(tokens, RightParenthesisToken)
            end
        else
            break
        end
    end
    isempty(lhs) && throw("Observable '$name_str': missing pattern, e.g. (p : Positive)")

    # ---- optional `within r`: pair observable over two unconnected nodes closer than r ----
    within = nothing
    j = peek_ix()
    if j <= length(tokens) && tokens[j] isa IdentifierToken &&
       String(tokens[j].position.value) == "within"
        for _ in 1:j; popfirst!(tokens); end              # newlines + `within`
        skip_eols!()
        r_tok = isempty(tokens) ? nothing : popfirst!(tokens)
        r_tok isa Union{FloatToken, IntegerToken} ||
            throw("Observable '$name_str': `within` expects a distance, e.g. `within 0.1`")
        within = tryparse(Float64, String(r_tok.position.value))
        (within === nothing || within <= 0) &&
            throw("Observable '$name_str': `within` expects a positive number, got '$(r_tok.position.value)'")
        (length(lhs) == 2 && all(!(el isa UndirectedTypeEdgeNode) for el in lhs)) ||
            throw("Observable '$name_str': `within r` needs a pattern of exactly two nodes, e.g. " *
                  "`(a : T) (b : T) within 0.1`")
    end

    # Connectivity: every node must be reachable from the first through the pattern's edges.
    # (A `within` pair observable is the exception: its two nodes are joined by distance.)
    node_name(v) = String(v.name.token.position.value)
    parent = Dict{String, String}()
    find(x) = (while parent[x] != x; parent[x] = parent[parent[x]]; x = parent[x]; end; x)
    for el in lhs
        if el isa UndirectedTypeEdgeNode
            for v in (el.left_vertex, el.right_vertex)
                haskey(parent, node_name(v)) || (parent[node_name(v)] = node_name(v))
            end
            parent[find(node_name(el.left_vertex))] = find(node_name(el.right_vertex))
        else
            haskey(parent, node_name(el)) || (parent[node_name(el)] = node_name(el))
        end
    end
    if within === nothing && length(unique(find(k) for k in keys(parent))) > 1
        throw("Observable '$name_str': the pattern must be connected (join its nodes with `--`)")
    end

    # ---- optional bindings << alias : node::Attribute, ... >> ----
    lhs_params = ParameterNode([])
    j = peek_ix()
    if j <= length(tokens) && tokens[j] isa LeftAngleBracketToken
        skip_eols!()
        lhs_params = parse_symbol_parameters!(tokens)
    end

    # ---- optional `return ( expr )` (`with` is reserved for rule propensities) ----
    expr = nothing
    j = peek_ix()
    if j <= length(tokens) && tokens[j] isa WithToken
        throw("Observable '$name_str': `with` is reserved for rule propensities; " *
              "write `return ( expression )` for the per-match value")
    end
    if j <= length(tokens) && tokens[j] isa ReturnToken
        skip_eols!()
        popfirst!(tokens)                                 # return
        skip_eols!()
        isa(lookahead(tokens), LeftParenthesisToken) ||
            throw("Observable '$name_str': expected `(` after `return`")
        expr = parse_expression!(tokens)
        if !(expr isa Union{FunctionNode, BinaryOpNode, IdentifierNode, FloatNode, IntegerNode, GroupNode})
            expr = GroupNode(expr)
        end
    elseif reduction == "count"
        expr = FloatNode(FloatToken(PositionToken("1.0", pos.line_no, pos.col_no)))
    else
        throw("Observable '$name_str': `$reduction` needs a value, add `return ( expression )`")
    end

    rule_name = IdentifierNode(IdentifierToken(PositionToken("obs_" * name_str, pos.line_no, pos.col_no)))
    rule = RuleNode(rule_name, lhs, lhs_params, copy(lhs), ParameterNode([]),
                    WithClauseNode(expr, WhereClauseNode(Node[])))
    return GlobalObservableNode(name, reduction, rule, within)
end

"""
Parse an observables section:

    observables <Name> {
        obs_name : Attribute << (p1 : TypeName << alias1: T1, ... >>) >> -> RetType := {
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
    global_observables = GlobalObservableNode[]

    while true
        while isa(lookahead(tokens), EndLineToken); popfirst!(tokens); end
        if isa(lookahead(tokens), RightBracketToken)
            popfirst!(tokens); break
        end

        # Two kinds of entries:
        #   name : Attribute  << (p1 : TypeName << fields... >>) >> -> Rehas ( expr )]  }   (per node)
        #   name : Observable := <reduction> <pattern> [<< bindings >>] [with ( expr )]      (global)
        # The former spelling `name : Observable << ... >>` is still accepted as an Attribute.
        obs_name = IdentifierNode(expect_token!(tokens, IdentifierToken))
        expect_token!(tokens, SingleColonToken)
        kind_tok = popfirst!(tokens)
        if kind_tok isa ObservableToken
            j = 1
            while j <= length(tokens) && tokens[j] isa EndLineToken; j += 1; end
            if j <= length(tokens) && tokens[j] isa DefineToken
                push!(global_observables, parse_global_observable!(tokens, obs_name))
                continue
            end
        elseif !(kind_tok isa AttributeToken)
            throw("Expected `Attribute` or `Observable` after `:` in observables section, got $kind_tok")
        end

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

    return ObservableSectionNode(section_name, definitions, global_observables)
end
