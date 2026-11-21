
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
    elseif lookahead(tokens) isa IdentifierToken
        # generator call, e.g. grid(Vertex, u, 0.0, [0.0, 1.0, 10])
        return SimulationGeneratorStateNode(parse_expression!(tokens))
    else
        throw(ErrorException("Expected LoadFileToken or a generator call got $(lookahead(tokens))"))
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

function expect_comma!(tokens, ctx)
    if !(lookahead(tokens) isa PunctuationToken)
        throw(ErrorException("$ctx: expected , got $(lookahead(tokens))"))
    end
    popfirst!(tokens)  # pop ,
    remove_endlines_func_args!(tokens)
end

function expect_ident!(tokens, ctx)
    if !(lookahead(tokens) isa IdentifierToken)
        throw(ErrorException("$ctx: expected an identifier got $(lookahead(tokens))"))
    end
    return IdentifierNode(popfirst!(tokens))
end

function parse_sim_functions!(tokens)
    """
    `f : SimulationFunctions := load("functions.fflow")` (after the type name was consumed).
    """
    if !(lookahead(tokens) isa DefineToken)
        throw(ErrorException("Expected := got $(lookahead(tokens))"))
    end
    popfirst!(tokens)  # pop :=
    if !(lookahead(tokens) isa LoadFileToken)
        throw(ErrorException("Expected load(...) got $(lookahead(tokens))"))
    end
    return SimulationFunctionsNode(parse_load_file!(tokens))
end

function parse_sim_grammar!(tokens)
    """
    `g : Grammar := Grammar(a, b, c [, d])` (after the type name was consumed). The arguments
    name SimulationFunctions / SimulationRules / SimulationTypes / SimulationObservables
    declarations; the IR classifies them by their declared type.
    """
    if !(lookahead(tokens) isa DefineToken)
        throw(ErrorException("Expected := got $(lookahead(tokens))"))
    end
    popfirst!(tokens)  # pop :=
    if !(lookahead(tokens) isa IdentifierToken && lookahead(tokens).position.value == "Grammar")
        throw(ErrorException("Expected Grammar(...) got $(lookahead(tokens))"))
    end
    popfirst!(tokens)  # pop Grammar
    if !(lookahead(tokens) isa LeftParenthesisToken)
        throw(ErrorException("Grammar: expected ( got $(lookahead(tokens))"))
    end
    popfirst!(tokens)  # pop (
    remove_endlines_func_args!(tokens)

    parts = IdentifierNode[]
    while true
        remove_endlines_func_args!(tokens)
        push!(parts, expect_ident!(tokens, "Grammar"))
        remove_endlines_func_args!(tokens)
        if lookahead(tokens) isa PunctuationToken
            popfirst!(tokens)  # pop ,
        else
            break
        end
    end
    if !(lookahead(tokens) isa RightParenthesisToken)
        throw(ErrorException("Grammar: expected ) got $(lookahead(tokens))"))
    end
    popfirst!(tokens)  # pop )
    return SimulationGrammarNode(parts)
end

function parse_simulation_config!(tokens)
    """
    `Simulation(["save.bin",] grammar, state, params, time [, checkpoint])`; the leading
    `Simulation` token is the current token.
    """
    popfirst!(tokens)  # pop Simulation
    if !(lookahead(tokens) isa LeftParenthesisToken)
        throw(ErrorException("Simulation: expected ( got $(lookahead(tokens))"))
    end
    popfirst!(tokens)  # pop (
    remove_endlines_func_args!(tokens)

    save_name = nothing
    if lookahead(tokens) isa StringToken
        save_name = StringNode(popfirst!(tokens))
        expect_comma!(tokens, "Simulation")
    end

    grammar = expect_ident!(tokens, "Simulation (grammar)")
    expect_comma!(tokens, "Simulation")
    state = expect_ident!(tokens, "Simulation (initial state)")
    expect_comma!(tokens, "Simulation")
    params = expect_ident!(tokens, "Simulation (parameters)")
    expect_comma!(tokens, "Simulation")

    steps = nothing
    if lookahead(tokens) isa IntegerToken
        steps = IntegerNode(popfirst!(tokens))
    elseif lookahead(tokens) isa IdentifierToken
        steps = IdentifierNode(popfirst!(tokens))
    else
        throw(ErrorException("Simulation: expected the simulation time got $(lookahead(tokens))"))
    end

    remove_endlines_func_args!(tokens)
    checkpoint = nothing
    dim = nothing
    while lookahead(tokens) isa PunctuationToken
        popfirst!(tokens)
        remove_endlines_func_args!(tokens)
        if length(tokens) >= 2 && lookahead(tokens) isa IdentifierToken && tokens[2] isa EqualToken
            # keyword option: `dim = 2|3`
            opt = expect_ident!(tokens, "Simulation (option)").token.position.value
            popfirst!(tokens)  # pop =
            if opt == "dim"
                if !(lookahead(tokens) isa IntegerToken)
                    throw(ErrorException("Simulation: dim must be 2 or 3, got $(lookahead(tokens))"))
                end
                dim = tryparse(Int, string(get_value(IntegerNode(popfirst!(tokens)))))
                if dim != 2 && dim != 3
                    throw(ErrorException("Simulation: dim must be 2 or 3, got $dim"))
                end
            else
                throw(ErrorException("Simulation: unknown option `$opt` (expected dim)"))
            end
        elseif checkpoint === nothing
            if lookahead(tokens) isa IntegerToken
                checkpoint = IntegerNode(popfirst!(tokens))
            elseif lookahead(tokens) isa IdentifierToken
                checkpoint = IdentifierNode(popfirst!(tokens))
            else
                throw(ErrorException("Simulation: expected the checkpoint frequency got $(lookahead(tokens))"))
            end
        else
            throw(ErrorException("Simulation: unexpected argument $(lookahead(tokens)) after the checkpoint frequency (options are written `dim = 2|3`)"))
        end
        remove_endlines_func_args!(tokens)
    end

    if !(lookahead(tokens) isa RightParenthesisToken)
        throw(ErrorException("Simulation: expected ) got $(lookahead(tokens))"))
    end
    popfirst!(tokens)  # pop )
    return SimulationConfigNode(save_name, grammar, state, params, steps, checkpoint, dim)
end

function parse_run_statement!(tokens)
    """
    `run(sim, N [, domain_parallel = true|false] [, jobs = K] [, seed = S] [, name = "exp"]
     [, listeners = true|false] [, log = true|false])`.
    """
    popfirst!(tokens)  # pop `run`
    popfirst!(tokens)  # pop (
    remove_endlines_func_args!(tokens)

    sim = expect_ident!(tokens, "run (simulation)")
    expect_comma!(tokens, "run")

    replicates = nothing
    if lookahead(tokens) isa IntegerToken
        replicates = IntegerNode(popfirst!(tokens))
    elseif lookahead(tokens) isa IdentifierToken
        replicates = IdentifierNode(popfirst!(tokens))
    else
        throw(ErrorException("run: expected the number of runs got $(lookahead(tokens))"))
    end

    domain_parallel = nothing
    jobs = nothing
    seed = nothing
    run_name = nothing
    listeners = nothing
    run_log = nothing
    while true
        remove_endlines_func_args!(tokens)
        lookahead(tokens) isa PunctuationToken || break
        popfirst!(tokens)  # pop ,
        remove_endlines_func_args!(tokens)

        opt = expect_ident!(tokens, "run (option)").token.position.value
        if !(lookahead(tokens) isa EqualToken)
            throw(ErrorException("run: expected `$opt = value`, got $(lookahead(tokens)) after $opt"))
        end
        popfirst!(tokens)  # pop =

        if opt == "domain_parallel" || opt == "listeners" || opt == "log"
            v = lookahead(tokens) isa IdentifierToken ? lookahead(tokens).position.value : ""
            if v != "true" && v != "false"
                throw(ErrorException("run: $opt must be true or false, got $(lookahead(tokens))"))
            end
            popfirst!(tokens)
            if opt == "domain_parallel"
                domain_parallel = (v == "true")
            elseif opt == "listeners"
                listeners = (v == "true")
            else
                run_log = (v == "true")
            end
        elseif opt == "jobs" || opt == "seed"
            val = nothing
            if lookahead(tokens) isa IntegerToken
                val = IntegerNode(popfirst!(tokens))
            elseif lookahead(tokens) isa IdentifierToken
                val = IdentifierNode(popfirst!(tokens))
            else
                throw(ErrorException("run: $opt must be an integer, got $(lookahead(tokens))"))
            end
            opt == "jobs" ? (jobs = val) : (seed = val)
        elseif opt == "name"
            lookahead(tokens) isa StringToken ||
                throw(ErrorException("run: name must be a string, e.g. name = \"exp_a\", got $(lookahead(tokens))"))
            run_name = StringNode(popfirst!(tokens))
        else
            throw(ErrorException("run: unknown option `$opt` (expected domain_parallel, jobs, seed, name, listeners or log)"))
        end
    end

    if !(lookahead(tokens) isa RightParenthesisToken)
        throw(ErrorException("run: expected ) got $(lookahead(tokens))"))
    end
    popfirst!(tokens)  # pop )
    while lookahead(tokens) isa EndLineToken
        popfirst!(tokens)
    end
    return RunNode(sim, replicates, domain_parallel, jobs, seed, run_name, listeners, run_log)
end

function parse_run_simulation!(tokens)

    if !(lookahead(tokens) isa DefineToken)
        throw(ErrorException("Expected Define Token got $(lookahead(tokens))"))
    end

    popfirst!(tokens)  # pop DefineToken

    if lookahead(tokens) isa SimulationToken
        return parse_simulation_config!(tokens)
    end

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

    # Optional first argument: save name for the last save-state bin file.
    # Provided as a string literal, e.g. RunSimulation("my_state.bin", ...).
    save_name = nothing
    if lookahead(tokens) isa StringToken
        save_name = StringNode(popfirst!(tokens))  # pop StringToken

        if lookahead(tokens) isa PunctuationToken
            popfirst!(tokens)  # pop CommaToken
        else
            throw(ErrorException("Expected CommaToken after save name got $(lookahead(tokens))"))
        end
        remove_endlines_func_args!(tokens)
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

    # Optional trailing arguments, in order: , checkpointFrequency , observablesVar
    #   checkpointFrequency: an integer literal (or an identifier bound to an integer)
    #                        = write output every N steps.
    #   observablesVar:      identifier bound to an observables file.
    # A single trailing identifier is kept as the observables argument for backwards
    # compatibility; the IR stage reclassifies it if it names an integer.
    remove_endlines_func_args!(tokens)
    sim_observables = nothing
    sim_checkpoint_frequency = nothing
    if isa(lookahead(tokens), PunctuationToken)
        popfirst!(tokens)  # pop CommaToken
        remove_endlines_func_args!(tokens)
        if lookahead(tokens) isa IntegerToken
            # Numeric literal (e.g. 60), not an `Integer := 60` declaration.
            sim_checkpoint_frequency = IntegerNode(popfirst!(tokens))
        elseif lookahead(tokens) isa IdentifierToken
            sim_observables = IdentifierNode(popfirst!(tokens))
        else
            throw(ErrorException("Expected IntegerToken or IdentifierToken after simulation time got $(lookahead(tokens))"))
        end

        remove_endlines_func_args!(tokens)
        if isa(lookahead(tokens), PunctuationToken)
            popfirst!(tokens)  # pop CommaToken
            remove_endlines_func_args!(tokens)
            if !(lookahead(tokens) isa IdentifierToken)
                throw(ErrorException("Expected IdentifierToken for observables arg got $(lookahead(tokens))"))
            end
            if sim_checkpoint_frequency === nothing
                # first extra arg was an identifier => it is the checkpoint frequency
                sim_checkpoint_frequency = sim_observables
            end
            sim_observables = IdentifierNode(popfirst!(tokens))
        end
    end

    if !(lookahead(tokens) isa RightParenthesisToken)
        throw(ErrorException("Expected RightParenthesisToken got $(lookahead(tokens))"))
    end
    popfirst!(tokens)  # pop RightParenthesisToken

    return RunSimulationNode(save_name, initial_state, sim_paramters, sim_rules, sim_types, sim_iterations,
                             sim_checkpoint_frequency, sim_observables)
end

function parse_sim_observables!(tokens)
    if !(lookahead(tokens) isa DefineToken)
        throw(ErrorException("Expected := got $(lookahead(tokens))"))
    end
    popfirst!(tokens)  # pop DefineToken
    if !(lookahead(tokens) isa LoadFileToken)
        throw(ErrorException("Expected LoadFileToken got $(lookahead(tokens))"))
    end
    load_node = parse_load_file!(tokens)
    return SimulationObservablesNode(load_node)
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
    elseif (lookahead(tokens) isa SimulationObservablesToken)
        popfirst!(tokens)
        node_value = parse_sim_observables!(tokens)
    elseif lookahead(tokens) isa SimulationToken
        popfirst!(tokens)  # pop SimulationToken
        node_value = parse_run_simulation!(tokens)
    elseif lookahead(tokens) isa IdentifierToken && lookahead(tokens).position.value == "Grammar"
        popfirst!(tokens)  # pop Grammar
        node_value = parse_sim_grammar!(tokens)
    elseif lookahead(tokens) isa IdentifierToken && lookahead(tokens).position.value == "SimulationFunctions"
        popfirst!(tokens)  # pop SimulationFunctions
        node_value = parse_sim_functions!(tokens)
    else
        throw(ErrorException("Expected SimulationParameterToken, StateToken, RulesToken, Grammar, SimulationFunctions, IntegerToken or FloatToken got $(lookahead(tokens))"))
    end

    while lookahead(tokens) isa EndLineToken
        popfirst!(tokens)  # pop EndlineToken
    end

    decl_node = SimDeclarationNode(node_name, node_value)
end

function parse_trace_statement!(tokens)
    """
    Parses `trace(simulationName, "rule", ...)`. Rule names are string literals
    ("*" selects every rule).
    """

    popfirst!(tokens)  # pop `trace`
    popfirst!(tokens)  # pop (
    remove_endlines_func_args!(tokens)

    if !(lookahead(tokens) isa IdentifierToken)
        throw(ErrorException("trace: expected a simulation name as first argument, got $(lookahead(tokens))"))
    end
    sim_node = IdentifierNode(popfirst!(tokens))

    rules = StringNode[]
    while true
        remove_endlines_func_args!(tokens)
        if lookahead(tokens) isa PunctuationToken
            popfirst!(tokens)  # pop ,
            remove_endlines_func_args!(tokens)
            if !(lookahead(tokens) isa StringToken)
                throw(ErrorException("trace: rule names must be string literals, e.g. trace(sim, \"crossover\"); got $(lookahead(tokens))"))
            end
            push!(rules, StringNode(popfirst!(tokens)))
        else
            break
        end
    end

    if !(lookahead(tokens) isa RightParenthesisToken)
        throw(ErrorException("trace: expected ) got $(lookahead(tokens))"))
    end
    popfirst!(tokens)  # pop )

    if isempty(rules)
        throw(ErrorException("trace($(sim_node.token.position.value)): give at least one rule name, e.g. trace($(sim_node.token.position.value), \"crossover\") or \"*\" for all rules"))
    end

    while lookahead(tokens) isa EndLineToken
        popfirst!(tokens)
    end
    return TraceNode(sim_node, rules)
end

function parse_measure_statement!(tokens)
    """
    Parses `measure(sim1, sim2, ..., "observable", ...)`. Simulation names are identifiers,
    observable names are string literals ("*" or none selects every observable).
    """

    popfirst!(tokens)  # pop `measure`
    popfirst!(tokens)  # pop (

    sims = IdentifierNode[]
    obs = StringNode[]
    while true
        remove_endlines_func_args!(tokens)
        tok = lookahead(tokens)
        if tok isa IdentifierToken
            isempty(obs) || throw(ErrorException("measure: list the simulations before the observable names, got $(tok.position.value) after a string"))
            push!(sims, IdentifierNode(popfirst!(tokens)))
        elseif tok isa StringToken
            push!(obs, StringNode(popfirst!(tokens)))
        else
            throw(ErrorException("measure: expected a simulation name or an observable name string, got $(tok)"))
        end
        remove_endlines_func_args!(tokens)
        if lookahead(tokens) isa PunctuationToken
            popfirst!(tokens)  # pop ,
        else
            break
        end
    end

    if !(lookahead(tokens) isa RightParenthesisToken)
        throw(ErrorException("measure: expected ) got $(lookahead(tokens))"))
    end
    popfirst!(tokens)  # pop )

    if isempty(sims)
        throw(ErrorException("measure: give at least one simulation, e.g. measure(mySim) or measure(simA, simB, \"n_tips\")"))
    end

    while lookahead(tokens) isa EndLineToken
        popfirst!(tokens)
    end
    return MeasureNode(sims, obs)
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

        # `trace(sim, "rule", ...)` statement (contextual keyword, not reserved).
        if lookahead(tokens) isa IdentifierToken && lookahead(tokens).position.value == "trace" &&
                length(tokens) >= 2 && tokens[2] isa LeftParenthesisToken
            push!(declarations, parse_trace_statement!(tokens))
            continue
        end

        # `measure(sim1, sim2, "observable", ...)` statement (contextual keyword, not reserved).
        if lookahead(tokens) isa IdentifierToken && lookahead(tokens).position.value == "measure" &&
                length(tokens) >= 2 && tokens[2] isa LeftParenthesisToken
            push!(declarations, parse_measure_statement!(tokens))
            continue
        end

        # `run(sim, N, domain_parallel = false, ...)` statement (contextual keyword, not reserved).
        if lookahead(tokens) isa IdentifierToken && lookahead(tokens).position.value == "run" &&
                length(tokens) >= 2 && tokens[2] isa LeftParenthesisToken
            push!(declarations, parse_run_statement!(tokens))
            continue
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
