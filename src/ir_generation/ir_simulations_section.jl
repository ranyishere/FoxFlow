
# Declared kind ("rules"/"types"/"observables"/"functions") of each file-loading declaration of the
# simulations section being processed; `Grammar(...)` arguments are classified with it.
const SIM_DECL_KINDS = Dict{String, String}()

function get_value_sim_node(node)
    if (node isa SimulationParametersNode 
        || node isa SimulationNode 
        || node isa SimulationRulesNode 
        || node isa SimulationTypesNode 
        || node isa SimulationObservablesNode
        || node isa SimulationStateNode)
        return node.value.filepath
    elseif node isa FloatNode
        return node.value
    elseif node isa IntegerNode
        return node.value
    elseif node isa StringNode
        return node.value
    else
        error("Unsupported node type for getting value: $(typeof(node))")
    end
end

function populate_sim_table(node_name, node_value, sim_table)
    if (node_value isa SimulationParametersNode  ||
          node_value isa SimulationNode 
          || node_value isa  SimulationRulesNode || node_value isa  SimulationTypesNode ||
          node_value isa SimulationObservablesNode || node_value isa SimulationFunctionsNode ||
          node_value isa SimulationStateNode)
        value = get_value(node_value.value.filepath)
        sim_table[node_name] = value
        # Remember what the declaration is so `Grammar(...)` can classify its arguments.
        kind = node_value isa SimulationRulesNode ? "rules" :
               node_value isa SimulationTypesNode ? "types" :
               node_value isa SimulationObservablesNode ? "observables" :
               node_value isa SimulationFunctionsNode ? "functions" : nothing
        kind === nothing || (SIM_DECL_KINDS[node_name] = kind)
    elseif node_value isa SimulationGeneratorStateNode
        value = node_value.value
        sim_table[node_name] = value
    elseif (node_value isa FloatNode || node_value isa IntegerNode)
        value = get_value(node_value)
        sim_table[node_name] = value
    elseif (node_value isa StringNode)
        value = get_value(node_value)
    else
        error("Unsupported simulation declaration node value type: $(typeof(node_value))")
    end
    return value
end

function ir_run_simulation_node!(ast, sim_table)
    """
    Validates the run simulation node and then 
    generates IR.
    """

    run_sim_node = ast.value

    initial_state = run_sim_node.initialState
    parameters = run_sim_node.parameters
    rules = run_sim_node.rules
    types = run_sim_node.types
    steps = run_sim_node.steps

    str_initial_state = get_value(initial_state)
    str_param_val = get_value(parameters)
    str_rules_val = get_value(rules)
    str_types_val = get_value(types)
    str_steps_val = get_value(steps)

    is_value = sim_table[str_initial_state]

    params_value = sim_table[str_param_val]
    rules_value = sim_table[str_rules_val]
    types_value = sim_table[str_types_val]
    steps_value = sim_table[str_steps_val]

    checkpoint_frequency_value = nothing
    if run_sim_node.checkpointFrequency !== nothing
        cf_node = run_sim_node.checkpointFrequency
        checkpoint_frequency_value = cf_node isa IdentifierNode ?
            sim_table[get_value(cf_node)] : get_value(cf_node)
    end

    observables_value = nothing
    if run_sim_node.observables !== nothing
        str_obs_val = get_value(run_sim_node.observables)
        obs_candidate = sim_table[str_obs_val]
        if obs_candidate isa Union{Real, AbstractString} && tryparse(Int, string(obs_candidate)) !== nothing &&
                checkpoint_frequency_value === nothing
            # `RunSimulation(..., time, every)` with `every` a named integer.
            checkpoint_frequency_value = obs_candidate
        else
            observables_value = obs_candidate
        end
    end

    if checkpoint_frequency_value !== nothing
        # Numeric literals/declarations are stored as their source text (e.g. "60").
        parsed_cf = tryparse(Int, string(checkpoint_frequency_value))
        if parsed_cf === nothing || parsed_cf < 1
            error("RunSimulation checkpoint frequency must be a positive integer, got $(checkpoint_frequency_value)")
        end
        checkpoint_frequency_value = parsed_cf
    end

    save_name_value = nothing
    if run_sim_node.saveName !== nothing
        save_name_value = get_value(run_sim_node.saveName)
    end

    ir_run_sim = new_stage(get_value(ast.name), is_value, params_value, rules_value, types_value,
                           steps_value, checkpoint_frequency_value, observables_value, save_name_value)

end


function new_stage(name, initial_state, parameters, rules, types, steps, checkpoint_frequency,
                   observables, save_name; functions=nothing, grammar=nothing, active=true)
    # One simulation stage. `RunSimulation(...)` stages always run (`active`); a `Simulation(...)`
    # stage only runs when a `run(...)` statement names it. `replicates`/`domain_parallel`/`jobs`/
    # `seed` are filled in by that statement; `domain_parallel` stays `nothing` for stages without
    # one, which run once, in this process, with the serial runner.
    return Dict{String, Any}(
        "name" => name,
        "initial_state" => initial_state,
        "parameters" => parameters,
        "rules" => rules,
        "types" => types,
        "steps" => steps,
        "checkpoint_frequency" => checkpoint_frequency,
        "observables" => observables,
        "save_name" => save_name,
        "functions" => functions,
        "grammar" => grammar,
        "active" => active,
        "replicates" => 1,
        "domain_parallel" => nothing,
        "jobs" => nothing,
        "seed" => nothing,
    )
end

function parse_positive_int(what, v; allow_zero=false)
    parsed = tryparse(Int, string(v))
    if parsed === nothing || parsed < (allow_zero ? 0 : 1)
        error("$what must be a $(allow_zero ? "non-negative" : "positive") integer, got $(v)")
    end
    return parsed
end

function sim_int_value(node, sim_table)
    # An integer literal, or an identifier bound to an integer declaration.
    return node isa IdentifierNode ? sim_table[get_value(node)] : get_value(node)
end

function ir_grammar_node!(ast, sim_table)
    """
    `g : Grammar := Grammar(functions, rules, types [, observables])`: groups already declared
    SimulationFunctions / SimulationRules / SimulationTypes / SimulationObservables entries.
    """
    gname = get_value(ast.name)
    g = Dict{String, Any}("__grammar__" => true, "name" => gname, "rules" => nothing,
                          "types" => nothing, "observables" => nothing, "functions" => nothing)
    for part in ast.value.parts
        pname = get_value(part)
        if !haskey(SIM_DECL_KINDS, pname)
            error("Grammar $gname: '$pname' is not a SimulationRules / SimulationTypes / " *
                  "SimulationObservables / SimulationFunctions declaration (declare it before the Grammar)")
        end
        kind = SIM_DECL_KINDS[pname]
        g[kind] === nothing ||
            error("Grammar $gname: more than one $kind ('$pname'); a grammar has exactly one of each")
        g[kind] = sim_table[pname]
    end
    g["rules"] === nothing && error("Grammar $gname: no SimulationRules among its arguments")
    g["types"] === nothing && error("Grammar $gname: no SimulationTypes among its arguments")
    sim_table[gname] = g
    return g
end

function ir_simulation_config_node!(ast, sim_table)
    """
    `s : Simulation := Simulation(["save.bin",] grammar, state, params, time [, checkpoint])`.
    """
    node = ast.value
    sname = get_value(ast.name)

    gname = get_value(node.grammar)
    g = get(sim_table, gname, nothing)
    (g isa Dict && get(g, "__grammar__", false)) ||
        error("Simulation $sname: '$gname' is not a Grammar declaration")

    state_name = get_value(node.initialState)
    param_name = get_value(node.parameters)
    haskey(sim_table, state_name) || error("Simulation $sname: unknown initial state '$state_name'")
    haskey(sim_table, param_name) || error("Simulation $sname: unknown parameters '$param_name'")

    steps_value = node.steps isa IdentifierNode ? sim_table[get_value(node.steps)] : get_value(node.steps)

    checkpoint_frequency_value = nothing
    if node.checkpointFrequency !== nothing
        checkpoint_frequency_value = parse_positive_int("Simulation $sname checkpoint frequency",
                                                        sim_int_value(node.checkpointFrequency, sim_table))
    end

    save_name_value = node.saveName === nothing ? nothing : get_value(node.saveName)

    return new_stage(sname, sim_table[state_name], sim_table[param_name], g["rules"], g["types"],
                     steps_value, checkpoint_frequency_value, g["observables"], save_name_value;
                     functions=g["functions"], grammar=gname, active=false)
end

function apply_run_statements!(ir_run_sims, ir_runs, sim_table)
    """
    `run(sim, N, domain_parallel = ..., jobs = ..., seed = ...)`: mark the simulation they name as
    one to execute and attach the run options. `domain_parallel = true` selects the OpenMP
    domain-decomposition runner (compiled in only for that simulation); the default, and `false`,
    is the serial runner.
    """
    seen = Set{String}()
    for r in ir_runs
        sim_name = get_value(r.simulation)
        ix = findfirst(s -> s["name"] == sim_name, ir_run_sims)
        if ix === nothing
            known = join([s["name"] for s in ir_run_sims], ", ")
            error("run($sim_name, ...): no Simulation named '$sim_name' in this simulations section" *
                  (isempty(known) ? "" : " (known: $known)"))
        end
        sim_name in seen && error("run($sim_name, ...): '$sim_name' is run more than once")
        push!(seen, sim_name)

        stage = ir_run_sims[ix]
        stage["active"] = true
        stage["replicates"] = parse_positive_int("run($sim_name, N): N", sim_int_value(r.replicates, sim_table))
        stage["domain_parallel"] = r.domainParallel === true
        r.jobs === nothing ||
            (stage["jobs"] = parse_positive_int("run($sim_name, jobs = ...)", sim_int_value(r.jobs, sim_table)))
        if r.seed !== nothing
            seed = parse_positive_int("run($sim_name, seed = ...)", sim_int_value(r.seed, sim_table); allow_zero=true)
            seed <= typemax(UInt32) || error("run($sim_name, seed = $seed): seed must fit in 32 bits")
            stage["seed"] = seed
        end
        if stage["replicates"] == 1 && r.jobs !== nothing
            println("Warning: run($sim_name, 1, jobs = ...): jobs is ignored for a single run")
        end
    end
end


function ir_sim_declarations!(ast, sim_table)
    """
    Generates the simulations section
    """

    ir_run_sims = []
    ir_traces = []
    ir_measures = []
    ir_runs = []
    empty!(SIM_DECL_KINDS)
    while length(ast) > 0
        sim_decl_node = popfirst!(ast)

        if sim_decl_node isa TraceNode
            push!(ir_traces, sim_decl_node)
            continue
        end

        if sim_decl_node isa MeasureNode
            push!(ir_measures, sim_decl_node)
            continue
        end

        if sim_decl_node isa RunNode
            push!(ir_runs, sim_decl_node)
            continue
        end

        node_name = get_value(sim_decl_node.name)
        node_value = sim_decl_node.value

        if (node_value isa RunSimulationNode)
            ir_run_sim = ir_run_simulation_node!(sim_decl_node, sim_table)
            push!(ir_run_sims, ir_run_sim)
        elseif node_value isa SimulationConfigNode
            push!(ir_run_sims, ir_simulation_config_node!(sim_decl_node, sim_table))
        elseif node_value isa SimulationGrammarNode
            ir_grammar_node!(sim_decl_node, sim_table)
        else
            populate_sim_table(node_name, node_value, sim_table)
        end

    end

    apply_run_statements!(ir_run_sims, ir_runs, sim_table)
    for s in ir_run_sims
        s["active"] || println("Warning: Simulation '$(s["name"])' is never run; add `run($(s["name"]), 1)` to execute it")
    end
    # Attach `trace(sim, "rule", ...)` statements to the RunSimulation they name.
    for ir_run_sim in ir_run_sims
        ir_run_sim["trace_rules"] = String[]
    end
    for t in ir_traces
        sim_name = get_value(t.simulation)
        ix = findfirst(s -> s["name"] == sim_name, ir_run_sims)
        if ix === nothing
            known = join([s["name"] for s in ir_run_sims], ", ")
            error("trace($sim_name, ...): no RunSimulation named '$sim_name' in this simulations section" *
                  (isempty(known) ? "" : " (known: $known)"))
        end
        for r in t.rules
            rule_name = get_value(r)
            rule_name in ir_run_sims[ix]["trace_rules"] || push!(ir_run_sims[ix]["trace_rules"], rule_name)
        end
    end

    # `measure(sim, ..., "obs", ...)`: every named simulation must be a RunSimulation. Those
    # simulations keep per-checkpoint snapshots (my_results/<sim>/state_<step>.bin) so the
    # observables can be evaluated on them after the run.
    for ir_run_sim in ir_run_sims
        ir_run_sim["measured"] = false
    end
    measure_specs = []
    for m in ir_measures
        sims = String[]
        for s in m.simulations
            sim_name = get_value(s)
            ix = findfirst(r -> r["name"] == sim_name, ir_run_sims)
            if ix === nothing
                known = join([r["name"] for r in ir_run_sims], ", ")
                error("measure($sim_name, ...): no RunSimulation named '$sim_name' in this simulations section" *
                      (isempty(known) ? "" : " (known: $known)"))
            end
            ir_run_sims[ix]["measured"] = true
            sim_name in sims || push!(sims, sim_name)
        end
        obs_names = String[get_value(o) for o in m.observables]
        "*" in obs_names && (obs_names = String[])
        push!(measure_specs, Dict{String, Any}("simulations" => sims, "observables" => unique(obs_names)))
    end
    return ir_run_sims, measure_specs
end


function ir_simulation_section!(ast, sim_table)
    """
    Generates the simulation
    """

    sim_section_name = get_value(ast.name)
    sim_declarations = ast.declarations

    # Return the list of simulation executions needed.
    ir_run_sims, measure_specs = ir_sim_declarations!(sim_declarations, sim_table)
    return sim_section_name, ir_run_sims, measure_specs
end
