
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
          node_value isa SimulationObservablesNode ||
          node_value isa SimulationStateNode)
        value = get_value(node_value.value.filepath)
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

    observables_value = nothing
    if run_sim_node.observables !== nothing
        str_obs_val = get_value(run_sim_node.observables)
        observables_value = sim_table[str_obs_val]
    end

    ir_run_sim = Dict(
        "initial_state" => is_value,
        "parameters" => params_value,
        "rules" => rules_value,
        "types" => types_value,
        "steps" => steps_value,
        "observables" => observables_value
    )

end


function ir_sim_declarations!(ast, sim_table)
    """
    Generates the simulations section
    """

    ir_run_sims = []
    while length(ast) > 0
        sim_decl_node = popfirst!(ast)
        node_name = get_value(sim_decl_node.name)
        node_value = sim_decl_node.value

        if (node_value isa RunSimulationNode)
            ir_run_sim = ir_run_simulation_node!(sim_decl_node, sim_table)
            push!(ir_run_sims, ir_run_sim)
        else
            populate_sim_table(node_name, node_value, sim_table)
        end

    end
    return ir_run_sims
end


function ir_simulation_section!(ast, sim_table)
    """
    Generates the simulation
    """

    sim_section_name = get_value(ast.name)
    sim_declarations = ast.declarations

    # Return the list of simulation executions needed.
    ir_run_sims = ir_sim_declarations!(sim_declarations, sim_table)
    return sim_section_name, ir_run_sims
end
