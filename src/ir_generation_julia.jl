"""
Generates Intermediate Code from AST.
"""

function initial_condition_to_state(initial_conditions, symbol_table)
    """
    Convert initial condition to state dictionary
    """

    initial_state = Dict()
    # FIXME: Assuming value is an integer
    for condition in initial_conditions.value
        initial_state[condition.name] = condition.value
    end
    # println("initial_state: ", initial_state)
    initial_state
end

function rules_to_rates(rules, symbol_table)
    """
    Converts rules to delta matrix and rate table
    """

    # Maybe ignore null?
    # for now no, we can keep track of death.

    # Note add +1 for null (Stochiometry matrix)
    delta_matrix = zeros(
                         Int, length(rules),
                         length(symbol_table)
                        )

    # Rate Matrix
    rate_matrix = zeros(
                        Int, length(rules)
                        # length(symbol_table)
                       )

    # Convert symbol to positon in array
    code_book = Dict()
    for (ix, symbol) in enumerate(collect(symbol_table))
        code_book[symbol[1]] = ix
    end

    # println("code_book: ", code_book)
    # println("delta_matrix before: ", delta_matrix)

    for (ix, rule) in enumerate(rules)

        println("ix: ", ix)

        for lhs_symbol in rule.lhs
            cb_loc = code_book[lhs_symbol.value]
            delta_matrix[ix,cb_loc] -= 1
        end

        for rhs_symbol in rule.rhs
            cb_loc = code_book[rhs_symbol.value]
            delta_matrix[ix,cb_loc] += 1
        end

	println("modify: ", rule.modify_clause)

    end

    # println("delta_matrix after: ", delta_matrix)
    println("rate_matrix: ", rate_matrix)

    delta_matrix, rate_matrix
end


function generate_im(grammars, symbol_table)
    """
    Generate Intermediate Code Representation
    """

    println("symbol_table: ", symbol_table)

    for grammar in grammars
        initial_state = initial_condition_to_state(grammar.initial_conditions, symbol_table)
        delta_matrix, rate_matrix = rules_to_rates(grammar.rules, symbol_table)
    end

    # println("grammar[1].InitialConditionListNode: ", 
            # grammar[1].initial_conditions)

    # for rule in grammar[1].rules
        # println("rule: ", rule)
    # end

end
