
"""
This parser should then take the tokens
that we created and use the metatheory library 
to then construct the correct semantics.
"""

# TODO: Add rule names
# TODO: Incorporate comments (#   #) //
# TODO: Line breaks in code?
# TODO: Incorporate/allow users to specify max time of simulation.

# Add types, Real, 
#


include("./tokens.jl")
include("./ast_nodes.jl")


function parse_grammar_signature!(tokens, symbol_table)
    """
    Parses Signature of Grammar
    """

    cur_token = popfirst!(tokens)

    grammar_signature = Token[]
    if isa(cur_token, LeftParenthesisToken)
        parenthesis_count = 1
    else
        println("Error expected grammar signature")
    end

    # If you see a token that is Left curly bracket break
    while parenthesis_count > 0

        cur_token = popfirst!(tokens)

        if isa(cur_token, LeftParenthesisToken)
            parenthesis_count += 1
        elseif isa(cur_token, RightParenthesisToken)
            parenthesis_count -= 1
        else
            push!(grammar_signature, cur_token)
        end

    end

    grammar_signature
end

function parse_with!(tokens, symbol_table)
    """
    Parses With Clause
    """

    cur_token = popfirst!(tokens)

    function_name = Token[]
    predicate = Token[]

    seen_where = false
    while isa(cur_token, PunctuationToken) == false
        # Parse function name
        if seen_where == false
            if isa(cur_token, WhereToken)
                seen_where = true
            else
                push!(function_name, cur_token)
            end
        else
            push!(predicate, cur_token)

            # Parse predicate
        end

        cur_token = popfirst!(tokens)

    end

    WithClauseNode(function_name, predicate)
end

function parse_solve!(tokens, symbol_table)
    """
    Parses Solve Clause
    """

    cur_token = popfirst!(tokens)

    # FIXME: Parsing predicate is not working.
    #        I don't think it can handle functions
    #        with multiple arguments.
    function_name = Token[]
    equation = Token[]

    seen_where = false
    while isa(cur_token, PunctuationToken)

        # Parse function name
        if seen_where == false
            if isa(cur_token, WhereToken)
                push!(function_name, cur_token)
            else
                seen_where = true
            end
        else
            # Parse predicate
            push!(equation, cur_token)
        end

        cur_token = popfirst!(tokens)
    end

    SolveClauseNode(function_name, equation)
end

function parse_rule!(passed_token, tokens, symbol_table)
    """
    Parse a single rule
    """

    cur_token = passed_token
    lhs = Token[]
    rhs = Token[]
    seen_arrow = false

    seen_solve = false
    seen_with = false

    # Actually should parse until seeing with or Solve Token
    while ((isa(cur_token, WithToken) != true) && (isa(cur_token, SolvingToken) != true))

        if isa(cur_token, RightArrowToken)
            seen_arrow = true
            cur_token = popfirst!(tokens)
            continue
        else

            if seen_arrow == false
                push!(lhs, cur_token)
            else
                push!(rhs, cur_token)
            end

        end
        cur_token = popfirst!(tokens)
    end

    if isa(cur_token, WithToken)
        #With Clause
        modify_clause = parse_with!(tokens, symbol_table)
    elseif isa(cur_token, SolvingToken)
        # Solve Clause
        modify_clause = parse_solve!(tokens, symbol_table)
    else
        modify_clause = ModifyClauseToken("", Token[])
    end

    lhs, rhs, modify_clause
end

function parse_initial_condition!(tokens, symbol_table)

    initial_condition_node = nothing

    cur_token = popfirst!(tokens)

    if !isa(cur_token, IdentifierToken)
        println("Expected an Identifier Token")
    end

    check_equal = popfirst!(tokens)

    if isa(check_equal, OperatorToken)

        if check_equal.value != "="
            println("Expected the equal operator")
        else
            condition_value = popfirst!(tokens)
            if isa(condition_value, LiteralToken)
                initial_condition_node = InitialConditionNode(cur_token,
                                                              condition_value)
                symbol_table[cur_token.value] = condition_value.value
            else
                println("Expected Literal Token received: $condition_value")
            end
        end
    end

    cur_token = popfirst!(tokens)

    if !isa(cur_token, PunctuationToken)
        println("Expected ; at end of initial condition declaration")
        println("Got: ", cur_token)
    end

    initial_condition_node
end

function parse_initial_condition_list!(tokens, symbol_table)

    # Begins Grammar Context
    cur_token = popfirst!(tokens)
    if !isa(cur_token, LeftBracketToken)
        println("Error expected LeftBracket to begin grammar context.")
    end

    initial_condition_list_nodes = []
    cur_token = popfirst!(tokens)

    if !isa(cur_token, InitialConditionToken)
        println("Error expected initial conditions first.")
    end

    cur_token = popfirst!(tokens)

    # Initial Condition Opening Bracket
    if !isa(cur_token, LeftBracketToken)
        println("Error initial condition expected opening left bracket.")
    end

    parenthesis_count = 1

    # TODO FIX parenthesis counting to exit initial_condition block
    while (parenthesis_count != 0 && length(tokens) != 0)

        look_ahead = tokens[1]

        if isa(look_ahead, RightBracketToken)
            # End initial conditions
            parenthesis_count -= 1
        else
            cur_initial_node = parse_initial_condition!(tokens, symbol_table)
            push!(initial_condition_list_nodes, cur_initial_node)
        end
    end

    if (parenthesis_count != 0)
        println("Error expected closing parenthesis for initial conditions")
    end

    cur_token = popfirst!(tokens)

    initial_condition_list_nodes
end

function parse_rules!(tokens, symbol_table)
    """
    Parse Rules
    """

    cur_token = popfirst!(tokens)

    bracket_counter = 1

        # rule_tokens = []
    rules = []
    while bracket_counter > 0

        cur_token = popfirst!(tokens)

        if isa(cur_token, RightBracketToken)
            bracket_counter -= 1
        elseif isa(cur_token, LeftBracketToken)
            bracket_counter += 1
        else

            lhs, rhs, modify_clause = parse_rule!(cur_token, tokens, symbol_table)

            rule_node = RuleNode(lhs,rhs, modify_clause)

            push!(rules, rule_node)
        end

    end
    # rule_tokens
    # RuleNode(rule_tokens)
    rules
end


function parse!(tokens, symbol_table)
    """
    Parser the Entire Grammar
    """

    i = 0
    grammars = GrammarNode[]

    # Keep runing until we run out of Grammars to parse
    while length(tokens) > 0
        cur_token = popfirst!(tokens)

        if isa(cur_token, GrammarToken)

            time_node = popfirst!(tokens)
            name_node = popfirst!(tokens)

            signature_node = parse_grammar_signature!(tokens,
                                                      symbol_table)

            # Parse initial condition
            initial_conditions = parse_initial_condition_list!(tokens,
                                                               symbol_table)

            # TODO: Need to move curly bracket outside of parse_rules.


            rules = parse_rules!(tokens, symbol_table)

            cur_grammar = GrammarNode(TimeTypeNode(time_node), IdentifierNode(name_node),
                                      GrammarSignatureNode(signature_node), 
                                      InitialConditionListNode(initial_conditions),
                                      rules)

            push!(grammars, cur_grammar)
        end

    end

    grammars
end

function parse_type()
    if isa(cur_token, IdentifierToken)
        type_name_inst = cur_token
    end
end

# if !isa(cur_token, 
function parse_type_class!(tokens, symbol_table)

    cur_token = popfirst!(tokens)
    class_name = nothing
    class_parameter = Token[]

    if isa(cur_token, IdentifierToken)
        class_name = IdentifierNode(cur_token)
        # Should be class name
    else
        println("Expected class name")
    end

    cur_token = popfirst!(tokens)

    if isa(cur_token, PunctuationToken)
        if (cur_token.value == ';')
        end

    # Class parameter
    elseif isa(cur_token, LeftAngleBracketToken)

        while !isa(cur_token, RightAngleBracketToken)

            cur_token = popfirst!(tokens)

            if isa(cur_token, IdentifierToken)
                # its a type.
            elseif isa(cur_token, PunctuationToken)

                tok_val = cur_token.value

                if cur_token.value != ','
                    println("Error expected , got $tok_val")
                elseif isa(cur_token, IdentifierToken)

                    tmp_class_param = IdentifierNode(IdentifierToken)
                    push!(class_parameter, tmp_class_param)

                else
                    continue
                end

            end
        end
    end

   TypeNode(class_name, class_parameter)
end

function parse_type!(tokens, symbol_table)

    # types = TypeNode[]
    type_name_inst = nothing
    parameter = Node[]
    # value = TypeNode[]

    cur_token = popfirst!(tokens)

    # Check for instance name
    if isa(cur_token, IdentifierToken)
        type_name_inst = IdentifierNode(cur_token)
    else
        println("Error expected Identifier Token")
    end

    # Check for Parameter
    
    lookahead = tokens[1]
    if isa(lookahead, LeftAngleBracketToken)
        cur_token = popfirst!(tokens)

        brack = 1
        while brack > 0 && length(tokens) > 0
            cur_token = popfirst!(tokens)
            if isa(cur_token, IdentiferToken)
                push!(
                      parameter,
                      IdentiferNode(cur_token)
                     )
            elseif isa(cur_token, RightAngleBracketToken)
                brack -= 1
            end
        end

        if brack > 0
            println("Error missing parameter closing bracket")
        end
    end

    # Check for type
    type = nothing
    type_name = nothing
    type_parameter = Node[]

    lookahead = tokens[1]
    if isa(lookahead, DoubleColonToken)
        popfirst!(tokens)

        cur_token = popfirst!(tokens)

        if isa(cur_token, IdentifierToken)
            type_name = IdentifierNode(cur_token)
        else
            println("Expected type name")
        end

        lookahead = tokens[1]

        # cur_token = popfirst!(tokens)

        # Type class parameter
        println("lookahead angle brack: $lookahead")
        if isa(lookahead, LeftAngleBracketToken)

            cur_token = popfirst!(tokens)
            brack = 1
            while brack > 0 && length(tokens) > 0

                cur_token = popfirst!(tokens)
                if isa(cur_token, IdentifierToken)
                    push!(
                          type_parameter,
                          IdentifierNode(cur_token)
                         )
                elseif isa(cur_token, RightAngleBracketToken)
                    brack -= 1
                end

            end

            if brack > 0
                println("Error missing closing right angle bracket")
            end

        end

    else
        println("Error expected type signature")
    end

    # Check for value assignment
    value = nothing
    # cur_token = popfirst!(tokens)
    lookahead = tokens[1]
    if isa(lookahead, OperatorToken)

        # println("cur_token.position.value: ", cur_token.position.value)
        cur_token = popfirst!(tokens)
        if (cur_token.position.value == "=")

            cur_token = popfirst!(tokens)

            if !isa(cur_token, LeftBracketToken)

                value = cur_token

            # NOTE:  It's assigning a list of types
            elseif isa(cur_token, LeftBracketToken)
                value = parse_list_types!(tokens, symbol_table)

            end

        else
            println("Error expected assignment =")
        end

    end

    lookahead = tokens[1]
    if isa(lookahead, PunctuationToken)
        cur_token = popfirst!(tokens)
        if (cur_token.position.value != ";")
            println("Error Expected end of type creation")
        end
    else
        println("Error Expected end of type creation else given: $lookahead")
    println("Tokens: $tokens")
    end
    
    type_class = TypeClassNode(type_name,
                               ParameterNode(type_parameter))

    TypeInstanceNode(type_name_inst,
                     ParameterNode(parameter),
                     type_class,
                     value)
end

function parse_list_types!(tokens, symbol_table)
    bracket_count = 1
    type_instances = TypeInstanceNode[]
    while bracket_count > 0
        lookahead = tokens[1]
        # cur_token = popfirst!(tokens)
        if isa(lookahead, RightBracketToken)
            popped = popfirst!(tokens)
            bracket_count -= 1;
        elseif isa(lookahead, LeftBracketToken)
            popfirst!(tokens)
            # type_instance = parse_list_type!(tokens, symbol_table)
            # push!(type_instances, type_instance)
            bracket_count += 1;
        else
            type_instance = parse_type!(tokens, symbol_table)
            push!(type_instances, type_instance)

        end
    end

    cur_token = popfirst!(tokens)
    if (cur_token.position.value != ";")
        println("expected semicolon")
    end

    type_instances
end

function parse_types_section!(tokens, symbol_table)

    section_name = nothing
    type_instances = TypeInstanceNode[]
    while length(tokens) > 0

        cur_token = popfirst!(tokens)

        # types keyword token
        if isa(cur_token, TypeSectionToken)

            cur_token = popfirst!(tokens)

            # Identifier
            if isa(cur_token, IdentifierToken)
                section_name = IdentifierNode(cur_token)
            else
                println("Error expected section name got $cur_token")
            end

            cur_token = popfirst!(tokens)

            if isa(cur_token, LeftBracketToken)
                type_instances = parse_list_types!(tokens, symbol_table)
            else
                println("Error expected curly brackets.")
            end

        end ## if end

        ##
    end # while end

    TypeSectionNode(section_name, type_instances)
end
