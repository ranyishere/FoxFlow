
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

# <grammar> ::= "grammar" <time-type> <symbol-name> "(" <grammar-signature>) "{" <initial-conditions> <rules> "}"
# <time-type> ::= "(discretetime)" | "(continuoustime)"
#
# <initial-conditions> ::= "initial_conditions" "{" <initial-condition-list> "}"
# <initial-condition-list> ::= <initial-condition> | <initial-condition> <initial-condition-list>
# <initial-condition> ::= <symbol-name> "=" <literal> ";"
#
# <rules> ::= <rules> | <rule> | <solve-rule>
#
# <grammar-signature> ::= <parameterized-types-left> "->" <parameterized-types-right>
#
# <rule> ::= <parameterized-types-left> "->" <parameterized-types-right> "with" <with-clause> ";"
# <solve-rule> ::= <parameterized-type-left> "->" <parameterized-type-right> "solving" <solve-clause> ";"
#
# <parameterized-types-left> ::= <parameterized-types> "," | <parameterized-type>
# <parameterized-types-right> ::= <parameterized-types> "," | <parameterized-type>
#
# <parameterized-types> ::= <parameterized-types> "," | <parameterized-type>
# <parameterized-type> ::= <symbol-name> "(" <parameters> ")" | "{" <parameterized-types> "|" <predicate> "}"
#
# <with-clause> ::= <function-type> "where" <predicate> ";" | <function-type> ";"
# <solve-clause> ::= <function-type> "where" <equations> ";" | <function-type> ";"
#
# <symbol-name> ::= [Aa-Zz] | [0-9] | _ | <symbol-name>
#
# <predicate> ::= <symbol-name> "!=" <symbol-name> 
#                 | <symbol-name> "<" <symbol-name> 
#                 | <symbol-name> ">" <symbol-name> 
#                 | <symbol-name> "==" <symbol-name>
#
# <function-type> ::= <symbol-name> "(" <symbol-name> ")"
#
# <literal> ::= [0-9] | [Aa-Zz] | <literal>


include("./tokens.jl")


abstract type Node end

struct TimeTypeNode <: Node
    token::Token
end

struct IdentifierNode <: Node
    token::Token
end

struct SymbolNode <: Node
    name::IdentifierNode
    parameters::Array{Node}
end

struct ParameterNode <: Node
    token::Array{Token}
end

struct GrammarSignatureNode <: Node
    token::Array{Token}
end

abstract type ModifyClauseNode <: Node
end

struct SolveClauseNode <: ModifyClauseNode
    name::Array{Token}
    clause::Array{Token}
end

struct WithClauseNode <: ModifyClauseNode
    name::Array{Token}
    clause::Array{Token}
end

struct RuleNode <: Node
    lhs::Array{Token}
    rhs::Array{Token}
    modify_clause::ModifyClauseNode
end

struct InitialConditionNode <: Node
    name::Token
    value::Token
end

struct InitialConditionListNode <: Node
    value::Array{InitialConditionNode}
end

struct GrammarNode <: Node
    time::TimeTypeNode
    name::IdentifierNode
    signature::GrammarSignatureNode
    initial_conditions::InitialConditionListNode
    rules::Array{RuleNode}
end


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
                                      InitialConditionListNode(initial_conditions), rules)

            push!(grammars, cur_grammar)
        end

    end

    grammars
end


function main()

    # TODO: Declare variable types
    source_code = """
        grammar discretetime symbol1 () {};
    """

    source_code_0 = """
        grammar discretetime symbol1 (Node(x) -> Node(y)) {
            Node(x) -> Node(y) with f(0);
        };
    """

    source_code_1 = """
        grammar discretetime symbol1 (Node(x) -> Node(y)) {
            Node(x) -> Node(y) with f(a) where a < 10;
            Node(a) -> Node(b) with g(b) where b < 10;
        };
    """

    source_code_2 = """
        grammar discretetime PredatorPrey () {
            rabbit -> null with alpha where alpha = 0.3;
            fox -> null with beta = 0.3;
            rabbit -> 2 * rabbit with gamma = 1.0;
            fox -> 2 * fox with f(fox, rabbit);
        };
    """

    tokens = tokenize(source_code_1)
    symbol_table = Dict()
    grammars = parse!(tokens, symbol_table)

end
# main()
