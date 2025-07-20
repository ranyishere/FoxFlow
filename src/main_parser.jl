include("./tokens.jl")
include("./lexer.jl")
include("./ast_nodes.jl")


# TODO: Support expression
function tokenize_file(file_name)
    lines_tokens = []
    open(file_name) do f
        line_no = 0
        while !eof(f)

            cur_line = readline(f)

            if cur_line != ""
                if cur_line[1] == '#'
                    line_no += 1
                    println("continuing")
                    continue
                end
            end

            line_token = tokenize(cur_line, line_no+1)
            # println("line_token :", line_token)

            line_no += 1

            if line_token != Token[]
                lines_tokens = [lines_tokens; line_token]
                push!(lines_tokens, EndLineToken())
            end

        end
    end
    lines_tokens
end


# Helper function: Peek at the next token without consuming it
function lookahead(tokens)
    if !isempty(tokens)
        return tokens[1]
    else
        throw("Lookahead requested on empty token list")
    end
end

# Parse a symbol name (identifier)
function parse_symbol_name!(tokens)
    cur_token = popfirst!(tokens)

    if isa(cur_token, IdentifierToken)
        return IdentifierNode(cur_token)
    else
        throw("Expected Identifier Token, got: $cur_token")
    end
end

# Parse symbol parameters `<<param>>` or empty
function parse_symbol_parameters!(tokens)
    params = []

    lookahead_token = lookahead(tokens)
    # println("oof $lookahead_token")
    
    if !isempty(tokens) && isa(lookahead(tokens), LeftAngleBracketToken)
        popfirst!(tokens) # Consume `<<`

        while !isempty(tokens) && !isa(lookahead(tokens), RightAngleBracketToken)

            # TODO: Handle left and right parenthesis '(' and ')' for grouping parameters
            # if isa(lookahead(tokens), IntegerToken)
                # push!(params, IntegerNode(popfirst!(tokens)))
            # elseif isa(lookahead(tokens), FloatToken)
                # push!(params, FloatNode(popfirst!(tokens)))
            # Handling nested parameters
            if isa(lookahead(tokens), LeftAngleBracketToken)
                params = [
                          params;
                          parse_symbol_parameters!(tokens)
                        ]
                exit(0)
            else

                cur_expr = parse_expression!(tokens)

                # cur_token = popfirst!(tokens)
                # throw("Error got unrecognized token: $cur_token")
            end

            if !isempty(tokens) && isa(lookahead(tokens), PunctuationToken)
                popfirst!(tokens) # Consume `,`
            end

        end

        if isempty(tokens) || !isa(popfirst!(tokens), RightAngleBracketToken)
            throw("Expected `>>` to close parameter list")
        end
    end
    
    return ParameterNode(params)
end

# Parse a type signature list (e.g., `A -> B -> C`)
function parse_type_signature_list!(tokens)

    typename = nothing
    if isa(lookahead(tokens), IdentifierToken)
        type_name = parse_symbol_name!(tokens)
    elseif typeof(lookahead(tokens)) in [IntegerToken, FloatToken]
        type_name = IdentifierNode(popfirst!(tokens))
    else
        cur_token = popfirst!(tokens)
        throw("Error couldnt determine type got: $cur_Token")
    end

    if !isempty(tokens) && isa(lookahead(tokens), RightArrowToken)
        popfirst!(tokens)  # Consume `->`
        return TypeClassNode(type_name, parse_type_signature_list!(tokens))
    else
        return TypeClassNode(type_name, ParameterNode([]))
    end
end

function parse_type_assignment!(tokens)
    """
    Parse Type Assignment
    """

    symbol_name = parse_symbol_name!(tokens)
    symbol_parameters = parse_symbol_parameters!(tokens)

    if !isa(popfirst!(tokens), SingleColonToken)
        throw("Expected `:` after symbol name")
    end

    type_signature_list = parse_type_signature_list!(tokens)

    if !isempty(tokens) && isa(lookahead(tokens), DefineToken)

        popfirst!(tokens)  # Consume `:=`

        if !isempty(tokens) && isa(lookahead(tokens), LeftBracketToken)
            popfirst!(tokens) # Consume `{`
            type_declarations = parse_type_declaration_list!(tokens)

            println("type_declarations: ", type_declarations)
            if isempty(tokens) || !isa(popfirst!(tokens), RightBracketToken)
                throw("Expected `}` to close type declaration list")
            end

            # println("type_declarations: ", type_declarations)
            return TypeInstanceNode(symbol_name, symbol_parameters, type_signature_list, type_declarations)
        else
            println("tokens: ", tokens)
            #TODO: Support expression

            # Assume it's a literal/expression
            tmp = parse_expression!(tokens)
            println("tmp: ", tmp)
            exit(0)
            literal = popfirst!(tokens)  

            return TypeInstanceNode(symbol_name, symbol_parameters, type_signature_list, literal)
        end
    end


    return TypeInstanceNode(symbol_name, symbol_parameters, type_signature_list, nothing)
end


function parse_type_update!(tokens)

    symbol_name = parse_symbol_name!(tokens)
    symbol_parameters = parse_symbol_parameters!(tokens)

    if !isa(popfirst!(tokens), SingleColonToken)
        throw("Expected `:` after symbol name")
    end

    type_signature_list = parse_type_signature_list!(tokens)

    if !isempty(tokens) && isa(lookahead(tokens), EqualToken)

        popfirst!(tokens)  # Consume `:=`

        if !isempty(tokens) && isa(lookahead(tokens), LeftBracketToken)
            popfirst!(tokens) # Consume `{`
            type_declarations = parse_type_declaration_list!(tokens)

            if isempty(tokens) || !isa(popfirst!(tokens), RightBracketToken)
                throw("Expected `}` to close type declaration list")
            end

            return TypeInstanceNode(symbol_name, symbol_parameters, type_signature_list, type_declarations)
        else
            #TODO: Support expression
            expression = []
            while !isa(lookahead(tokens), EndLineToken) && !isa(lookahead(tokens), RightBracketToken)
                literal = popfirst!(tokens)  # Assume it's a literal/expression
                push!(expression,literal)
            end

            return TypeInstanceNode(symbol_name, symbol_parameters, type_signature_list, expression)

        end
    end


    return TypeInstanceNode(symbol_name, symbol_parameters, type_signature_list, nothing)
end



# Parse a type declaration (non-assignment case)
function parse_type_declaration!(tokens)

    type_name = parse_symbol_name!(tokens)
    symbol_parameters = parse_symbol_parameters!(tokens)

    if !isa(popfirst!(tokens), SingleColonToken)
        throw("Expected `:` in type declaration")
    end

    type_signature_list = parse_type_signature_list!(tokens)

    return TypeInstanceNode(type_name, symbol_parameters, type_signature_list, nothing)
end

# Parse a list of type declarations recursively
function parse_type_declaration_list!(tokens)
    declarations = []

    while !isempty(tokens) && isa(lookahead(tokens), IdentifierToken)
        push!(declarations, parse_type_declaration!(tokens))

        if !isempty(tokens) && isa(lookahead(tokens), EndLineToken)
            popfirst!(tokens)  # Consume newline
        end
    end

    return declarations
end

# Parse a type section `{ ... }` into an AST node
function parse_types_section!(tokens)

    if !isa(lookahead(tokens), TypeSectionToken)
        throw("Expected `types` section")
    end

    popfirst!(tokens)  # Consume `types`

    section_name = parse_symbol_name!(tokens)

    if !isa(popfirst!(tokens), LeftBracketToken)
        throw("Expected `{` to open type section")
    end

    if isa(lookahead(tokens), EndLineToken)
        popfirst!(tokens)
    end
    type_declarations = parse_type_declaration_list!(tokens)

    if !isa(popfirst!(tokens), RightBracketToken)
        throw("Expected `}` to close type section")
    end

    return TypeSectionNode(section_name, type_declarations)
end

function parse_factor!(tokens)
    # <factor> ::= <symbol>
    # | <literal>
    # | "(" <expression> ")"
    lookahead_token = lookahead(tokens)
    # exit(0)
   

    if lookahead(tokens) isa LeftParenthesisToken
        popfirst!(tokens)

        expr = parse_expression!(tokens)

        cur_look = lookahead(tokens)

        if !isa(cur_look, RightParenthesisToken)
            throw("Expected RightParenthesisToken got: $cur_look")
        end
        return expr
    end

    if lookahead(tokens) isa IdentifierToken
        return parse_symbol_name!(tokens)
    elseif lookahead(tokens) isa IdentifierToken
        return parse_symbol_name!(tokens)
        # push!(params, parse_symbol_name!(tokens))
    elseif isa(lookahead(tokens), LiteralToken)
        return LiteralNode(popfirst!(tokens))
        # push!(params, LiteralNode(popfirst!(tokens)))
    elseif isa(lookahead(tokens),
                       LeftParenthesisToken) || isa(lookahead(tokens),
                                                    RightParenthesisToken)
        popfirst!(tokens)
    elseif isa(lookahead(tokens), IdentifierToken)
        return IdentifierNode(IdentifierToken(popfirst!(tokens)))
        # push!(params, IdentifierToken(popfirst!(tokens)))
    end

end

function parse_term!(tokens)
    # <term> ::= <factor>
            # | <term> "*" <factor>
            # | <term> "/" <factor>
    factor = parse_factor!(tokens)

    # println("factor: ", factor)
    # println("lookahead(tokens): ", lookahead(tokens))
    # exit(0)

    if lookahead(tokens) isa AsteriskToken
        popfirst!(tokens)
        return MultiplyNode(factor, parse_term!(tokens))
    elseif lookahead(tokens) isa SlashToken
        popfirst!(tokens)
        return DivideNode(factor, parse_term!(tokens))
    elseif lookahead(tokens) isa RightAngleBracketToken
        return factor
    else
        # throw("Expected operator got: $(lookahead(tokens))")
        return factor
    end

end


# <expression> ::= <term>
                # | <expression> "+" <term> 
                # | <expression> "-" <term>

function parse_expression!(tokens)

    # <expression> ::= <term>
    # | <expression> "+" <term> 
    # | <expression> "-" <term>

    term = parse_term!(tokens)

    if term == nothing
        println("tokens left: ", tokens)
        throw("Error parsing term")
    end

    print("term: ", term)
    # println(lookahead_token.

    if lookahead(tokens) isa PlusToken
        popfirst!(tokens)
        tmp = PlusNode(term, parse_expression!(tokens))
        print("tmp: ", tmp)
        return tmp

    elseif lookahead(tokens) isa MinusToken
        tmp = parse_expression!(tokens)
        popfirst!(tokens)
        print("tmp: ", tmp)
        println("tokens minus: ", tokens)
        exit(0)
        return MinusNode(term, tmp)
    else
        # println("returning term")
        # throw("PArsing expression got: $(lookahead(tokens))")
        return term
    end

    # println("term: ", term)
    # exit(0)
    # term
end



function parse_function_type!(tokens)
    # println("--- in func tokens: ", tokens)

    function_name = parse_symbol_name!(tokens)

    cur_token = popfirst!(tokens)
    if !isa(cur_token, LeftParenthesisToken)
        throw("Expected LeftParenthesisToken got: $cur_token")
    end

    function_args = []
    while !isa(lookahead(tokens), RightParenthesisToken)
        # println("lookahead(tokens): ", lookahead(tokens))

        # TODO: Have this handle expression
        cur_arg = popfirst!(tokens)
        # cur_arg = parse_symbol_name!(tokens)
        push!(function_args, cur_arg)
        if isa(lookahead(tokens), CommaToken)
            popfirst!(tokens)
        end
    end

    cur_token = popfirst!(tokens)
    if !isa(cur_token, RightParenthesisToken)
        throw("Expected RightParenthesisToken got: $cur_token")
    end

    FunctionNode(function_name, function_args)
end

function parse_type_assignment_list!(tokens)
    """
    Parse Type Assignment List
    """

    type_assignments = []

    while isa(lookahead(tokens), EndLineToken)
        popfirst!(tokens)
    end

    cur_type_assignment = parse_type_assignment!(tokens)

    push!(type_assignments, cur_type_assignment)

    if isa(lookahead(tokens), EndLineToken)
        tmp = parse_type_assignment_list!(tokens)
        type_assignments = [type_assignments;tmp]
    end

    return type_assignments
end

function parse_type_update_list!(tokens)

    type_updates = []
    if isa(lookahead(tokens), RightBracketToken)
        type_updates
    else
        while isa(lookahead(tokens), EndLineToken)
            popfirst!(tokens)
        end

        cur_type_update = parse_type_update!(tokens)

        push!(type_updates, cur_type_update)

        if isa(lookahead(tokens), EndLineToken)
            popfirst!(tokens)
            tmp = parse_type_update_list!(tokens)
            type_updates = [type_updates;tmp]
        end

        return type_updates
    end
end


function parse_where_clause!(tokens)
    """
    Parse Where Clause
    """

    WhereClauseNode(
        parse_type_update_list!(tokens)
   )
end

function parse_with_clause!(tokens)
    """
    Parse With Clause
    """

    cur_token = popfirst!(tokens)
    if !isa(cur_token, LeftParenthesisToken)
        println("Expected LeftParenthesis got $cur_token")
    end
    # println(
    # "----------- parsing with clause -----------: ",
    # tokens
    # )
    function_node = parse_function_type!(tokens)

    cur_token = popfirst!(tokens)
    if !isa(cur_token, RightParenthesisToken)
        println("Expected RightParenthesis got $cur_token")
    end

    cur_token = popfirst!(tokens)
    if !isa(cur_token, WhereToken)
        throw("Expected WhereToken got $cur_token")
    end

    cur_token = popfirst!(tokens)

    if !isa(cur_token, LeftBracketToken)
        throw("Expected LeftBracketToken got $cur_token")
    end

    where_clause = parse_where_clause!(tokens)

    while isa(lookahead(tokens), EndLineToken)
        popfirst!(tokens)
    end

    # if !isa(cur_token, RightBracketToken)
        # throw("Expected RightBracketToken got $cur_token")
    # end

    WithClauseNode(function_node, where_clause)
end

function parse_parameterized_type!(tokens)
    popfirst!(tokens)
end

function parse_rule!(tokens)
    """
    Parse Rule
    """

    println("parsing rule")
    rule_name = parse_symbol_name!(tokens)

    cur_token = popfirst!(tokens)

    if !isa(cur_token, DefineToken)
        throw("Expected assignment Token ':=' got $cur_token")
    end

    lhs = []
    while !isempty(tokens) && !isa(lookahead(tokens), LeftAngleBracketToken) && !isa(lookahead(tokens), RightArrowToken)
        push!(lhs, popfirst!(tokens))
    end

    left_rule_parameter_node = ParameterNode([])
    if isa(lookahead(tokens), LeftAngleBracketToken)
        left_rule_parameter_node = parse_symbol_parameters!(tokens)
    end

    cur_token = popfirst!(tokens)

    # Gets rid of new lines
    while isa(lookahead(tokens), EndLineToken)
        popfirst!(tokens)
    end

    cur_token = popfirst!(tokens)

    if isempty(tokens) || !isa(cur_token, RightArrowToken)
        throw("Expected `->` in rule definition got $cur_token")
    end

    rhs = []
    # while !isempty(tokens) && !isa(lookahead(tokens), WithToken) && !isa(lookahead(tokens), SolvingToken)
        # push!(rhs, popfirst!(tokens))
    # end

    while !isempty(tokens) && !isa(lookahead(tokens), LeftAngleBracketToken)
        push!(rhs, popfirst!(tokens))
    end

    right_rule_parameter_node = ParameterNode([])
    if isa(lookahead(tokens), LeftAngleBracketToken)
        # println("Parsing symb param: ", tokens)
        # exit(0)
        right_rule_parameter_node = parse_symbol_parameters!(tokens)
    end
    println("done parsing right param node")

    # TODO: check parameters
    modify_clause = nothing
    if !isempty(tokens) && (isa(lookahead(tokens), WithToken) || isa(lookahead(tokens), SolvingToken))

        # Consume `with` or `solving`
        clause_token = popfirst!(tokens) 
        clause_content = []

        
        # paren_count = 0
        # while !isempty(tokens) && !isa(lookahead(tokens), EndLineToken)
            # push!(clause_content, popfirst!(tokens))
        # end

        if isa(clause_token, WithToken)
            println("clause_token: $clause_token \n clause_content: $clause_content")
            with_clause = parse_with_clause!(tokens)
            println("with_clause: ", with_clause)
            modify_clause = with_clause
            # modify_clause = WithClauseNode(clause_token, clause_content)

        else
            modify_clause = SolveClauseNode(clause_token, clause_content)
        end
    end

    # Find where token
    return RuleNode(
                    rule_name, 
                    lhs, left_rule_parameter_node,
                    rhs, right_rule_parameter_node,
                    modify_clause
                   )
end

# Parse the entire file into AST nodes
function parse_file!(tokens)
    ast_nodes = []

    while !isempty(tokens)

        cur_token = lookahead(tokens)
        # println("cur_token: ", cur_token)
        # println("typeof(cur_token): ", typeof(cur_token))

        if isa(cur_token, TypeSectionToken)
            push!(ast_nodes, parse_types_section!(tokens))
        elseif isa(cur_token, ParameterSectionToken)
            push!(ast_nodes, parse_parameters_section!(tokens))
        elseif isa(cur_token, RuleSectionToken)
            # push!(ast_nodes, parse_rule!(tokens))
            push!(
                  ast_nodes,
                  parse_rules_section!(tokens)
              )

        elseif isa(cur_token, EndLineToken)
            popfirst!(tokens)
        else
            throw("Unexpected token: $cur_token")
        end
    end


    return ast_nodes
end

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

    if isa(lookahead(tokens), EndLineToken)
        popfirst!(tokens)
    elseif isa(lookahead(tokens), RightBracketToken)
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

function parse_rules_list!(tokens)
    """
    parse rules list
    """

    rules = []
    while isa(lookahead(tokens),EndLineToken)
        popfirst!(tokens)
    end

    rule = parse_rule!(tokens)

    print("rule: ", rule)
    push!(rules, rule)

    while isa(lookahead(tokens),EndLineToken)
        popfirst!(tokens)
    end

    println("====== tokens ====== : ", tokens)

    # Still another rule
    if isa(lookahead(tokens), IdentifierToken)
        rules = [rules;parse_rules_list!(tokens)]
    end

    println("returning rule")

    return rules
end

function parse_rules_section!(tokens)
    """
    Rule Section
    """

    cur_token = popfirst!(tokens)
    section_name = nothing
    rules_list = []
    if isa(cur_token, RuleSectionToken)
        section_name = parse_symbol_name!(tokens)

        cur_token = popfirst!(tokens)
        if !isa(cur_token, LeftBracketToken)
            throw("Expected Left Bracket Token got $cur_token ")
        end

        rules_list = parse_rules_list!(tokens)
        cur_token = popfirst!(tokens)

        if !isa(cur_token, RightBracketToken)
            throw("Expected Right Bracket Token got: $cur_token")
        end

    else
        throw("Error expected Rule Section Token")
    end

    while isa(lookahead(tokens), EndLineToken)
        popfirst!(tokens)
    end

    if isa(lookahead(tokens), RightBracketToken)
        popfirst!(tokens)
    end

    RuleSectionNode(section_name, rules_list)
end

# Entry point for parsing
function main()

    # tokens = tokenize_file("tests/test_types_3.fflow")

    # tokens = tokenize_file("tests/rules/test_rule_1.fflow")

    tokens = tokenize_file("tests/expressions/expressions.fflow")

    # for token in tokens
        # println(token)
    # end
    # exit(0)
    ast = parse_file!(tokens)
end

main()
