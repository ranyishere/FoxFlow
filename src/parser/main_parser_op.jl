import ..Tokens:
    Token, IdentifierToken, IntegerToken, FloatToken, EndLineToken,
    LeftParenthesisToken, RightParenthesisToken, LeftAngleBracketToken,
    RightAngleBracketToken, SingleColonToken, DefineToken, LeftBracketToken,
    RightBracketToken, PunctuationToken, TypeSectionToken, ParameterSectionToken,
    RuleSectionToken, SlashToken, AsteriskToken, PlusToken, MinusToken, RightArrowToken, WithToken, WhereToken, SolvingToken, EqualToken, CommaToken, EdgeToken

import ..AstNodes:
    Node, IdentifierNode, ParameterNode, TypeInstanceNode, TypeClassNode,
    TypeSectionNode, ParameterSectionNode, RuleSectionNode, RuleNode,
    WhereClauseNode, WithClauseNode, SolveClauseNode, CallNode,
    GroupNode, BinaryOpNode, FunctionNode, FloatNode, IntegerNode, UndirectedTypeEdgeNode

const BUILT_IN_FUNC = [
    "heaviside", "sqrt", "normal_distr",
    "uniform_distr", "cos", "sin", "inverse", "pow"
]

const PRECEDENCE = Dict(
    "+" => 1,
    "-" => 1,
    "*" => 2,
    "/" => 2,
)

function expect_token!(tokens, ::Type{T}) where T <: Token

    token = popfirst!(tokens)

    if !(token isa T)
        error("Expected $(T), got $(typeof(token))")
    end
    return token
end


# TODO: Support expression
function tokenize_file(file_name)

    println("Tokenizing file: $file_name")

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

function tokenize_string(code)
    check = tokenize(code, 0)
    # println("check: ", check)
    check
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
    
    if !isempty(tokens) && isa(lookahead(tokens), LeftAngleBracketToken)

        popfirst!(tokens) # Consume `<<`

        while !isempty(tokens) && !isa(lookahead(tokens), RightAngleBracketToken)

            # Handling nested parameters
            if isa(lookahead(tokens), LeftAngleBracketToken)
                nested_params = parse_symbol_parameters!(tokens)
                push!(params, nested_params)  # Add nested params to params array

            elseif isa(lookahead(tokens), LeftParenthesisToken) || isa(lookahead(tokens), RightParenthesisToken)
                popfirst!(tokens)
            else
                cur_expr = parse_expression!(tokens)
                push!(params, cur_expr)  # Add parsed expression to params array
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

            if isempty(tokens) || !isa(popfirst!(tokens), RightBracketToken)
                throw("Expected `}` to close type declaration list")
            end

            # println("type_declarations: ", type_declarations)
            return TypeInstanceNode(symbol_name, symbol_parameters, type_signature_list, type_declarations)
        else

            #TODO: Support expression

            # Assume it's a literal/expression
            tmp = parse_expression!(tokens)
            # literal = popfirst!(tokens)  
            #

            return TypeInstanceNode(symbol_name, symbol_parameters, type_signature_list, tmp)
        end
    end


    return TypeInstanceNode(symbol_name,
                            symbol_parameters,
                            type_signature_list,
                            nothing)

end

function parse_type!(tokens)
    """
    Parse Type
    """

    symbol_name = parse_symbol_name!(tokens)

    if !isa(popfirst!(tokens), SingleColonToken)
        throw("Expected `:` in type declaration")
    end

    type_signature_list = parse_type_signature_list!(tokens)

    return TypeInstanceNode(symbol_name, nothing, type_signature_list, nothing)
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
    """
    Parse Factor
    """

    if lookahead(tokens) isa LeftParenthesisToken
        popfirst!(tokens)
        expr = parse_expression!(tokens)

        lookahead_token = lookahead(tokens)
        if lookahead_token isa RightParenthesisToken
            popfirst!(tokens)
            return GroupNode(expr)
        else
            return expr
        end

        # expect_token!(tokens, RightParenthesisToken)


    elseif lookahead(tokens) isa IdentifierToken
        id_token = popfirst!(tokens)
        if !isempty(tokens) && lookahead(tokens) isa LeftParenthesisToken
            popfirst!(tokens)  # consume '('
            args = Node[]
            while !(lookahead(tokens) isa RightParenthesisToken)
                push!(args, parse_expression!(tokens))
                if lookahead(tokens) isa CommaToken
                    popfirst!(tokens)
                end
            end
            popfirst!(tokens)  # consume ')'
            return CallNode(IdentifierNode(id_token), args)
        else
            return IdentifierNode(id_token)
        end

    elseif lookahead(tokens) isa IntegerToken
        return IntegerNode(popfirst!(tokens))

    elseif lookahead(tokens) isa FloatToken
        return FloatNode(popfirst!(tokens))

    elseif lookahead(tokens) isa RightParenthesisToken
        return popfirst!(tokens)
    else
        error("Unexpected token in factor: $(lookahead(tokens))")
    end
end


function parse_term!(tokens)
    left = parse_factor!(tokens)

    while !isempty(tokens) && (lookahead(tokens) isa AsteriskToken || lookahead(tokens) isa SlashToken)
        op = popfirst!(tokens)
        right = parse_factor!(tokens)
        left = BinaryOpNode(op, left, right)
    end

    return left
end



# <expression> ::= <term>
                # | <expression> "+" <term> 
                # | <expression> "-" <term>
function parse_primary!(tokens)
    token = popfirst!(tokens)

    if token isa IntegerToken
        return IntegerNode(token)
    elseif token isa FloatToken
        return FloatNode(token)
    elseif token isa IdentifierToken
        if !isempty(tokens) && tokens[1] isa LeftParenthesisToken
            popfirst!(tokens)  # consume '('
            args = Node[]
            while !(tokens[1] isa RightParenthesisToken)
                push!(args, parse_expression!(tokens))
                if tokens[1] isa CommaToken
                    popfirst!(tokens)
                end
            end
            popfirst!(tokens)  # consume ')'
            return CallNode(IdentifierNode(token), args)
        else
            return IdentifierNode(token)
        end
    elseif token isa LeftParenthesisToken
        expr = parse_expression!(tokens)
        expect_token!(tokens, RightParenthesisToken)
        return GroupNode(expr)
    else
        error("Unexpected token in expression: $token")
    end
end

function parse_binary_op!(tokens, min_prec)
    left = parse_primary!(tokens)

    while !isempty(tokens)
        op_token = tokens[1]
        op_str = string(typeof(op_token))
        op_str = replace(op_str, r".*\\.(\\w+)Token" => s"\\1")

        if !haskey(PRECEDENCE, op_str)
            break
        end

        prec = PRECEDENCE[op_str]
        if prec < min_prec
            break
        end

        popfirst!(tokens)  # consume operator
        right = parse_binary_op!(tokens, prec + 1)
        left = BinaryOpNode(op_token, left, right)
    end

    return left
end

function parse_expression!(tokens)
    left = parse_term!(tokens)

    while !isempty(tokens) && (lookahead(tokens) isa PlusToken || lookahead(tokens) isa MinusToken)
        op = popfirst!(tokens)
        right = parse_term!(tokens)
        left = BinaryOpNode(op, left, right)
    end

    return left
end


function parse_function_type!(tokens)

    function_name = parse_symbol_name!(tokens)

    cur_token = popfirst!(tokens)
    if !isa(cur_token, LeftParenthesisToken)
        throw("Expected LeftParenthesisToken got: $cur_token")
    end

    function_args = []
    while !isa(lookahead(tokens), RightParenthesisToken)

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

    # check lookahead
    # if isa(lookahead(tokens), LeftParenthesisToken)
    if isa(lookahead(tokens), IdentifierToken) && tokens[2] isa LeftParenthesisToken
        propensity = parse_function_type!(tokens)
    else
        # Probably some expression
        propensity = parse_expression!(tokens)
    end

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

    cur_token = popfirst!(tokens)
    if !isa(cur_token, RightBracketToken)
        throw("Expected RightBracketToken got $cur_token")
    end

    # if !isa(cur_token, RightBracketToken)
        # throw("Expected RightBracketToken got $cur_token")
    # end

    WithClauseNode(propensity, where_clause)
end

function parse_parameterized_type!(tokens)
    popfirst!(tokens)
end

function parse_rule!(tokens)
    """
    Parse Rule
    """

    rule_name = parse_symbol_name!(tokens)
    cur_token = popfirst!(tokens)

    if !isa(cur_token, DefineToken)
        throw("Expected assignment Token ':=' got $cur_token")
    end

    lhs = []
    while !isempty(tokens) && !isa(lookahead(tokens), LeftAngleBracketToken) && !isa(lookahead(tokens), RightArrowToken)

        println("typeof: $(typeof(lookahead(tokens)))")

        if isa(lookahead(tokens), EndLineToken)
            popfirst!(tokens)

        elseif isa(lookahead(tokens), LeftParenthesisToken) || isa(lookahead(tokens), RightParenthesisToken)
            popfirst!(tokens)

        # TODO: Handle directed edges
        elseif isa(lookahead(tokens), EdgeToken)

            # TODO: Handle things like Node -- Node -- Node
            vert_0 = pop!(lhs)

            if vert_0 isa UndirectedTypeEdgeNode
                # If the last node was an edge node, we need to pop it
                # tmp = vert_0.left_vertex
                tmp = vert_0.right_vertex
                push!(lhs, vert_0)
                vert_0 = tmp
            end

            # Popping edge token
            popfirst!(tokens)

            # Popping left parenthesis token
            popfirst!(tokens)

            vert_1 = parse_type_assignment!(tokens)

            # Popping Right parenthesis token
            popfirst!(tokens)

            edge_node = UndirectedTypeEdgeNode(vert_0, vert_1)
            # edge_node.left_vertex = vert_0
            # edge_node.right_vertex = vert_1
            push!(lhs, edge_node)

        else
            tmp = parse_type_assignment!(tokens)
            push!(lhs, tmp)
            # push!(lhs, popfirst!(tokens))
        end
    end

    left_rule_parameter_node = ParameterNode([])
    if isa(lookahead(tokens), LeftAngleBracketToken)
        left_rule_parameter_node = parse_symbol_parameters!(tokens)
    end

    # println("tokens left: $tokens")
    # cur_token = popfirst!(tokens)

    # Gets rid of new lines
    while isa(lookahead(tokens), EndLineToken)
        popfirst!(tokens)
    end

    cur_token = popfirst!(tokens)

    if isempty(tokens) || !isa(cur_token, RightArrowToken)
        throw("Expected `->` in rule definition got $cur_token")
    end

    rhs = []
    while !isempty(tokens) && !isa(lookahead(tokens), LeftAngleBracketToken)
        if isa(lookahead(tokens), EndLineToken)
            popfirst!(tokens)
        elseif isa(lookahead(tokens), LeftParenthesisToken) || isa(lookahead(tokens), RightParenthesisToken)
            popfirst!(tokens)

        # TODO: Handle directed edges
        elseif isa(lookahead(tokens), EdgeToken)

            # TODO: Handle things like Node -- Node -- Node
            vert_0 = pop!(rhs)

            # println("check vert_0 ---------> ", vert_0)
# 
            if vert_0 isa UndirectedTypeEdgeNode
                # If the last node was an edge node, we need to pop it
                tmp = vert_0.right_vertex
                push!(rhs, vert_0)
                vert_0 = tmp
            end

            # Popping edge token
            popfirst!(tokens)

            # Popping left parenthesis token
            popfirst!(tokens)

            vert_1 = parse_type_assignment!(tokens)

            # Popping Right parenthesis token
            popfirst!(tokens)

            edge_node = UndirectedTypeEdgeNode(vert_0, vert_1)
            push!(rhs, edge_node)

        else
            tmp = parse_type_assignment!(tokens)
            push!(rhs, tmp)
            # push!(lhs, popfirst!(tokens))
        end
        # push!(rhs, popfirst!(tokens))
    end

    right_rule_parameter_node = ParameterNode([])
    if isa(lookahead(tokens), LeftAngleBracketToken)
        # println("Parsing symb param: ", tokens)
        right_rule_parameter_node = parse_symbol_parameters!(tokens)
    end

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

            with_clause = parse_with_clause!(tokens)

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

        if isa(cur_token, TypeSectionToken)
            push!(ast_nodes, parse_types_section!(tokens))
        elseif isa(cur_token, ParameterSectionToken)
            push!(ast_nodes, parse_parameters_section!(tokens))

        elseif isa(cur_token, RuleSectionToken)
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

function parse_rules_list!(tokens)
    """
    parse rules list
    """

    rules = []

    # Remove extra lines
    while isa(lookahead(tokens),EndLineToken)
        popfirst!(tokens)
    end

    # Pares rule
    rule = parse_rule!(tokens)
    push!(rules, rule)

    # Remove nedline tokens and return the rules
    # if there are no more tokens
    while isa(lookahead(tokens),EndLineToken)
        popfirst!(tokens)
        if isempty(tokens)
            return rules
        end
    end

    # If there are no more endline tokens
    # There are still rules left.
    if length(tokens) != 0
        # Still another rule
        if isa(lookahead(tokens), IdentifierToken)
        rules = [rules;parse_rules_list!(tokens)]
        end
    end

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

    # Check if tokens are empty
    if isempty(tokens)
        return RuleSectionNode(section_name, rules_list)
    else

        while !isempty(tokens)
            cur_token = popfirst!(tokens)
            if isa(cur_token, EndLineToken)
                continue
            end
            throw("Expected empty token list, got: $(tokens)")
        end
        return RuleSectionNode(section_name, rules_list)
    end
end

# Entry point for parsing
function main()
    # tokens = tokenize_file("tests/rules/test_edge_rule.fflow")
    # ast = parse_file!(tokens)
    # println("ast: ", ast)

    check = "
    create_edge_2 := (p1 : ParticleNode) (p2: ParticleNode) << (x_0, y_0), (x_1, y_1) >> ->
    (p1 : ParticleNode) -- (p2 : ParticleNode) -- (p3 : ParticleNode) << (x_0, y_0), (x_1, y_1), (x_2, y_2) >> with (heaviside(10, 1)) where { 

        x_2 : Float = x_2 + 0.3 
        y_2 : Float = y_2 + 0.3 

    } 
    "

    res = tokenize_string(check)

    oof = parse_rule!(res)

    rhs_check = oof.rhs

    rhs_param = oof.rhs_parameter
    println("rhs_param: $rhs_param")

end

# main()
