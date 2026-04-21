import ..Tokens:
    Token, IdentifierToken, IntegerToken, FloatToken, EndLineToken,
    LeftParenthesisToken, RightParenthesisToken, LeftAngleBracketToken,
    RightAngleBracketToken, SingleColonToken, DefineToken, LeftBracketToken,
    RightBracketToken, PunctuationToken, TypeSectionToken, ParameterSectionToken,
    RuleSectionToken, SlashToken, AsteriskToken, PlusToken, MinusToken, RightArrowToken, WithToken, WhereToken, SolvingToken,
    EqualToken, CommaToken, EdgeToken, ODEToken, SampleToken, BackslashToken,
    LtToken, GtToken, LtEqToken, GtEqToken, EqEqToken, NotEqToken,
    NotToken, AndAndToken, OrOrToken, TypeToken, FunctionSectionToken, FunctionToken, CaretToken, ReturnToken,
    StateToken, SimulationTypesToken, RulesToken, DotToken, ParameterToken, SimulationToken, RunSimulationToken, QuoteToken, SimulationSectionToken, LoadFileToken, StringToken, SimulationParametersToken, DoubleColonToken,
    LeftSquareBracketToken, RightSquareBracketToken
import ..AstNodes:
    Node, IdentifierNode, ParameterNode, TypeInstanceNode, TypeInstanceUpdateNode, TypeClassNode,
    TypeSectionNode, ParameterSectionNode, RuleSectionNode, RuleNode,
    WhereClauseNode, WithClauseNode, SolveClauseNode, CallNode,
    GroupNode, BinaryOpNode, FunctionNode, FloatNode, IntegerNode, UndirectedTypeEdgeNode,
    BindingVariableNode, ODENode, UnaryOpNode, DefinitionNode, FunctionArgNode, FunctionSignatureNode,
    FunctionBodyExpressionNode, FunctionDefinitionExpressionNode, FunctionDefinitionNode,
    ReturnNode, FunctionSectionNode, SimulationSectionNode, StringNode, LoadNode, SimDeclarationNode,
    SimulationNode, SimulationRulesNode, SimulationTypesNode, SimulationStateNode, SimulationParametersNode,
    RunSimulationNode, NamedParameterNode, IndexAccessNode, ArrayLiteralNode, ModelLoadNode, SliceNode

include("expression_parser.jl")
include("function_parser.jl")
include("parser_utils.jl")
include("sim_parser.jl")

const PRECEDENCE = Dict(
    "+" => 1,
    "-" => 1,
    "*" => 2,
    "/" => 2,
    "^" => 3
)

function expect_token!(tokens, ::Type{T}) where T <: Token

    token = popfirst!(tokens)

    if !(token isa T)
        error("Expected $(T), got $(typeof(token))")
    end
    return token
end

function tokenize_file(file_name)

    lines_tokens = []
    open(file_name) do f
        line_no = 0
        while !eof(f)

            cur_line = readline(f)
            cur_line = String(lstrip(cur_line))
            if cur_line != ""
                if cur_line[1] == '#'
                    line_no += 1
                    continue
                end
            end

            line_token = tokenize(cur_line, line_no+1)

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
    check
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

            println("lookahead tokens: $(lookahead(tokens))")
            # Handling nested parameters
            if isa(lookahead(tokens), LeftAngleBracketToken)
                nested_params = parse_symbol_parameters!(tokens)
                push!(params, nested_params)  # Add nested params to params array

            elseif isa(lookahead(tokens), LeftParenthesisToken) || isa(lookahead(tokens), RightParenthesisToken)
                popfirst!(tokens)
            elseif isa(lookahead(tokens), EndLineToken)
                # Pop EndLineToken and ignore
                popfirst!(tokens)

            elseif isa(lookahead(tokens), IdentifierToken) && length(tokens) > 1 && isa(tokens[2], SingleColonToken)
                # Handle cases like `param: Type`
                param_name = parse_symbol_name!(tokens)
                popfirst!(tokens) # Consume `:`
                param_type = parse_type_signature_list!(tokens)

                # println("param_name: $param_name, param_type: $param_type")

                push!(params, NamedParameterNode(param_name, param_type))
                # push!(params, (param_name, param_type))  # Add parameter name and type as a tuple to params array

            else

                # It could jsut be a regular identifer token. Is this valid?
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

function parse_type_symbol_parameters!(tokens)
    params = []

    lookahead_token = lookahead(tokens)
    
    if !isempty(tokens) && isa(lookahead(tokens), LeftAngleBracketToken)

        popfirst!(tokens) # Consume `<<`

        while !isempty(tokens) && !isa(lookahead(tokens), RightAngleBracketToken)

            if !(lookahead(tokens) isa IdentifierToken)
                throw("Expected IdentifierToken in type parameter list, got $(lookahead(tokens))")
            end

            param_name = parse_symbol_name!(tokens)
            if !(lookahead(tokens) isa SingleColonToken)
                throw("Expected SingleColonToken after parameter name in type parameter list, got $(lookahead(tokens))")
            end
            popfirst!(tokens) # Consume `:`

            # Handling nested parameters
            if isa(lookahead(tokens), LeftAngleBracketToken)
                nested_params = parse_symbol_parameters!(tokens)
                push!(params, nested_params)  # Add nested params to params array

            elseif isa(lookahead(tokens), LeftParenthesisToken) || isa(lookahead(tokens), RightParenthesisToken)
                popfirst!(tokens)
            elseif isa(lookahead(tokens), EndLineToken)
                # Pop EndLineToken and ignore
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


function parse_type_signature_list!(tokens)

    typename = nothing
    if isa(lookahead(tokens), IdentifierToken)
        type_name = parse_symbol_name!(tokens)
    elseif typeof(lookahead(tokens)) in [IntegerToken, FloatToken]
        type_name = IdentifierNode(popfirst!(tokens))
    elseif isa(lookahead(tokens), TypeToken)
        type_name = IdentifierNode(popfirst!(tokens))
    else
        cur_token = popfirst!(tokens)
        throw("Error couldnt determine type got: $cur_token")
    end

    if !isempty(tokens) && isa(lookahead(tokens), RightArrowToken)
        popfirst!(tokens)  # Consume `->`
        return TypeClassNode(type_name, parse_type_signature_list!(tokens))
    elseif !isempty(tokens) && isa(lookahead(tokens), LeftAngleBracketToken)
        # popfirst!(tokens)  # Consume `->`
        param_node = parse_symbol_parameters!(tokens)
        return TypeClassNode(type_name, param_node)
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
        throw("Expected `:` after symbol name got $(lookahead(tokens))")
    end

    # Parse Type Signature
    type_signature_list = parse_type_signature_list!(tokens)

    # println("type_signature_list: $type_signature_list")

    if !isempty(tokens) && isa(lookahead(tokens), DefineToken)

        popfirst!(tokens)  # Consume `=`

        if !isempty(tokens) && isa(lookahead(tokens), LeftBracketToken)
            popfirst!(tokens) # Consume `{`
            type_declarations = parse_type_declaration_list!(tokens)

            if isempty(tokens) || !isa(popfirst!(tokens), RightBracketToken)
                throw("Expected `}` to close type declaration list")
            end

            return TypeInstanceNode(symbol_name, symbol_parameters, type_signature_list, type_declarations)
        else

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

"""
    parse_index_element!(tokens)

Parse a single index element inside [...]. This can be:
  - A plain expression:        tensor[i]
  - A slice with start:stop:    tensor[0:3]
  - A slice with start:stop:step: tensor[1:5:2]
  - A bare colon (select all):  tensor[:]
  - A colon with step (::step):  tensor[::2]

Returns either a normal expression Node or a SliceNode.
"""
function parse_index_element!(tokens)
    # Check for bare colon first: [:] or [::step] or [:stop] or [:stop:step]
    if !isempty(tokens) && lookahead(tokens) isa SingleColonToken
        popfirst!(tokens)  # consume ':'
        start_expr = nothing
        # Check for second colon (::step pattern) or stop expression
        if !isempty(tokens) && lookahead(tokens) isa SingleColonToken
            # ::step
            popfirst!(tokens)  # consume second ':'
            if !isempty(tokens) && !(lookahead(tokens) isa RightSquareBracketToken) && !(lookahead(tokens) isa PunctuationToken)
                step_expr = parse_expression!(tokens)
                return SliceNode(nothing, nothing, step_expr)
            else
                return SliceNode(nothing, nothing, nothing)
            end
        elseif !isempty(tokens) && !(lookahead(tokens) isa RightSquareBracketToken) && !(lookahead(tokens) isa PunctuationToken)
            # :stop or :stop:step
            stop_expr = parse_expression!(tokens)
            if !isempty(tokens) && lookahead(tokens) isa SingleColonToken
                popfirst!(tokens)  # consume ':'
                step_expr = parse_expression!(tokens)
                return SliceNode(nothing, stop_expr, step_expr)
            end
            return SliceNode(nothing, stop_expr, nothing)
        else
            # bare : (select all)
            return SliceNode(nothing, nothing, nothing)
        end
    end

    # Parse the first expression (could be start of a slice or a plain index)
    expr = parse_expression!(tokens)

    # Check if followed by colon -> slice
    if !isempty(tokens) && lookahead(tokens) isa SingleColonToken
        popfirst!(tokens)  # consume ':'
        start_expr = expr
        # Check for second colon right away (start::step)
        if !isempty(tokens) && lookahead(tokens) isa SingleColonToken
            popfirst!(tokens)  # consume second ':'
            if !isempty(tokens) && !(lookahead(tokens) isa RightSquareBracketToken) && !(lookahead(tokens) isa PunctuationToken)
                step_expr = parse_expression!(tokens)
                return SliceNode(start_expr, nothing, step_expr)
            else
                return SliceNode(start_expr, nothing, nothing)
            end
        elseif !isempty(tokens) && !(lookahead(tokens) isa RightSquareBracketToken) && !(lookahead(tokens) isa PunctuationToken)
            # start:stop possibly followed by :step
            stop_expr = parse_expression!(tokens)
            if !isempty(tokens) && lookahead(tokens) isa SingleColonToken
                popfirst!(tokens)  # consume ':'
                step_expr = parse_expression!(tokens)
                return SliceNode(start_expr, stop_expr, step_expr)
            end
            return SliceNode(start_expr, stop_expr, nothing)
        else
            # start: (from start to end)
            return SliceNode(start_expr, nothing, nothing)
        end
    end

    return expr
end

function parse_type_update!(tokens)
    # Updates to existing types done by a rule or sees a definition of a tmp variable

    symbol_name = parse_symbol_name!(tokens)

    # Check for indexed LHS: name[idx] = expr
    lhs_node = symbol_name
    if !isempty(tokens) && lookahead(tokens) isa LeftSquareBracketToken
        while !isempty(tokens) && lookahead(tokens) isa LeftSquareBracketToken
            popfirst!(tokens)  # consume '['
            indices = Node[]
            push!(indices, parse_index_element!(tokens))
            while !isempty(tokens) && lookahead(tokens) isa PunctuationToken
                popfirst!(tokens)  # consume ','
                push!(indices, parse_index_element!(tokens))
            end
            if isempty(tokens) || !(lookahead(tokens) isa RightSquareBracketToken)
                error("Expected ']' after index expression in LHS")
            end
            popfirst!(tokens)  # consume ']'
            lhs_node = IndexAccessNode(lhs_node, indices)
        end
    end

    symbol_parameters = parse_symbol_parameters!(tokens)
    type_signature_list = nothing

    # Its a definition node
    if lookahead(tokens) isa SingleColonToken
        popfirst!(tokens)  # Consume `:`

        type_signature_list = parse_type_signature_list!(tokens)

        # This is an intermediate value and not a type update
        if isa(lookahead(tokens), DefineToken)
            popfirst!(tokens)  # Consume `:=`
            expression = []
            while !isa(lookahead(tokens), EndLineToken)
                literal = popfirst!(tokens)  # Assume it's a literal/expression
                push!(expression,literal)
            end

            expression_nodes = parse_expression!(expression)
            return DefinitionNode(symbol_name, type_signature_list, expression_nodes)
        end


    end
    
    if !isempty(tokens) && isa(lookahead(tokens), EqualToken)

        popfirst!(tokens)  # Consume `=`

        if !isempty(tokens) && isa(lookahead(tokens), LeftBracketToken)

            popfirst!(tokens) # Consume `{`
            type_declarations = parse_type_declaration_list!(tokens)

            if isempty(tokens) || !isa(popfirst!(tokens), RightBracketToken)
                throw("Expected `}` to close type declaration list")
            end


            # return TypeInstanceNode(symbol_name, symbol_parameters, type_signature_list, type_declarations)
            return TypeInstanceUpdateNode(lhs_node, type_declarations)
        else

            expression = []
            while !isa(lookahead(tokens), EndLineToken) && !isa(lookahead(tokens), RightBracketToken)
                literal = popfirst!(tokens)  # Assume it's a literal/expression
                push!(expression,literal)
            end

            expression_nodes = parse_expression!(expression)
            # return TypeInstanceNode(symbol_name, symbol_parameters, type_signature_list, expression_nodes)
            return TypeInstanceUpdateNode(lhs_node, expression_nodes)

        end
    end

    return TypeInstanceNode(symbol_name, symbol_parameters, type_signature_list, nothing)
end


# Parse a type declaration (non-assignment case)
function parse_type_declaration!(tokens)

    # Nucleator
    type_name = parse_symbol_name!(tokens)

    # << param >>
    symbol_parameters = parse_symbol_parameters!(tokens)

    if !isa(lookahead(tokens), SingleColonToken)
        throw("Expected `:` in type declaration got $(lookahead(tokens))")
    else
        popfirst!(tokens)  # Consume `:`
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

    # Remove EndLineTokens
    while isa(lookahead(tokens), EndLineToken)
        popfirst!(tokens)
    end

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
        namespace = nothing

        if (length(tokens) > 0)
            if (lookahead(tokens) isa DoubleColonToken)
                # popping Double Colon
                popfirst!(tokens)
                namespace = id_token
                if lookahead(tokens) isa IdentifierToken
                    id_token = popfirst!(tokens)
                else
                    throw("Error expected an IdentifierToken after '::'")
                end
            end
        end

        if !isempty(tokens) && lookahead(tokens) isa LeftParenthesisToken

            
            popfirst!(tokens)  # consume '('
            args = Node[]

            while !(lookahead(tokens) isa RightParenthesisToken)
                push!(args, parse_expression!(tokens))
                if lookahead(tokens) isa PunctuationToken
                    popfirst!(tokens)
                end

            end
            popfirst!(tokens)  # consume ')'

            # FIXME: Namespace not appearing with variable.
            result = CallNode(IdentifierNode(id_token, namespace), args)
        else
            result = IdentifierNode(id_token)
        end

        # Postfix index access: expr[i] or expr[i, j] or expr[0:3] etc.
        while !isempty(tokens) && lookahead(tokens) isa LeftSquareBracketToken
            popfirst!(tokens)  # consume '['
            indices = Node[]
            push!(indices, parse_index_element!(tokens))
            while !isempty(tokens) && lookahead(tokens) isa PunctuationToken
                popfirst!(tokens)  # consume ','
                push!(indices, parse_index_element!(tokens))
            end
            if isempty(tokens) || !(lookahead(tokens) isa RightSquareBracketToken)
                error("Expected ']' after index expression")
            end
            popfirst!(tokens)  # consume ']'
            result = IndexAccessNode(result, indices)
        end

        return result

    elseif lookahead(tokens) isa IntegerToken
        return IntegerNode(popfirst!(tokens))

    elseif lookahead(tokens) isa FloatToken
        return FloatNode(popfirst!(tokens))

    elseif lookahead(tokens) isa RightParenthesisToken
        return popfirst!(tokens)

    elseif lookahead(tokens) isa MinusToken
        minus_token = popfirst!(tokens)  # consume '-'
        right = parse_factor!(tokens)
        return UnaryOpNode(minus_token, right)
    elseif lookahead(tokens) isa PlusToken
        plus_token = popfirst!(tokens)  # consume '+'
        right = parse_factor!(tokens)
        return UnaryOpNode(plus_token, right)
    elseif lookahead(tokens) isa NotToken
        not_token = popfirst!(tokens)  # consume '!'
        right = parse_factor!(tokens)
        return UnaryOpNode(not_token, right)
    elseif lookahead(tokens) isa SampleToken
        sample_token = popfirst!(tokens)  # consume '~'
        right = parse_factor!(tokens)
        return UnaryOpNode(sample_token, right)

    elseif lookahead(tokens) isa LeftSquareBracketToken
        # Array/tensor literal: [expr, expr, ...] or [[expr,...], [expr,...], ...]
        popfirst!(tokens)  # consume '['
        elements = Node[]
        if !(lookahead(tokens) isa RightSquareBracketToken)
            push!(elements, parse_expression!(tokens))
            while !isempty(tokens) && lookahead(tokens) isa PunctuationToken
                popfirst!(tokens)  # consume ','
                push!(elements, parse_expression!(tokens))
            end
        end
        if isempty(tokens) || !(lookahead(tokens) isa RightSquareBracketToken)
            error("Expected ']' after array literal")
        end
        popfirst!(tokens)  # consume ']'
        return ArrayLiteralNode(elements)

    # elseif lookahead(tokens) isa RightAngleBracketToken
        # right_angle_token = popfirst!(tokens)  # consume '>'
        # right = parse_factor!(tokens)
        # return UnaryOpNode(right_angle_token, right)
    else
        println("tokens: $tokens")
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

function skip_eol!(tokens)
    while !isempty(tokens) && isa(lookahead(tokens), EndLineToken)
        popfirst!(tokens)
    end
    return nothing
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

# function parse_type_assignment_list!(tokens)
#     """
#     Parse Type Assignment List
#     """
# 
#     type_assignments = []
# 
#     while isa(lookahead(tokens), EndLineToken)
#         popfirst!(tokens)
#     end
# 
#     cur_type_assignment = parse_type_assignment!(tokens)
# 
#     push!(type_assignments, cur_type_assignment)
# 
#     if isa(lookahead(tokens), EndLineToken)
#         tmp = parse_type_assignment_list!(tokens)
#         type_assignments = [type_assignments;tmp]
#     end
# 
#     return type_assignments
# end

function parse_type_update_list!(tokens)

    type_updates = []

    while isa(lookahead(tokens), EndLineToken)
        popfirst!(tokens)
    end

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

    # cur_token = popfirst!(tokens)
    if !isa(lookahead(tokens), LeftParenthesisToken)
        println("Expected LeftParenthesis got $cur_token")
    end

    propensity = parse_expression!(tokens)

    # remove endline tokens
    while isa(lookahead(tokens), EndLineToken)
        popfirst!(tokens)
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

function parse_solving_content(tokens)
    """
    Parse Solving Content
    """

    content = []

    while !isempty(tokens) && !isa(lookahead(tokens), EndLineToken)
        if isa(lookahead(tokens), EndLineToken)
            popfirst!(tokens)
        elseif isa(lookahead(tokens), LeftParenthesisToken) || isa(lookahead(tokens), RightParenthesisToken)
            popfirst!(tokens)
        else
            push!(content, parse_expression!(tokens))
        end
    end

    return content
end

function parse_binding_variable!(tokens)
    """
    Parse Binding Variable
    """

    # Should be identifier name, possibly with index access like dpos[0]
    var_name = popfirst!(tokens)
    if !isa(var_name, IdentifierToken)
        throw("Expected IdentifierToken got $var_name")
    end

    # Check for index access: dpos[0]
    var_node = IdentifierNode(var_name)
    if !isempty(tokens) && lookahead(tokens) isa LeftSquareBracketToken
        popfirst!(tokens)  # consume '['
        indices = Node[]
        push!(indices, parse_expression!(tokens))
        while !isempty(tokens) && lookahead(tokens) isa PunctuationToken
            popfirst!(tokens)  # consume ','
            push!(indices, parse_expression!(tokens))
        end
        if isempty(tokens) || !(lookahead(tokens) isa RightSquareBracketToken)
            error("Expected ']' after index in binding variable")
        end
        popfirst!(tokens)  # consume ']'
        var_node = IndexAccessNode(var_node, indices)
    end

    # Should pop := define symbol
    cur_token = popfirst!(tokens)
    if !(cur_token isa DefineToken)
        throw("Expected DefineToken got $cur_token")
    end

    D_token = popfirst!(tokens)
    if D_token.position.value != "D"
        throw("Expected Derivative Operator D(..) got $D_token")
    end

    left_paren = popfirst!(tokens)
    if !isa(left_paren, LeftParenthesisToken)
        throw("Expected LeftParenthesisToken got $left_paren")
    end

    ode_variables = []
    while !(lookahead(tokens) isa RightParenthesisToken)

        if lookahead(tokens) isa PunctuationToken
            popfirst!(tokens)
        else
            ode_variable = parse_symbol_name!(tokens)
            # Check for index access on ODE variable: D(im_pos[0], t)
            if !isempty(tokens) && lookahead(tokens) isa LeftSquareBracketToken
                popfirst!(tokens)  # consume '['
                indices = Node[]
                push!(indices, parse_expression!(tokens))
                while !isempty(tokens) && lookahead(tokens) isa PunctuationToken
                    popfirst!(tokens)  # consume ','
                    push!(indices, parse_expression!(tokens))
                end
                if isempty(tokens) || !(lookahead(tokens) isa RightSquareBracketToken)
                    error("Expected ']' after index in ODE variable")
                end
                popfirst!(tokens)  # consume ']'
                ode_variable = IndexAccessNode(ode_variable, indices)
            end
            push!(ode_variables, ode_variable)
        end

    end

    right_paren = popfirst!(tokens)
    if !isa(right_paren, RightParenthesisToken)
        throw("Expected RightParenthesisToken got $right_paren")
    end

    BindingVariableNode(var_node, ode_variables)
end

# TODO: Add support for DefinitionNode
function parse_solve_clause!(tokens)
    """
    Parses Solve Clause and returns an ODENode.
    Handles both simple names (dx : ODE = expr) and
    indexed names (dx[0] : ODE = expr) for tensor component ODEs.
    """

    var_name = popfirst!(tokens)

    if !(var_name isa IdentifierToken)
        throw("Expected Identifier token got $var_name")
    end

    var_name = IdentifierNode(var_name)

    # Check for index access: dx[i] : ODE = expr or dx[0:3] : ODE = expr
    if !isempty(tokens) && lookahead(tokens) isa LeftSquareBracketToken
        popfirst!(tokens)  # consume '['
        indices = Node[]
        push!(indices, parse_index_element!(tokens))
        while !isempty(tokens) && lookahead(tokens) isa PunctuationToken
            popfirst!(tokens)  # consume ','
            push!(indices, parse_index_element!(tokens))
        end
        if isempty(tokens) || !(lookahead(tokens) isa RightSquareBracketToken)
            error("Expected ']' after index expression in ODE name")
        end
        popfirst!(tokens)  # consume ']'
        var_name = IndexAccessNode(var_name, indices)
    end

    colon_token = popfirst!(tokens)
    if !(colon_token isa SingleColonToken)
        throw("Error expected single colon token")
    end

    # Should be an ode token (Type)
    ode_token = popfirst!(tokens)
    if !((ode_token isa ODEToken) || (ode_token isa FloatToken) || (ode_token isa IntegerToken))
        throw("Error expected an ODE type or Float or Integer. Got: $ode_token.")
    end

    eq_token = popfirst!(tokens)
    if !(eq_token isa EqualToken || eq_token isa DefineToken)
        throw("Error expected '=' or ':=' got $eq_token")
    end

    expression = parse_expression!(tokens)

    if eq_token isa EqualToken
        return ODENode(var_name, expression)
    elseif eq_token isa DefineToken
        return DefinitionNode(var_name, TypeClassNode(IdentifierNode(ode_token), ParameterNode([])), expression)
    else
        throw("Unexpected token in solve clause: $eq_token")
    end
end

function parse_rule_solve!(tokens)
    """
    Parses the Binding Variables
    """

    # Pops left parenthesis
    cur_token = popfirst!(tokens)
    if !isa(cur_token, LeftParenthesisToken)
        throw("Expected LeftParenthesisToken got $cur_token")
    end

    total_binding_variables = []
    while !(lookahead(tokens) isa RightParenthesisToken)

        if (isa(lookahead(tokens), EndLineToken) 
            || isa(lookahead(tokens), PunctuationToken))
            popfirst!(tokens)
        else
            push!(total_binding_variables, parse_binding_variable!(tokens))
        end

    end

    # Pops right parenthesis
    cur_token = popfirst!(tokens)
    if !isa(cur_token, RightParenthesisToken)
        throw("Expected RightParenthesisToken got $cur_token")
    end

    # Need to parse the  { ... } at the end

    # { token
    cur_token = popfirst!(tokens)
    if !isa(cur_token, LeftBracketToken)
        throw("Expected LeftBracketToken got $cur_token")
    end

    solve_clause = []
    while !isempty(tokens) && !isa(lookahead(tokens), RightBracketToken)

        if isa(lookahead(tokens), EndLineToken)
            popfirst!(tokens)
        elseif isa(lookahead(tokens), LeftParenthesisToken) || isa(lookahead(tokens), RightParenthesisToken)
            popfirst!(tokens)
            # Check its a definition token
        else
            ode_node = parse_solve_clause!(tokens)
            # push!(solve_clause, parse_expression!(tokens))
            push!(solve_clause, ode_node)
        end
    end

    cur_token = popfirst!(tokens)
    if !isa(cur_token, RightBracketToken)
        throw("Expected RightBracketToken got $cur_token")
    end

    SolveClauseNode(
        total_binding_variables,
        solve_clause
    )
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

        if isa(lookahead(tokens), EndLineToken)
            popfirst!(tokens)

        elseif isa(lookahead(tokens), LeftParenthesisToken) || isa(lookahead(tokens), RightParenthesisToken)
            popfirst!(tokens)

        # TODO: Handle directed edges
        elseif isa(lookahead(tokens), EdgeToken)

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

    # println("left rule parameter node: $left_rule_parameter_node")

    # Gets rid of new lines
    while isa(lookahead(tokens), EndLineToken)
        popfirst!(tokens)
    end

    cur_token = popfirst!(tokens)

    if isempty(tokens) || !isa(cur_token, RightArrowToken)
        throw("Expected `->` in rule definition got $cur_token")
    end

    # Parsing RHS
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
        right_rule_parameter_node = parse_symbol_parameters!(tokens)
    end

    # TODO: check parameters
    modify_clause = nothing

    # Removes Endlines Before Solving
    while(isa(lookahead(tokens), EndLineToken))
        popfirst!(tokens)  # Consume EndLineToken
    end

    if (!isempty(tokens) 
        && (isa(lookahead(tokens), WithToken) 
            || isa(lookahead(tokens), SolvingToken))
       )

        # Consume `with` or `solving`
        clause_token = popfirst!(tokens) 
        clause_content = []
        
        if isa(clause_token, WithToken)

            with_clause = parse_with_clause!(tokens)

            modify_clause = with_clause
            # modify_clause = WithClauseNode(clause_token, clause_content)
        else

            modify_clause = parse_rule_solve!(tokens)
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

function parse_function_section!(tokens)
    popfirst!(tokens)
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
        elseif isa(cur_token, FunctionSectionToken)
            push!(ast_nodes, parse_functions_section!(tokens))

        elseif isa(cur_token, SimulationSectionToken)
            push!(ast_nodes, parse_simulations_section!(tokens))

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

    # Parse rule
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

