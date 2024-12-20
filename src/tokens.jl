
"""
Defines Tokens
"""

# TODO: Add --

abstract type Token end

struct KeywordToken <: Token
    value::String
end

struct GrammarToken <: Token
    value::String
end

struct GrammarTimeToken <: Token
    value::String
end

struct WithToken <: Token
    value::String
end

struct SolvingToken <: Token
    value::String
end

struct WhereToken <: Token
    value::String
end

struct IdentifierToken <: Token
    value::String
end

struct OperatorToken <: Token
    value::String
end

struct RightArrowToken <: Token
    value::String
end

struct PunctuationToken <: Token
    value::String
end

struct LeftParenthesisToken <: Token
    value::String
end

struct RightParenthesisToken <: Token
    value::String
end

struct LeftBracketToken <: Token
    value::String
end

struct RightBracketToken <: Token
    value::String
end

struct LiteralToken <: Token
    value::String
end

struct ErrorToken <: Token
    message::String
end

struct InitialConditionToken <: Token
    value::String
end

function isalpha(input::Char)::Bool
    isletter(input)
end

function tokenize(source::String)
    """
    Tokenize input
    """

    tokens = Token[]
    i = 1

    while i <= length(source)
        c = source[i]

        # Ignore whitespace
        if isspace(c)
            i += 1
            continue

        # New starting point when a character is found
        elseif isalpha(c)
            start = i

        # Keep on counting 
	    while i <= length(source) && (isalpha(source[i]) || isdigit(source[i]) || source[i] == '_')
                i += 1
            end

            value = source[start:i-1]
            if value in [
			 "grammar", "discretetime", "continuoustime",
			 "with", "solving", "where", "initial_conditions",
			 ]
                if string(value) == "grammar"
                    push!(tokens, GrammarToken(value))
                elseif string(value) == "discretetime" || string(value) == "continuoustime"
                    push!(tokens, GrammarTimeToken(value))
                elseif string(value) == "with"
                    push!(tokens, WithToken(value))
                elseif string(value) == "solving"
                    push!(tokens, SolvingToken(value))
                elseif string(value) == "where"
                    push!(tokens, WhereToken(value))
                elseif string(value) == "initial_conditions"
                    push!(tokens, InitialConditionToken(value))
                else
                    push!(tokens, KeywordToken(value))
                end
            else
                push!(tokens, IdentifierToken(value))
            end

        elseif isdigit(c)
            start = i
            while i <= length(source) && (isdigit(source[i]) || source[i] == '.')
                i += 1
            end
            value = source[start:i-1]
            push!(tokens, LiteralToken(value))
        elseif c == '-'
            if i + 1 <= length(source) && source[i + 1] == '>'
                push!(tokens, RightArrowToken("->"))
                i += 2
            else
                push!(tokens, ErrorToken("Unexpected character: $c"))
                i += 1
            end
        elseif c in ['!', '=', '<', '>']
            if i + 1 <= length(source) && source[i + 1] == '='
                push!(tokens, OperatorToken(string(c, '=')))
                i += 2
            else
                push!(tokens, OperatorToken(string(c)))
                i += 1
            end
        elseif c in ['(', ')', '{', '}', ',', ';']

            if c == '('
                push!(tokens, LeftParenthesisToken(string(c)))
            elseif c == ')'
                push!(tokens, RightParenthesisToken(string(c)))
            elseif c == '{'
                push!(tokens, LeftBracketToken(string(c)))
            elseif c == '}'
                push!(tokens, RightBracketToken(string(c)))
            else
                push!(tokens, PunctuationToken(string(c)))
            end

            i += 1
        else
            push!(tokens, ErrorToken("Unexpected character: $c"))
            i += 1
        end
    end
    return tokens
end

# source_code = """
# grammar discretetime symbol1 (symbol2 -> symbola  with function(symbol4) where symbol5 == symbol6;
# """

# source_code2 = """
	# grammar discretetime test_grammar(){
	# }
# """

# source_code3 = """
   # nodeset(x) -> node(x), {child(x) | 1 <= i <= n} with q(n) where n >= 0;
# """

# tokens = tokenize(source_code2)

# for token in tokens
    # println(token)
# end
