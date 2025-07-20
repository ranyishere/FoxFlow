include("./tokens.jl")


function isalpha(input::Char)::Bool
    isletter(input)
end

function tokenize(source::String, line_no::Integer)
    """
    Tokenize input
    """

    tokens = Token[]
    i = 1
    while i <= length(source)

        c = source[i]

        # println("c: $c")

        # Ignore whitespace
        if isspace(c)
            # println("space: $i")
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
                 "with", "solving", "where"
             ]
            if string(value) == "grammar"
                push!(tokens, GrammarToken(PositionToken(value, line_no, i-1)))
            elseif string(value) == "discretetime" || string(value) == "continuoustime"
                push!(tokens, GrammarTimeToken(PositionToken(value, line_no, i-1)))
            elseif string(value) == "with"
                push!(tokens, WithToken(PositionToken(value, line_no, i-1)))
            elseif string(value) == "solving"
                push!(tokens, SolvingToken(PositionToken(value, line_no, i-1)))
            elseif string(value) == "where"
                push!(tokens, WhereToken(PositionToken(value, line_no, i-1)))
            else
                push!(tokens, KeywordToken(PositionToken(value, line_no, i-1)))
            end

            elseif value in [
                     "types", "parameters",
                     "functions", "rules",
                     "observables", "states",
                     "simulations"
                ]

                if value == "types"
                    push!(tokens, TypeSectionToken(PositionToken(value, line_no, i-1)))
                elseif value == "parameters"
                    push!(tokens, ParameterSectionToken(PositionToken(value, line_no, i-1)))
                elseif value == "functions"
                    push!(tokens, FunctionSectionToken(PositionToken(value, line_no, i-1)))
                elseif value == "rules"
                    push!(tokens, RuleSectionToken(PositionToken(value, line_no, i-1)))
                elseif value == "observables"
                    push!(tokens, ObservableSectionToken(PositionToken(value, line_no, i-1)))
                elseif value == "states"
                    push!(tokens, StateSectionToken(PositionToken(value, line_no, i-1)))
                elseif value == "simulations"
                    push!(tokens, SimulationSectionToken(PositionToken(value, line_no, i-1)))
            end

            elseif value in [
                     "type", "parameter",
                     "function", "rule",
                     "observable", "state",
                     "simulation"
                ]

                if value == "type"
                    push!(tokens, TypeToken(PositionToken(value, line_no, i-1)))
                elseif value == "parameter"
                    push!(tokens, ParameterToken(PositionToken(value, line_no, i-1)))
                elseif value == "function"
                    push!(tokens, FunctionToken(PositionToken(value, line_no, i-1)))
                elseif value == "rule"
                    push!(tokens, RuleToken(PositionToken(value, line_no, i-1)))
                elseif value == "observable"
                    push!(tokens, ObservableToken(PositionToken(value, line_no, i-1)))
                elseif value == "state"
                    push!(tokens, StateToken(PositionToken(value, line_no, i-1)))
                elseif value == "simulation"
                    push!(tokens, SimulationToken(PositionToken(value, line_no, i-1)))
            end

            elseif value in ["Float", "Integer"]
                if value == "Float"
                    push!(tokens, FloatToken(PositionToken(value, line_no, i-1)))
                else
                    push!(tokens, IntegerToken(PositionToken(value, line_no, i-1)))
            end

        else


            push!(tokens, IdentifierToken(PositionToken(value, line_no, i-1)))
        end

        elseif isdigit(c)
            start = i
            while i <= length(source) && (isdigit(source[i]) || source[i] == '.')
                i += 1
            end
            value = source[start:i-1]
            push!(tokens, IntegerToken(PositionToken(value, line_no, i)))
            # push!(tokens, LiteralToken(PositionToken(value, line_no, i)))

        # TODO: What about +,-,/, and * symbols?
        elseif c == '-'
            if i + 1 <= length(source) && source[i + 1] == '>'
                push!(tokens, RightArrowToken(PositionToken("->", line_no, i)))
                i += 2
            elseif i + 1 <= length(source) && source[i + 1] == '-'
                push!(tokens, EdgeToken(PositionToken("--", line_no, i)))
                i += 2
            else
                # push!(tokens, ErrorToken(PositionToken("Unexpected character: $c", line_no, i)))
                push!(
                  tokens,
                  MinusToken(
                         PositionToken("-", line_no, i)
                     )
                  )
                i += 1
            end
        elseif c in ['!', '=', '<', '>', '*', '+', '/']

            if c == '<'
                if i + 1 <= length(source) && source[i + 1] == '<'
                    push!(tokens, LeftAngleBracketToken(PositionToken("<<", line_no, i)))
                    i += 2
                elseif i + 1 <= length(source) && source[i + 1] == '='
                    push!(tokens, LeftAngleBracketToken(PositionToken("<<", line_no, i)))
                    i += 2

                else
                    push!(tokens, OperatorToken(PositionToken("<=", line_no, i)))
                    i += 1
                end
            elseif c == '>'
                if i + 1 <= length(source) && source[i + 1] == '>'
                    push!(tokens, RightAngleBracketToken(PositionToken(">>", line_no, i)))
                    i += 2
                elseif i + 1 <= length(source) && source[i + 1] == '='
                    push!(tokens, OperatorToken(PositionToken(">=", line_no, i)))
                    i += 2
                else
                    check = string(c, '>')
                    push!(tokens, OperatorToken(PositionToken(string(c), line_no, i)))
                    i += 1
                end

            elseif i + 1 <= length(source) && source[i + 1] == '='
                push!(tokens, OperatorToken(PositionToken(string(c, '='), line_no, i)))
                i += 2
            elseif c == '='
                push!(tokens, EqualToken(PositionToken("=", line_no, i)))
                i += 2
            elseif c == '/'
                push!(tokens, SlashToken(PositionToken("/", line_no, i)))
                i += 1
            elseif c == '*'
                push!(tokens, AsteriskToken(PositionToken("*", line_no, i)))
                i += 1
            elseif c == '+'
                push!(tokens, PlusToken(PositionToken("+", line_no, i)))
                i += 1
            else
                push!(tokens, OperatorToken(PositionToken(string(c), line_no, i)))
                i += 1
            end

        elseif c in ['(', ')', '{', '}', ',', ';']

            if c == '('
                push!(tokens, LeftParenthesisToken(PositionToken(string(c), line_no, i)))
            elseif c == ')'
                push!(tokens, RightParenthesisToken(PositionToken(string(c), line_no, i)))
            elseif c == '{'
                push!(tokens, LeftBracketToken(PositionToken(string(c), line_no, i)))
            elseif c == '}'
                push!(tokens, RightBracketToken(PositionToken(string(c), line_no, i)))
            else
                push!(tokens, PunctuationToken(PositionToken(string(c), line_no, i)))
            end
            i += 1

        elseif c == '#'
            push!(tokens, CommentToken(PositionToken("#", line_no, i)))
            i += 1

        elseif c == '.'
            push!(tokens, DotToken(PositionToken(".", line_no, i)))
            i += 1

        elseif c == ':'

            if (source[i+1] == ':') && (i + 1 <= length(source))
                push!(tokens, DoubleColonToken(PositionToken("::", line_no, i)))
                i += 2
            elseif (source[i+1] == '=') && (i + 1 <= length(source))
                push!(tokens, DefineToken(PositionToken(":=", line_no, i)))
                i += 2
            else
                push!(tokens, SingleColonToken(PositionToken(":", line_no, i)))
                i += 1
            end

        else
            push!(tokens, ErrorToken(PositionToken("Unexpected character: $c", line_no, i)))
            i += 1
        end

    end
    return tokens
end
