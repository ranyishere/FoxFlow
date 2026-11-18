# Shared parser helpers: operator precedence, token utilities, tokenizer wrappers.

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
        error("Expected $(T), got $(typeof(token)), $token")
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

function skip_eol!(tokens)
    while !isempty(tokens) && isa(lookahead(tokens), EndLineToken)
        popfirst!(tokens)
    end
    return nothing
end
