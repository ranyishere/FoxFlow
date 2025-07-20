include("./tokens.jl")
include("./lexer.jl")

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
            line_no += 1

            if line_token != Token[]
                push!(lines_tokens, line_token)
                push!(lines_tokens, EndLineToken())
            end

        end
    end
    lines_tokens
end


function lookahead(tokens)
    if length(tokens) > 1
        println("tokens[1]: ", tokens[1])
        return tokens[1]
    else
        println("Lookahead doesn't exist")
    end

end

function parse_symbol_name!(tokens)
    cur_token = popfirst!(tokens)

    if isa(cur_token,IdentifierToken)
        return cur_token
    else
        throw("Expected Identifer Token got: $cur_token")
    end
end

function parse_symbol_parameters!(tokens)
end

function parse_type_signature_list!(tokens)

    symbol_name = parse_symbol_name!(tokens)

    if isa(lookahead(tokens), RightArrowToken)

        popfirst!(tokens)

        cur_token = parse_type_signature_list!(tokens)

    elseif isa(lookahead(tokens), FloatToken) || isa(lookahead(tokens), FloatToken)
        return [popfirst!(tokens)]
    end


end


function parse_type_declaration!(tokens)

    type_name = parse_symbol_name!(tokens)
    symbol_parameters = parse_symbol_parameters!(tokens)

    if !isa(popfirst!(tokens), SingleColonToken)
        throw("Expected a SingleColonToken")
    end

    type_signature_list = parse_type_signature_list!(tokens)

end

function parse_type_declaration_list!(tokens)

    type_declaration = parse_type_declaration!(tokens)
    if isempty(tokens)
        return [type_declaration]
    end

    cur_token = popfirst!(tokens)
    if isa(cur_token, NewLineToken)
        return [type_declaration]
    else
        type_declaration_list = parse_type_declaration_list!(tokens)
        return [type_declaration;type_declaration_list]
    end

end

function parse_types_section!(tokens)

    if isa(lookahead(tokens), TypeSectionToken)

        popfirst!(tokens)

        type_section_name = parse_symbol_name!(tokens)

        cur_token = popfirst!(tokens)

        if isa(cur_token, LeftBracketToken)

            type_declarations = parse_type_declaration_list!(tokens)

            cur_token = popfirst!(tokens)

            if !isa(cur_token, RightBracketToken)
                throw("Error expected right bracket token got: $cur_token")
            end

        else
            throw("Error expected left bracket token got: $cur_token")
        end
    end

end

function main()

    # tokens = tokenize_file("examples/random_particle.fflow")
    # tokens = tokenize_file("tests/test_types_2.fflow")
    tokens = tokenize_file("tests/rules/test_rule_1.fflow")
    # parse_types_section!(tokens[1])
    # println("tokens: ", tokens)
end

main()
