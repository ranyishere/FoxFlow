include("../lexer.jl")
include("../parser.jl")
include("../utils/dot_generation.jl")


function tokenize_file_parser(file_name)
    lines_tokens = []
    open(file_name) do f
        line_no = 0
        while !eof(f)
            cur_line = readline(f)
            line_token = tokenize(cur_line, line_no+1)
            line_no += 1

            if line_token != Token[]
                push!(lines_tokens, line_token)
            end

        end
    end
    total_tokens = collect(Iterators.flatten(lines_tokens))
end

function test_parser_runner()

    lines_tokens = tokenize_file_parser("test_types_2.fflow")
    symbol_table = Dict()

    cur_ast = parse_types_section!(
                                   lines_tokens,
                                   symbol_table
                                  )
end

test_parser_runner()
