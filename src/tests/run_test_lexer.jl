include("../lexer.jl")


function tokenize_file(file_name)
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
    lines_tokens
end

function test_comment_lexer()
    lines_tokens = tokenize_file("test_comment_lexer.fflow")
    correct = [
        LeftParenthesisToken(PositionToken("(", 1, 1))
        CommentToken(PositionToken("#", 1, 2))
        IdentifierToken(PositionToken("Test", 1, 7))
        CommentToken(PositionToken("#", 1, 9))
        RightParenthesisToken(PositionToken(")", 1, 10))
     ]

    each_line_token = lines_tokens[1]
    for (ix, each_token) in enumerate(each_line_token)
        check = each_token == correct[ix]
        if check == false
            println("Expected $(correct[ix]) got $each_token")
        end
    end

end

function test_section_header_lexer()

    lines_tokens = tokenize_file(
        "test_sections_lexer.fflow"
    )

    correct = [
        TypeSectionToken(PositionToken("types", 1, 5)),
        LeftBracketToken(PositionToken("{", 1, 7)),
        RightBracketToken(PositionToken("}", 2, 1)),
        ParameterSectionToken(PositionToken("parameters", 4, 10)),
        LeftBracketToken(PositionToken("{", 4, 12)),
        RightBracketToken(PositionToken("}", 5, 1)),
        FunctionSectionToken(PositionToken("functions", 7, 9)),
        LeftBracketToken(PositionToken("{", 7, 11)),
        RightBracketToken(PositionToken("}", 8, 1)),
        RuleSectionToken(PositionToken("rules", 10, 5)),
        LeftBracketToken(PositionToken("{", 10, 7)),
        RightBracketToken(PositionToken("}", 11, 1)),
        ObservableSectionToken(PositionToken("observables", 13, 11)),
        LeftBracketToken(PositionToken("{", 13, 13)),
        RightBracketToken(PositionToken("}", 14, 1)),
        StateSectionToken(PositionToken("states", 16, 6)),
        LeftBracketToken(PositionToken("{", 16, 8)),
        RightBracketToken(PositionToken("}", 17, 1)),
        SimulationSectionToken(PositionToken("simulations", 19, 11)),
        LeftBracketToken(PositionToken("{", 19, 13)),
        RightBracketToken(PositionToken("}", 20, 1)) 
    ]

    total_tokens = collect(Iterators.flatten(lines_tokens))

    for (ix, each_token) in enumerate(total_tokens)
        check = each_token == correct[ix]
        if check == false
            println("Expected $(correct[ix]) got $each_token")
        end
    end

end

function test_lexer_types()

    lines_tokens = tokenize_file("test_types.fflow")

    correct = [
               Token[TypeSectionToken(PositionToken("types", 1, 5)), IdentifierToken(PositionToken("Plant", 1, 11)), LeftBracketToken(PositionToken("{", 1, 13))],
Token[IdentifierToken(PositionToken("Node", 2, 8)), DoubleColonToken(PositionToken("::", 2, 10)), IdentifierToken(PositionToken("Type", 2, 16)), PunctuationToken(PositionToken(";", 2, 17))],
Token[IdentifierToken(PositionToken("dim", 3, 7)), DoubleColonToken(PositionToken("::", 3, 9)), IdentifierToken(PositionToken("Integer", 3, 18)), OperatorToken(PositionToken("=", 3, 20)), LiteralToken(PositionToken("3", 3, 23)), PunctuationToken(PositionToken(";", 3, 23))],
Token[IdentifierToken(PositionToken("UnitNode", 4, 12)), DoubleColonToken(PositionToken("::", 4, 14)), IdentifierToken(PositionToken("Node", 4, 20)), OperatorToken(PositionToken("=", 4, 22)), LeftBracketToken(PositionToken("{", 4, 24))],
Token[IdentifierToken(PositionToken("unit_vec", 5, 16)), DoubleColonToken(PositionToken("::", 5, 17)), IdentifierToken(PositionToken("FixedList", 5, 27)), LeftAngleBracketToken(PositionToken("<<", 5, 28)), IdentifierToken(PositionToken("dim", 5, 32)), PunctuationToken(PositionToken(",", 5, 33)), IdentifierToken(PositionToken("Float", 5, 39)), RightAngleBracketToken(PositionToken(">>", 5, 40)), PunctuationToken(PositionToken(";", 5, 42))],
Token[RightBracketToken(PositionToken("}", 6, 5)), PunctuationToken(PositionToken(";", 6, 6))],
Token[IdentifierToken(PositionToken("zipper", 7, 10)), DoubleColonToken(PositionToken("::", 7, 12)), IdentifierToken(PositionToken("UnitNode", 7, 22)), PunctuationToken(PositionToken(";", 7, 23))],
Token[IdentifierToken(PositionToken("I", 8, 5)), DoubleColonToken(PositionToken("::", 8, 7)), IdentifierToken(PositionToken("Float", 8, 14)), PunctuationToken(PositionToken(";", 8, 15))],
Token[RightBracketToken(PositionToken("}", 9, 1)), PunctuationToken(PositionToken(";", 9, 2))]
      ]

    for (ix, each_token) in enumerate(lines_tokens)
	# println("each_token: $each_token")
        check = each_token == correct[ix]
        if check == false
            println("Expected $(correct[ix]) got $each_token")
        end
    end

end

function test_lexer_runner()
    println("Running Tests")
    println("Comment Lexer Test")
    test_comment_lexer()
    println("Section Header Lexer Test")
    test_section_header_lexer()
    println("Types Section Lexer Test")
    test_lexer_types()
end

test_lexer_runner()
