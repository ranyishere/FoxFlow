
# Helper function: Peek at the next token without consuming it
function lookahead(tokens)
    if !isempty(tokens)
        return tokens[1]
    else
        throw("Lookahead requested on empty token list")
    end
end


