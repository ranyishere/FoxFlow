function parse_expression!(tokens)
    return parse_assignment!(tokens)
end

function parse_assignment!(tokens)  # lowest precedence, right-assoc
    left = parse_logical_or!(tokens)
    if !isempty(tokens) && (lookahead(tokens) isa EqualToken || lookahead(tokens) isa DefineToken)
        skip_eol!(tokens)
        op = popfirst!(tokens)               # '=' or compound like PlusEqToken if you add them
        right = parse_assignment!(tokens)    # right-associative
        skip_eol!(tokens)
        return AssignNode(op, left, right)
    end
    return left
end

function parse_logical_or!(tokens)
    left = parse_logical_and!(tokens)
    while !isempty(tokens) && (lookahead(tokens) isa OrOrToken)
        skip_eol!(tokens)
        op = popfirst!(tokens)
        right = parse_logical_and!(tokens)

        skip_eol!(tokens)
        left = BinaryOpNode(op, left, right)  # keep a distinct node type if you need short-circuit codegen
    end
    return left
end

function parse_logical_and!(tokens)
    left = parse_equality!(tokens)
    while !isempty(tokens) && (lookahead(tokens) isa AndAndToken)
        skip_eol!(tokens)
        op = popfirst!(tokens)
        right = parse_equality!(tokens)
        skip_eol!(tokens)
        left = BinaryOpNode(op, left, right)
    end
    return left
end

function parse_equality!(tokens)  # ==, !=
    left = parse_relational!(tokens)
    while !isempty(tokens) && (lookahead(tokens) isa EqEqToken || lookahead(tokens) isa NotEqToken)
        skip_eol!(tokens)
        op = popfirst!(tokens)
        right = parse_relational!(tokens)
        skip_eol!(tokens)
        left = BinaryOpNode(op, left, right)
    end
    return left
end

function parse_relational!(tokens)  # <, <=, >, >=
    left = parse_additive!(tokens)
    while !isempty(tokens) && (
        lookahead(tokens) isa LtToken || lookahead(tokens) isa LtEqToken ||
        lookahead(tokens) isa GtToken || lookahead(tokens) isa GtEqToken
    )
        skip_eol!(tokens)
        op = popfirst!(tokens)
        right = parse_additive!(tokens)
        skip_eol!(tokens)
        left = BinaryOpNode(op, left, right)
    end
    return left
end

# rename your existing parse_expression!/parse_term! to these:
function parse_additive!(tokens)    # +, -
    left = parse_multiplicative!(tokens)
    while !isempty(tokens) && (lookahead(tokens) isa PlusToken || lookahead(tokens) isa MinusToken)
        skip_eol!(tokens)
        op = popfirst!(tokens)
        right = parse_multiplicative!(tokens)
        skip_eol!(tokens)
        left = BinaryOpNode(op, left, right)
    end
    return left
end

function parse_exponential!(tokens)
    left = parse_factor!(tokens)
    while !isempty(tokens) && lookahead(tokens) isa CaretToken
        skip_eol!(tokens)
        op = popfirst!(tokens)
        # Right-associative, so we recursively call parse_exponential!
        right = parse_exponential!(tokens)
        skip_eol!(tokens)
        left = BinaryOpNode(op, left, right)
    end
    return left
end


function parse_multiplicative!(tokens)   # *, /
    left = parse_exponential!(tokens)
    while !isempty(tokens) && (lookahead(tokens) isa AsteriskToken || lookahead(tokens) isa SlashToken)
        skip_eol!(tokens)
        op = popfirst!(tokens)
        right = parse_factor!(tokens)
        skip_eol!(tokens)
        left = BinaryOpNode(op, left, right)
    end
    return left
end
