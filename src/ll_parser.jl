include("./lexer.jl")

mutable struct SymbolTable
    parent::Union{SymbolTable, Nothing}
    symbols::Dict{String, Any}
end

function create_symbol_table(parent::Union{SymbolTable, Nothing}=nothing)
    return SymbolTable(parent, Dict{String, Any}())
end

function add_symbol(table::SymbolTable, name::String, attributes::Dict)
    table.symbols[name] = attributes
end

function lookup_symbol(table::SymbolTable, name::String)
    while table !== nothing
        if haskey(table.symbols, name)
            return table.symbols[name]
        end
        table = table.parent
    end
    return nothing  # Not found
end



function llparse(tokens)

    while length(tokens) > 0
        cur_token = popfirst!(tokens)
        println("$cur_token")
        exit(0)
    end

end

function main()
    types_test = """types Plant {
        Node :: Type;
        dim :: Integer = 3;
        UnitNode :: Node = {
            unit_vec::FixedList<<dim, Float>>;
        };
        zipper :: UnitNode;
        I :: Float;
    };"""

    tokens = tokenize(types_test, 0)
    llparse(tokens)
end

main()
