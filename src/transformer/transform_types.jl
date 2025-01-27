include("../ast_nodes.jl")


abstract type TransformedToken end

struct TransformedHistoryToken <: TransformedToken
    foxtoken :: Node
end

struct TransformedType <: TransformedToken
end

struct TransformedTypeToken <: TransformedToken
    history :: TransformedHistoryToken
    name :: String
    type :: String
    value :: Union{Nothing, String}
end

struct StructTransformedToken <: TransformedToken
    name :: String
    type :: String
end

struct_tmp = """
    struct {{name}} 
    {
        {{#:types}} {{.}} {{/:types}}
    };
"""

type_tmp = """
    {{type}} 
"""

function transform_type_section(type_section)
    println("type_section: $type_section")
end
