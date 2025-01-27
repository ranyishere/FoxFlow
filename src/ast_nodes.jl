include("./tokens.jl")


abstract type Node end

struct TimeTypeNode <: Node
    token::Token
end

struct IdentifierNode <: Node
    token::Token
end

struct SymbolNode <: Node
    name::IdentifierNode
    parameters::Array{Node}
end

struct ParameterNode <: Node
    token::Array{Node}
end

struct GrammarSignatureNode <: Node
    token::Array{Token}
end

abstract type ModifyClauseNode <: Node
end

struct SolveClauseNode <: ModifyClauseNode
    name::Array{Token}
    clause::Array{Token}
end

struct WithClauseNode <: ModifyClauseNode
    name::Array{Token}
    clause::Array{Token}
end

struct RuleNode <: Node
    lhs::Array{Token}
    rhs::Array{Token}
    modify_clause::ModifyClauseNode
end

struct InitialConditionNode <: Node
    name::Token
    value::Token
end

struct InitialConditionListNode <: Node
    value::Array{InitialConditionNode}
end

struct TypeClassNode <: Node
    name :: IdentifierNode
    parameter :: ParameterNode
end

struct GrammarNode <: Node
    time::TimeTypeNode
    name::IdentifierNode
    signature::GrammarSignatureNode
    initial_conditions::InitialConditionListNode
    rules::Array{RuleNode}
end

struct TypeInstanceNode <: Node
    name  :: IdentifierNode
    parameter :: ParameterNode
    type  :: TypeClassNode
     # Value can be a list of types or a single value.
    value :: Union{
                    Token, Nothing,
                    Array
                   }
end

struct TypeSectionNode <: Node
    name :: IdentifierNode
    types::Array{TypeInstanceNode}
end
