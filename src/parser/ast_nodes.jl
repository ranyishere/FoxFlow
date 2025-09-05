"""
This file defines the abstract syntax tree (AST) nodes for the FoxFlow language.
"""

module AstNodes

    import ..Tokens: Token, OperatorToken, MinusToken,
                    PlusToken, AsteriskToken,
                    SlashToken, NotToken, SampleToken


    abstract type Node end
    abstract type LiteralNode <: Node end
    abstract type ModifyClauseNode <: Node end

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

    struct FixedListParameterNode <: Node
        name::IdentifierNode
        parameters::Array{Node}
    end

    struct ParameterNode <: Node
        token::Array{Node} # Can be a ParameterNode or a IdentifierNode
    end

    struct FloatNode <: LiteralNode
        token::Token
    end

    struct IntegerNode <: LiteralNode
        token::Token
    end

    struct BinaryOpNode <: Node
        expression :: OperatorToken
        lhs :: Node
        rhs :: Node
    end

    struct GroupNode <: Node
        expression :: Node
    end

    struct UnaryOpNode <: Node
        expression :: Union{MinusToken, PlusToken, NotToken, SampleToken}
        operand :: Node
    end

    struct FunctionNode <: Node
        name::IdentifierNode
        args::Array{Union{
                          Token, IntegerNode, FloatNode,
                          IdentifierNode, BinaryOpNode,
                          GroupNode, UnaryOpNode
                         }}
    end

    struct GrammarSignatureNode <: Node
        token::Array{Token}
    end

    
    struct BindingVariableNode <: Node
        name::IdentifierNode
        # What the derivative is in respect to.
        value::Array{Any}
    end

    struct ODENode <: Node
        name :: IdentifierNode
        value :: Union{BinaryOpNode,
                       IdentifierNode,
                        IntegerNode,
                        FloatNode
                      }
    end

    struct SolveClauseNode <: ModifyClauseNode
        variables::Array{BindingVariableNode}
        clause::Array{ODENode}
    end

    struct WhereClauseNode <: ModifyClauseNode
        clause::Array{Node}
    end

    struct WithClauseNode <: ModifyClauseNode
        # name::Token
        function_node:: Union{
                              FunctionNode,
                              BinaryOpNode,
                              IdentifierNode,
                              FloatNode,
                              IntegerNode, GroupNode
                             }
        # clause::Array{Token}
        where_clause::WhereClauseNode
    end

    struct TypeClassNode <: Node
        name :: IdentifierNode
        parameter :: ParameterNode
    end

    struct DefinitionNode <: Node
        name :: IdentifierNode
        type  :: TypeClassNode
        value :: Node
    end

    struct TypeInstanceNode <: Node
        name  :: IdentifierNode
        parameter :: ParameterNode
        type  :: TypeClassNode
         # Value can be a list of types or a single value.
        value :: Union{
                        Token, Nothing,
                        Array, Node
                       }
    end

    abstract type EdgeNode <: Node
    end

    struct UndirectedTypeEdgeNode <: EdgeNode
        left_vertex :: TypeInstanceNode
        right_vertex :: TypeInstanceNode
    end

    struct DirectedTypeEdgeNode <: EdgeNode
        left_vertex :: TypeInstanceNode
        right_vertex :: TypeInstanceNode
        direction :: String
    end

    struct RuleNode <: Node

        # Identifier Node
        name::IdentifierNode
        # lhs::Array{Token}
        lhs::Array{
                   Union{TypeInstanceNode,
                        UndirectedTypeEdgeNode,
                        DirectedTypeEdgeNode}
                  }
        lhs_parameter::ParameterNode

        # rhs::Array{Token}
        rhs::Array{
                   Union{TypeInstanceNode,
                        UndirectedTypeEdgeNode,
                        DirectedTypeEdgeNode}
                  }
        rhs_parameter::ParameterNode

        modify_clause::ModifyClauseNode
    end

    struct InitialConditionNode <: Node
        name::Token
        value::Token
    end

    struct InitialConditionListNode <: Node
        value::Array{InitialConditionNode}
    end

    struct GrammarNode <: Node
        time::TimeTypeNode
        name::IdentifierNode
        signature::GrammarSignatureNode
        initial_conditions::InitialConditionListNode
        rules::Array{RuleNode}
    end

    struct TypeSectionNode <: Node
        name :: IdentifierNode
        types::Array{TypeInstanceNode}
    end

    struct ParameterSectionNode <: Node
        name :: IdentifierNode
        parameter_list :: Array{Union{ParameterNode, TypeInstanceNode}}
    end

    struct FunctionSectionNode <: Node
        name :: IdentifierNode
        # types ::Array{TypeInstanceNode}
    end

    struct RuleSectionNode <: Node
        name :: IdentifierNode
        # rules_list ::Array{TypeInstanceNode}
        rules_list ::Array{RuleNode}
    end

    struct ObservableSectionNode <: Node
        name :: IdentifierNode
        # types ::Array{TypeInstanceNode}
    end

    struct GrammarSectionNode <: Node
        name :: IdentifierNode
        # types ::Array{TypeInstanceNode}
    end

    struct ExpressionNode <: Node
        expression :: OperatorToken
        lhs :: Token
        rhs :: Token
    end

    struct MinusNode <: Node
        expression :: OperatorToken
        lhs :: Token
    end

    struct PlusNode <: Node
        expression :: OperatorToken
        lhs :: Token
    end

    struct CallNode <: Node
        function_node :: FunctionNode
        args :: Array{
                      Union{Token, IntegerNode, FloatNode,
                            IdentifierNode, BinaryOpNode, GroupNode,
                            UnaryOpNode, CallNode
                           }
                     }
    end

    struct MultiplyNode <: Node
        expression :: OperatorToken
        lhs :: Token
        rhs :: Token
    end

    struct DivideNode <: Node
        expression :: OperatorToken
        lhs :: Token
        rhs :: Token
    end

    
    struct FoxFlowNode <: Node
        type_section :: TypeSectionNode
        parameter_section :: ParameterSectionNode
        function_section :: FunctionSectionNode
        rule_section :: RuleSectionNode
        observable_section :: ObservableSectionNode
        grammar_section :: GrammarSectionNode
    end

end
