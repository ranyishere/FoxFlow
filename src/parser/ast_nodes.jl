"""
This file defines the abstract syntax tree (AST) nodes for the FoxFlow language.
"""

module AstNodes

    import ..Tokens: Token, OperatorToken, MinusToken,
                    PlusToken, AsteriskToken,
                    SlashToken, NotToken, SampleToken, StringToken

    abstract type Node end
    abstract type LiteralNode <: Node end
    abstract type ModifyClauseNode <: Node end

    struct FloatNode <: LiteralNode
        token::Token
    end

    struct TimeTypeNode <: Node
        token::Token
    end

    Base.@kwdef struct IdentifierNode <: Node
        token::Token
        namespace::Union{Token, Nothing} = nothing
    end
    IdentifierNode(token::Token) = IdentifierNode(token, nothing)

    struct StringNode <: Node
        token :: StringToken
    end

    struct SymbolNode <: Node
        name::IdentifierNode
        parameters::Array{Node}
    end

    struct FixedListParameterNode <: Node
        name::IdentifierNode
        parameters::Array{Node}
    end

    struct NamedParameterNode <: Node
        name::IdentifierNode
        parameter::Node
    end

    struct ParameterNode <: Node
        token:: Union{Array{Node}, Nothing} # Can be a ParameterNode or a IdentifierNode
    end
    
    struct IntegerNode <: LiteralNode
        token::Token
    end

    struct TypeClassNode <: Node
        name :: Union{IdentifierNode, FloatNode, IntegerNode}
        parameter :: ParameterNode
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

    # Represents an inline array/tensor literal, e.g. [1,2,3] or [[1,2],[3,4]]
    struct ArrayLiteralNode <: Node
        elements :: Array{Node}
    end

    # Represents an index access operation, e.g. tensor[i] or tensor[i, j]
    struct IndexAccessNode <: Node
        object :: Node
        indices :: Array{Node}
    end

    # Represents a slice expression inside an index, e.g. 0:3, ::2, :, 1:5:2
    # Any of start, stop, step can be `nothing` to indicate omission.
    # Examples:
    #   :       -> SliceNode(nothing, nothing, nothing)  -- select all
    #   0:3     -> SliceNode(0, 3, nothing)               -- range [0, 3)
    #   ::2     -> SliceNode(nothing, nothing, 2)          -- every 2nd element
    #   1:5:2   -> SliceNode(1, 5, 2)                      -- range [1,5) step 2
    struct SliceNode <: Node
        start :: Union{Node, Nothing}
        stop  :: Union{Node, Nothing}
        step  :: Union{Node, Nothing}
    end

    struct CallNode <: Node
        function_node :: IdentifierNode
        args :: Array{
                      Union{Token, IntegerNode, FloatNode,
                            IdentifierNode, BinaryOpNode, GroupNode,
                            UnaryOpNode, CallNode, IndexAccessNode,
                            ArrayLiteralNode, StringNode
                           }
                     }
    end

    # Encapsulates any node in a namespace.
    # This is used to allow for nested namespaces and to allow for namespaces to contain any type of node.
    struct NamespaceNode <: Node
        name :: IdentifierNode
        value :: Node
    end

    struct FunctionNode <: Node
        name::IdentifierNode
        args::Array{Union{
                          Token, IntegerNode, FloatNode,
                          IdentifierNode, BinaryOpNode,
                          GroupNode, UnaryOpNode, IndexAccessNode,
                          ArrayLiteralNode
                         }}
    end

    struct GrammarSignatureNode <: Node
        token::Array{Token}
    end

    struct BindingVariableNode <: Node
        name::Union{IdentifierNode, IndexAccessNode}
        # What the derivative is in respect to.
        value::Array{Any}
    end

    struct ODENode <: Node
        name :: Union{IdentifierNode, IndexAccessNode}
        value :: Union{BinaryOpNode,
                       IdentifierNode,
                        IntegerNode,
                        FloatNode,
                        GroupNode, CallNode,
                        ArrayLiteralNode
                      }
    end

    struct DefinitionNode <: Node
        name :: IdentifierNode
        type  :: TypeClassNode
        value :: Node
    end

    struct SolveClauseNode <: ModifyClauseNode
        variables::Array{BindingVariableNode}
        clause::Array{Union{ODENode, DefinitionNode}}
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

    struct TypeInstanceNode <: Node
        name  :: IdentifierNode
        parameter :: ParameterNode
        type  :: Union{TypeClassNode, Nothing}
         # Value can be a list of types or a single value.
        value :: Union{
                        Token, Nothing,
                        Array, Node
                       }
    end

    struct TypeInstanceUpdateNode <: Node
        name :: Union{IdentifierNode, IndexAccessNode}
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

    # struct GrammarNode <: Node
        # time::TimeTypeNode
        # name::IdentifierNode
        # signature::GrammarSignatureNode
        # initial_conditions::InitialConditionListNode
        # rules::Array{RuleNode}
    # end

    struct TypeSectionNode <: Node
        name :: IdentifierNode
        types::Array{TypeInstanceNode}
    end

    struct ParameterSectionNode <: Node
        name :: IdentifierNode
        parameter_list :: Array{Union{ParameterNode, TypeInstanceNode}}
    end

    struct RuleSectionNode <: Node
        name :: IdentifierNode
        # rules_list ::Array{TypeInstanceNode}
        rules_list ::Array{RuleNode}
    end

    struct FunctionArgNode <: Node
        name :: IdentifierNode
        type :: TypeClassNode
    end

    # A destructured param:  name : TypeName << alias1 : T1, alias2 : T2, ... >>
    # The field aliases are in scope inside the body without dot access.
    # Reuses FunctionArgNode for each alias binding.
    struct ObservableDestructuredParamNode <: Node
        name      :: IdentifierNode                # e.g. p1
        type_name :: IdentifierNode                # e.g. Layer
        fields    :: Vector{FunctionArgNode}       # alias : TypeSig bindings
    end

    struct FunctionDefinitionExpressionNode <: Node
        name :: IdentifierNode
        type :: Union{IdentifierNode, FloatNode}
        value :: Node
    end

    struct FunctionBodyExpressionNode <: Node
        expressions :: Array{FunctionDefinitionExpressionNode}
    end

    struct ReturnNode <: Node
        value :: Union{Node, IdentifierNode, FloatNode, CallNode, BinaryOpNode}
    end

    # A local Function inside an Observable body:
    #   name : Function << (p1 : TypeName << alias1: T1, ... >>) >> -> RetType := { ... }
    struct ObservableLocalFunctionNode <: Node
        name        :: IdentifierNode
        param       :: ObservableDestructuredParamNode
        return_type :: IdentifierNode
        body        :: Union{FunctionBodyExpressionNode, Nothing}
        body_return :: ReturnNode
    end

    # A single observable definition:
    #   name : Observable << (p1 : TypeName << alias1: T1, ... >>) >> -> RetType := { ... }
    struct ObservableDefinitionNode <: Node
        name        :: IdentifierNode
        param       :: ObservableDestructuredParamNode
        return_type :: IdentifierNode
        local_fns   :: Vector{ObservableLocalFunctionNode}
        body_return :: ReturnNode
    end

    struct ObservableSectionNode <: Node
        name        :: IdentifierNode
        definitions :: Vector{ObservableDefinitionNode}
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

    
    struct AssignNode <: Node
        name :: IdentifierNode
        value :: Node
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



    struct FunctionSignatureNode <: Node
        args :: Array{FunctionArgNode}
        output :: TypeClassNode
    end


    # Represents a model loaded from a file path, e.g. load("policy_net.pt")
    struct ModelLoadNode <: Node
        filepath :: StringNode
    end

    struct FunctionDefinitionNode <: Node
        name :: IdentifierNode
        signature :: FunctionSignatureNode
        body :: Union{FunctionBodyExpressionNode, Nothing}
        fun_return :: Union{Node, Nothing}
        model_load :: Union{ModelLoadNode, Nothing}
    end

    struct FunctionSectionNode <: Node
        name :: IdentifierNode
        functions :: Array{FunctionDefinitionNode}
    end

    struct GrammarNode <: Node
        type_section :: TypeSectionNode
        parameter_section :: ParameterSectionNode
        function_section :: FunctionSectionNode
        rule_section :: RuleSectionNode
        observable_section :: ObservableSectionNode
    end

    struct FoxFlowNode <: Node
        grammars :: Array{GrammarNode}
        parameters :: ParameterSectionNode
    end

    struct LoadNode <: Node
        filepath :: StringNode
    end

    struct SimulationNode <: Node
        name :: IdentifierNode
        # rules :: Array{RuleNode}
        # types :: Array{TypeInstanceNode}
    end

    struct SimulationParametersNode <: Node
        value :: LoadNode
    end

    struct SimulationRulesNode <: Node
        value :: LoadNode
    end

    struct SimulationTypesNode <: Node
        value :: LoadNode
    end

    struct SimulationObservablesNode <: Node
        value :: LoadNode
    end

    struct SimulationStateNode <: Node
        value :: LoadNode
    end

    struct RunSimulationNode <: Node
        initialState :: IdentifierNode
        parameters :: IdentifierNode
        rules :: IdentifierNode
        types :: IdentifierNode
        steps :: Union{IntegerNode, FloatNode, IdentifierNode}
        observables :: Union{IdentifierNode, Nothing}
    end

    struct SimDeclarationNode <: Node 
        name :: IdentifierNode
        value :: Union{FloatNode, IntegerNode, SimulationParametersNode,
              SimulationNode, SimulationRulesNode, SimulationTypesNode,
              SimulationObservablesNode, SimulationStateNode, StringNode, RunSimulationNode}
    end

    struct SimulationSectionNode <: Node
        name :: IdentifierNode
        declarations :: Array{SimDeclarationNode}
    end
    
end
