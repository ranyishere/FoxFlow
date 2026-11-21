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

    # One entry of a PDE `solving` block, e.g.
    #   du   : PDE       = kappa * laplacian(u)
    #   mesh : Mesh      = interval(L, n)
    #   left : Dirichlet = boundary(left, u_left)
    #   flux : Neumann   = boundary(top, g)
    #   init : Initial   = sin(pi * x / L)
    # `kind` is one of PDE_ENTRY_KINDS (contextual identifiers, not lexer keywords).
    struct PDEEntryNode <: Node
        name :: IdentifierNode
        kind :: String
        value :: Node
    end

    const PDE_ENTRY_KINDS = ("PDE", "Element", "Dirichlet", "Neumann", "Mesh", "Initial")

    struct SolveClauseNode <: ModifyClauseNode
        variables::Array{BindingVariableNode}
        clause::Array{Union{ODENode, DefinitionNode, PDEEntryNode}}
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

    # Clocked rule: `every (n_r, p_r, T_clock) [additive] where { ... }`.
    # Any of the three schedule expressions may be omitted (`nothing`).
    struct EveryClauseNode <: ModifyClauseNode
        n_r     :: Union{Node, Nothing}
        p_r     :: Union{Node, Nothing}
        t_clock :: Union{Node, Nothing}
        additive :: Bool
        where_clause::WhereClauseNode
    end

    # Sub-grammar call: `with (propensity) via g(t, unpack: <graph> [<<params>>] [where {..}], pack: <graph> [<<params>>] [where {..}])`.
    # A stochastic rule (fires with `function_node` as its propensity): `unpack` rewrites the rule's
    # LHS match into the starting graph of grammar `g`; `g` runs for `t` time units; `pack` rewrites
    # the final graph back onto the rule's RHS.
    struct ViaClauseNode <: ModifyClauseNode
        function_node    :: Union{
                              FunctionNode,
                              BinaryOpNode,
                              IdentifierNode,
                              FloatNode,
                              IntegerNode, GroupNode
                             }
        grammar          :: IdentifierNode
        duration         :: Node
        unpack_graph     :: Array
        unpack_parameter :: ParameterNode
        unpack_where     :: WhereClauseNode
        pack_graph       :: Array
        pack_parameter   :: ParameterNode
        pack_where       :: WhereClauseNode
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

    # A global observable: a functional of the whole state,
    #   name : Observable := <reduction> <pattern> [<< bindings >>] [with ( expr )]
    # The pattern is stored as a synthesized RuleNode (rhs == lhs, the `with`
    # expression as the propensity) so it is generated through the rule machinery:
    # each match of the pattern contributes the value of the expression, and
    # `reduction` (sum, mean, count, min, max) combines the contributions.
    struct GlobalObservableNode <: Node
        name      :: IdentifierNode
        reduction :: String
        rule      :: RuleNode
        # `within r`: a pair observable over two unconnected nodes closer than r (nothing otherwise)
        within    :: Union{Nothing, Float64}
    end
    GlobalObservableNode(name, reduction, rule) = GlobalObservableNode(name, reduction, rule, nothing)

    const OBSERVABLE_REDUCTIONS = ("sum", "mean", "count", "min", "max")

    # `definitions` are the per-node `Attribute` definitions (ParaView columns);
    # `global_observables` are the whole-state `Observable` definitions.
    struct ObservableSectionNode <: Node
        name        :: IdentifierNode
        definitions :: Vector{ObservableDefinitionNode}
        global_observables :: Vector{GlobalObservableNode}
    end
    ObservableSectionNode(name, definitions) =
        ObservableSectionNode(name, definitions, GlobalObservableNode[])

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

    # `f : SimulationFunctions := load("functions.fflow")`
    struct SimulationFunctionsNode <: Node
        value :: LoadNode
    end

    # `g : Grammar := Grammar(functions, rules, types [, observables])`: the compile-time unit of a
    # simulation. Each argument names a SimulationFunctions / SimulationRules / SimulationTypes /
    # SimulationObservables declaration (in any order).
    struct SimulationGrammarNode <: Node
        parts :: Vector{IdentifierNode}
    end

    # `s : Simulation := Simulation(["save.bin",] grammar, state, params, time [, checkpoint] [, dim = 2|3])`
    struct SimulationConfigNode <: Node
        saveName :: Union{StringNode, Nothing}
        grammar :: IdentifierNode
        initialState :: IdentifierNode
        parameters :: IdentifierNode
        steps :: Union{IntegerNode, FloatNode, IdentifierNode}
        checkpointFrequency :: Union{IntegerNode, IdentifierNode, Nothing}
        # Dimension of the expanded cell complex (2 or 3); `nothing` = the default (3)
        dim :: Union{Int, Nothing}
    end

    # `State := grid(...)`: an initial graph produced by a generator call (see ir_mesh_generator)
    struct SimulationGeneratorStateNode <: Node
        value :: Node
    end

    struct RunSimulationNode <: Node
        saveName :: Union{StringNode, Nothing}
        initialState :: IdentifierNode
        parameters :: IdentifierNode
        rules :: IdentifierNode
        types :: IdentifierNode
        steps :: Union{IntegerNode, FloatNode, IdentifierNode}
        # Optional: write output/state every N steps (and always on the final step).
        checkpointFrequency :: Union{IntegerNode, IdentifierNode, Nothing}
        observables :: Union{IdentifierNode, Nothing}
    end

    struct SimDeclarationNode <: Node 
        name :: IdentifierNode
        value :: Union{FloatNode, IntegerNode, SimulationParametersNode,
              SimulationNode, SimulationRulesNode, SimulationTypesNode,
              SimulationObservablesNode, SimulationStateNode, SimulationGeneratorStateNode,
              StringNode, RunSimulationNode, SimulationFunctionsNode, SimulationGrammarNode,
              SimulationConfigNode}
    end

    # `trace(simulationName, "rule1", "rule2", ...)`: record every fire of the named rules
    # while `simulationName` (a RunSimulation declaration) runs, and render the before/after
    # report. A rule name of "*" selects every rule.
    struct TraceNode <: Node
        simulation :: IdentifierNode
        rules :: Vector{StringNode}
    end

    # `measure(sim1, sim2, ..., "observable", ...)`: after the run, re-evaluate the global
    # observables on every saved checkpoint of each named simulation. No observable names (or
    # "*") means all observables.
    struct MeasureNode <: Node
        simulations :: Vector{IdentifierNode}
        observables :: Vector{StringNode}
    end

    # `run(sim, N [, domain_parallel = true|false] [, jobs = K] [, seed = S] [, name = "exp"]
    #  [, listeners = true|false] [, log = true|false])`: execute the
    # Simulation `sim` N times (a Simulation only runs when a run statement names it).
    # `domain_parallel = true` compiles in the OpenMP geocell domain-decomposition runner
    # (default: serial runner); with it off the replicates run concurrently, `jobs` at a time.
    # `name` names the experiment: its runs are written to runs/<name>/run_<i>/ (also for N = 1).
    struct RunNode <: Node
        simulation :: IdentifierNode
        replicates :: Union{IntegerNode, IdentifierNode}
        domainParallel :: Union{Bool, Nothing}
        jobs :: Union{IntegerNode, IdentifierNode, Nothing}
        seed :: Union{IntegerNode, IdentifierNode, Nothing}
        name :: Union{StringNode, Nothing}
        listeners :: Union{Bool, Nothing}   # listeners = true|false: liveliness/cell listener files
        log :: Union{Bool, Nothing}         # log = true|false: standard output of the run
    end

    struct SimulationSectionNode <: Node
        name :: IdentifierNode
        declarations :: Array{Union{SimDeclarationNode, TraceNode, MeasureNode, RunNode}}
    end
    
end
