import ..Tokens:
    Token, IdentifierToken, IntegerToken, FloatToken, EndLineToken,
    LeftParenthesisToken, RightParenthesisToken, LeftAngleBracketToken,
    RightAngleBracketToken, SingleColonToken, DefineToken, LeftBracketToken,
    RightBracketToken, PunctuationToken, TypeSectionToken, ParameterSectionToken,
    RuleSectionToken, SlashToken, AsteriskToken, PlusToken, MinusToken, RightArrowToken, WithToken, WhereToken, SolvingToken,
    EqualToken, CommaToken, EdgeToken, ODEToken, SampleToken, BackslashToken,
    LtToken, GtToken, LtEqToken, GtEqToken, EqEqToken, NotEqToken,
    NotToken, AndAndToken, OrOrToken, TypeToken, FunctionSectionToken, FunctionToken, CaretToken, ReturnToken,
    OperatorToken,
    StateToken, SimulationTypesToken, RulesToken, DotToken, ParameterToken, SimulationToken, RunSimulationToken, QuoteToken, SimulationSectionToken, LoadFileToken, StringToken, SimulationParametersToken, DoubleColonToken,
    LeftSquareBracketToken, RightSquareBracketToken, ObservableSectionToken, ObservableToken, SimulationObservablesToken
import ..AstNodes:
    Node, IdentifierNode, ParameterNode, TypeInstanceNode, TypeInstanceUpdateNode, TypeClassNode,
    TypeSectionNode, ParameterSectionNode, RuleSectionNode, RuleNode,
    WhereClauseNode, WithClauseNode, EveryClauseNode, SolveClauseNode, CallNode,
    GroupNode, BinaryOpNode, FunctionNode, FloatNode, IntegerNode, UndirectedTypeEdgeNode,
    BindingVariableNode, ODENode, UnaryOpNode, DefinitionNode, FunctionArgNode, FunctionSignatureNode,
    FunctionBodyExpressionNode, FunctionDefinitionExpressionNode, FunctionDefinitionNode,
    ReturnNode, FunctionSectionNode, SimulationSectionNode, StringNode, LoadNode, SimDeclarationNode,
    SimulationNode, SimulationRulesNode, SimulationTypesNode, SimulationObservablesNode, SimulationStateNode, SimulationParametersNode,
    RunSimulationNode, NamedParameterNode, IndexAccessNode, ArrayLiteralNode, ModelLoadNode, SliceNode,
    ObservableDestructuredParamNode, ObservableLocalFunctionNode, ObservableDefinitionNode, ObservableSectionNode

include("expression_parser.jl")
include("function_parser.jl")
include("parser_utils.jl")
include("sim_parser.jl")
include("parser_common.jl")
include("type_parser.jl")
include("expression_core.jl")
include("rule_parser.jl")
include("observable_parser.jl")
include("parameter_parser.jl")

function parse_function_section!(tokens)
    popfirst!(tokens)
end

# Parse the entire file into AST nodes
function parse_file!(tokens)

    ast_nodes = []
    while !isempty(tokens)

        cur_token = lookahead(tokens)

        if isa(cur_token, TypeSectionToken)
            push!(ast_nodes, parse_types_section!(tokens))
        elseif isa(cur_token, ParameterSectionToken)
            push!(ast_nodes, parse_parameters_section!(tokens))

        elseif isa(cur_token, RuleSectionToken)
            push!(
                  ast_nodes,
                  parse_rules_section!(tokens)
            )
        elseif isa(cur_token, FunctionSectionToken)
            push!(ast_nodes, parse_functions_section!(tokens))

        elseif isa(cur_token, SimulationSectionToken)
            push!(ast_nodes, parse_simulations_section!(tokens))

        elseif isa(cur_token, ObservableSectionToken)
            push!(ast_nodes, parse_observables_section!(tokens))

        elseif isa(cur_token, EndLineToken)
            popfirst!(tokens)
        else
            throw("Unexpected token: $cur_token")
        end
    end
    return ast_nodes
end
