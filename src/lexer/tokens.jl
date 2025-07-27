"""
Defines Tokens
"""

module Tokens

    abstract type Token end

    struct PositionToken <: Token
        value::String
        line_no::Integer
        col_no::Integer
    end

    struct KeywordToken <: Token
        position::PositionToken
    end

    struct GrammarToken <: Token
        position::PositionToken
    end

    struct GrammarTimeToken <: Token
        position::PositionToken
    end

    struct WithToken <: Token
        position::PositionToken
    end
    struct SolvingToken <: Token
        position::PositionToken
    end

    struct WhereToken <: Token
        position::PositionToken
    end

    struct IdentifierToken <: Token
        position::PositionToken
    end

    struct FunctionToken <: Token
        position::PositionToken
    end

    struct OperatorToken <: Token
        position::PositionToken
    end

    struct PlusToken <: Token
        position::PositionToken
    end

    struct MinusToken <: Token
        position::PositionToken
    end

    struct AsteriskToken <: Token
        position::PositionToken
    end

    struct SlashToken <: Token
        position::PositionToken
    end

    struct RightArrowToken <: Token
        position::PositionToken
    end

    struct EdgeToken <: Token
        position::PositionToken
    end

    struct PunctuationToken <: Token
        position::PositionToken
    end

    struct LeftParenthesisToken <: Token
        position::PositionToken
    end

    struct LeftAngleBracketToken <: Token
        position::PositionToken
    end

    struct RightAngleBracketToken <: Token
        position::PositionToken
    end

    struct RightParenthesisToken <: Token
        position::PositionToken
    end

    struct LeftBracketToken <: Token
        position::PositionToken
    end

    struct RightBracketToken <: Token
        position::PositionToken
    end

    struct LiteralToken <: Token
        position::PositionToken
    end

    struct ErrorToken <: Token
        position::PositionToken
    end

    struct InitialConditionToken <: Token
        position::PositionToken
    end

    struct CommentToken <: Token
        position::PositionToken
    end

    struct TypeSectionToken <: Token
        position::PositionToken
    end

    struct TypeToken <: Token
        position::PositionToken
    end

    struct DoubleColonToken <: Token
        position::PositionToken
    end

    struct SingleColonToken <: Token
        position::PositionToken
    end

    struct AssignToken <: Token
        position::PositionToken
    end

    struct DefineToken <: Token
        position::PositionToken
    end

    struct ParameterSectionToken <: Token
        position::PositionToken
    end

    struct ParameterToken <: Token
        position::PositionToken
    end

    struct FunctionSectionToken <: Token
        position::PositionToken
    end

    struct FunctionToken <: Token
        position::PositionToken
    end

    struct RuleSectionToken <: Token
        position::PositionToken
    end

    struct RuleToken <: Token
        position::PositionToken
    end

    struct ObservableSectionToken <: Token
        position::PositionToken
    end

    struct ObservableToken <: Token
        position::PositionToken
    end

    struct StateSectionToken <: Token
        position::PositionToken
    end

    struct StateToken <: Token
        position::PositionToken
    end

    struct SimulationSectionToken <: Token
        position::PositionToken
    end

    struct SimulationToken <: Token
        position::PositionToken
    end

    struct DotToken <: Token
        position::PositionToken
    end

    struct FloatToken <: Token
        position::PositionToken
    end

    struct IntegerToken <: Token
        position::PositionToken
    end

    struct CommaToken <: Token
        position::PositionToken
    end

    struct EqualToken <: Token
        position::PositionToken
    end

    struct EndLineToken <: Token
    end
end
