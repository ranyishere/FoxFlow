#!/usr/bin/env julia
"""
Tests for ArrayLiteralNode parsing.
Tests that inline tensor literals like [1,2,3] and [[1,2],[3,4]] are correctly parsed.
"""

# Load the full FoxFlow pipeline
include("../ir_generation/ir.jl")

using .AstNodes: ArrayLiteralNode, IntegerNode, FloatNode, BinaryOpNode, 
                 IdentifierNode, UnaryOpNode

function tokenize_line(line)
    tokens = tokenize(line, 1)
    return tokens
end

function test_1d_array_literal()
    tokens = tokenize_line("[1, 2, 3]")
    result = parse_expression!(tokens)
    @assert result isa ArrayLiteralNode "Expected ArrayLiteralNode, got $(typeof(result))"
    @assert length(result.elements) == 3 "Expected 3 elements, got $(length(result.elements))"
    @assert result.elements[1] isa IntegerNode "Expected IntegerNode"
    @assert result.elements[2] isa IntegerNode "Expected IntegerNode"
    @assert result.elements[3] isa IntegerNode "Expected IntegerNode"
    println("✓ Test 1D array literal [1, 2, 3] passed")
end

function test_2d_array_literal()
    tokens = tokenize_line("[[1, 2], [3, 4]]")
    result = parse_expression!(tokens)
    @assert result isa ArrayLiteralNode "Expected ArrayLiteralNode, got $(typeof(result))"
    @assert length(result.elements) == 2 "Expected 2 elements, got $(length(result.elements))"
    @assert result.elements[1] isa ArrayLiteralNode "Expected nested ArrayLiteralNode"
    @assert result.elements[2] isa ArrayLiteralNode "Expected nested ArrayLiteralNode"
    @assert length(result.elements[1].elements) == 2 "Expected 2 inner elements"
    @assert length(result.elements[2].elements) == 2 "Expected 2 inner elements"
    println("✓ Test 2D array literal [[1,2],[3,4]] passed")
end

function test_array_with_floats()
    tokens = tokenize_string("[1.5, 2.0, 3.14]")
    result = parse_expression!(tokens)
    @assert result isa ArrayLiteralNode "Expected ArrayLiteralNode, got $(typeof(result))"
    @assert length(result.elements) == 3 "Expected 3 elements, got $(length(result.elements))"
    # Note: lexer tokenizes decimal numbers as IntegerToken (with value containing '.')
    @assert result.elements[1] isa IntegerNode "Expected IntegerNode (lexer treats decimals as IntegerToken)"
    println("✓ Test array literal with floats [1.5, 2.0, 3.14] passed")
end

function test_array_with_expressions()
    tokens = tokenize_line("[1 + 2, 3 * 4]")
    result = parse_expression!(tokens)
    @assert result isa ArrayLiteralNode "Expected ArrayLiteralNode, got $(typeof(result))"
    @assert length(result.elements) == 2 "Expected 2 elements"
    @assert result.elements[1] isa BinaryOpNode "Expected BinaryOpNode, got $(typeof(result.elements[1]))"
    @assert result.elements[2] isa BinaryOpNode "Expected BinaryOpNode, got $(typeof(result.elements[2]))"
    println("✓ Test array literal with expressions [1+2, 3*4] passed")
end

function test_array_with_identifiers()
    tokens = tokenize_line("[x, y, z]")
    result = parse_expression!(tokens)
    @assert result isa ArrayLiteralNode "Expected ArrayLiteralNode, got $(typeof(result))"
    @assert length(result.elements) == 3 "Expected 3 elements"
    @assert result.elements[1] isa IdentifierNode "Expected IdentifierNode"
    println("✓ Test array literal with identifiers [x, y, z] passed")
end

function test_array_with_negatives()
    tokens = tokenize_line("[-1, 2, -3]")
    result = parse_expression!(tokens)
    @assert result isa ArrayLiteralNode "Expected ArrayLiteralNode, got $(typeof(result))"
    @assert length(result.elements) == 3 "Expected 3 elements"
    @assert result.elements[1] isa UnaryOpNode "Expected UnaryOpNode for -1"
    @assert result.elements[3] isa UnaryOpNode "Expected UnaryOpNode for -3"
    println("✓ Test array literal with negatives [-1, 2, -3] passed")
end

function test_empty_array()
    tokens = tokenize_line("[]")
    result = parse_expression!(tokens)
    @assert result isa ArrayLiteralNode "Expected ArrayLiteralNode, got $(typeof(result))"
    @assert length(result.elements) == 0 "Expected 0 elements, got $(length(result.elements))"
    println("✓ Test empty array literal [] passed")
end

function test_3d_nested_array()
    tokens = tokenize_line("[[[1, 2], [3, 4]], [[5, 6], [7, 8]]]")
    result = parse_expression!(tokens)
    @assert result isa ArrayLiteralNode "Expected ArrayLiteralNode"
    @assert length(result.elements) == 2 "Expected 2 top-level elements"
    @assert result.elements[1] isa ArrayLiteralNode "Expected nested ArrayLiteralNode"
    @assert result.elements[1].elements[1] isa ArrayLiteralNode "Expected doubly nested ArrayLiteralNode"
    @assert length(result.elements[1].elements[1].elements) == 2 "Expected 2 inner elements"
    println("✓ Test 3D nested array literal [[[1,2],[3,4]],[[5,6],[7,8]]] passed")
end

function run_all_tests()
    println("Running ArrayLiteralNode parser tests...")
    println()
    test_1d_array_literal()
    test_2d_array_literal()
    test_array_with_floats()
    test_array_with_expressions()
    test_array_with_identifiers()
    test_array_with_negatives()
    test_empty_array()
    test_3d_nested_array()
    println()
    println("All 8 ArrayLiteralNode parser tests passed! ✓")
end

run_all_tests()
