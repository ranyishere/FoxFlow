#!/usr/bin/env julia
"""
Tests for ArrayLiteralNode codegen.
Tests that inline tensor literals generate correct C++ torch::tensor() calls.
"""

include("../ir_generation/ir.jl")

using .AstNodes: ArrayLiteralNode, IntegerNode, FloatNode, BinaryOpNode, 
                 IdentifierNode, UnaryOpNode, GroupNode
using .IRRuleGeneration: ir_array_literal

function test_1d_codegen()
    tokens = tokenize_string("[1, 2, 3]")
    node = parse_expression!(tokens)
    result = ir_array_literal(node)
    expected = "torch::tensor({1, 2, 3}, torch::kFloat64)"
    @assert result == expected "Expected '$expected', got '$result'"
    println("✓ Test 1D codegen: $result")
end

function test_2d_codegen()
    tokens = tokenize_string("[[1, 2], [3, 4]]")
    node = parse_expression!(tokens)
    result = ir_array_literal(node)
    expected = "torch::tensor({{1, 2}, {3, 4}}, torch::kFloat64)"
    @assert result == expected "Expected '$expected', got '$result'"
    println("✓ Test 2D codegen: $result")
end

function test_3d_codegen()
    tokens = tokenize_string("[[[1, 2], [3, 4]], [[5, 6], [7, 8]]]")
    node = parse_expression!(tokens)
    result = ir_array_literal(node)
    expected = "torch::tensor({{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}}, torch::kFloat64)"
    @assert result == expected "Expected '$expected', got '$result'"
    println("✓ Test 3D codegen: $result")
end

function test_negative_values_codegen()
    tokens = tokenize_string("[-1, 2, -3]")
    node = parse_expression!(tokens)
    result = ir_array_literal(node)
    expected = "torch::tensor({-1, 2, -3}, torch::kFloat64)"
    @assert result == expected "Expected '$expected', got '$result'"
    println("✓ Test negative values codegen: $result")
end

function test_empty_codegen()
    tokens = tokenize_string("[]")
    node = parse_expression!(tokens)
    result = ir_array_literal(node)
    expected = "torch::tensor({}, torch::kFloat64)"
    @assert result == expected "Expected '$expected', got '$result'"
    println("✓ Test empty codegen: $result")
end

function test_identifiers_codegen()
    tokens = tokenize_string("[x, y, z]")
    node = parse_expression!(tokens)
    result = ir_array_literal(node)
    expected = "torch::tensor({x, y, z}, torch::kFloat64)"
    @assert result == expected "Expected '$expected', got '$result'"
    println("✓ Test identifiers codegen: $result")
end

function test_generated_rules_h()
    # Verify the full pipeline generates correct code for count[0] = 1, count[1] = 0
    rules_path = joinpath(@__DIR__, "generated_tests/generated_2/rules.h")
    if isfile(rules_path)
        content = read(rules_path, String)
        @assert occursin(".Direction[(0)] = (1);", content) "Expected indexed tensor assignment .Direction[(0)] = (1); in generated rules.h"
        @assert occursin(".Direction[(1)] = (0);", content) "Expected indexed tensor assignment .Direction[(1)] = (0); in generated rules.h"
        println("✓ Test generated rules.h contains correct indexed tensor assignments")
    else
        println("⚠ Skipping generated rules.h test (file not found)")
    end
end

function run_all_tests()
    println("Running ArrayLiteralNode codegen tests...")
    println()
    test_1d_codegen()
    test_2d_codegen()
    test_3d_codegen()
    test_negative_values_codegen()
    test_empty_codegen()
    test_identifiers_codegen()
    test_generated_rules_h()
    println()
    println("All 7 ArrayLiteralNode codegen tests passed! ✓")
end

run_all_tests()
