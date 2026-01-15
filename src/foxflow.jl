#!/usr/bin/env julia
include("ir_generation/ir.jl")

function compile_foxflow(file_name::String)

    # Step 1: Lexical Analysis
    tokens = tokenize_file(file_name)
    println("Tokens:", tokens)
    exit(0)

    # Step 2: Parsing
    ast = parse(tokens)

    # Step 3: IR Generation
    ir = generate_ir(ast)
    return ir
end

function main()
    if length(ARGS) != 1
	println("Usage: julia foxflow.jl <source_file.fox>")
	exit(1)
    end

    source_file = ARGS[1]
    ir = compile_foxflow(source_file)
    # For demonstration, print the generated IR
    println("Generated Intermediate Representation:")
    println(ir)
end

main()
