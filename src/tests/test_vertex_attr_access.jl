#
# Focused test for vertex attribute access in rule parameter blocks:
#   << (rho1 : c1::rho, rho2 : c2::rho) >>
#
# Validates the lexer (`::`), parser (NamedParameterNode with a namespaced
# IdentifierNode value), and the IR generation of a `with` rule that uses
# these bindings in the where body.
#
include("../ir_generation/ir.jl")

using .Tokens: Token, EndLineToken
using .AstNodes: NamedParameterNode, IdentifierNode
using .IRUtils: get_value
using .IRRuleGeneration: ir_rules_section!

function tokenize_lines(code)
    lines_tokens = []
    line_no = 0
    for cur_line in split(code, '\n')
        cur_line = String(lstrip(cur_line))
        if cur_line == "" || cur_line[1] == '#'
            line_no += 1
            continue
        end
        line_token = tokenize(cur_line, line_no + 1)
        line_no += 1
        if line_token != Token[]
            lines_tokens = [lines_tokens; line_token]
            push!(lines_tokens, EndLineToken())
        end
    end
    lines_tokens
end

# --- 1. Parser test ------------------------------------------------------
rules_src = """
rules Fluids {
    equalize_density := (c1 : Fluid) -- (c2 : Fluid)
        << (rho1 : c1::rho, rho2 : c2::rho) >>
        ->
        (c1 : Fluid) -- (c2 : Fluid)
        << (rho1 : c1::rho, rho2 : c2::rho) >>
        with (1.0) where {
            rho1 = (rho1 + rho2) / 2.0
        }
}
"""

tokens = tokenize_lines(rules_src)
ast = parse_rules_section!(tokens)
rule = ast.rules_list[1]

println("Rule name: ", get_value(rule.name))

lhs_params = rule.lhs_parameter.token
@assert length(lhs_params) == 2 "Expected 2 LHS bindings, got $(length(lhs_params))"

for p in lhs_params
    @assert p isa NamedParameterNode "Expected NamedParameterNode, got $(typeof(p))"
    @assert p.parameter isa IdentifierNode "Expected IdentifierNode value, got $(typeof(p.parameter))"
    @assert p.parameter.namespace !== nothing "Expected a namespace on the value"
    alias = get_value(p.name)
    node = p.parameter.namespace.position.value
    attr = get_value(p.parameter)
    println("  binding: $alias -> $node::$attr")
end

# Verify exact bindings
b1, b2 = lhs_params
@assert get_value(b1.name) == "rho1"
@assert b1.parameter.namespace.position.value == "c1"
@assert get_value(b1.parameter) == "rho"
@assert get_value(b2.name) == "rho2"
@assert b2.parameter.namespace.position.value == "c2"
@assert get_value(b2.parameter) == "rho"

println("PARSER OK: vertex attribute access bindings parsed correctly.")

# --- 2. IR generation test ----------------------------------------------
# Minimal symbol table: Fluid type with a single scalar attribute `rho`.
symbol_tables = Dict()
symbol_tables["Fluid"] = Dict()
symbol_tables["Fluid"][1] = ("Float", "rho", (false, nothing))

rules_table = Dict()
propensity_table = Dict()
propensity_table["parameter_table"] = Dict()
propensity_table["function_table"] = Dict()

ir_code = ir_rules_section!(ast, rules_table, symbol_tables, propensity_table, "Fluids")

println("\n===== Generated Rule IR =====")
println(ir_code)

@assert occursin("rho1", ir_code) "Generated IR should reference binding rho1"
@assert occursin("Fluid", ir_code) "Generated IR should reference node type Fluid"
@assert occursin(".rho", ir_code) "Generated IR should access the rho attribute"

println("\nIR OK: vertex attribute access generated C++ successfully.")
