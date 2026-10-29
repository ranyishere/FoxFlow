"""
    IRRuleGeneration:
        Module for generating Intermediate Representation (IR) for rules in a grammar.
"""

module IRRuleGeneration

    import ..IRUtils: get_value, convert_type_name, convert_type_node, is_list_type_tensor,
            is_tf_array_type, is_torch_tensor_type, tf_array_size, fixed_list_dims
    import ..IRBuildUtils: emit, build, IRBuilder, build_sameline
    import ..IRTypeInheritance: expand_rule_inheritance
    import ...AstNodes: RuleNode, IdentifierNode, TypeInstanceNode,
            UndirectedTypeEdgeNode, BinaryOpNode, ExpressionNode, ParameterNode,
            WithClauseNode, EveryClauseNode, SolveClauseNode, UnaryOpNode, FunctionNode, GroupNode,
            CallNode, LiteralNode, DefinitionNode, ODENode, BindingVariableNode, WhereClauseNode, NamedParameterNode,
            IndexAccessNode, ArrayLiteralNode, SliceNode,
            IntegerNode, FloatNode, StringNode,
            FunctionDefinitionNode, FunctionSignatureNode, FunctionArgNode,
            FunctionBodyExpressionNode, FunctionDefinitionExpressionNode, ReturnNode

    export ir_rules_section!

    using OrderedCollections

    const BUILT_IN_FUNC = [
        "heaviside", "sqrt", "normal_distr",
        "uniform_distr", "cos", "sin", "inverse", "pow", "indicator",
        # Matrix / tensor operations (FixedList<<M,N,Float>> treated as 2-D tensor)
        "zeros_matrix", "ones_matrix", "rand_matrix", "eye_matrix",
        "mat_add", "mat_mul", "mat_dot", "transpose",
        # N-D tensor operations
        "tensordot", "einsum", "permute",
        # Autograd
        "backward", "grad",
        # Single-call automatic differentiation: autodiff(expr, var)
        "autodiff"
    ]

    # Module-level type namespace, set by ir_rules_section! at entry
    _type_namespace = Ref("")

    # Global counter for generating unique array literal temp variable names
    _arr_tmp_counter = Ref(0)

    function ir_rules_section!(
            ast, rules_table, symbol_tables,
            propensity_table, type_namespace;
            stage_index::Union{Int, Nothing}=nothing
        )
        """
        Generate intermediate rules
        for the section.
        """

        rule_section_name = get_value(ast.name)
        global rule_namespace = rule_section_name
        _type_namespace[] = type_namespace

        # Generate a unique header guard per stage
        guard_suffix = stage_index !== nothing ? "_STAGE_$(stage_index)" : ""
        guard_name = "DGGML_RULES$(guard_suffix)_HPP"

        rule_section_header = [
            "#ifndef $guard_name",
            "#define $guard_name",
            "#include \"types.h\"",
            "#include \"parameters.h\"",
            "#include \"functions.h\"",
            "namespace $rule_section_name {",
            "using GT = $type_namespace::graph_type;"
        ]

        rules_ir = IRBuilder([])

        map(
            (hdr) -> begin
                emit(rules_ir, hdr)
            end,
            rule_section_header
        )


        map((rule) -> begin
                build_rules_table(rule,
                                  rules_table,
                                  type_namespace)

            # Clear rule-local binding variables leaked into the global
            # propensity_table by a previous rule's where-clause. ir_where_left_clause!
            # writes LHS binding names to the top level of propensity_table, and those
            # entries are otherwise never cleared. Without this, a binding name reused
            # across rules (e.g. `pos_pos` that is an LHS node in one rule but an
            # RHS-only node in another) resolves to the previous rule's lhs[m1[...]]
            # read and indexes a non-existent match vertex -> segfault at runtime.
            reserved_keys = ("parameter_table", "function_table", "var_local_table", "random_device")
            for k in collect(keys(propensity_table))
                if !(k in reserved_keys)
                    delete!(propensity_table, k)
                end
            end

            # Reset the random_device each time and the local scope
            var_local_table = Dict{String, Any}()

            var_local_table["random_device"] = false
            var_local_table["rule_lhs"] = Dict()
            var_local_table["rule_lhs"]["declared"] = []
            var_local_table["propensity"] = Dict()
            var_local_table["propensity"]["declared"] = []
            var_local_table["rule_rhs"] = Dict()
            var_local_table["rule_rhs"]["declared"] = []

            propensity_table["var_local_table"] = var_local_table

            ir_rule!(rules_ir, rule,
                type_namespace,
                symbol_tables,
                propensity_table)

                end, expand_rule_inheritance(ast.rules_list)
        )

        emit(rules_ir, "}\n#endif")
        build(rules_ir)
    end

    # Rule-generation pieces (textually included into this module).
    include("rule_generation/array_index.jl")
    include("rule_generation/bindings.jl")
    include("rule_generation/expressions.jl")
    include("rule_generation/propensity.jl")
    include("rule_generation/ode.jl")
    include("rule_generation/every.jl")
    include("rule_generation/where_clause.jl")
    include("rule_generation/where_assignment.jl")
    include("rule_generation/rule.jl")

end
