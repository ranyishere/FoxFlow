"""
    IRTypeInheritance:
        Vertex type inheritance for FoxFlow.

        A type declared as `Shock << extra : Float >> : Fluid` is a subtype of `Fluid`.
        It inherits every attribute of `Fluid` (same attribute indices) and appends its own.

        The C++ graph stores vertices in a `std::variant`, so a `Shock` vertex can never be
        fetched as a `Fluid`. A rule written over a base type is therefore expanded at
        transpile time into one concrete rule per combination of types in the base type's
        family (the type itself plus all descendants), like the hand-written
        `add_center_variants_link` in the EveryLBM example.
"""
module IRTypeInheritance

    import ..IRUtils: get_value
    import ...AstNodes: RuleNode, TypeInstanceNode, TypeClassNode, IdentifierNode,
            UndirectedTypeEdgeNode, DirectedTypeEdgeNode
    import ...Tokens: IdentifierToken, PositionToken

    export reset_type_hierarchy!, register_type!, has_subtypes, type_family,
           expand_rule_inheritance

    # child type name => parent type name
    const TYPE_PARENT = Dict{String, String}()
    # All declared type names in declaration order
    const TYPE_ORDER = String[]

    function reset_type_hierarchy!()
        empty!(TYPE_PARENT)
        empty!(TYPE_ORDER)
    end

    function register_type!(name::String, parent::Union{String, Nothing}=nothing)
        push!(TYPE_ORDER, name)
        if parent !== nothing
            TYPE_PARENT[name] = parent
        end
    end

    function _is_descendant(name::String, ancestor::String)
        cur = get(TYPE_PARENT, name, nothing)
        while cur !== nothing
            cur == ancestor && return true
            cur = get(TYPE_PARENT, cur, nothing)
        end
        return false
    end

    has_subtypes(name::String) = any(t -> _is_descendant(t, name), TYPE_ORDER)

    """
    The type itself followed by all of its descendants (declaration order).
    """
    function type_family(name::String)
        family = String[name]
        for t in TYPE_ORDER
            if _is_descendant(t, name)
                push!(family, t)
            end
        end
        family
    end

    function _vertex_nodes(items)
        out = TypeInstanceNode[]
        for it in items
            if it isa TypeInstanceNode
                push!(out, it)
            elseif it isa UndirectedTypeEdgeNode || it isa DirectedTypeEdgeNode
                push!(out, it.left_vertex)
                push!(out, it.right_vertex)
            end
        end
        out
    end

    function _retype(node::TypeInstanceNode, new_type::String)
        old = node.type.name
        pos = old.token.position
        new_token = IdentifierToken(PositionToken(new_type, pos.line_no, pos.col_no))
        new_name = IdentifierNode(new_token, old.namespace)
        TypeInstanceNode(node.name, node.parameter,
                         TypeClassNode(new_name, node.type.parameter),
                         node.value)
    end

    function _substitute(item, assign, alias_type)
        if item isa TypeInstanceNode
            item.type === nothing && return item
            alias = get_value(item.name)
            # Only substitute nodes that were declared with the same base type as the
            # LHS node of the same alias (this keeps preserved RHS nodes in sync).
            if haskey(assign, alias) && get_value(item.type.name) == alias_type[alias]
                return _retype(item, assign[alias])
            end
            return item
        elseif item isa UndirectedTypeEdgeNode
            return UndirectedTypeEdgeNode(
                _substitute(item.left_vertex, assign, alias_type),
                _substitute(item.right_vertex, assign, alias_type))
        elseif item isa DirectedTypeEdgeNode
            return DirectedTypeEdgeNode(
                _substitute(item.left_vertex, assign, alias_type),
                _substitute(item.right_vertex, assign, alias_type),
                item.direction)
        end
        item
    end

    """
        expand_rule_inheritance(rules) -> Vector{RuleNode}

    Expands every rule that matches a vertex of a type with subtypes into one concrete rule
    per combination of types. Rules that match no such vertex are returned unchanged.

    - LHS vertices with a base type expand over the type's family.
    - An RHS vertex that reuses an LHS alias (declared with the same type) keeps the matched
      concrete type.
    - RHS-only vertices are created with exactly the type written.
    """
    function expand_rule_inheritance(rules)
        expanded = RuleNode[]
        for rule in rules
            append!(expanded, _expand_rule(rule))
        end
        expanded
    end

    function _expand_rule(rule::RuleNode)
        isempty(TYPE_PARENT) && return [rule]

        aliases = String[]
        alias_type = Dict{String, String}()
        for node in _vertex_nodes(rule.lhs)
            node.type === nothing && continue
            alias = get_value(node.name)
            haskey(alias_type, alias) && continue
            type_name = get_value(node.type.name)
            if has_subtypes(type_name)
                push!(aliases, alias)
                alias_type[alias] = type_name
            end
        end

        isempty(aliases) && return [rule]

        families = [type_family(alias_type[a]) for a in aliases]
        base_name = get_value(rule.name)
        result = RuleNode[]

        for combo in Iterators.product(families...)
            assign = Dict{String, String}(zip(aliases, combo))

            lhs = [_substitute(x, assign, alias_type) for x in rule.lhs]
            rhs = [_substitute(x, assign, alias_type) for x in rule.rhs]

            old = rule.name.token.position
            new_name = base_name * "_" * join(combo, "_")
            name_node = IdentifierNode(
                IdentifierToken(PositionToken(new_name, old.line_no, old.col_no)),
                rule.name.namespace)

            push!(result, RuleNode(name_node, lhs, rule.lhs_parameter,
                                   rhs, rule.rhs_parameter, rule.modify_clause))
        end

        println("Rule $base_name expanded over type inheritance into $(length(result)) rules")
        result
    end

end
