# Clocked (`every`) rule validation and IR helpers.
# Included into module IRRuleGeneration (see ../ir_rule_generation.jl).

    # ── Clocked (`every`) rules ──────────────────────────────────────────────
    # An `every` rule is a data-only rewrite: LHS and RHS must describe the same
    # nodes (same aliases/types, same first-appearance order) and the same edges.

    function _every_edge_set(items)
        edges = Set{Tuple{String, String}}()
        for it in items
            if it isa UndirectedTypeEdgeNode
                a = get_value(it.left_vertex.name)
                b = get_value(it.right_vertex.name)
                push!(edges, a < b ? (a, b) : (b, a))
            end
        end
        return edges
    end

    function validate_every_rule(rule_name, rule_node, lhs_symbol, rhs_symbol)
        if lhs_symbol != rhs_symbol
            throw("Clocked rule '$rule_name': the LHS and RHS must contain the same " *
                  "nodes (same alias and type) listed in the same order. " *
                  "LHS: $(sort(collect(lhs_symbol), by = p -> p[2])), " *
                  "RHS: $(sort(collect(rhs_symbol), by = p -> p[2])).")
        end
        if _every_edge_set(rule_node.lhs) != _every_edge_set(rule_node.rhs)
            throw("Clocked rule '$rule_name': the LHS and RHS must contain the same " *
                  "edges; an `every` rule only updates node data.")
        end
    end

    function _every_schedule_expr(expr, propensity_table, rule_name, what)
        expr === nothing && return nothing
        ctx = IRBuilder([])
        val = IRBuilder([])
        ir_definition!(expr, val, propensity_table, ctx)
        if !isempty(strip(build(ctx)))
            throw("Clocked rule '$rule_name': `$what` must be built from numbers and " *
                  "parameters; it cannot reference rule variables.")
        end
        return build_sameline(val)
    end

    # Additive combine: live += (after - before) for every attribute of the
    # principal node bound on the RHS.
    function ir_every_combine(type_namespace, principal_type, rhs_assgn_to_node)
        bound = Dict{String, CPPVariable}()
        for (_, rp) in rhs_assgn_to_node
            rp.index == 1 || continue
            bound[rp.cpp_var.name] = rp.cpp_var
        end
        lines = String[]
        numeric = ("double", "float", "int", "int64_t", "long", "size_t")
        for var in sort(collect(values(bound)), by = v -> v.index)
            if is_list_type_tensor(var.type)
                push!(lines, "L.$(var.name) = L.$(var.name) + (A.$(var.name) - B.$(var.name));")
                if var.index == 1
                    push!(lines, "for (int d = 0; d < 3; ++d) live.position[d] += after.position[d] - before.position[d];")
                end
            elseif var.type in numeric
                push!(lines, "L.$(var.name) += A.$(var.name) - B.$(var.name);")
                if var.index <= 3
                    push!(lines, "live.position[$(var.index - 1)] += after.position[$(var.index - 1)] - before.position[$(var.index - 1)];")
                end
            end
        end
        body = join(lines, "\n")
        return "[](auto &live, auto &before, auto &after) {\n" *
               "auto &L = std::get<$(type_namespace)::$(principal_type)>(live.data);\n" *
               "auto &B = std::get<$(type_namespace)::$(principal_type)>(before.data);\n" *
               "auto &A = std::get<$(type_namespace)::$(principal_type)>(after.data);\n" *
               body * "\n}"
    end
