# PDE rules (MFEM-backed finite elements).
# Included into module IRRuleGeneration (see ../ir_rule_generation.jl).
#
# A PDE `solving` rule is written in the SAME differential form as an ODE `solving` rule,
# and its mesh is driven by the graph exactly like an ODE system is:
#
#   heat := (a : Vertex) -- (b : Vertex)  << (ua : a::u),  (ub : b::u)  >>
#        -> (a : Vertex) -- (b : Vertex)  << (dua : a::u), (dub : b::u) >>
#   solving (dua := D(ua, t), dub := D(ub, t)) {
#       cell  : Element   = segment(a, b)
#       dua   : PDE       = kappa * laplacian(ua)
#       left  : Dirichlet = boundary(x < 2.000001, u_left)
#       right : Dirichlet = boundary(x > 3.999999, u_right)
#   }
#
# * Every match of the LHS pattern is one finite element (`Element` names the element shape and
#   which pattern nodes are its vertices, in the shape's local order). The union of the matches
#   is the mesh, so growing, cutting or moving the graph changes the discretisation: the
#   operator rebuilds its mesh from the matches whenever the solver is (re)built, exactly as
#   the ODE solver rebuilds its equation system.
# * The degrees of freedom ARE the bound attribute of the element nodes (P1: one per vertex),
#   so they sit in the same SUNDIALS state vector as every other solved variable and the PDE
#   composes with ODE solving rules and with ordinary rules on the same nodes.
# * Dirichlet / Neumann conditions select boundary vertices / faces with a position predicate.
#   A boundary face is a face that belongs to exactly one element.
#
# `du/dt = F(u)` with F a linear spatial operator. MFEM works with the weak form, so the
# compiler converts the differential form to a Galerkin system
#
#       M du/dt = L u + b          (M: mass matrix)
#
# term by term (v is the test function, Gamma the boundary):
#
#   c * laplacian(u),  div(c * grad(u))  ->  -K(c) u       K(c)_ij = int c grad(phi_j).grad(phi_i)
#                                              (integration by parts; the boundary term is the
#                                               natural/Neumann boundary condition)
#   dot([bx, by, bz], grad(u)), dx(u) ... ->  +C(b) u      C(b)_ij = int phi_i (b . grad(phi_j))
#   c * u                                 ->  +R(c) u      R(c)_ij = int c phi_j phi_i
#   u-free term f                         ->  +b_i = int f phi_i
#   Neumann boundary(pred, g)             ->  +b_i = int_Gamma g phi_i   (conormal flux g)
#
# The generated operator plugs into the DGGML Solver (graph_driven contract): Setup builds the
# mesh from the rule matches; Mult is the implicit term; PrecSetup/PrecSolve precondition it.

    const PDE_FUNCS = Dict(
        "sin" => "std::sin", "cos" => "std::cos", "tan" => "std::tan",
        "exp" => "std::exp", "log" => "std::log", "sqrt" => "std::sqrt",
        "abs" => "std::abs", "tanh" => "std::tanh", "pow" => "std::pow",
        "min" => "std::min", "max" => "std::max",
    )

    # Coordinates inside generated coefficient lambdas are named px_/py_/pz_.
    const PDE_XYZ_RE = r"\bp[xyz]_(?!\w)"

    # Element shapes (MFEM vertex numbering). `faces` are the facets (lists of local vertices)
    # used to find the boundary; `edges` are the connections the LHS pattern must contain so
    # that a graph match enumerates the vertices in a valid cyclic order.
    const PDE_ELEMENTS = Dict(
        "segment" => (dim = 1, n = 2,
            faces = [[0], [1]],
            edges = [(0, 1)],
            add_el = "mesh->AddSegment(e[0], e[1]);",
            add_bdr = "mesh->AddBdrPoint(f[0]);",
            orient = ""),
        "triangle" => (dim = 2, n = 3,
            faces = [[0, 1], [1, 2], [2, 0]],
            edges = [(0, 1), (1, 2), (2, 0)],
            add_el = "mesh->AddTriangle(e[0], e[1], e[2]);",
            add_bdr = "mesh->AddBdrSegment(f[0], f[1]);",
            orient = "{ const auto& p0 = X[e[0]]; const auto& p1 = X[e[1]]; const auto& p2 = X[e[2]];\n" *
                     "  const double a = (p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]);\n" *
                     "  if (a < 0.0) std::swap(e[1], e[2]); }"),
        "quad" => (dim = 2, n = 4,
            faces = [[0, 1], [1, 2], [2, 3], [3, 0]],
            edges = [(0, 1), (1, 2), (2, 3), (3, 0)],
            add_el = "mesh->AddQuad(e[0], e[1], e[2], e[3]);",
            add_bdr = "mesh->AddBdrSegment(f[0], f[1]);",
            orient = "{ double a = 0.0;\n" *
                     "  for (int j = 0; j < 4; ++j) { const auto& p = X[e[j]]; const auto& q = X[e[(j + 1) % 4]]; a += p[0] * q[1] - q[0] * p[1]; }\n" *
                     "  if (a < 0.0) std::swap(e[1], e[3]); }"),
        "tet" => (dim = 3, n = 4,
            faces = [[1, 2, 3], [0, 3, 2], [0, 1, 3], [0, 2, 1]],
            edges = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)],
            add_el = "mesh->AddTet(e[0], e[1], e[2], e[3]);",
            add_bdr = "mesh->AddBdrTriangle(f[0], f[1], f[2]);",
            orient = "{ const auto& p0 = X[e[0]]; const auto& p1 = X[e[1]]; const auto& p2 = X[e[2]]; const auto& p3 = X[e[3]];\n" *
                     "  const double a[3] = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};\n" *
                     "  const double b[3] = {p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]};\n" *
                     "  const double c[3] = {p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2]};\n" *
                     "  const double d = a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]);\n" *
                     "  if (d < 0.0) std::swap(e[2], e[3]); }"),
        "hex" => (dim = 3, n = 8,
            faces = [[3, 2, 1, 0], [0, 1, 5, 4], [1, 2, 6, 5], [2, 3, 7, 6], [3, 0, 4, 7], [4, 5, 6, 7]],
            edges = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7)],
            add_el = "mesh->AddHex(e[0], e[1], e[2], e[3], e[4], e[5], e[6], e[7]);",
            add_bdr = "mesh->AddBdrQuad(f[0], f[1], f[2], f[3]);",
            orient = "{ const auto& p0 = X[e[0]]; const auto& p1 = X[e[1]]; const auto& p3 = X[e[3]]; const auto& p4 = X[e[4]];\n" *
                     "  const double a[3] = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};\n" *
                     "  const double b[3] = {p3[0] - p0[0], p3[1] - p0[1], p3[2] - p0[2]};\n" *
                     "  const double c[3] = {p4[0] - p0[0], p4[1] - p0[1], p4[2] - p0[2]};\n" *
                     "  const double d = a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]);\n" *
                     "  if (d < 0.0) { std::swap(e[1], e[3]); std::swap(e[5], e[7]); } }"),
    )
    PDE_ELEMENTS["hexahedron"] = PDE_ELEMENTS["hex"]

    is_pde_clause(::Any) = false
    is_pde_clause(sc::SolveClauseNode) = any(e -> e isa PDEEntryNode, sc.clause)

    struct PDECtx
        rule::String
        field::String
        params::Set{String}
        dim::Int
    end

    struct PDETerm
        kind::Symbol            # :diffusion | :convection | :reaction | :source
        scale::String           # C++ scalar expression (sign and u-free factors included)
        vec::Vector{String}     # convection direction components (before scaling)
    end

    _pde_err(ctx::PDECtx, msg) = throw("PDE rule '$(ctx.rule)': $msg")

    _pde_name(node) = string(node.token.position.value)

    function _pde_depends(node, field)
        node isa IdentifierNode && return _pde_name(node) == field
        node isa GroupNode && return _pde_depends(node.expression, field)
        node isa UnaryOpNode && return _pde_depends(node.operand, field)
        node isa BinaryOpNode && return _pde_depends(node.lhs, field) || _pde_depends(node.rhs, field)
        node isa CallNode && return any(a -> _pde_depends(a, field), node.args)
        node isa ArrayLiteralNode && return any(a -> _pde_depends(a, field), node.elements)
        return false
    end

    # ── Scalar expression -> C++ (double) ──────────────────────────────────────
    function _pde_expr(node, ctx::PDECtx)
        if node isa IntegerNode || node isa FloatNode
            s = _pde_name(node)
            return occursin('.', s) ? s : s * ".0"
        elseif node isa IdentifierNode
            n = _pde_name(node)
            n == ctx.field && _pde_err(ctx, "the field '$n' may only appear inside a linear operator term")
            n in ("x", "y", "z") && return "p$(n)_"
            n == "pi" && return "M_PI"
            n in ctx.params && return "static_cast<double>(settings.$n)"
            _pde_err(ctx, "unknown identifier '$n' (expected a parameter, x, y, z or pi)")
        elseif node isa GroupNode
            return "(" * _pde_expr(node.expression, ctx) * ")"
        elseif node isa UnaryOpNode
            op = string(node.expression.position.value)
            op in ("-", "+") || _pde_err(ctx, "unsupported unary operator '$op'")
            return "($(op)$(_pde_expr(node.operand, ctx)))"
        elseif node isa BinaryOpNode
            op = string(node.expression.position.value)
            l = _pde_expr(node.lhs, ctx)
            r = _pde_expr(node.rhs, ctx)
            op in ("+", "-", "*", "/") && return "($l $op $r)"
            op == "^" && return "std::pow($l, $r)"
            _pde_err(ctx, "unsupported operator '$op' in a coefficient expression")
        elseif node isa CallNode
            fname = string(get_value(node.function_node))
            haskey(PDE_FUNCS, fname) ||
                _pde_err(ctx, "unsupported function '$fname' in a coefficient expression " *
                              "(supported: $(join(sort(collect(keys(PDE_FUNCS))), ", ")))")
            args = join([_pde_expr(a, ctx) for a in node.args], ", ")
            return "$(PDE_FUNCS[fname])($args)"
        end
        _pde_err(ctx, "unsupported expression $(typeof(node)) in a coefficient")
    end

    # ── Position predicate -> C++ (bool) ───────────────────────────────────────
    function _pde_pred(node, ctx::PDECtx)
        node = _pde_unwrap(node)
        if node isa BinaryOpNode
            op = string(node.expression.position.value)
            if op in ("&&", "||")
                return "(" * _pde_pred(node.lhs, ctx) * " $op " * _pde_pred(node.rhs, ctx) * ")"
            elseif op in ("<", "<=", ">", ">=", "==", "!=")
                return "(" * _pde_expr(node.lhs, ctx) * " $op " * _pde_expr(node.rhs, ctx) * ")"
            end
        end
        _pde_err(ctx, "a boundary predicate must be a comparison of x, y, z and parameters " *
                      "(<, <=, >, >=, ==, !=), optionally combined with && and ||")
    end

    # ── Differential form -> weak-form terms ───────────────────────────────────
    _pde_unwrap(node) = node isa GroupNode ? _pde_unwrap(node.expression) : node

    function _pde_is_call(node, name)
        node = _pde_unwrap(node)
        return node isa CallNode && string(get_value(node.function_node)) == name
    end

    function _pde_is_grad_u(node, ctx)
        node = _pde_unwrap(node)
        _pde_is_call(node, "grad") || return false
        length(node.args) == 1 || return false
        a = _pde_unwrap(node.args[1])
        return a isa IdentifierNode && _pde_name(a) == ctx.field
    end

    function _pde_operator_call!(terms, node::CallNode, scale, ctx::PDECtx)
        fname = string(get_value(node.function_node))
        args = node.args

        if fname == "laplacian"
            length(args) == 1 && _pde_unwrap(args[1]) isa IdentifierNode &&
                _pde_name(_pde_unwrap(args[1])) == ctx.field ||
                _pde_err(ctx, "laplacian() takes the field, e.g. laplacian($(ctx.field))")
            push!(terms, PDETerm(:diffusion, scale, String[]))

        elseif fname == "div"
            length(args) == 1 || _pde_err(ctx, "div() takes one argument, e.g. div(c * grad($(ctx.field)))")
            a = _pde_unwrap(args[1])
            if _pde_is_grad_u(a, ctx)
                push!(terms, PDETerm(:diffusion, scale, String[]))
            elseif a isa BinaryOpNode && string(a.expression.position.value) == "*" &&
                   (_pde_is_grad_u(a.lhs, ctx) != _pde_is_grad_u(a.rhs, ctx))
                coef = _pde_is_grad_u(a.lhs, ctx) ? a.rhs : a.lhs
                _pde_depends(coef, ctx.field) && _pde_err(ctx, "the diffusion coefficient cannot depend on '$(ctx.field)' (linear PDEs only)")
                push!(terms, PDETerm(:diffusion, "($scale * $(_pde_expr(coef, ctx)))", String[]))
            elseif a isa BinaryOpNode && string(a.expression.position.value) == "/" && _pde_is_grad_u(a.lhs, ctx)
                _pde_depends(a.rhs, ctx.field) && _pde_err(ctx, "the diffusion coefficient cannot depend on '$(ctx.field)' (linear PDEs only)")
                push!(terms, PDETerm(:diffusion, "($scale / $(_pde_expr(a.rhs, ctx)))", String[]))
            else
                _pde_err(ctx, "div() must contain grad($(ctx.field)) optionally multiplied by a coefficient")
            end

        elseif fname == "dot"
            length(args) == 2 || _pde_err(ctx, "dot() takes a vector and grad($(ctx.field))")
            vec_node, grad_node = _pde_is_grad_u(args[2], ctx) ? (args[1], args[2]) : (args[2], args[1])
            _pde_is_grad_u(grad_node, ctx) || _pde_err(ctx, "dot() needs grad($(ctx.field)) as one argument")
            vec_node = _pde_unwrap(vec_node)
            vec_node isa ArrayLiteralNode ||
                _pde_err(ctx, "dot() needs an array literal vector, e.g. dot([bx, by], grad($(ctx.field)))")
            length(vec_node.elements) == ctx.dim ||
                _pde_err(ctx, "the convection vector has $(length(vec_node.elements)) components but the elements are $(ctx.dim)D")
            comps = String[]
            for e in vec_node.elements
                _pde_depends(e, ctx.field) && _pde_err(ctx, "the convection velocity cannot depend on '$(ctx.field)' (linear PDEs only)")
                push!(comps, _pde_expr(e, ctx))
            end
            push!(terms, PDETerm(:convection, scale, comps))

        elseif fname in ("dx", "dy", "dz")
            axis = Dict("dx" => 1, "dy" => 2, "dz" => 3)[fname]
            axis <= ctx.dim || _pde_err(ctx, "$fname($(ctx.field)) used with $(ctx.dim)D elements")
            length(args) == 1 && _pde_unwrap(args[1]) isa IdentifierNode &&
                _pde_name(_pde_unwrap(args[1])) == ctx.field ||
                _pde_err(ctx, "$fname() takes the field, e.g. $fname($(ctx.field))")
            comps = [i == axis ? "1.0" : "0.0" for i in 1:ctx.dim]
            push!(terms, PDETerm(:convection, scale, comps))

        else
            _pde_err(ctx, "operator '$fname' applied to '$(ctx.field)' is not supported " *
                          "(use laplacian, div(c*grad(u)), dot([..], grad(u)), dx/dy/dz, or c*u)")
        end
    end

    function _pde_decompose!(terms, node, scale, ctx::PDECtx)
        node = _pde_unwrap(node)
        f = ctx.field

        if node isa UnaryOpNode
            op = string(node.expression.position.value)
            if op == "-"
                return _pde_decompose!(terms, node.operand, "(-$scale)", ctx)
            elseif op == "+"
                return _pde_decompose!(terms, node.operand, scale, ctx)
            end
            _pde_err(ctx, "unsupported unary operator '$op'")
        elseif node isa BinaryOpNode
            op = string(node.expression.position.value)
            if op == "+"
                _pde_decompose!(terms, node.lhs, scale, ctx)
                return _pde_decompose!(terms, node.rhs, scale, ctx)
            elseif op == "-"
                _pde_decompose!(terms, node.lhs, scale, ctx)
                return _pde_decompose!(terms, node.rhs, "(-$scale)", ctx)
            end
            if !_pde_depends(node, f)
                return push!(terms, PDETerm(:source, "($scale * $(_pde_expr(node, ctx)))", String[]))
            end
            if op == "*"
                dl, dr = _pde_depends(node.lhs, f), _pde_depends(node.rhs, f)
                dl && dr && _pde_err(ctx, "nonlinear term: '$f' appears in both factors of a product (only linear PDEs are supported)")
                return dr ? _pde_decompose!(terms, node.rhs, "($scale * $(_pde_expr(node.lhs, ctx)))", ctx) :
                            _pde_decompose!(terms, node.lhs, "($scale * $(_pde_expr(node.rhs, ctx)))", ctx)
            elseif op == "/"
                _pde_depends(node.rhs, f) && _pde_err(ctx, "nonlinear term: '$f' appears in a denominator")
                return _pde_decompose!(terms, node.lhs, "($scale / $(_pde_expr(node.rhs, ctx)))", ctx)
            end
            _pde_err(ctx, "operator '$op' cannot be applied to an expression containing '$f' (only linear PDEs are supported)")
        elseif !_pde_depends(node, f)
            return push!(terms, PDETerm(:source, "($scale * $(_pde_expr(node, ctx)))", String[]))
        elseif node isa IdentifierNode
            return push!(terms, PDETerm(:reaction, scale, String[]))
        elseif node isa CallNode
            return _pde_operator_call!(terms, node, scale, ctx)
        end
        _pde_err(ctx, "unsupported expression $(typeof(node)) containing '$f'")
    end

    # ── Coefficient declarations ───────────────────────────────────────────────
    function _pde_scalar_coef_decl(nm, expr)
        if occursin(PDE_XYZ_RE, expr)
            return "mfem::FunctionCoefficient $nm([&](const mfem::Vector& X_) { " *
                   "const double px_ = X_.Size() > 0 ? X_(0) : 0.0; " *
                   "const double py_ = X_.Size() > 1 ? X_(1) : 0.0; " *
                   "const double pz_ = X_.Size() > 2 ? X_(2) : 0.0; " *
                   "(void)px_; (void)py_; (void)pz_; return $expr; });"
        end
        return "mfem::ConstantCoefficient $nm($expr);"
    end

    function _pde_vector_coef_decl(nm, comps, dim)
        if any(c -> occursin(PDE_XYZ_RE, c), comps)
            sets = join(["V_($(i - 1)) = $(c);" for (i, c) in enumerate(comps)], " ")
            return "mfem::VectorFunctionCoefficient $nm($dim, [&](const mfem::Vector& X_, mfem::Vector& V_) { " *
                   "const double px_ = X_.Size() > 0 ? X_(0) : 0.0; " *
                   "const double py_ = X_.Size() > 1 ? X_(1) : 0.0; " *
                   "const double pz_ = X_.Size() > 2 ? X_(2) : 0.0; " *
                   "(void)px_; (void)py_; (void)pz_; V_.SetSize($dim); $sets });"
        end
        sets = join(["$(nm)v_($(i - 1)) = $(c);" for (i, c) in enumerate(comps)], " ")
        return "mfem::Vector $(nm)v_($dim); $sets mfem::VectorConstantCoefficient $nm($(nm)v_);"
    end

    # C++ lambda over a point: [&](double px_, double py_, double pz_) -> ret { return expr; }
    function _pde_point_lambda(nm, ret, expr)
        return "auto $nm = [&](double px_, double py_, double pz_) -> $ret { " *
               "(void)px_; (void)py_; (void)pz_; return $expr; };"
    end

    # `boundary(<predicate>, <value>)` -> (predicate C++, value C++)
    function _pde_boundary_entry(entry, ctx::PDECtx)
        call = entry.value
        (call isa CallNode && string(get_value(call.function_node)) == "boundary" && length(call.args) == 2) ||
            _pde_err(ctx, "'$(_pde_name(entry.name))' must be written as boundary(<position predicate>, <value>), " *
                          "e.g. boundary(x < 0.001, 1.0)")
        _pde_depends(call.args[2], ctx.field) &&
            _pde_err(ctx, "boundary values cannot depend on '$(ctx.field)'")
        return _pde_pred(call.args[1], ctx), _pde_expr(call.args[2], ctx)
    end

    # ── The rule ───────────────────────────────────────────────────────────────
    function ir_pde_rule!(ir_builder, rule_node, type_namespace, symbol_tables, propensity_table,
                          visit_rule_node!)
        rule_name = string(get_value(rule_node.name))
        sc = rule_node.modify_clause
        pre = "PDE rule '$rule_name': "

        for e in sc.clause
            e isa PDEEntryNode ||
                throw(pre * "a PDE `solving` block cannot mix PDE entries with ODE/definition entries")
        end

        entries = Dict{String, Vector{PDEEntryNode}}(k => PDEEntryNode[] for k in PDE_ENTRY_KINDS)
        for e in sc.clause
            push!(entries[e.kind], e)
        end
        isempty(entries["Mesh"]) ||
            throw(pre * "`Mesh` entries are gone: the mesh is the set of graph matches of this rule. " *
                        "Name the element shape with `cell : Element = segment(a, b)` and generate the initial " *
                        "graph in simulation.fflow")
        isempty(entries["Initial"]) ||
            throw(pre * "`Initial` entries are gone: the field is a node attribute; " *
                        "set its initial value in the initial-state generator in simulation.fflow")
        length(entries["Element"]) == 1 ||
            throw(pre * "exactly one `Element` entry is required, e.g. `cell : Element = segment(a, b)`")
        length(entries["PDE"]) == 1 ||
            throw(pre * "exactly one `PDE` entry is required, e.g. `dua : PDE = kappa * laplacian(ua)`")

        params = Set{String}(string(k) for k in keys(get(propensity_table, "parameter_table", Dict())))

        # ---- LHS pattern: node ids (1-based by first appearance), types, bound attributes ----
        lhs_names, lhs_types, lhs_order = build_node_pos_to_type(rule_node.lhs)
        _, lhs_assgn = link_param_to_nodes(rule_node.lhs_parameter, lhs_types, symbol_tables, lhs_names)

        # ---- Element ----
        el_call = entries["Element"][1].value
        (el_call isa CallNode && haskey(PDE_ELEMENTS, string(get_value(el_call.function_node)))) ||
            throw(pre * "Element must be one of segment(a, b), triangle(a, b, c), quad(a, b, c, d), " *
                        "tet(a, b, c, d), hex(a, b, c, d, e, f, g, h) with the pattern node names as arguments")
        el_kind = string(get_value(el_call.function_node))
        el = PDE_ELEMENTS[el_kind]
        length(el_call.args) == el.n || throw(pre * "$el_kind() takes $(el.n) pattern nodes")
        el_nodes = String[]
        for a in el_call.args
            a isa IdentifierNode || throw(pre * "$el_kind() arguments must be pattern node names")
            nm = _pde_name(a)
            haskey(lhs_order, nm) || throw(pre * "'$nm' is not a node of the LHS pattern")
            push!(el_nodes, nm)
        end
        length(unique(el_nodes)) == el.n || throw(pre * "$el_kind() needs $(el.n) distinct nodes")
        Set(keys(lhs_order)) == Set(el_nodes) ||
            throw(pre * "every LHS node must be a vertex of the element (LHS nodes: " *
                        "$(join(collect(keys(lhs_order)), ", ")))")

        # The pattern must contain the element's connections (so that matches are valid elements)
        pattern_edges = Set{Tuple{String, String}}()
        for n in rule_node.lhs
            if n isa UndirectedTypeEdgeNode
                a, b = get_value(n.left_vertex.name), get_value(n.right_vertex.name)
                push!(pattern_edges, a < b ? (a, b) : (b, a))
            end
        end
        for (i, j) in el.edges
            a, b = el_nodes[i + 1], el_nodes[j + 1]
            (a < b ? (a, b) : (b, a)) in pattern_edges ||
                throw(pre * "the LHS pattern must connect '$a' -- '$b' ($el_kind vertices $i and $j are adjacent)")
        end

        # ---- Field binding on each element node ----
        length(sc.variables) == el.n ||
            throw(pre * "expected one binding variable per element node ($(el.n)), e.g. " *
                        "`solving (dua := D(ua, t), dub := D(ub, t))`")
        node_field = Dict{String, Tuple{String, String, String}}()   # node => (alias, attr, type)
        pde_bv = nothing
        pde_name = _pde_name(entries["PDE"][1].name)
        for bv in sc.variables
            (bv.name isa IdentifierNode && bv.value[1] isa IdentifierNode) ||
                throw(pre * "binding variables and fields must be plain names (no indexing)")
            alias = _pde_name(bv.value[1])
            haskey(lhs_assgn, alias) ||
                throw(pre * "'$alias' is not bound on the LHS, e.g. `(a : Vertex) << ($alias : a::u) >>`")
            rp = lhs_assgn[alias]
            nodename = first(k for (k, v) in lhs_order if v == rp.index)
            haskey(node_field, nodename) &&
                throw(pre * "node '$nodename' has more than one binding variable")
            rp.cpp_var.type in ("double", "float") ||
                throw(pre * "the field attribute '$(rp.cpp_var.name)' must be a Float")
            rp.cpp_var.index > 1 ||
                throw(pre * "the field attribute '$(rp.cpp_var.name)' cannot be the position attribute")
            node_field[nodename] = (alias, rp.cpp_var.name, rp.type)
            _pde_name(bv.name) == pde_name && (pde_bv = bv)
        end
        pde_bv === nothing &&
            throw(pre * "the PDE entry must be named after one of the binding variables ($(join([_pde_name(b.name) for b in sc.variables], ", ")))")
        field = _pde_name(pde_bv.value[1])
        field_attr = node_field[first(k for (k, v) in node_field if v[1] == field)][2]
        node_type = node_field[el_nodes[1]][3]
        for nm in el_nodes
            (_, attr, ty) = node_field[nm]
            (attr == field_attr && ty == node_type) ||
                throw(pre * "all element nodes must have the same type and the same field attribute " *
                            "('$node_type::$field_attr'), but '$nm' has '$ty::$attr'")
        end

        ctx = PDECtx(rule_name, field, params, el.dim)

        # ---- PDE terms ----
        terms = PDETerm[]
        _pde_decompose!(terms, entries["PDE"][1].value, "1.0", ctx)
        has_operator = any(t -> t.kind in (:diffusion, :convection, :reaction), terms)
        has_operator || throw(pre * "the PDE has no term containing '$field' (nothing to solve)")

        decls = String[]
        ncoef = Ref(0)
        function new_scalar!(expr)
            ncoef[] += 1
            nm = "c$(ncoef[])_"
            push!(decls, _pde_scalar_coef_decl(nm, expr))
            return nm
        end
        function new_vector!(comps)
            ncoef[] += 1
            nm = "c$(ncoef[])_"
            push!(decls, _pde_vector_coef_decl(nm, comps, el.dim))
            return nm
        end

        l_integrators = String[]
        b_integrators = String[]
        for t in terms
            if t.kind == :diffusion
                push!(l_integrators, "L_form.AddDomainIntegrator(new mfem::DiffusionIntegrator($(new_scalar!("(-1.0 * $(t.scale))"))));")
            elseif t.kind == :reaction
                push!(l_integrators, "L_form.AddDomainIntegrator(new mfem::MassIntegrator($(new_scalar!(t.scale))));")
            elseif t.kind == :convection
                comps = ["($(t.scale) * $c)" for c in t.vec]
                push!(l_integrators, "L_form.AddDomainIntegrator(new mfem::ConvectionIntegrator($(new_vector!(comps)), 1.0));")
            else
                push!(b_integrators, "b_form.AddDomainIntegrator(new mfem::DomainLFIntegrator($(new_scalar!(t.scale))));")
            end
        end

        # Neumann: flux g on the boundary faces whose points satisfy the predicate
        for e in entries["Neumann"]
            pred, val = _pde_boundary_entry(e, ctx)
            ncoef[] += 1
            k = ncoef[]
            push!(decls, _pde_point_lambda("nq$(k)_", "bool", pred))
            push!(decls, _pde_point_lambda("nf$(k)_", "double", val))
            push!(decls, "mfem::FunctionCoefficient nc$(k)_([&](const mfem::Vector& X_) { " *
                         "const double px_ = X_.Size() > 0 ? X_(0) : 0.0; " *
                         "const double py_ = X_.Size() > 1 ? X_(1) : 0.0; " *
                         "const double pz_ = X_.Size() > 2 ? X_(2) : 0.0; " *
                         "return nq$(k)_(px_, py_, pz_) ? nf$(k)_(px_, py_, pz_) : 0.0; });")
            push!(b_integrators, "b_form.AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(nc$(k)_));")
        end
        has_b = !isempty(b_integrators)

        # Dirichlet: vertices satisfying the predicate are held at the value
        dirichlet_ir = String[]
        for e in entries["Dirichlet"]
            pred, val = _pde_boundary_entry(e, ctx)
            push!(dirichlet_ir, """
        {
            $(_pde_point_lambda("dq_", "bool", pred))
            $(_pde_point_lambda("df_", "double", val))
            for (int i = 0; i < n; ++i)
                if (dq_(X[i][0], X[i][1], X[i][2])) { ess_set.insert(i); ess_value[i] = df_(X[i][0], X[i][1], X[i][2]); }
        }""")
        end

        op_name = "$(rule_name)_Operator"
        nnode = el.n
        el_lhs = join([string(lhs_order[nm]) for nm in el_nodes], ", ")
        faces_init = join(["{" * join(f, ", ") * "}" for f in el.faces], ", ")
        facen = length(el.faces[1])
        nface = length(el.faces)
        field_ptr(g) = "&std::get<$type_namespace::$node_type>($g.data).$field_attr"

        implicit_ir = has_b ? "" : """
    // ---- implicit integration (homogeneous linear system: Mult is linear in y) ----
    // Newton matrix on the element dofs: I - gamma J = I - gamma M^-1 L  =>  z = (M - gamma L)^-1 M r.
    void PrecSetup(double gamma)
    {
        if (!active || gamma == prec_gamma) return;
        // Install the new preconditioner before the old one/matrix are destroyed:
        // the solver keeps raw pointers to both.
        std::unique_ptr<mfem::SparseMatrix> A_new(mfem::Add(1.0, *M, -gamma, *L_dir));
        auto P_new = std::make_unique<mfem::GSSmoother>(*A_new);
        A_solver->SetPreconditioner(*P_new);
        A_solver->SetOperator(*A_new);
        A = std::move(A_new);
        A_prec = std::move(P_new);
        prec_gamma = gamma;
    }

    void PrecSolve(const double* r, double* z, double gamma)
    {
        if (!active) return;
        PrecSetup(gamma);
        mfem::Vector rv(nv), zv(nv), bb(nv);
        for (int i = 0; i < nv; ++i) rv(i) = r[yidx[i]];
        M->Mult(rv, bb);
        zv = 0.0;
        A_solver->Mult(bb, zv);
        for (int i = 0; i < nv; ++i) z[yidx[i]] = zv(i);
    }
"""
        implicit_members = has_b ? "" : """
    std::unique_ptr<mfem::SparseMatrix> L_dir, A;
    std::unique_ptr<mfem::GSSmoother> A_prec;
    std::unique_ptr<mfem::GMRESSolver> A_solver;
    double prec_gamma = -1.0;
"""
        implicit_teardown = has_b ? "" : """
        A_solver.reset(); A_prec.reset(); A.reset(); L_dir.reset();
        prec_gamma = -1.0;"""
        implicit_setup = has_b ? "" : """
        L_dir = std::make_unique<mfem::SparseMatrix>(*L);
        for (int d : ess) L_dir->EliminateRowCol(d, mfem::Operator::DIAG_ZERO);
        A_solver = std::make_unique<mfem::GMRESSolver>();
        A_solver->SetRelTol(1e-10);
        A_solver->SetMaxIter(500);
        A_solver->SetPrintLevel(0);"""

        b_ir = !has_b ? "" : """
        b.SetSize(nv);
        {
            mfem::LinearForm b_form(fespace.get());
            $(join(b_integrators, "\n            "))
            b_form.Assemble();
            b = b_form;
        }"""
        b_apply = has_b ? "rhs += b;" : ""

        class_ir = """
// ----------------------------------------------------------------------------
// MFEM operator generated from PDE rule '$rule_name':  M du/dt = L u$(has_b ? " + b" : "")
// (weak form of the differential equation, P1 H1 elements). The mesh is built from the
// graph matches of the rule: one $el_kind element per match.
// ----------------------------------------------------------------------------
class $op_name {
public:
    static constexpr bool graph_driven = true;
    int ndofs = 0;                 // no separate FEM segment: the dofs are grammar variables
    int fem_offset = 0;

    explicit $op_name(const Parameters& settings_) : settings(settings_) {}

    // Build the mesh from the matches of this rule: one element per distinct node set.
    // `matches`: rule name -> list of (LHS node id -> graph key); `varmap`: field pointer -> y index.
    template <typename Graph, typename Matches, typename VarMap>
    bool Setup(Graph& graph, const Matches& matches, const VarMap& varmap)
    {
        const auto mit = matches.find("$rule_name");
        if (mit == matches.end() || mit->second.empty()) {
            active = false;
            return false;
        }

        std::map<std::vector<std::size_t>, std::array<std::size_t, NNODE>> distinct;
        for (const auto& m : mit->second) {
            std::array<std::size_t, NNODE> keys;
            for (int j = 0; j < NNODE; ++j) keys[j] = m.at(ELEM_LHS[j]);
            std::vector<std::size_t> sorted(keys.begin(), keys.end());
            std::sort(sorted.begin(), sorted.end());
            distinct.emplace(std::move(sorted), keys);
        }

        std::map<std::size_t, int> vid;
        for (const auto& kv : distinct)
            for (std::size_t k : kv.second) vid.emplace(k, 0);
        int count = 0;
        for (auto& kv : vid) kv.second = count++;

        std::vector<std::size_t> vkey(count);
        std::vector<std::array<double, 3>> X(count);
        yidx.assign(count, 0);
        for (const auto& kv : vid) {
            auto& vertex = graph[kv.first];
            vkey[kv.second] = kv.first;
            X[kv.second] = {vertex.position[0], vertex.position[1], vertex.position[2]};
            double* fp = $(field_ptr("vertex"));
            yidx[kv.second] = static_cast<int>(varmap.at(fp));
        }

        std::vector<std::array<int, NNODE>> elems;
        elems.reserve(distinct.size());
        for (const auto& kv : distinct) {
            std::array<int, NNODE> e;
            for (int j = 0; j < NNODE; ++j) e[j] = vid.at(kv.second[j]);
$(el.orient)
            elems.push_back(e);
        }

        // Same discretisation as before (typical for every SSA window): keep the matrices
        if (!(built && vkey == sig_keys && X == sig_x && elems == sig_elems)) {
            build(X, elems);
            sig_keys = std::move(vkey);
            sig_x = std::move(X);
            sig_elems = std::move(elems);
            built = true;
        }
        active = true;
        return true;
    }

    // ---- interface required by DGGML::Solver ----

    void Mult(double /*t*/, N_Vector y, N_Vector ydot)
    {
        if (!active) return;
        const double* yp = N_VGetArrayPointer(y);
        mfem::Vector uy(nv), rhs(nv), du(nv);
        for (int i = 0; i < nv; ++i) uy(i) = yp[yidx[i]];
        L->Mult(uy, rhs);
        $b_apply
        for (int d : ess) rhs(d) = 0.0;
        du = 0.0;
        M_solver->Mult(rhs, du);

        double* ydp = N_VGetArrayPointer(ydot);
        for (int i = 0; i < nv; ++i) ydp[yidx[i]] += du(i);
    }

    // Held (Dirichlet) vertices start at their prescribed value
    void Load(N_Vector y)
    {
        if (!active) return;
        double* yp = N_VGetArrayPointer(y);
        for (std::size_t k = 0; k < ess.size(); ++k) yp[yidx[ess[k]]] = ess_val[k];
    }

    void Sync(N_Vector) {}   // the field lives in the graph; the Solver writes it back
$implicit_ir
private:
    static constexpr int NNODE = $nnode;
    static constexpr int NFACE = $nface;
    static constexpr int FACEN = $facen;
    static constexpr int ELEM_LHS[NNODE] = {$el_lhs};
    static constexpr int FACES[NFACE][FACEN] = {$faces_init};

    const Parameters& settings;
    bool active = false, built = false;
    int nv = 0;
    std::vector<int> yidx;                    // dof i -> index in the SUNDIALS state
    std::vector<std::size_t> sig_keys;
    std::vector<std::array<double, 3>> sig_x;
    std::vector<std::array<int, NNODE>> sig_elems;

    std::unique_ptr<mfem::Mesh> mesh;
    std::unique_ptr<mfem::H1_FECollection> fec;
    std::unique_ptr<mfem::FiniteElementSpace> fespace;
    std::unique_ptr<mfem::SparseMatrix> M, L;
    std::vector<int> ess;
    std::vector<double> ess_val;
    mfem::Vector b;
    std::unique_ptr<mfem::DSmoother> M_prec;
    std::unique_ptr<mfem::CGSolver> M_solver;
$implicit_members
    void teardown()
    {
        M_solver.reset(); M_prec.reset();$implicit_teardown
        M.reset(); L.reset();
        fespace.reset(); fec.reset(); mesh.reset();
    }

    void build(const std::vector<std::array<double, 3>>& X,
               const std::vector<std::array<int, NNODE>>& elems)
    {
        teardown();
        const int n = static_cast<int>(X.size());

        // Boundary faces: faces that belong to exactly one element
        std::map<std::vector<int>, int> face_count;
        std::map<std::vector<int>, std::vector<int>> face_order;
        for (const auto& e : elems)
            for (int fi = 0; fi < NFACE; ++fi) {
                std::vector<int> fv(FACEN);
                for (int j = 0; j < FACEN; ++j) fv[j] = e[FACES[fi][j]];
                std::vector<int> key(fv);
                std::sort(key.begin(), key.end());
                if (++face_count[key] == 1) face_order[key] = fv;
            }
        int nbdr = 0;
        for (const auto& kv : face_count) if (kv.second == 1) ++nbdr;

        mesh = std::make_unique<mfem::Mesh>($(el.dim), n, static_cast<int>(elems.size()), nbdr, $(el.dim));
        for (int i = 0; i < n; ++i) mesh->AddVertex(X[i][0], X[i][1], X[i][2]);
        for (const auto& e : elems) { $(el.add_el) }
        for (const auto& kv : face_count)
            if (kv.second == 1) { const auto& f = face_order[kv.first]; $(el.add_bdr) }
        mesh->FinalizeMesh(1);

        fec = std::make_unique<mfem::H1_FECollection>(1, $(el.dim));
        fespace = std::make_unique<mfem::FiniteElementSpace>(mesh.get(), fec.get());
        nv = fespace->GetNDofs();
        if (nv != n)
            throw std::runtime_error("$rule_name: expected one P1 dof per mesh vertex");

        $(join(decls, "\n        "))

        {
            mfem::BilinearForm M_form(fespace.get());
            M_form.AddDomainIntegrator(new mfem::MassIntegrator());
            M_form.Assemble();
            M_form.Finalize();
            M = std::make_unique<mfem::SparseMatrix>(M_form.SpMat());
        }
        {
            mfem::BilinearForm L_form(fespace.get());
            $(join(l_integrators, "\n            "))
            L_form.Assemble();
            L_form.Finalize();
            L = std::make_unique<mfem::SparseMatrix>(L_form.SpMat());
        }
$b_ir

        // Dirichlet vertices are held at their value
        std::set<int> ess_set;
        std::map<int, double> ess_value;
$(join(dirichlet_ir, "\n"))
        ess.assign(ess_set.begin(), ess_set.end());
        ess_val.clear();
        for (int d : ess) ess_val.push_back(ess_value[d]);

        // Held dofs do not move: remove them from the mass solve
        for (int d : ess) M->EliminateRowCol(d);
        M_prec = std::make_unique<mfem::DSmoother>(*M);
        M_solver = std::make_unique<mfem::CGSolver>();
        M_solver->SetOperator(*M);
        M_solver->SetPreconditioner(*M_prec);
        M_solver->SetRelTol(1e-12);
        M_solver->SetMaxIter(500);
        M_solver->SetPrintLevel(0);
$implicit_setup
    }
};
"""
        emit(ir_builder, class_ir)

        # ---- the rule: the LHS pattern defines the elements; the dofs are the solved variables ----
        emit(ir_builder, """
void $rule_name(DGGML::Grammar<$type_namespace::graph_type> &gamma,
           $type_namespace::graph_type &system_graph,
           Parameters &settings) {
""")
        lhs_name = "$(rule_name)_lhs"
        emit(ir_builder, "GT $lhs_name;")
        lhs_symbol = Dict()
        for each_node in rule_node.lhs
            visit_rule_node!(each_node, lhs_symbol, ir_builder, lhs_name, type_namespace)
        end
        var_inserts = join(["varset.insert(&std::get<$type_namespace::$node_type>(lhs[m1[$(lhs_order[nm])]].data).$field_attr);"
                            for nm in el_nodes], "\n        ")
        emit(ir_builder, """
    DGGML::SolvingRule<GT> $rule_name("$rule_name", $lhs_name, $lhs_name,
    $nnode,
    [](auto &lhs, auto &m1, auto &varset) {
        $var_inserts
    },
    // The PDE term is the operator's (implicit Mult); nothing to add per match.
    [](auto &lhs, auto &m1, auto y, auto ydot, auto &varmap) {}
    );
    gamma.addRule($rule_name);
};""")
        return nothing
    end
