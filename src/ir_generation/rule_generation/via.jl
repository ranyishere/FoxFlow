# Sub-grammar (`via`) rule IR.
# Included into module IRRuleGeneration (see ../ir_rule_generation.jl).

    # ── `via` rules ──────────────────────────────────────────────────────────
    #
    #   r := (c : Cell) << (a : c::A) >> -> (c : Cell) << (b : c::A) >>
    #        with (k) via inner(t) {
    #            unpack = (p : P) << (a : p::X) >>
    #            pack   = (p : P) << (b : p::X) >>
    #        }
    #
    # A `via` rule is a stochastic (`with`) rule: it fires with propensity `k` like any other, and
    # its update is a sub-grammar call. Each firing *unpacks* the matched LHS nodes into a fresh graph
    # (the `unpack:` pattern, filled in by the same shared-name copies / `where` assignments as any
    # rule), runs grammar `inner` on that graph for `t` time units (a nested simulation, see
    # `via_run_<grammar>` in model.h), and *packs* the final graph back: the `pack:` pattern is matched
    # in it (first match, nodes of equal type, pattern edges present) and its bindings are written
    # onto the rule's RHS. LHS nodes missing from the RHS are destroyed, new RHS nodes are created
    # (at the position of the first LHS node unless the pack `where` sets it).

    const FFLOW_VIA_RUNTIME = raw"""
#ifndef FFLOW_VIA_RUNTIME_H
#define FFLOW_VIA_RUNTIME_H
#include <algorithm>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <set>
#include <streambuf>
#include <string>
#include <vector>
namespace fflow_via {

// Maps every node of `pattern` to a distinct node of `target` of the same type such that each
// pattern edge is also an edge of `target` (the pattern may be disconnected). Nodes are tried in
// ascending key order and the first complete assignment wins. out: pattern key -> target key.
template <class G, class Map>
bool match(G &pattern, G &target, Map &out) {
    using K = typename G::key_type;
    std::vector<K> pkeys, tkeys;
    for (auto &kv : pattern.getNodeSetRef()) pkeys.push_back(kv.first);
    for (auto &kv : target.getNodeSetRef()) tkeys.push_back(kv.first);
    std::sort(pkeys.begin(), pkeys.end());
    std::sort(tkeys.begin(), tkeys.end());
    if (pkeys.size() > tkeys.size()) return false;
    std::vector<std::size_t> pick(pkeys.size(), 0);
    std::vector<bool> used(tkeys.size(), false);
    std::function<bool(std::size_t)> extend = [&](std::size_t i) -> bool {
        if (i == pkeys.size()) return true;
        auto ptype = pattern.findNode(pkeys[i])->second.getData().type;
        for (std::size_t j = 0; j < tkeys.size(); ++j) {
            if (used[j]) continue;
            if (!(target.findNode(tkeys[j])->second.getData().type == ptype)) continue;
            bool ok = true;
            for (std::size_t k = 0; k < i && ok; ++k) {
                if (pattern.adjacent(pkeys[i], pkeys[k]) || pattern.adjacent(pkeys[k], pkeys[i]))
                    ok = target.adjacent(tkeys[pick[k]], tkeys[j]) || target.adjacent(tkeys[j], tkeys[pick[k]]);
            }
            if (!ok) continue;
            used[j] = true;
            pick[i] = j;
            if (extend(i + 1)) return true;
            used[j] = false;
        }
        return false;
    };
    if (!extend(0)) return false;
    out.clear();
    for (std::size_t i = 0; i < pkeys.size(); ++i) out[pkeys[i]] = tkeys[pick[i]];
    return true;
}

inline void warn_no_match(const char *rule) {
    static std::set<std::string> seen;
    if (seen.insert(rule).second)
        std::cerr << "via rule '" << rule << "': the pack pattern matched nothing in the sub-grammar's "
                     "final graph; the rule's RHS keeps its default values (reported once)" << std::endl;
}

struct NullBuf : std::streambuf {
    int overflow(int c) override { return c; }
};

// Silences std::cout and the listener files for the duration of a nested sub-grammar run.
class QuietScope {
    NullBuf null_;
    std::streambuf *saved_;
    bool had_env_ = false;
    std::string old_env_;
public:
    QuietScope() : saved_(std::cout.rdbuf(&null_)) {
        if (const char *e = std::getenv("FOXFLOW_LISTENERS")) { had_env_ = true; old_env_ = e; }
        setenv("FOXFLOW_LISTENERS", "0", 1);
    }
    ~QuietScope() {
        std::cout.rdbuf(saved_);
        if (had_env_) setenv("FOXFLOW_LISTENERS", old_env_.c_str(), 1);
        else unsetenv("FOXFLOW_LISTENERS");
    }
    QuietScope(const QuietScope &) = delete;
    QuietScope &operator=(const QuietScope &) = delete;
};

} // namespace fflow_via
#endif
"""

    # Forget the where-clause locals left behind by the previous lambda, so the next one declares
    # its own variables (a `via` rule emits several lambdas for one rule).
    function _via_reset_local_state!(propensity_table)
        reserved_keys = ("parameter_table", "function_table", "var_local_table", "random_device")
        for k in collect(keys(propensity_table))
            k in reserved_keys || delete!(propensity_table, k)
        end
        var_local_table = Dict{String, Any}()
        var_local_table["random_device"] = false
        var_local_table["rule_lhs"] = Dict()
        var_local_table["rule_lhs"]["declared"] = []
        var_local_table["propensity"] = Dict()
        var_local_table["propensity"]["declared"] = []
        var_local_table["rule_rhs"] = Dict()
        var_local_table["rule_rhs"]["declared"] = []
        propensity_table["var_local_table"] = var_local_table
    end

    function validate_via_rule(rule_name, via)
        for (what, graph) in (("unpack", via.unpack_graph), ("pack", via.pack_graph))
            for it in graph
                (it isa TypeInstanceNode || it isa UndirectedTypeEdgeNode) ||
                    throw("Via rule '$rule_name': unsupported item in `$what:` pattern ($(typeof(it)))")
            end
        end
    end

    """
    Emits the update lambda of a `via` rule (a C++ `[&](auto &lhs, auto &rhs, auto &m1, auto &m2)`):
    unpack -> `via_run_<grammar>` -> pack. `lhs_*` / `rhs_*` describe the rule's own LHS / RHS.
    """
    function ir_via_update(rule_name, via, visit_rule_node!, type_namespace, symbol_tables,
                           propensity_table, lhs_names, rhs_names,
                           lhs_param_to_node, lhs_assgn_to_node,
                           rhs_param_to_node, rhs_assgn_to_node, rhs_node_locs)
        validate_via_rule(rule_name, via)
        gname = get_value(via.grammar)

        unpack_names, unpack_types, unpack_order = build_node_pos_to_type(via.unpack_graph)
        pack_names, pack_types, pack_order = build_node_pos_to_type(via.pack_graph)

        # Unpack / pack aliases must not collide with the rule's own aliases: the shared-alias
        # carry-over of where clauses would otherwise copy attributes between unrelated nodes.
        rule_aliases = Set(vcat(collect(lhs_names), collect(rhs_names)))
        for (what, order) in (("unpack", unpack_order), ("pack", pack_order))
            clash = intersect(rule_aliases, Set(keys(order)))
            isempty(clash) ||
                throw("Via rule '$rule_name': the `$what:` pattern reuses the node alias(es) " *
                      "$(join(sort(collect(clash)), ", ")) of the rule's LHS/RHS; use different aliases.")
        end

        # unpack: RHS-like side is the sub-grammar graph (m2 = identity map), LHS is the rule's LHS.
        u_param_to_node, u_assgn_to_node = link_param_to_nodes(via.unpack_parameter, unpack_types,
                                                symbol_tables, unpack_names; ordinals=unpack_order)
        u_node_locs = Dict()
        for (ix, nm) in enumerate(unpack_names)
            haskey(u_node_locs, nm) || (u_node_locs[nm] = (unpack_order[nm], unpack_types[ix]))
        end

        # pack: LHS-like side is the sub-grammar graph (m1 = pack match), RHS is the rule's RHS.
        p_param_to_node, p_assgn_to_node = link_param_to_nodes(via.pack_parameter, pack_types,
                                                symbol_tables, pack_names; ordinals=pack_order)

        b = IRBuilder([])
        emit(b, "[&](auto &lhs, auto &rhs, auto &m1, auto &m2) {")

        # 1. the sub-grammar's starting graph
        emit(b, "GT via_sub;")
        u_symbol = Dict()
        for item in via.unpack_graph
            visit_rule_node!(item, u_symbol, b, "via_sub", type_namespace; ordinals=unpack_order)
        end
        n_unpack = length(unpack_order)
        emit(b, "std::decay_t<decltype(m1)> via_um;")
        emit(b, "for (std::size_t i = 1; i <= $n_unpack; ++i) via_um[i] = i;")
        # every unpacked node starts at the position of the rule's principal (first LHS) node
        emit(b, "for (std::size_t i = 1; i <= $n_unpack; ++i) std::copy(std::begin(lhs[m1[1]].position), " *
                "std::end(lhs[m1[1]].position), std::begin(via_sub[i].position));")

        _via_reset_local_state!(propensity_table)
        unpack_ir = IRBuilder([])
        ir_where_clause!(via.unpack_where.clause, unpack_ir,
                         lhs_param_to_node, lhs_assgn_to_node,
                         u_param_to_node, u_assgn_to_node,
                         symbol_tables, type_namespace, propensity_table, u_node_locs)
        emit(b, "auto via_unpack = " * build(unpack_ir) * ";")
        emit(b, "via_unpack(lhs, via_sub, m1, via_um);")

        # 2. run the sub-grammar for `t`
        t_expr = _every_schedule_expr(via.duration, propensity_table, rule_name, "t")
        t_expr === nothing && throw("Via rule '$rule_name': missing duration")
        emit(b, "const double via_origin[3] = {lhs[m1[1]].position[0], lhs[m1[1]].position[1], lhs[m1[1]].position[2]};")
        emit(b, "via_run_$gname(via_sub, settings, static_cast<double>($t_expr), via_origin);")

        # 3. pack the final graph back onto the RHS
        emit(b, "GT via_pat;")
        p_symbol = Dict()
        for item in via.pack_graph
            visit_rule_node!(item, p_symbol, b, "via_pat", type_namespace; ordinals=pack_order)
        end
        emit(b, "std::decay_t<decltype(m1)> via_pm;")
        # new RHS nodes start where the rule's principal LHS node is (the pack `where` may move them)
        for (nm, (ord, _)) in sort(collect(rhs_node_locs); by=kv -> kv[2][1])
            nm in lhs_names && continue
            emit(b, "std::copy(std::begin(via_origin), std::end(via_origin), std::begin(rhs[m2[$ord]].position));")
        end
        emit(b, "if (fflow_via::match(via_pat, via_sub, via_pm)) {")
        _via_reset_local_state!(propensity_table)
        pack_ir = IRBuilder([])
        ir_where_clause!(via.pack_where.clause, pack_ir,
                         p_param_to_node, p_assgn_to_node,
                         rhs_param_to_node, rhs_assgn_to_node,
                         symbol_tables, type_namespace, propensity_table, rhs_node_locs)
        emit(b, "auto via_pack = " * build(pack_ir) * ";")
        emit(b, "via_pack(via_sub, rhs, via_pm, m2);")
        emit(b, "} else { fflow_via::warn_no_match(\"$rule_name\"); }")
        emit(b, "}")

        return gname, t_expr, build(b)
    end
