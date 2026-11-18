"""
IR Models Section
"""

function ir_load_graph_node(types, type_namespace)
    """
    Load Graph Node IR.
    Given the types file generated from the DSL,
    generate the code to load each type.
    """

    ir = IRBuilder([])

    emit(ir, "std::visit([&](auto &alt) {")
    emit(ir, "using T = std::decay_t<decltype(alt)>;")

    for (ix, type) in enumerate(types)

        if ix == 1
            emit(ir, "if constexpr (std::is_same_v<T, $type_namespace::$type>) {")
        else
            emit(ir, "} else if constexpr (std::is_same_v<T, $type_namespace::$type>) {")
        end
        copy_ir = "$type_namespace::$type copy = alt;"
        add_node_ir = "system_graph.addNode({new_node_key, {copy,
        node_data.position[0], node_data.position[1], node_data.position[2]} });"
        emit(ir, copy_ir)
        emit(ir, add_node_ir)
    end

    emit(ir, "} else {")
    emit(ir, "std::cout << \"  Type: Unknown\\n\";exit(0);")
    emit(ir, "}")
    emit(ir, "}, node_data.data);")

    return build(ir)
end

function ir_load_graph(symbol_table, type_namespace)
    """
    Loads the intermediate representation for the graph
    """

    visit_code = ir_load_graph_node(collect(keys(symbol_table)), type_namespace)

    types = collect(keys(symbol_table))
    # Prepend with namespace
    types = map(t -> "$type_namespace::$t", types)

    # join by comma now
    type_list = join(types, ", ")

    load_graph_ir = "
void load_graph(DGGML::Grammar<$type_namespace::graph_type> &grammar,
               $type_namespace::graph_type &system_graph,
               Parameters &settings, DGGML::KeyGenerator<key_type> &gen,
               std::string filename = \"simulation_state.bin\"
               ) {

    std::ifstream is(filename, std::ios::binary);
    if (!is) {
        std::cerr << \"Error opening file: \" << filename << std::endl;
        return;
    }

    // Create the input archive
    cereal::BinaryInputArchive archive(is);

    size_t num_nodes;
    archive(num_nodes);
    std::cout << \"Number of nodes to load: \" << num_nodes << std::endl;

    // Check to see if key already exists in the graph
    std::unordered_map<key_type, key_type> existing_keys;

    std::unordered_map<key_type, std::unordered_set<key_type>> out_edges_map;
    std::unordered_map<key_type, std::unordered_set<key_type>> in_edges_map;

    for (unsigned int i=0; i < num_nodes; i++) {

        key_type old_key;

        SpatialNode3D<$type_list> node_data;
        std::unordered_set<key_type> out_neighbors;
        std::unordered_set<key_type> in_neighbors;

        archive(old_key, node_data);
        archive(out_neighbors);
        archive(in_neighbors);

        key_type new_node_key;
        new_node_key = gen.get_key();
        existing_keys[old_key] = new_node_key;

        out_edges_map[new_node_key] = out_neighbors;
        in_edges_map[new_node_key] = in_neighbors;

        $visit_code

    }

    // Adding outgoing edges
    for(const auto& [old_key, out_neighbors] : out_edges_map) {
        for (const auto& old_neighbor_key : out_neighbors) {
            key_type new_neighbor_key = existing_keys[old_neighbor_key];
            system_graph.addEdge(old_key, new_neighbor_key);
        }
    }

    // Add incoming edges
    for(const auto& [old_key, in_neighbors] : in_edges_map) {
        for (const auto& old_neighbor_key : in_neighbors) {
            key_type new_neighbor_key = existing_keys[old_neighbor_key];
            system_graph.addEdge(new_neighbor_key, old_key);
        }
    }

};
"

load_graph_ir
end

function ir_save_graph(type_namespace)
    """
    Generate Intermediate Representation
    for the Save Graph
    """

    models_ir = IRBuilder([])

    tmp_ir = "void save_graph(DGGML::Grammar<$type_namespace::graph_type> &grammar,
            $type_namespace::graph_type &system_graph,
            Parameters &settings, std::string filename = \"simulation_state.bin\"
        ) {

        std::ofstream os(filename, std::ios::binary);
        cereal::BinaryOutputArchive archive(os);

        size_t num_nodes = system_graph.numNodes();
        archive(num_nodes);

        auto check = system_graph.getNodeSetRef();

        int count_nodes = 0;
        // Loop through system graph
        for (auto i=system_graph.node_list_begin(); i != system_graph.node_list_end(); ++i) {
            count_nodes++;
            auto& key =  i->first;
            auto& node_data =  i->second.getData();

            archive(key, node_data);
            archive(system_graph.out_neighbors(i->second));
            archive(system_graph.in_neighbors(i->second));

        }

    };
    "
    emit(models_ir, tmp_ir)

    build(models_ir)
end

function ir_grid_file_Writer(ast, type_namespace)
    """
    Generate Intermediate Representation
    for the Grid File Writer
    """

    models_ir = IRBuilder([])
end


function ir_vtk_file_writer(ast, type_namespace, symbol_tables)
    """
    Generate Intermediate Representation
    for the VTK File Writer.

    Should return code that takes all attributes
    for each vertex and collect them.
    """

    models_ir = IRBuilder([])
end

# ---- Observable expression → C++ string ----
# Built-in function names that map to torch tensor methods returning a scalar
const OBS_BUILTIN_METHODS = Dict(
    "mean" => "mean()",
    "norm" => "norm()",
    "min"  => "min()",
    "max"  => "max()",
    "sum"  => "sum()",
)

"""
Build a Dict{alias => cpp_expr} from an ObservableDestructuredParamNode.
`alias_map` maps the user's field alias name → C++ access expression.
Fields are matched by position to the ordered attributes in symbol_tables.
e.g. `p1 : Layer << position: FixedList, cur_weights: FixedList, layer_id: Integer >>`
  → { "position" => "alt.Position", "cur_weights" => "alt.Weights", "layer_id" => "alt.LayerID" }
The param_var_name is what the user calls the node in the body (e.g. "p1"), which maps to "alt"
when we are in the top-level observable, or to the lambda arg name in a local function.
"""
function _build_alias_map(param, symbol_tables, cpp_var)
    type_name = get_value(param.type_name)
    attrs = get(symbol_tables, type_name, Dict())
    # attrs is ordered dict: index → (cpp_type, field_name, attr_info)
    ordered_fields = sort(collect(attrs), by = x -> x[1])

    alias_map = Dict{String, String}()
    param_name = get_value(param.name)
    # The param name itself maps to the C++ variable
    alias_map[param_name] = cpp_var

    for (i, binding) in enumerate(param.fields)
        alias = get_value(binding.name)
        if i <= length(ordered_fields)
            _, attr_info = ordered_fields[i]
            field_name = attr_info[2]
            alias_map[alias] = "$cpp_var.$field_name"
        end
    end
    return alias_map
end

"""
Recursively emit a C++ expression from an observable body AST node.
alias_map maps FoxFlow identifier names → C++ expressions.
"""
function _emit_obs_expr(node, alias_map)
    if node isa CallNode
        fname = get_value(node.function_node)
        args_cpp = [_emit_obs_expr(a, alias_map) for a in node.args]
        if haskey(OBS_BUILTIN_METHODS, fname)
            # e.g. mean(cur_weights) → cur_weights.mean().template item<double>()
            return "$(args_cpp[1]).$(OBS_BUILTIN_METHODS[fname]).template item<double>()"
        else
            # local function call: mean_layer(p1) → mean_layer(alt)
            return "$fname($(join(args_cpp, ", ")))"
        end
    elseif node isa IdentifierNode
        name = get_value(node)
        return get(alias_map, name, name)
    elseif node isa BinaryOpNode
        lhs = _emit_obs_expr(node.lhs, alias_map)
        rhs = _emit_obs_expr(node.rhs, alias_map)
        op  = node.expression.position.value
        return "($lhs $op $rhs)"
    elseif node isa GroupNode
        return "(" * _emit_obs_expr(node.expression, alias_map) * ")"
    elseif node isa FloatNode || node isa IntegerNode
        return get_value(node)
    else
        return "/* unhandled $(typeof(node)) */"
    end
end

function _emit_obs_local_fn(fn_node, type_namespace, symbol_tables)
    """Emit a C++ auto lambda for an ObservableLocalFunctionNode.
    The lambda takes the node by reference and binds field aliases as local refs."""
    fn_name      = get_value(fn_node.name)
    param        = fn_node.param
    ctx_type     = get_value(param.type_name)
    ctx_var      = get_value(param.name)
    ret_type_str = get_value(fn_node.return_type)
    cpp_ret      = ret_type_str == "Float" ? "double" :
                   ret_type_str == "Integer" ? "int64_t" : "double"
    cpp_ctx_type = "$type_namespace::$ctx_type"

    # Build alias map: alias → ctx_var.FieldName
    alias_map = _build_alias_map(param, symbol_tables, ctx_var)

    lines = String[]
    push!(lines, "auto $fn_name = [&]($cpp_ctx_type& $ctx_var) -> $cpp_ret {")

    # Emit local refs for each field alias
    for (alias, cpp_expr) in alias_map
        if alias != ctx_var  # skip the node itself
            push!(lines, "    auto& $alias = $cpp_expr;")
        end
    end

    # Pre-return body expressions
    if fn_node.body !== nothing
        for expr_node in fn_node.body.expressions
            expr_name = get_value(expr_node.name)
            expr_cpp  = _emit_obs_expr(expr_node.value, alias_map)
            push!(lines, "    auto $expr_name = $expr_cpp;")
        end
    end

    ret_cpp = _emit_obs_expr(fn_node.body_return.value, alias_map)
    push!(lines, "    return $ret_cpp;")
    push!(lines, "};")
    return join(lines, "\n")
end

function ir_get_type_attr(ast, type_namespace, symbol_tables, observable_section=nothing)
    """
    IR Get Type Attribute

    Generates a C++ function that collects all node attributes into
    a vector of (name, values) pairs for VTK output.

    Scalar (Float/Integer) attributes are always emitted as columns.

    If an ObservableSectionNode is provided, each ObservableDefinitionNode
    generates one additional column (named after the observable) computed
    by evaluating the user's expression over each matching node.
    """

    # ---- Build type_to_attrs: list of (attr_name, is_tensor) per type ----
    type_to_attrs = Dict{String, Vector{Tuple{String, Bool}}}()
    for (key, value) in symbol_tables
        attrs = Tuple{String, Bool}[]
        for (_, attr) in value
            attr_type = attr[1]
            is_tensor = (is_list_type_tensor(attr_type) && attr[3][1] == true)
            push!(attrs, (attr[2], is_tensor))
        end
        type_to_attrs[String(key)] = attrs
    end

    # ---- Fill loop: scalar attrs + observable-computed columns ----
    function ir_loop_type_attr(type_to_attrs, observable_section, type_namespace)
        ir_loop = IRBuilder([])
        emit(ir_loop, "std::size_t row = 0;")
        emit(ir_loop, "for (auto it = system_graph.node_list_begin(); it != system_graph.node_list_end(); ++it, ++row) {")
        emit(ir_loop, "auto &n = it->second.getData();")
        emit(ir_loop, "std::visit([&](auto &alt) {")
        emit(ir_loop, "using T = std::decay_t<decltype(alt)>;")

        for (key, attrs) in type_to_attrs
            emit(ir_loop, "if constexpr (std::is_same_v<T, $type_namespace::$key>) {")
            # Scalar attrs
            for (attr_name, is_tensor) in attrs
                if !is_tensor
                    emit(ir_loop, "set(\"$key.$attr_name\", row, static_cast<double>(alt.$attr_name));")
                end
            end
            # Observable-defined columns targeting this type
            if observable_section !== nothing
                for def in observable_section.definitions
                    if get_value(def.param.type_name) != key; continue; end
                    obs_name  = get_value(def.name)
                    alias_map = _build_alias_map(def.param, symbol_tables, "alt")
                    emit(ir_loop, "{")
                    for fn in def.local_fns
                        emit(ir_loop, _emit_obs_local_fn(fn, type_namespace, symbol_tables))
                    end
                    ret_cpp = _emit_obs_expr(def.body_return.value, alias_map)
                    emit(ir_loop, "set(\"$obs_name\", row, static_cast<double>($ret_cpp));")
                    emit(ir_loop, "}")
                end
            end
            emit(ir_loop, "}")
        end

        emit(ir_loop, "}, n.data);")
        emit(ir_loop, "}")
        return build(ir_loop)
    end

    ir_get_type_attr_fn = IRBuilder([])
    emit(ir_get_type_attr_fn, "template <typename GraphType>")
    emit(ir_get_type_attr_fn, "auto get_type_attributes(GraphType &system_graph) {")
    emit(ir_get_type_attr_fn, "std::vector<std::pair<std::string, std::vector<double>>> type_attributes;")
    emit(ir_get_type_attr_fn, "double nan_value = std::numeric_limits<double>::quiet_NaN();")
    emit(ir_get_type_attr_fn, "std::vector<std::string> col_names;")
    emit(ir_get_type_attr_fn, "std::size_t N = system_graph.numNodes();")
    emit(ir_get_type_attr_fn, "// Register column names")

    # Scalar attrs (always emitted)
    for (key, attrs) in type_to_attrs
        for (attr_name, is_tensor) in attrs
            if !is_tensor
                emit(ir_get_type_attr_fn, "col_names.push_back(\"$key.$attr_name\");")
            end
        end
    end

    # Observable-defined columns
    if observable_section !== nothing
        for def in observable_section.definitions
            obs_name = get_value(def.name)
            emit(ir_get_type_attr_fn, "col_names.push_back(\"$obs_name\");")
        end
    end

    emit(ir_get_type_attr_fn, "std::vector<std::vector<double>> columns(col_names.size());")
    emit(ir_get_type_attr_fn, "for (auto &col : columns) { col.resize(N, nan_value); }")
    emit(ir_get_type_attr_fn, "std::unordered_map<std::string, std::size_t> col_idx;")
    emit(ir_get_type_attr_fn, "for (std::size_t i = 0; i < col_names.size(); ++i) { col_idx[col_names[i]] = i; }")
    emit(ir_get_type_attr_fn, """
        auto set = [&](const std::string &col_name, std::size_t row, double value) {
            columns[col_idx.at(col_name)][row] = value;
        };
    """)

    emit(ir_get_type_attr_fn, ir_loop_type_attr(type_to_attrs, observable_section, type_namespace))

    emit(ir_get_type_attr_fn, "for (std::size_t i = 0; i < col_names.size(); ++i) {
        type_attributes.emplace_back(col_names[i], std::move(columns[i]));
    }")
    emit(ir_get_type_attr_fn, "return type_attributes;}")
    return build(ir_get_type_attr_fn)
end


"""
Generate the C++ members that evaluate the global observables (`name : Observable :=
reduction pattern return (expr)`) on the model's graph; they are only called by
measure(...). Each observable was synthesized into a hidden rule
`obs_<name>` (see parse_global_observable!); its propensity lambda is the per-match value.
"""
function ir_global_observables(global_observables, observable_section)
    isempty(global_observables) && return ""
    ns = get_value(observable_section.name)

    blocks = String[]
    for obs in global_observables
        name = get_value(obs.name)
        rule = get_value(obs.rule.name)
        accumulate, final = if obs.reduction == "sum"
            ("acc += v;", "acc")
        elseif obs.reduction == "mean"
            ("acc += v;", "(n > 0 ? acc / static_cast<double>(n) : nan_value)")
        elseif obs.reduction == "count"
            ("if (v != 0.0) acc += 1.0;", "acc")
        elseif obs.reduction == "min"
            ("acc = (n == 0) ? v : std::min(acc, v);", "(n > 0 ? acc : nan_value)")
        else  # max
            ("acc = (n == 0) ? v : std::max(acc, v);", "(n > 0 ? acc : nan_value)")
        end
        if obs.within !== nothing
            rc = obs.within
            push!(blocks, """
        if (wanted.empty() || std::find(wanted.begin(), wanted.end(), "$name") != wanted.end())
        {
            // $name := $(obs.reduction) over the ordered pairs of nodes within $rc of each other
            // (pattern nodes 1 and 2 of $rule), found with a uniform grid.
            graph_grammar_t obs_gamma;
            $ns::$rule(obs_gamma, this->system_graph, settings);
            auto &r = obs_gamma.stochastic_rules.at("$rule");
            using match_t = YAGL::Mtype<std::remove_reference_t<decltype(this->system_graph)>>;
            using key_t = typename match_t::key_type;
            const auto type_a = r.lhs_graph.findNode(1)->second.getData().type;
            const auto type_b = r.lhs_graph.findNode(2)->second.getData().type;
            const double rc = $rc;
            const double rc2 = rc * rc;
            auto cell_of = [&](double c) { return static_cast<long long>(std::floor(c / rc)); };
            auto cell_key = [](long long ix, long long iy, long long iz) {
                return ((ix + 100000) * 200003LL + (iy + 100000)) * 200003LL + (iz + 100000);
            };
            struct Item { key_t key; double x, y, z; };
            std::vector<Item> items_a, items_b;
            for (auto it = this->system_graph.node_list_begin(); it != this->system_graph.node_list_end(); ++it) {
                auto &nd = it->second.getData();
                Item item{it->first, static_cast<double>(nd.position[0]),
                          static_cast<double>(nd.position[1]), static_cast<double>(nd.position[2])};
                if (nd.type == type_a) items_a.push_back(item);
                if (nd.type == type_b) items_b.push_back(item);
            }
            std::unordered_map<long long, std::vector<std::size_t>> grid;
            for (std::size_t k = 0; k < items_b.size(); ++k)
                grid[cell_key(cell_of(items_b[k].x), cell_of(items_b[k].y), cell_of(items_b[k].z))].push_back(k);
            double acc = 0.0;
            std::size_t n = 0;
            for (auto &a : items_a) {
                const long long ix = cell_of(a.x), iy = cell_of(a.y), iz = cell_of(a.z);
                for (long long dx = -1; dx <= 1; ++dx)
                for (long long dy = -1; dy <= 1; ++dy)
                for (long long dz = -1; dz <= 1; ++dz) {
                    auto cell = grid.find(cell_key(ix + dx, iy + dy, iz + dz));
                    if (cell == grid.end()) continue;
                    for (auto k : cell->second) {
                        auto &b = items_b[k];
                        if (b.key == a.key) continue;
                        const double ddx = a.x - b.x, ddy = a.y - b.y, ddz = a.z - b.z;
                        if (ddx * ddx + ddy * ddy + ddz * ddz > rc2) continue;
                        match_t m;
                        m[1] = a.key;
                        m[2] = b.key;
                        double v = r.propensity(this->system_graph, m);
                        $accumulate
                        ++n;
                    }
                }
            }
            result.emplace_back("$name", $final);
        }
        """)
            continue
        end
        push!(blocks, """
        if (wanted.empty() || std::find(wanted.begin(), wanted.end(), "$name") != wanted.end())
        {
            // $name := $(obs.reduction) over the matches of the pattern of $rule
            graph_grammar_t obs_gamma;
            $ns::$rule(obs_gamma, this->system_graph, settings);
            auto &r = obs_gamma.stochastic_rules.at("$rule");
            auto matches = YAGL::subgraph_isomorphism2(r.lhs_graph, this->system_graph);
            double acc = 0.0;
            std::size_t n = 0;
            for (auto &m : matches) {
                double v = r.propensity(this->system_graph, m);
                $accumulate
                ++n;
            }
            result.emplace_back("$name", $final);
        }
        """)
    end

    return """
    std::vector<std::pair<std::string, double>> compute_observables(
            const std::vector<std::string> &wanted = {}) {
        std::vector<std::pair<std::string, double>> result;
        const double nan_value = std::numeric_limits<double>::quiet_NaN();
        $(join(blocks, "\n"))
        return result;
    }
    """
end

function ir_serialize_torch_tensor()
    ir = "#include <cereal/archives/binary.hpp>
    #include <cereal/types/vector.hpp>
    #include <torch/torch.h>

    namespace cereal {
        // SAVE function
        template<class Archive>
        void save(Archive& ar, const torch::Tensor& tensor) {
            // 1. Save metadata (dimensions and type)
            std::vector<int64_t> size = tensor.sizes().vec();
            int64_t type = static_cast<int64_t>(tensor.scalar_type());
            ar(size, type);

            // 2. Save raw data
            auto contiguous_tensor = tensor.contiguous();
            size_t bytes = contiguous_tensor.nbytes();
            ar(binary_data(contiguous_tensor.data_ptr(), bytes));
        }

        // LOAD function
        template<class Archive>
        void load(Archive& ar, torch::Tensor& tensor) {
            // 1. Load metadata
            std::vector<int64_t> size;
            int64_t type_int;
            ar(size, type_int);

            // 2. Prepare tensor and load raw data
            tensor = torch::empty(size, torch::TensorOptions().dtype(static_cast<torch::ScalarType>(type_int)));
            ar(binary_data(tensor.data_ptr(), tensor.nbytes()));
        }
    }
"
end

function ir_serialize_spatialnode()
    ir = "namespace cereal {
        // This function tells cereal how to handle SpatialNode3D without
        // modifying the DGGML source code.
        template <class Archive, typename... Ts>
        void serialize(Archive& ar, SpatialNode3D<Ts...>& node) {
            // Map the internal members of the class to the archive
            ar(node.position);
            ar(node.data);
        }
    }
    "
end

function ir_models_section(ast, type_namespace, symbol_tables)
    # Single-stage backwards-compatible wrapper
    ir_models_section_multistage(
        [ast], type_namespace, symbol_tables,
        [Main.rules_table], ["rules_0.h"]
    )
end


function ir_models_section_multistage(all_stages, type_namespace, symbol_tables,
                                       stage_rules_tables, rules_includes;
                                       observable_section=nothing,
                                       with_global_observables=false)
    """
    Generate model.h with one Model_N class per simulation stage.
    Each Model_N only registers the rules for its stage.
    Shared infrastructure (serializers, load_graph, save_graph, etc.) is emitted once.
    """

    println("ir_models_section_multistage: Generating Models Section IR for $(length(all_stages)) stage(s)")

    # Use first stage's type info for the shared helpers
    first_stage = all_stages[1]

    ir_load_graph_ir = ir_load_graph(symbol_tables, type_namespace)
    ir_get_type_fn = ir_get_type_attr(first_stage, type_namespace, symbol_tables, observable_section)
    # Global observables are only evaluated by measure(...) (compute_observables), never during a run.
    global_observables = (observable_section === nothing || !with_global_observables) ? [] :
                         observable_section.global_observables
    ir_global_obs_fn = ir_global_observables(global_observables, observable_section)

    serializer_spatialnode_ir = ir_serialize_spatialnode()
    serializer_torch_tensor_ir = ir_serialize_torch_tensor()

    model_section_name = type_namespace

    # ===== Common header (includes, namespace, shared functions) =====
    model_header = [
        "#ifndef DGGML_MODELS_HPP\n",
        "#define DGGML_MODELS_HPP\n",
        "#include <fstream>\n",
        "#include <filesystem>\n",
        "#include <algorithm>\n",
        "#include <iomanip>\n",
        "#include <limits>\n",
        "#include \"DGGML.h\"\n",
    ]

    # Include all per-stage rule headers
    for rfile in rules_includes
        push!(model_header, "#include \"$rfile\"\n")
    end

    push!(model_header, "#include \"simdjson.h\"\n")
    if any(rt -> any(r -> r["activated"] == true && is_pde_clause(r["modify_clause"]), values(rt)), stage_rules_tables)
        push!(model_header, "#include \"mfem.hpp\"\n")
    end
    push!(model_header, "#include \"ExpandedComplex2D.hpp\"\n")
    push!(model_header, "#include \"ExpandedComplex3D.hpp\"\n")
    push!(model_header, "#include \"YAGL_Algorithms.hpp\"\n")
    push!(model_header, "#include \"parameters.h\"\n")
    push!(model_header, "#include <cereal/archives/binary.hpp>\n")
    push!(model_header, "#include <cereal/types/variant.hpp>\n")
    push!(model_header, "#include <cereal/types/vector.hpp>\n")
    push!(model_header, "#include <cereal/types/unordered_set.hpp>\n")
    push!(model_header, serializer_spatialnode_ir)
    push!(model_header, serializer_torch_tensor_ir)
    push!(model_header, "namespace $model_section_name {\n")
    push!(model_header, "using graph_grammar_t = DGGML::Grammar<$type_namespace::graph_type>;\n")

    # ===== Generate one Model_N class per stage =====
    model_classes = []

    for (stage_ix, stage_info) in enumerate(all_stages)
        stage_idx = stage_ix - 1  # 0-based

        stage_rt = stage_rules_tables[stage_ix]

        initial_state_file = stage_info["initial_state"]

        generator_spec = initial_state_file isa AbstractString ? nothing : initial_state_file
        load_initial_state = nothing
        if generator_spec !== nothing || initial_state_file == "default"
            load_initial_state = "false"
            initial_state_file = "\"simulation_state.bin\""
        else
            load_initial_state = "true"
            initial_state_file = "\"$initial_state_file\""
        end

        model_class_name = "Model_$stage_idx"

        # Create the add_default_type template function (only emitted in the class that needs it)
        model_create_initial_type = 
"template<typename GraphType, typename CplexType, typename ParamType, typename GenType>
    void add_default_type(
                    GraphType &graph,
                    CplexType &cplex,
                    ParamType &settings, GenType &gen
                    ) {

            using node_type = typename GraphType::node_type;
            const double min_x = 0.0;
            const double min_y = 0.0;
            const double min_z = 0.0;

            const double max_x = settings.CELL_NX-1;
            const double max_y = settings.CELL_NY-1;

            const std::size_t max_nx = 16;
            const std::size_t max_ny = 16;
            const std::size_t max_nz = 16;

            DGGML::CartesianGrid3D &reaction_grid = cplex.reaction_grid;
            std::size_t max_particles = reaction_grid.totalNumCells();

            std::vector<std::size_t> selected_cells;
            for (int i = max_particles/2; i < max_particles; i++) {
                selected_cells.push_back(i);
            }

            auto centerCardinal = reaction_grid.cardinalCellIndex(
                        reaction_grid._nx/2, reaction_grid._ny/2, reaction_grid._nz/2
                    );

            double center_x, center_y, center_z;
            reaction_grid.cardinalCellToPoint(center_x, center_y, center_z, centerCardinal);
            $type_namespace::graph_type tg;
            $type_namespace::StartType tmp = $type_namespace::StartType{};

            tmp.start_location[0] = center_x;
            tmp.start_location[1] = center_y;
            tmp.start_location[2] = center_z;

            node_type oof_node = {gen.get_key(),
                                {tmp, 
                                 center_x, center_y, center_z}
                        };

            graph.addNode(oof_node);
            };"

        # Build rule registration calls for this stage only
        call_rule_ir = map(
            (rule_name) -> begin
                if stage_rt[rule_name]["activated"] == true
                    tmp_ir = "\t$(stage_rt[rule_name]["namespace"])::$rule_name(gamma, this->system_graph, settings);\n"
                    return tmp_ir
                end
            end,
            collect(keys(stage_rt))
        )

        initial_rules_ir = join(call_rule_ir)

        # PDE rules: each one owns an MFEM operator (created when its rule is registered).
        pde_rule_names = String[]
        for rule_name in sort(collect(keys(stage_rt)))
            if stage_rt[rule_name]["activated"] == true && is_pde_clause(stage_rt[rule_name]["modify_clause"])
                push!(pde_rule_names, rule_name)
            end
        end
        has_non_pde_rules = any(
            r -> stage_rt[r]["activated"] == true && !is_pde_clause(stage_rt[r]["modify_clause"]),
            keys(stage_rt))
        pde_members_ir = ""
        if !isempty(pde_rule_names)
            op_types = ["$(stage_rt[r]["namespace"])::$(r)_Operator" for r in pde_rule_names]
            fem_type = length(op_types) == 1 ? op_types[1] :
                       "DGGML::FemOperatorGroup<$(join(op_types, ", "))>"
            make_ops = length(op_types) == 1 ? "std::make_shared<$(op_types[1])>(settings)" :
                "std::make_shared<fem_operator_type>($(join(["std::make_shared<$t>(settings)" for t in op_types], ", ")))"
            # Graph-driven operators: every geocell solver gets its own operator instance, which
            # builds its mesh from that geocell's rule matches (see FemOperatorGroup / Solver).
            pde_members_ir = "using fem_operator_type = $fem_type;\n" *
                "std::map<std::size_t, std::shared_ptr<fem_operator_type>> fem_ops;\n" *
                "std::shared_ptr<std::mutex> fem_ops_mutex = std::make_shared<std::mutex>();\n" *
                "fem_operator_type* fem_operator(std::size_t k, bool has_rules) {\n" *
                "    if (!has_rules) return nullptr;\n" *
                "    std::lock_guard<std::mutex> lock(*fem_ops_mutex);\n" *
                "    auto& slot = fem_ops[k];\n" *
                "    if (!slot) slot = $make_ops;\n" *
                "    return slot.get();\n}\n"
            # A pure-PDE stage with a generator needs no StartType: its graph is the generated
            # lattice. When other (ODE / propensity) rules are active they are seeded from the
            # StartType node as usual, so keep it alongside the lattice.
            if generator_spec !== nothing && !has_non_pde_rules
                model_create_initial_type = ""
            end
        end
        default_type_ir = "this->add_default_type(this->system_graph,
                        geoplex2D,
                        settings,
                        this->gen);"
        initial_state_ir = default_type_ir
        if generator_spec !== nothing
            generator_ir = ir_mesh_generator(generator_spec, type_namespace, symbol_tables,
                                             keys(propensity_table["parameter_table"]))
            initial_state_ir = has_non_pde_rules ? default_type_ir * "\n" * generator_ir : generator_ir
        end

        save_graph_ir = ir_save_graph(type_namespace)

        # Checkpoint function

        # Name of the "latest" save-state bin file. Provided as the first
        # argument to RunSimulation, e.g. RunSimulation("my_state.bin", ...).
        # Falls back to the default when not specified.
        stage_save_name = get(stage_info, "save_name", nothing)
        if stage_save_name === nothing || stage_save_name == ""
            stage_save_name = "simulation_state_latest.bin"
        end

        # `measure(sim, ...)` evaluates observables on this simulation's saved checkpoints, so a
        # measured simulation keeps them in my_results/<sim>/state_<step>.bin (cleared at step 0)
        # instead of the shared my_results/simulation_state_<step>.bin.
        stage_measured = get(stage_info, "measured", false)
        stage_sim_name = string(get(stage_info, "name", "stage$(stage_ix - 1)"))
        snapshot_prefix = stage_measured ? "/$(stage_sim_name)/state_" : "/simulation_state_"
        snapshot_dir_ir = stage_measured ?
            "std::filesystem::create_directories(results_dir_name + \"/$(stage_sim_name)\");
                        for (const auto& snap : std::filesystem::directory_iterator(results_dir_name + \"/$(stage_sim_name)\"))
                        {
                            if (snap.is_regular_file() && snap.path().extension() == \".bin\")
                                std::filesystem::remove(snap.path());
                        }" : ""

        # Clear stale transient output (VTK/CSV) from previous runs at step 0,
        # but preserve saved-state (.bin) files. A blanket remove_all would
        # delete the state file that a later stage — or a standalone rerun of a
        # later stage — needs to load (e.g. "my_results/boundary_latest.bin"),
        # which is what caused an empty simulation when running a stage on its
        # own. Selective cleanup keeps output tidy while protecting state.
        dir_setup_ir =
            "// Make sure the results directory exists.
                        if (!std::filesystem::exists(results_dir_name))
                            std::filesystem::create_directory(results_dir_name);
                        // Clear stale transient output (VTK/CSV) from a previous
                        // run, but keep saved-state (.bin) files so earlier
                        // stages' checkpoints remain loadable.
                        for (const auto& entry : std::filesystem::directory_iterator(results_dir_name))
                        {
                            if (entry.is_regular_file() && entry.path().extension() != \".bin\")
                                std::filesystem::remove(entry.path());
                        }
                        $snapshot_dir_ir"

        # Write VTK/CSV/state output every `checkpoint_freq` steps (and always on the
        # final step so the stage's "latest" state file is complete). Default: every step.
        checkpoint_freq = get(stage_info, "checkpoint_frequency", nothing)
        if checkpoint_freq === nothing
            checkpoint_freq = 1
        end

        check_point = 
        "void checkpoint(std::size_t step) override {"*
                    "std::string results_dir_name = \"my_results\";"*
                    "const std::size_t checkpoint_every = $checkpoint_freq;
                    const bool checkpoint_due = (step != 0) &&
                        (step % checkpoint_every == 0 || step == static_cast<std::size_t>(settings.NUM_STEPS));"*
                    "if(step == 0)
                    {
                        // Create the local save directory
                        $dir_setup_ir

                        DGGML::GridFileWriter3D grid_writer;
                        grid_writer.save({geoplex2D.reaction_grid,geoplex2D.dim_label},
                                         results_dir_name+\"/expanded_cell_complex\");

                        std::string title = results_dir_name+\"/simulation_step_\";
                        DGGML::VtkFileWriterComplete<graph_type> vtk_writer;
                        auto attr = this->get_type_attributes(this->system_graph);
                        vtk_writer.set_extra_point_data(attr);

                        vtk_writer.save(system_graph, title+std::to_string(step));
                        collect(step);
                        { std::ofstream ts(results_dir_name+\"/timesteps.csv\", std::ios::trunc);
                          ts << \"step,time\\n\" << step << \",\" << (static_cast<double>(step) * settings.DELTA) << \"\\n\"; }
                    }
                    if( checkpoint_due )
                    {
                        std::string title = results_dir_name+\"/simulation_step_\";
                        DGGML::VtkFileWriterComplete<graph_type> vtk_writer;

                        auto attr = this->get_type_attributes(this->system_graph);
                        vtk_writer.set_extra_point_data(attr);

                        vtk_writer.save(system_graph, title+std::to_string(step));
                        collect(step);
                        { std::ofstream ts(results_dir_name+\"/timesteps.csv\", std::ios::app);
                          ts << step << \",\" << (static_cast<double>(step) * settings.DELTA) << \"\\n\"; }"*"}\n "*"if (save_system_graph && (step == 0 || checkpoint_due)) {save_graph(gamma, this->system_graph, settings, results_dir_name+\"$snapshot_prefix\"+std::to_string(step)+\".bin\");"*
                        "\nsave_graph(gamma, this->system_graph, settings, results_dir_name+\"/$stage_save_name\");}\n }\n"

        # Per-stage simulation time override
        stage_steps = stage_info["steps"]
        stage_time_override = ""
        if stage_steps !== nothing && stage_steps != ""
            stage_time_override = "
             // Per-stage simulation time override
             settings.TOTAL_TIME = $stage_steps;
             settings.NUM_STEPS = static_cast<int>(settings.TOTAL_TIME / settings.DELTA);
            "
        end

        # Rule-fire tracing requested with `trace(sim, "rule", ...)`. Always emitted so a stage
        # without a trace statement switches off tracing left on by an earlier stage.
        stage_trace_rules = get(stage_info, "trace_rules", String[])
        stage_trace_file = "trace_" * string(get(stage_info, "name", "stage$(stage_ix - 1)")) * ".jsonl"
        stage_trace_ir = "DGGML::Trace::begin_stage(std::vector<std::string>{" *
            join(["\"$(r)\"" for r in stage_trace_rules], ", ") *
            "}, \"$(stage_trace_file)\");\n"

        # Assemble the class
        class_ir = [
            "class $model_class_name : public DGGML::Model3D<graph_grammar_t> {\n",
            "\tpublic:\n",
            "\tParameters settings;\n",
            "bool save_system_graph = true;\n",
            "bool load_initial_state = $load_initial_state;\n",
            "std::string initial_state_filename = fflow_input_path($initial_state_file);\n",
            pde_members_ir,
            ""*model_create_initial_type,
            ""*ir_get_type_fn,
            ""*ir_global_obs_fn,
            ""*ir_load_graph_ir,
            "void initialize() override {\n",
            stage_time_override,
            stage_trace_ir,
            "int geoplex_size = settings.CELL_NX;",
            "
             int cell_nx = geoplex_size;
             int cell_ny = geoplex_size;
             int cell_nz = geoplex_size;

             double cell_dx = settings.CELL_DX;
             double cell_dy = settings.CELL_DY;
             double cell_dz = settings.CELL_DZ;

             double geoplex_epsilon = settings.varepsilon;
            ",
            "geoplex2D.init(
                            cell_nx,
                            cell_ny,
                            cell_nz,
                            cell_dx,
                            cell_dy,
                            cell_dz,
                            false,
                            geoplex_epsilon
                    );
            "*initial_rules_ir,
            "",
            "\n\t\t if (load_initial_state == false) {$initial_state_ir} else {load_graph(gamma, this->system_graph, settings, this->gen, 
                        this->initial_state_filename);}\n",
                "}\n",
            save_graph_ir,
            check_point,
            "}; // end $model_class_name\n\n",
        ]

        push!(model_classes, join(class_ir))
    end

    # ===== Post-run measurement (`measure(sim, ...)`) =====
    # Evaluates the observables on every snapshot my_results/<sim>/state_<step>.bin of a simulation
    # (of every replicate runs/<sim>/run_<i>/my_results/<sim>/state_<step>.bin when the simulation
    # was run N>1 times) and appends one row per snapshot to `out`
    # (columns: simulation,run,step,time,<observables>; run is 0 for a single run).
    measure_helper_ir = isempty(global_observables) ? "" : """
    template <class ModelT>
    void measure_simulation(const std::string &sim_name, const std::vector<std::string> &wanted,
                            std::ostream &out, bool &header_written, int replicates = 1) {
        namespace fs = std::filesystem;
        const int n_runs = replicates < 1 ? 1 : replicates;
        for (int run = 0; run < n_runs; ++run) {
        const std::string dir = replicates > 1
            ? "runs/" + sim_name + "/run_" + std::to_string(run) + "/my_results/" + sim_name
            : "my_results/" + sim_name;
        std::vector<std::pair<std::size_t, std::string>> snapshots;
        if (fs::exists(dir)) {
            for (const auto &e : fs::directory_iterator(dir)) {
                const std::string fname = e.path().filename().string();
                if (e.path().extension() != ".bin" || fname.rfind("state_", 0) != 0) continue;
                snapshots.emplace_back(std::stoull(fname.substr(6, fname.size() - 10)), e.path().string());
            }
        }
        if (snapshots.empty()) {
            std::cerr << "measure: no snapshots found in " << dir << " (was '" << sim_name << "' run?)\\n";
            continue;
        }
        std::sort(snapshots.begin(), snapshots.end());
        for (const auto &snap : snapshots) {
            ModelT m;
            m.load_graph(m.gamma, m.system_graph, m.settings, m.gen, snap.second);
            auto values = m.compute_observables(wanted);
            if (!header_written) {
                out << "simulation,run,step,time";
                for (const auto &v : values) out << "," << v.first;
                out << "\\n";
                header_written = true;
            }
            out << std::setprecision(10) << sim_name << "," << run << "," << snap.first << ","
                << (static_cast<double>(snap.first) * m.settings.DELTA);
            for (const auto &v : values) out << "," << v.second;
            out << "\\n";
        }
        }
    }
    """

    # ===== Footer =====
    model_footer = [measure_helper_ir, "};\n", "#endif\n"]

    join([model_header; model_classes; model_footer])
end

