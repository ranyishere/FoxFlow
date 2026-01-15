"""
IR Models Section
"""

function ir_load_graph()
    """
    Loads the intermediate representation for the graph
    """
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

function ir_get_type_attr(ast, type_namespace, symbol_tables)
    """
    IR Get Type Attribute
    """

    function ir_loop_type_attr(type_to_attrs, type_namespace)
        """
        Loop through type attributes and generate code
        that fetches attributes for all types
        """

        ir_loop = IRBuilder([])
        emit(ir_loop, "std::size_t row = 0;")
        emit(ir_loop, "for (auto it = system_graph.node_list_begin(); it != system_graph.node_list_end(); ++it, ++row) {")

        emit(ir_loop, "auto &n = it->second.getData();")


        emit(ir_loop, "std::visit([&](auto &alt) {")

        emit(ir_loop, "using T = std::decay_t<decltype(alt)>;")
        for (key, value) in type_to_attrs
            emit(ir_loop, "if constexpr (std::is_same_v<T, $type_namespace::$key>) {")
            for (col_name, attr) in value
                # Generate code to fetch attribute
                ir_fetch = IRBuilder([])
                emit(ir_fetch, "set(\"$col_name\", row, alt.$attr);")
                ir_fetch = build(ir_fetch)
                emit(ir_loop, ir_fetch)
            end
            emit(ir_loop, "}")
        end

        emit(ir_loop, "}, n.data);")

        emit(ir_loop, "}")

        return build(ir_loop)
    end

    ir_get_type_attr_fn = IRBuilder([])
    ir_hdr = "template <typename GraphType>"
    emit(ir_get_type_attr_fn, ir_hdr)

    ir_func = "auto get_type_attributes(GraphType &system_graph) {"
    emit(ir_get_type_attr_fn, ir_func)

    emit(ir_get_type_attr_fn, "std::vector<std::pair<std::string, std::vector<double>>> type_attributes;")

    ir_nan_value = "double nan_value = std::numeric_limits<double>::quiet_NaN();"
    emit(ir_get_type_attr_fn, ir_nan_value)

    ir_col_names = IRBuilder([])
    # Collect all type attributes
    col_names = []

    type_to_attrs = Dict{Any, Any}()
    for (key, value) in symbol_tables
        pos_to_attr = Dict{String, String}()
        for (i, attr) in value
            push!(col_names, "$key.$i")
            pos_to_attr["$key.$i"] = attr[2]
        end
        type_to_attrs[key] = pos_to_attr
    end

    emit(ir_col_names, "const std::array<std::string, $(length(col_names))> col_names = {")
    map(name -> emit(ir_col_names, "\"$name\", "), col_names)
    emit(ir_col_names, "};\n")
    ir_col_names = build(ir_col_names)
    emit(ir_get_type_attr_fn, ir_col_names)

    emit(ir_get_type_attr_fn, "std::array<std::vector<double>, $(length(col_names))> columns;")
    emit(ir_get_type_attr_fn, "std::size_t N = system_graph.numNodes();")

    # Fill NaN values
    emit(ir_get_type_attr_fn, "for (auto &col : columns) {
        col.resize(N, nan_value);
        }")

    # col_idx
    emit(ir_get_type_attr_fn, "std::unordered_map<std::string, std::size_t> col_idx;")
    emit(ir_get_type_attr_fn, "for (std::size_t i = 0; i < col_names.size(); ++i) {
        col_idx[col_names[i]] = i;
        }")

    # Setter helper function
    ir_setter = IRBuilder([])
    emit(ir_setter, """
        auto set = [&](const std::string &col_name, std::size_t row, double value) {
            columns[col_idx.at(col_name)][row] = value;
        };
    """)

    ir_setter = build(ir_setter)
    emit(ir_get_type_attr_fn, ir_setter)

    # Loop through graph nodes
    ir_loop_type_attr_code = ir_loop_type_attr(type_to_attrs, type_namespace)
    emit(ir_get_type_attr_fn, ir_loop_type_attr_code)

    # Emplace back
    emit(ir_get_type_attr_fn, "for (std::size_t i = 0; i < col_names.size(); ++i) {
        type_attributes.emplace_back(col_names[i], std::move(columns[i]));
    }")
    
    emit(ir_get_type_attr_fn, "return type_attributes;}")
    return build(ir_get_type_attr_fn)
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
    """
    Generate Intermediate Models
    for the Section
    """

    models_ir = IRBuilder([])

    ir_get_type_fn = ir_get_type_attr(ast, type_namespace, symbol_tables)

    # emit(models_ir, ir_get_type_fn)
    # model_boundary = ir_boundary!(models_ir)

    model_section_name = type_namespace

    # model_section_name = get_value(ast.name)
    # Create boundary
    
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

            // const double max_x = 4.98;
            const double max_x = settings.CELL_NX-1;
            // const double max_y = 4.98;
            const double max_y = settings.CELL_NY-1;

            const std::size_t max_nx = 16;
            const std::size_t max_ny = 16;
        const std::size_t max_nz = 16;

            // DGGML::CartesianGrid2D &reaction_grid = cplex.reaction_grid;
            DGGML::CartesianGrid3D &reaction_grid = cplex.reaction_grid;
            std::size_t max_particles = reaction_grid.totalNumCells();

            //step 1 create an ordered number of selectable cells
            std::vector<std::size_t> selected_cells;
            for (int i = max_particles/2; i < max_particles; i++) {
                selected_cells.push_back(i);
            }

        auto centerCardinal = reaction_grid.cardinalCellIndex(
                        reaction_grid._nx/2, reaction_grid._ny/2, reaction_grid._nz/2
                    );

            double center_x, center_y, center_z;
            // reaction_grid.cardinalCellToPoint(center_x, center_y, centerCardinal);
            reaction_grid.cardinalCellToPoint(center_x, center_y, center_z, centerCardinal);
            $type_namespace::graph_type tg;
            $type_namespace::StartType tmp = $type_namespace::StartType{};

            // Initial placement of particle node creator
            tmp.start_location[0] = center_x; // X
            tmp.start_location[1] = center_y; // Y
            tmp.start_location[2] = center_z; // Y
                                                       //
            node_type oof_node = {gen.get_key(),//i*segments,
                                {tmp, 
                                 center_x, center_y, center_z}
                        };

            graph.addNode(oof_node);
            };"                               

    # Create rules
    call_rule_ir = map(
        (rule_name) -> begin
            if rules_table[rule_name]["activated"] == true
                tmp_ir = "\t$(rules_table[rule_name]["namespace"])::$rule_name(gamma, this->system_graph, settings);\n"
                return tmp_ir
            end
        end,
        collect(keys(rules_table))
    )

    serializer_spatialnode_ir = ir_serialize_spatialnode()

    initial_rules_ir = join(call_rule_ir)


    model_section_header = [
        "#ifndef DGGML_MODELS_HPP\n",
        "#define DGGML_MODELS_HPP\n",
        "#include <fstream>\n",
        "#include \"DGGML.h\"\n",
        "#include \"rules.h\"\n",
        "#include \"simdjson.h\"\n",
        "#include \"ExpandedComplex2D.hpp\"\n",
        "#include \"ExpandedComplex3D.hpp\"\n",
        "#include \"YAGL_Algorithms.hpp\"\n",
        "#include \"parameters.h\"\n",
        "#include <cereal/archives/binary.hpp>\n",
        "#include <cereal/types/variant.hpp>\n",
        "#include <cereal/types/vector.hpp>\n",
        "#include <cereal/types/unordered_set.hpp>\n",
        serializer_spatialnode_ir,
        "namespace $model_section_name {\n",
        "using graph_grammar_t = DGGML::Grammar<$type_namespace::graph_type>;\n",
        "class Model : public DGGML::Model3D<graph_grammar_t> {\n",
        "\tpublic:\n",
        "\tParameters settings;\n",
        ""*model_create_initial_type,
        # ""*model_boundary,
        ""*ir_get_type_fn,
        "void initialize() override {\n",
        "int geoplex_size = 3;
         int cell_nx = geoplex_size;
         int cell_ny = geoplex_size;
         int cell_nz = geoplex_size;

         double cell_dx = 1.0;
         double cell_dy = 1.0;
         double cell_dz = 1.0;

         double maximal_rx_radius = 0.05;
        ",
        "geoplex2D.init(
                        cell_nx,
                        cell_ny,
                        cell_nz,
                        cell_dx,
                        cell_dy,
                        cell_dz,
                        false,
                        maximal_rx_radius
                );
        "*initial_rules_ir,
        "\n\t\tthis->add_default_type(this->system_graph,
                    geoplex2D,
                    settings,
                    this->gen);",

                    # "\n\t\tthis->add_boundary(
                        # geoplex2D.reaction_grid,
                        # this->system_graph,
                        # geoplex2D,
                        # settings,
                        # this->gen);",
            "}\n"
    ]

    model_section_body = []

    save_graph_ir = ir_save_graph(type_namespace)

    # Add the file writer
    check_point = 
    "void checkpoint(std::size_t step) override {"*
                "std::string results_dir_name = \"my_results\";"*
                "if(step == 0)
                {
                    // Create the local save directory
                    std::filesystem::remove_all(results_dir_name);
                    std::filesystem::create_directory(results_dir_name);

                    DGGML::GridFileWriter3D grid_writer;
                    grid_writer.save({geoplex2D.reaction_grid,geoplex2D.dim_label},
                                     results_dir_name+\"/expanded_cell_complex\");

                    std::string title = results_dir_name+\"/simulation_step_\";
                    DGGML::VtkFileWriterComplete<graph_type> vtk_writer;
                    auto attr = this->get_type_attributes(this->system_graph);
                    vtk_writer.set_extra_point_data(attr);

                    vtk_writer.save(system_graph, title+std::to_string(step));
                    collect(step);
                }
                if( step != 0 & step % 1 == 0)
                {
                    std::string title = results_dir_name+\"/simulation_step_\";
                    DGGML::VtkFileWriterComplete<graph_type> vtk_writer;

                    auto attr = this->get_type_attributes(this->system_graph);
                    vtk_writer.set_extra_point_data(attr);

                    vtk_writer.save(system_graph, title+std::to_string(step));
                    collect(step);"*"}\n }"
            # "save_graph(gamma, this->system_graph, settings, results_dir_name+\"/simulation_state_\"+std::to_string(step)+\".bin\");"

    model_section_footer = [save_graph_ir, check_point, "};\n};\n", "#endif"]
    model_section_data = [model_section_header; model_section_footer]

    # emit(models_ir, model_section_data)
    # build(models_ir)
    join(model_section_data)
end

