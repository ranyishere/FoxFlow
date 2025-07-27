"""
IR Models Section
"""

function ir_models_section(ast)
    """
    Generate Intermediate Models
    for the Section
    """

    models_ir = IRBuilder([])

    model_boundary = ir_boundary!(models_ir)
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
            // const double max_x = 4.98;
            const double max_x = settings.CELL_NX-1;
            // const double max_y = 4.98;
            const double max_y = settings.CELL_NY-1;
            const std::size_t max_nx = 16;
            const std::size_t max_ny = 16;

            DGGML::CartesianGrid2D &reaction_grid = cplex.reaction_grid;
            std::size_t max_particles = reaction_grid.totalNumCells();

            //step 1 create an ordered number of selectable cells
            std::vector<std::size_t> selected_cells;
            for (int i = max_particles/2; i < max_particles; i++) {
                selected_cells.push_back(i);
            }

        auto centerCardinal = reaction_grid.cardinalCellIndex(
                        reaction_grid._nx/2, reaction_grid._ny/2
                    );

            double center_x, center_y;
            reaction_grid.cardinalCellToPoint(center_x, center_y, centerCardinal);
            Particles::graph_type tg;
            Particles::StartType tmp = Particles::StartType{};

            // Initial placement of particle node creator
            tmp.start_location[0] = center_x; // X
            tmp.start_location[1] = center_y; // Y
                                                       //
            node_type oof_node = {gen.get_key(),//i*segments,
                                {tmp, 
                                 center_x, center_y, 0.0}
                        };

            graph.addNode(oof_node);
            };"                               

    call_rule_ir = map(
        (rule_name) -> begin
            if rules_table[rule_name]["activated"] == true
                tmp_ir = "\t$(rules_table[rule_name]["namespace"])::$rule_name(gamma, this->system_graph, settings);\n"
                return tmp_ir
            end
        end,
        collect(keys(rules_table))
       )

    initial_rules_ir = join(call_rule_ir)

    model_section_header = [
        "#ifndef DGGML_MODELS_HPP\n",
        "#define DGGML_MODELS_HPP\n",
        "#include <fstream>\n",
        "#include \"DGGML.h\"\n",
        "#include \"rules.h\"\n",
        "#include \"simdjson.h\"\n",
        "#include \"ExpandedComplex2D.hpp\"\n",
        "#include \"YAGL_Algorithms.hpp\"\n",
        "#include \"parameters.h\"\n",
        "namespace $model_section_name {\n",
        "using graph_grammar_t = DGGML::Grammar<Particles::graph_type>;\n",
        "class Model : public DGGML::Model<graph_grammar_t> {\n",
        "\tpublic:\n",
        "\tParameters settings;\n",
        ""*model_create_initial_type,
    ""*model_boundary,
        "void initialize() override {\n",
        "geoplex2D.init(
                        settings.CELL_NX,
                        settings.CELL_NY,
                        settings.CELL_DX,
                        settings.CELL_DY,
                        settings.MAXIMAL_REACTION_RADIUS
                );
"*initial_rules_ir,
        "\n\t\tthis->add_default_type(this->system_graph,
                    geoplex2D,
                    settings,
                    this->gen);",

                    "\n\t\tthis->add_boundary(
                        geoplex2D.reaction_grid,
                        this->system_graph,
                        geoplex2D,
                        settings,
                        this->gen);",
            "}\n"
    ]

    model_section_body = []

    # Add the file writer
    check_point = 
    "void checkpoint(std::size_t step) override {"*
                "std::string results_dir_name = \"my_results\";"*
                "if(step == 0)
                {
                    // Create the local save directory
                    std::filesystem::remove_all(results_dir_name);
                    std::filesystem::create_directory(results_dir_name);

                        DGGML::GridFileWriter grid_writer;
                    grid_writer.save({geoplex2D.reaction_grid,geoplex2D.dim_label},
                                     results_dir_name+\"/expanded_cell_complex\");

                    std::string title = results_dir_name+\"/simulation_step_\";
                    DGGML::VtkFileWriter<graph_type> vtk_writer;
                    vtk_writer.save(system_graph, title+std::to_string(step));
                    collect(step);
                }
                if( step != 0 & step % 5 == 0)
                {
                    DGGML::GridFileWriter grid_writer;
                    grid_writer.save({geoplex2D.reaction_grid,geoplex2D.dim_label},
                                     results_dir_name+\"/expanded_cell_complex\");

                    std::string title = results_dir_name+\"/simulation_step_\";
                    DGGML::VtkFileWriter<graph_type> vtk_writer;
                    vtk_writer.save(system_graph, title+std::to_string(step));
                    collect(step);
                }
                }"

    model_section_footer = [check_point, "};\n};\n", "#endif"]
    model_section_data = [model_section_header; model_section_footer]

    # emit(models_ir, model_section_data)
    # build(models_ir)
    join(model_section_data)
end

