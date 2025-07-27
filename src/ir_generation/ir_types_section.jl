

function ir_types_section(ast)

    section_name = get_value(ast.name)

    global type_namespace = section_name

    section_data = [
        "#ifndef DGGML_$(section_name)_TYPES_HPP\n",
        "#define DGGML_$(section_name)_TYPES_HPP\n",
        "#include \"YAGL_Graph.hpp\" \n",
        "#include \"YAGL_Node.hpp\" \n",
        "#include \"SpatialData3D.hpp\" \n",
        "namespace $(section_name) {\n",
        "\tstruct Type {};\n",
        "\tstruct Boundary {};\n",
        "\tstruct StartType {
            float start_location[2];\n
       };\n"
    ]

    symbol_tables["StartType"] = Dict()
    start_type_info = (
                   true, 
                   IntegerNode(IntegerToken(PositionToken("2", -1, -1))), 
                   FloatNode(FloatToken(PositionToken("Float", -1, -1)))
               )

    symbol_tables["StartType"][1] = ("float", "start_location[0]", start_type_info)
    symbol_tables["StartType"][2] = ("float", "start_location[1]", start_type_info)

    type_instances_ast = ast.types
    type_instance_names = []
    for type_inst_ast in type_instances_ast

        cur_type_inst = ir_type_instance(type_inst_ast, true)
        cur_type_inst = cur_type_inst*"};"
        type_name = get_value(type_inst_ast.name)
        push!(type_instance_names, type_name)
        push!(section_data, cur_type_inst)
    end

    spatial_node_3d = create_spatial_node_3d(type_instance_names)

    bottom_section_data = ["\tusing graph_type = YAGL::Graph<std::size_t,",
                               spatial_node_3d, ">;\n",
                           "};\n", "#endif"
                        ]

    section_data = [section_data; bottom_section_data]

    final_str = join(section_data)
    final_str
end

