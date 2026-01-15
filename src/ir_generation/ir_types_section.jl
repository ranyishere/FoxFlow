

function ir_serialize(fields)

    temp_ir = "\ttemplate <class Archive>\n"
    hdr_ir = "\tvoid serialize(Archive& archive) {\n"
    temp_ir *= hdr_ir
    for field in fields
        field_name = field
        temp_ir *= "\t\tarchive($(field_name));\n"
    end
    # hdr_ir *= "\t}\n"
    temp_ir *= "\t}\n"
end


function generate_ir_type_struct(name, fields)
    """
    Generates the IR for a type struct given its name and fields.
    """

    struct_ir = "struct $(name) {\n"
    if length(fields) == 0
        field_names = [last(x) for x in fields]
        ir_serialize_ir = ir_serialize(field_names)
        struct_ir *= ir_serialize_ir
    else
        for field in fields

            field_type = field[1]
            field_name = last(field)

            if field_type == "list"
                size = field[2]
                elem_type = field[3]
                struct_ir *= "\t$(elem_type) $(field_name)[$(size)];\n"
            else 
                struct_ir *= "\t$(field_type) $(field_name);\n"
            end
        end
        field_names = [last(x) for x in fields]
        ir_serialize_ir = ir_serialize(field_names)
        struct_ir *= ir_serialize_ir
    end
    struct_ir *= "};\n"
    struct_ir
end

function ir_type_struct(ast, symbol_tables)
    """
    Generates the IR for a type struct.
    """
end

function ir_types_section(ast, symbol_tables)

    section_name = get_value(ast.name)

    global type_namespace = section_name

    section_data = [
        "#ifndef DGGML_$(section_name)_TYPES_HPP\n",
        "#define DGGML_$(section_name)_TYPES_HPP\n",
        "#include \"YAGL_Graph.hpp\" \n",
        "#include \"YAGL_Node.hpp\" \n",
        "#include \"SpatialData3D.hpp\" \n",
        "namespace $(section_name) {\n",

        # "\tstruct Type {};\n",
        # "\tstruct Boundary {
            # float boundary_location[3];\n
        # };\n",
        # "\tstruct StartType {
            # float start_location[3];\n
       # };\n"
    ]

    default_type = ["Type", []]

    generated_types = Dict(
        "Type" => default_type,
        "StartType" => ["StartType", [("list", "3", "float", "start_location")]],
        "Boundary" => ["Boundary", [("list", "3", "float", "boundary_location")]]
    )

    for (type_name, (struct_name, fields)) in generated_types
        type_struct_ir = generate_ir_type_struct(struct_name, fields)
        println(type_struct_ir)
        push!(section_data, type_struct_ir)
    end

    symbol_tables["StartType"] = Dict()
    symbol_tables["Boundary"] = Dict()
    start_type_info = (
                   true, 
                   IntegerNode(IntegerToken(PositionToken("3", -1, -1))), 
                   FloatNode(FloatToken(PositionToken("Float", -1, -1)))
               )
    boundary_type_info = (
                   true, 
                   IntegerNode(IntegerToken(PositionToken("3", -1, -1))), 
                   FloatNode(FloatToken(PositionToken("Float", -1, -1)))
               )


    symbol_tables["StartType"][1] = ("float", "start_location[0]", start_type_info)
    symbol_tables["StartType"][2] = ("float", "start_location[1]", start_type_info)
    symbol_tables["StartType"][3] = ("float", "start_location[2]", start_type_info)

    symbol_tables["Boundary"][1] = ("float", "boundary_location[0]", boundary_type_info)
    symbol_tables["Boundary"][2] = ("float", "boundary_location[1]", boundary_type_info)
    symbol_tables["Boundary"][3] = ("float", "boundary_location[2]", boundary_type_info)

    type_instances_ast = ast.types
    type_instance_names = []
    for type_inst_ast in type_instances_ast

        cur_type_inst, param_names = ir_type_instance(type_inst_ast,
                                                      true,
                                                      symbol_tables)

        ir_serialize_ir = ir_serialize(param_names)

        # cur_type_inst = cur_type_inst*ir_serialize_ir*"};\n"

        cur_type_inst = cur_type_inst*ir_serialize_ir*"\n"

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

