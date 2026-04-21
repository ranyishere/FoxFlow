

function ir_type_declaration(ast, define_type=false, symbol_tables=nothing)
    """
    IR Type Instance
    """

    type_name = get_value(ast.name)

    type_parameters = []

    # Track for each type, what parameters does it belong to?
    if define_type == true
        symbol_tables[type_name] = Dict()
    end

    # Associate each parameter with the type and parameter position.
    type_parameters, param_names, param_types, is_list = ir_parameter(ast.parameter)

    println("type_name: ", type_name)
    if define_type == true
        param_count = 1
        println("param_names: ", param_names)
        for (ix, name) in enumerate(param_names)

            println("is_list: ", is_list[ix])
            if is_list[ix][1] == true

                list_size = is_list[ix][2]
                
                # println("Tensor size: ", tensor_size)
                # exit(0)
                # list_size = is_list[ix][2].token.position.value
                # list_size = tensor_size

                # list_type = is_list[ix][3].token.position.value
                list_type = "torch::Tensor"
                # println("list_type: ", list_type)
                # println("is_list[ix]: ", is_list[ix])
                # exit(0)

                symbol_tables[type_name][param_count] = (list_type, name, is_list[ix])
                param_count += 1

                # list_size_int = parse(Int64, list_size)
                # for i in 1:list_size_int
                    # symbol_tables[type_name][param_count] = (list_type, name*"[$(i-1)]", is_list[ix])
                    # param_count += 1
                # end

            else
                println("param_types[ix] ========> : ", param_types[ix])
                # It's not a list
                symbol_tables[type_name][param_count] = (param_types[ix], name, is_list[ix])
                param_count += 1
            end

        end
    end

    type_class_name = ir_type_class(ast.type)
    begin_struct = []

    # TODO: Handle identifier names
    if ast.value == nothing
        begin_struct = ["\tstruct $(type_name) : $(type_class_name) {\n"]
        begin_struct = [begin_struct;type_parameters]

        # operator_like_array = ir_struct_like_array(param_names)

        # begin_struct = [begin_struct;operator_like_array]
    elseif isa(ast.value, LiteralToken)

        ir_type = convert_type_name(type_class_name)
        begin_struct = ["\t"*ir_type*" $(type_name) "*"= $(ast.value.position.value);"]
    elseif isa(ast.value, IntegerNode)

        ir_type = convert_type_name(get_value(ast.type.name))
        begin_struct = ["\t"*ir_type*" $(type_name) "*"= $(ast.value.token.position.value);"]
    elseif isa(ast.value, BinaryOpNode)
        # It's a binary operation
        expression = ir_value(ast.value)
        ir_type = convert_type_name(type_class_name)
        begin_struct = ["\t"*ir_type*" $(type_name) "*"= $(expression);"]

    elseif isa(ast.value, IdentifierNode)
        ir_type = convert_type_name(type_class_name)
        begin_struct = ["\t"*ir_type*" $(type_name) "*"= $(ast.value.token.position.value);"]
    end

    push!(begin_struct,"\n")
    join(begin_struct), param_names
end



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
                # row of zeros
                struct_ir *= "\t$(elem_type) $(field_name) = torch::zeros($(size), torch::kFloat64);\n"
            elseif field_type == "double"
                struct_ir *= "\t$(field_type) $(field_name) = 0.0;\n"
            elseif field_type == "int"
                struct_ir *= "\t$(field_type) $(field_name) = 0;\n"
            else 
                struct_ir *= "\t$(field_type) $(field_name){};\n"
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
        "#include \"torch/torch.h\"\n",
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
        "StartType" => ["StartType", [("list", "3", "torch::Tensor", "start_location")]],
        "Boundary" => ["Boundary", [("list", "3", "torch::Tensor", "boundary_location")]]
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


    # symbol_tables["StartType"][1] = ("double", "start_location[0]", start_type_info)
    # symbol_tables["StartType"][2] = ("double", "start_location[1]", start_type_info)
    # symbol_tables["StartType"][3] = ("double", "start_location[2]", start_type_info)

    # symbol_tables["Boundary"][1] = ("double", "boundary_location[0]", boundary_type_info)
    # symbol_tables["Boundary"][2] = ("double", "boundary_location[1]", boundary_type_info)
    # symbol_tables["Boundary"][3] = ("double", "boundary_location[2]", boundary_type_info)

    # TODO: check if we need to convert kFloat64 when exporting to csv? 
    symbol_tables["StartType"][1] = ("torch::Tensor", "start_location", start_type_info)
    symbol_tables["Boundary"][1] = ("torch::Tensor", "boundary_location", boundary_type_info)

    type_instances_ast = ast.types
    type_instance_names = []
    for type_inst_ast in type_instances_ast

        # println("type_inst_ast: ", type_inst_ast)
        # for param in type_inst_ast.parameter.token
            # println("param: ", param)
        # end
        # exit(0)
        cur_type_inst, param_names = ir_type_declaration(type_inst_ast,
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

