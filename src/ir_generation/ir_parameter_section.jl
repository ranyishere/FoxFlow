
function ir_grammar_parameter(ast)
    """
    IR Parameter
    """

    if isa(ast, TypeInstanceNode)
        grammar_parameter = ir_type_instance(ast)
        return grammar_parameter
    else
        println("did not expect this $(ast)")
    end
end

function ir_parameter_section(ast)
    parameter_section_data = [
                            "#ifndef DGGML_PARAMETERS_HPP\n",
                            "#define DGGML_PARAMETERS_HPP\n",
                            "#include <string>\n",
                            "#include \"simdjson.h\"\n"
                        ]

    section_name = get_value(ast.name)

    # TODO: Reference Parameter section name
    push!(parameter_section_data, "struct Parameters {\n")
    for (ix, param) in enumerate(ast.parameter_list)
        push!(parameter_section_data, ir_grammar_parameter(param))
    end

    push!(parameter_section_data, "};\n")

    push!(parameter_section_data, "#endif")

    join(parameter_section_data)
end
