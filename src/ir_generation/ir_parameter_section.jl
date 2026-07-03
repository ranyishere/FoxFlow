
import ..AstNodes: CallNode, StringNode, TypeInstanceNode, TypeClassNode, IdentifierNode, IntegerNode, FloatNode

function ir_grammar_parameter!(ast, parameter_table)
    """
    IR Parameter
    """

    if isa(ast, TypeInstanceNode)

        param_name = get_value(ast.name)
        type_class_name = get_value(ast.type.name)

        # FixedList<<dims..., ElemType>>  emit torch::Tensor
        if type_class_name == "FixedList"
            inner = ast.type.parameter.token  # [IntegerNode/FloatNode, ...]
            dims = String[]
            elem_type = "Float"
            for p in inner
                val = p.token.position.value
                if val == "Float" || val == "Integer"
                    elem_type = val
                else
                    push!(dims, val)
                end
            end
            cpp_type = elem_type == "Float" ? "torch::kFloat64" : "torch::kInt64"
            dim_str  = join(dims, ", ")

            if isa(ast.value, CallNode) && get_value(ast.value.function_node) == "from_file"
                path = ast.value.args[1].token.position.value
                parameter_table[param_name] = ast
                # torch::from_file memory-maps a raw binary file — no pickle overhead
                total_elems = prod(parse.(Int, dims))
                return "\ttorch::Tensor $(param_name) = torch::from_file(\"$(path)\", false, $(total_elems), torch::TensorOptions().dtype($(cpp_type))).reshape({$(dim_str)});
"
            else
                parameter_table[param_name] = ast
                return "\ttorch::Tensor $(param_name) = torch::zeros({$(dim_str)}, $(cpp_type));
"
            end
        end

        grammar_parameter, param_names = ir_type_instance(ast)

        type_name = get_value(ast.name)
        parameter_table[type_name] = ast

        return grammar_parameter
    else
        println("did not expect this $(ast)")
    end

end

function ir_parameter_section(ast, propensity_table)

    parameter_section_data = [
                            "#ifndef DGGML_PARAMETERS_HPP\n",
                            "#define DGGML_PARAMETERS_HPP\n",
                            "#include <string>\n",
                            "#include \"simdjson.h\"\n"
                        ]

    section_name = get_value(ast.name)

    parameter_table = Dict()

    # TODO: Reference Parameter section name
    push!(parameter_section_data, "struct Parameters {\n")
    for (ix, param) in enumerate(ast.parameter_list)
        push!(parameter_section_data, ir_grammar_parameter!(param, parameter_table))
    end

    push!(parameter_section_data, "};\n")
    push!(parameter_section_data, "#endif")

    propensity_table["parameter_table"] = parameter_table
    join(parameter_section_data)
end
