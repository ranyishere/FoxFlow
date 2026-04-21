module IRUtils

    export
        get_value,
        convert_type_name,
        write_file, is_list_type_tensor


    function get_value(ast)
            ast.token.position.value
    end

    function convert_type_name(type_name, list=false)
        """
        Convert type name to a proper c++ type
        """

        if type_name == "Float"
            # return "double"
            if list == true
                return "torch::kFloat64"
            else
                return "double"
            end

        elseif type_name == "Integer"
            if list == true
                # return "std::vector<int>"
                return "torch::kInt64"
            else
                return "int"
            end
            # return "int"
        elseif type_name == "FixedList"
            return "torch::Tensor"
        else
            return type_name
        end

    end

    function is_list_type_tensor(type_name)
        if type_name == "torch::Tensor"
            return true
        else
            return false
        end
    end

    function write_file(file_name, ans)

        open("$(file_name)", "w+") do file
            write(file, ans)
        end

    end

end
