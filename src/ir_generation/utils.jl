

module IRUtils

    export
        get_value,
        convert_type_name,
        write_file


    function get_value(ast)
            ast.token.position.value
    end

    function convert_type_name(type_name)
        """
        Convert type name to a proper c++ type
        """

        if type_name == "Float"
            return "double"
        elseif type_name == "Integer"
            return "int"
        else
            return type_name
        end
    end


    function write_file(file_name, ans)
        open("$(file_name)", "w+") do file
            write(file, ans)
        end
    end

end
