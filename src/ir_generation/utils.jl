module IRUtils

    export
        get_value,
        convert_type_name,
        convert_type_node,
        write_file,
        is_list_type_tensor,
        is_tf_array_type,
        is_torch_tensor_type,
        tf_array_size,
        fixed_list_dims,
        TF_RUNTIME_HEADER

    # Name of the generated runtime header that defines TFArray<N>.
    const TF_RUNTIME_HEADER = "tf_array.h"

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
        elseif type_name == "TFTensor"
            return "torch::Tensor"
        elseif type_name == "Tensor"
            return "torch::Tensor"
        else
            return type_name
        end

    end

    # True for the C++ type of a torch tensor attribute/variable.
    function is_torch_tensor_type(type_name)
        return type_name == "torch::Tensor"
    end

    # True for the C++ type of a fixed-size stack array (TFArray<N>).
    function is_tf_array_type(type_name)
        return type_name isa AbstractString && startswith(type_name, "TFArray<")
    end

    # True for any list-like attribute type: a torch tensor or a TFArray<N>.
    function is_list_type_tensor(type_name)
        return is_torch_tensor_type(type_name) || is_tf_array_type(type_name)
    end

    # Number of elements of a "TFArray<N>" type string.
    function tf_array_size(type_name)
        m = match(r"^TFArray<(\d+)>$", type_name)
        m === nothing && throw("Not a TFArray type: $type_name")
        return parse(Int, m.captures[1])
    end

    function _param_value(node)
        try
            return string(node.token.position.value)
        catch
            return nothing
        end
    end

    """
        fixed_list_dims(type_node) -> (dims::Vector{String}, elem_type::String)

    Reads the `<<d1, d2, ..., Elem>>` parameters of a FixedList/TFTensor type node.
    """
    function fixed_list_dims(type_node)
        dims = String[]
        elem = "Float"
        params = type_node.parameter.token
        params === nothing && return dims, elem
        for p in params
            val = _param_value(p)
            val === nothing && continue
            if val == "Float" || val == "Integer"
                elem = val
            else
                push!(dims, val)
            end
        end
        return dims, elem
    end

    """
        convert_type_node(type_node; dimless="auto") -> String

    C++ type of a FoxFlow type node.

    - `FixedList<<N, Float>>`  -> `TFArray<N>` (stack array)
    - `FixedList` (no dims)    -> `dimless` (default `auto`, deduced from the value)
    - `TFTensor`               -> `torch::Tensor`
    - anything else            -> `convert_type_name`
    """
    function convert_type_node(type_node; dimless="auto")
        name = get_value(type_node.name)
        if name == "FixedList"
            dims, elem = fixed_list_dims(type_node)
            if isempty(dims)
                return dimless
            end
            if length(dims) != 1 || elem != "Float"
                throw("FixedList supports only one-dimensional Float arrays " *
                      "(got <<$(join(dims, ", ")), $elem>>). " *
                      "Use TFTensor for multi-dimensional or Integer data.")
            end
            return "TFArray<$(dims[1])>"
        end
        return convert_type_name(name)
    end

    function write_file(file_name, ans)

        open("$(file_name)", "w+") do file
            write(file, ans)
        end

    end

end
