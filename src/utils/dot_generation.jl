include("../ast_nodes.jl")
include("../tokens.jl")

function save_file(file_name, inp)
    open(file_name, "w") do file
        write(file, inp)
    end

end


function generate_type_section_edges(ast, name, ix)
    edges = ""

    # TODO: Add that it creates the 
    # branch of the AST also of these one
    for (jx, each_ast) in enumerate(ast)
        inst_value = replace(string(each_ast.name), "\"" => "\'")
        new_ix = jx + ix
        edges *= " \"$name\" -> \"$inst_value $new_ix\""
    end
    return edges
end

function generate_type_section_dot(ast)
    stack = []

    name = string(ast.name)
    name = replace(name, "\"" => "\'")

    type_verts = []

    dot_graph = """

    digraph {
        node [shape=rect];
    """
    for (ix, type_inst) in enumerate(ast.types)
        inst_name = string(type_inst.name)
        inst_name = replace(inst_name, "\"" => "\'")

        inst_param = string(type_inst.parameter)

        inst_type = string(type_inst.type)

        val_edge = nothing
        if isa(type_inst.value, Array)
            println("=== inside ast ===")
            val_edge = generate_type_section_edges(type_inst.value, inst_name, ix)
        else
            inst_value = string(type_inst.value)
            inst_value = replace(inst_value, "\"" => "\'")
            val_edge = "\"$inst_name\" -> \"$inst_value $ix \""
        end

        inst_param = replace(inst_param, "\"" => "\'")

        inst_type = replace(inst_type, "\"" => "\'")


        tmp_edges = """
        \"Type Section: $name\" -> \"$inst_name\"
        \"$inst_name\" -> \"Parameter: $inst_param $ix\"
        \"$inst_name\" -> \"Type: $inst_type $ix \"
        """

        dot_graph *= tmp_edges
        dot_graph *= val_edge

    end

    dot_graph *= "}"
    println("dot_graph: ", dot_graph)

    dot_graph
end

