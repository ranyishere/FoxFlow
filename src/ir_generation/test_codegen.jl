include("../lexer/tokens.jl")
include("../parser/ast_nodes.jl")
include("../lexer/lexer.jl")
include("../parser/main_parser_op.jl")
include("./utils.jl")
include("./ir_builder.jl")
include("./codegen.jl")
include("./visualize_codegen.jl")

import .CodeGen: CGCtx
import .IRBuildUtils: IRBuilder, emit, build
import ..AstNodes: Node, LiteralNode, ModifyClauseNode
import .IRUtils: get_value, convert_type_name
import ..Tokens: IntegerToken, FloatToken, PositionToken, LiteralToken, ErrorToken, OperatorToken
import .PrintCodegen: print_node

function _eval(node::LiteralNode, ctx::CGCtx)
    """
    Evaluate a node
    """
    value = get_value(node)
    ctx.ir = emit(ctx.ir, value)
end


abstract type FileGenerator end

struct CPPFileGenerator <: FileGenerator
    input_file_name::String
    output_file_name::String
end

compose(f,g) = x -> f(g(x))

write_file(file::CPPFileGenerator) = ctx :: CGCtx  -> begin
    file_name = file.output_file_name
    ans = build(ctx.ir)
    open("$(file_name)", "w+") do file
        write(file, ans)
    end
end

struct ParameterFile <: FileGenerator 
    base :: CPPFileGenerator
end
ParameterFile(input_file_name::String) = ParameterFile(CPPFileGenerator(input_file_name, "parameters.h"))

eval!(node:: Node) = ctx :: CGCtx -> _eval!(node, ctx)

_eval!(node :: IdentifierNode, ctx :: CGCtx) = begin

    println("IdentifierNode: ", node)
    return ctx
end


_eval!(ast::ParameterSectionNode, ctx :: CGCtx) = begin
    # Compose
    # eval!(ast.name)(ctx)
    map((x) -> print_node(x)(ctx), ast.parameter_list)
    exit(0)
    # return compose ( (x) -> foldr(_print_node, x, ast.parameters), eval!(ast.name))(ctx)
end

eval!(p::ParameterFile) = ctx :: CGCtx -> begin
    tokens = tokenize_file(p.base.input_file_name)
    ast = parse_file!(tokens)[1]

    compose(write_file(p.base), eval!(ast))(ctx)
end

function main()
    """
    Main
    """

    param_ir = IRBuilder([])
    symbol_table = Dict{Symbol,String}()
    param_count = 0
    param_namespace = ""
    param_context = CGCtx(param_ir, symbol_table, 0, param_namespace)

    test_folder = "microtubules"
    name_space = "Microtubule"
    base = "../tests/generated_tests/generated_2/"

    test_base = "../tests/"
    grammar_param = ParameterFile(test_base*"$test_folder/params.fflow")
    eval!(grammar_param)(param_context)

    exit(0)

end

main()
