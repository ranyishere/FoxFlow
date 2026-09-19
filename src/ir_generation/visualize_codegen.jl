
module PrintCodegen

    import ..CodeGen: CGCtx
    import ..IRBuildUtils: IRBuilder, emit, build

    import ..AstNodes: Node, LiteralNode, ModifyClauseNode, IdentifierNode,
                        IntegerNode, ParameterNode, TypeClassNode, TypeInstanceNode,
                        BinaryOpNode

    import ..IRUtils: get_value, convert_type_name
    import ..Tokens: IntegerToken, FloatToken,
                    PositionToken, LiteralToken, ErrorToken, OperatorToken


    print_node(node :: Node) = ctx :: CGCtx -> _print_node(node, ctx)

    # Fallback for unimplemented nodes
    _print_node(node, ctx :: CGCtx) = begin println("Node: ", node) end

    _print_node(node :: IdentifierNode, ctx :: CGCtx) = println("IdentifierNode: ",
                                                                node.token.position.value)
    _print_node(node :: IntegerNode, ctx :: CGCtx) = println("IntegerNode: ",
                                                             node.token.position.value)

    _print_node(node :: ParameterNode, ctx :: CGCtx) = begin
        println("ParameterNode: ")
        map(_print_node, node.token)
    end

    _print_node(node :: TypeClassNode, ctx :: CGCtx) = begin
        println("=====")
        println("TYPECLASS")
        println("------------")
        ctx |> print_node(node.name)
        ctx |> print_node(node.parameter)
        println("=====")
    end

    _print_node(node :: TypeInstanceNode, ctx :: CGCtx) = begin
        println("=====")
        println("ASSIGNMENT")
        println("------------")
        ctx |> print_node(node.name)
        ctx |> print_node(node.parameter)
        ctx |> print_node(node.type)
        ctx |> print_node(node.value)
        println("=====")
    end

    _print_node(node :: BinaryOpNode, ctx :: CGCtx) = begin
        println("BinaryOpNode: ", node.expression.position.value)
        ctx |> print_node(node.lhs)
        ctx |> print_node(node.rhs)
    end

end
