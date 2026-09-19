

module CodeGen

    export CGCtx

    import ..IRBuildUtils: IRBuilder

    struct CGCtx
        ir::IRBuilder
        env::Dict{Symbol,String}
        tmp::Int                           # counter for fresh temps
        type_namespace::String 
    end

end
