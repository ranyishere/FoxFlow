
module IRBuildUtils

    export IRBuilder, emit, build

    mutable struct IRBuilder
        instructions::Array{Any}
    end

    function emit(ir_builder::IRBuilder, instruction::String)
        push!(ir_builder.instructions, instruction)
    end

    function build(ir_builder::IRBuilder)
        return join(ir_builder.instructions, "\n")
    end


end
