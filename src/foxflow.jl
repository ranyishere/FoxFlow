include("./parser.jl")
include("./lexer.jl")
include("./semantics.jl")
include("./ir_generation.jl")


function save_files(code)
    """
    Saving File
    """

    cur_wd = pwd()
    println("before cur_wd: $cur_wd")
    cd("/home/rany/Work/research/FoxFlow/src")
    cur_wd = pwd()
    println("after cur_wd: $cur_wd")

    header_files = []

    for (key, value) in code

        println("Saving: $key")

        if (key == "settings")
            open("../generated/$key.json", "w") do f
                write(f, value)
                flush(f) # Ensures the content is written to disk
            end
        elseif (key == "main")
            open("../generated/$key.cpp", "w") do f
                write(f, value)
                flush(f)
            end
        else
            open("../generated/$key.h", "w") do f
                write(f, value)
                flush(f)
            end
        end
    end

end

function compiles()
    """
    Compiles the simulation
    """

    # TODO: Only run build if things have changed drastically?
    build_folder = dirname(pwd())*"/build/"
    compile_cmd = Cmd(`./build.sh`, dir=build_folder)

    run(compile_cmd)
    build_folder_gen = dirname(pwd())*"/build/generated"
    compile_cmd = Cmd(`make`, dir=build_folder_gen)
    run(compile_cmd)

end

function execute_sim()
    """
    Executes the Simulation
    """

    build_folder_gen = dirname(pwd())*"/build/generated"
    exec_sim = Cmd(`./main settings.json`, dir=build_folder_gen)
    # exec_sim = Cmd(`cd generated/\; ./main settings.json`)
    run(exec_sim)
end

function foxflow_run(source_code)

    # println("source_code: ", source_code)
    tokens = tokenize(source_code, 0)
    println("running thingamabob")
    println("tokens: ", tokens)
    symbol_table = Dict()
    symbol_table["null"] = nothing
    grammars = parse!(tokens, symbol_table)

    println("Showing grammars")
    # println("grammars: ", grammars)
    generated_code = generate_im(grammars, symbol_table)

    println("Done generating code")
    # TODO: Give return code 0 or 1 depending on how the run is.
    save_files(generated_code)
    println("done Saving file")
    compiles()
    println("Done compiling")
    execute_sim()
    println("done executing sim")
end

if length(ARGS) != 0
    foxflow_run(ARGS[1])
else
    println("Nothing was sent")
end

