include("./parser.jl")
include("./tokens.jl")
include("./semantics.jl")
include("./ir_generation.jl")

function save_files(code)
    cur_wd = pwd()
    # println("cur_wd: ", cur_wd)

    header_files = []
    for (key, value) in code
        if (key == "settings")
            write("../generated/$key.json", value)
        elseif (key == "main")
            write("../generated/$key.cpp", value)
        else
            write("../generated/$key.h", value)
        end
        # write("../generated/main_simulator.cpp", code["main"])
    end
end

function compiles()
    """
    Compiles the simulation
    """

    # println("pwd: ", pwd())
    # cd(pwd())

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
    tokens = tokenize(source_code)
    println("tokens: ", tokens)

    symbol_table = Dict()
    symbol_table["null"] = nothing
    grammars = parse!(tokens, symbol_table)

    # println("grammars: ", grammars)
    generated_code = generate_im(grammars, symbol_table)
    # TODO: Give return code 0 or 1 depending on how the run is.
    save_files(generated_code)
    compiles()
    execute_sim()
end

if length(ARGS) != 0
    foxflow_run(ARGS[1])
else
    println("Nothing was sent")
end
