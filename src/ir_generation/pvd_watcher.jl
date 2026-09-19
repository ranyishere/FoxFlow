"""
PvdWatcher — dependency-free ParaView `.pvd` collection generator/watcher.

The generated DGGML model writes a VTK time series into `my_results/`:
    my_results/simulation_step_<step>.vtu   (one per output step)

ParaView groups those into a series, but only by integer frame index. This
module writes/refreshes a `.pvd` collection file so ParaView sees a proper
time series (timestep == step). Run it alongside a simulation and it keeps the
`.pvd` up to date; open the `.pvd` in ParaView and hit "Reload Files" to watch
new frames appear.

Keeps the visualization pipeline entirely in Julia + C++ (no Python/pvpython).
"""
module PvdWatcher

export watch_pvd, write_pvd, run_cli

# Match `..._step_<N>.vtu` and capture the integer step.
const STEP_RE = r"_step_(\d+)\.vtu$"

"""
    read_time_map(dir; name) -> Dict{Int,Float64}

Read the optional `timesteps.csv` sidecar (written by the generated
`checkpoint()`), mapping each integer step to its real simulation time.
Returns an empty dict if the file is missing/unreadable so callers can fall
back to the step index.
"""
function read_time_map(dir::AbstractString; name::AbstractString="timesteps.csv")
    m = Dict{Int,Float64}()
    path = joinpath(dir, name)
    isfile(path) || return m
    try
        for (n, line) in enumerate(eachline(path))
            line = strip(line)
            isempty(line) && continue
            n == 1 && startswith(lowercase(line), "step") && continue  # header
            parts = split(line, ',')
            length(parts) >= 2 || continue
            step = tryparse(Int, strip(parts[1]))
            t = tryparse(Float64, strip(parts[2]))
            (step === nothing || t === nothing) && continue
            m[step] = t
        end
    catch
        # A concurrent write can truncate mid-read; just use what we parsed.
    end
    return m
end

"""
    collect_frames(dir; pattern, settle) -> Vector{Tuple{Int,String}}

Return `(step, filename)` pairs for VTK frames in `dir` whose name contains
`pattern` and ends in `_step_<N>.vtu`, sorted by step. Files modified within
`settle` seconds are skipped so a frame that is still being written is not
indexed mid-write.
"""
function collect_frames(dir::AbstractString; pattern::AbstractString="simulation_step_",
                        settle::Real=0.5)
    frames = Tuple{Int,String}[]
    isdir(dir) || return frames
    now_t = time()
    for f in readdir(dir)
        occursin(pattern, f) || continue
        m = match(STEP_RE, f)
        m === nothing && continue
        full = joinpath(dir, f)
        # Skip partially-written files (recently modified).
        if settle > 0 && (now_t - mtime(full)) < settle
            continue
        end
        push!(frames, (parse(Int, m.captures[1]), f))
    end
    sort!(frames; by = x -> x[1])
    return frames
end

"""
    write_pvd(dir, frames; out, time_map) -> path

Write a ParaView collection (`.pvd`) referencing `frames` (from
`collect_frames`) by relative filename. If `time_map` supplies a real
simulation time for a frame's step it is used as the `timestep`; otherwise the
integer step index is used. ParaView shows this value as "Time".
"""
function write_pvd(dir::AbstractString, frames::Vector{Tuple{Int,String}};
                   out::AbstractString="simulation.pvd",
                   time_map::Dict{Int,Float64}=Dict{Int,Float64}())
    out_path = isabspath(out) ? out : joinpath(dir, out)
    open(out_path, "w") do io
        println(io, "<?xml version=\"1.0\"?>")
        println(io, "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">")
        println(io, "  <Collection>")
        for (step, file) in frames
            t = get(time_map, step, Float64(step))
            println(io,
                "    <DataSet timestep=\"$(t)\" group=\"\" part=\"0\" file=\"$(file)\"/>")
        end
        println(io, "  </Collection>")
        println(io, "</VTKFile>")
    end
    return out_path
end

"""
    watch_pvd(dir; pattern, out, interval, settle, once)

Poll `dir` and rewrite the `.pvd` whenever the set of frames changes.
Set `once=true` to write a single snapshot and return (no polling loop).

`stop` is an optional zero-argument predicate; when it returns `true` the
watcher performs one final flush (picking up any last frames) and returns. This
lets a caller run the watcher alongside a simulation and stop it automatically
when the process exits.
"""
function watch_pvd(dir::AbstractString; pattern::AbstractString="simulation_step_",
                   out::AbstractString="simulation.pvd", interval::Real=1.0,
                   settle::Real=0.5, once::Bool=false,
                   stop::Union{Function,Nothing}=nothing)
    if once
        frames = collect_frames(dir; pattern=pattern, settle=0.0)
        tmap = read_time_map(dir)
        path = write_pvd(dir, frames; out=out, time_map=tmap)
        println("[pvd] Wrote $(length(frames)) frame(s) -> $path")
        return path
    end

    println("[pvd] Watching '$dir' (pattern='$pattern', interval=$(interval)s)")
    println("[pvd] Open the .pvd in ParaView and use 'Reload Files' to advance. Ctrl-C to stop.")
    last_count = -1
    out_path = isabspath(out) ? out : joinpath(dir, out)
    try
        while true
            frames = collect_frames(dir; pattern=pattern, settle=settle)
            if length(frames) != last_count
                tmap = read_time_map(dir)
                write_pvd(dir, frames; out=out, time_map=tmap)
                last_count = length(frames)
                if isempty(frames)
                    latest = "-"
                else
                    s = frames[end][1]
                    latest = haskey(tmap, s) ? "t=$(tmap[s])" : "step $s"
                end
                println("[pvd] $(length(frames)) frame(s) (latest $latest) -> $out_path")
            end
            if stop !== nothing && stop()
                # Final flush: no settle window so the very last frame is caught.
                frames = collect_frames(dir; pattern=pattern, settle=0.0)
                tmap = read_time_map(dir)
                write_pvd(dir, frames; out=out, time_map=tmap)
                println("[pvd] Simulation finished — $(length(frames)) frame(s) -> $out_path")
                break
            end
            sleep(interval)
        end
    catch e
        e isa InterruptException || rethrow(e)
        println("\n[pvd] Stopped.")
    end
    return out_path
end

"""
    run_cli(args) -> Int

Standalone CLI:
  julia pvd_watcher.jl <results_dir> [--out FILE] [--pattern P]
                       [--interval SECS] [--settle SECS] [--once]
"""
function run_cli(args)
    if isempty(args) || args[1] in ("-h", "--help", "help")
        println("""
        PvdWatcher — keep a ParaView .pvd collection in sync with a running sim.

        Usage:
          julia pvd_watcher.jl <results_dir> [options]

        Arguments:
          results_dir     Directory the model writes .vtu frames into (e.g. my_results).

        Options:
          --out FILE       Output .pvd name/path (default: simulation.pvd in results_dir).
          --pattern P      Frame filename substring (default: simulation_step_).
          --interval SECS  Poll interval (default: 1.0).
          --settle SECS    Ignore files modified within this window (default: 0.5).
          --once           Write a single snapshot and exit (no watching).
        """)
        return 0
    end

    dir = args[1]
    out = "simulation.pvd"
    pattern = "simulation_step_"
    interval = 1.0
    settle = 0.5
    once = false

    i = 2
    while i <= length(args)
        a = args[i]
        if a == "--out";           out = args[i+1];               i += 2
        elseif a == "--pattern";   pattern = args[i+1];           i += 2
        elseif a == "--interval";  interval = parse(Float64, args[i+1]); i += 2
        elseif a == "--settle";    settle = parse(Float64, args[i+1]);   i += 2
        elseif a == "--once";      once = true;                   i += 1
        else
            println("Unknown option: $a"); return 1
        end
    end

    watch_pvd(dir; pattern=pattern, out=out, interval=interval, settle=settle, once=once)
    return 0
end

end # module PvdWatcher

# Allow running this file directly: `julia pvd_watcher.jl <results_dir> ...`
if abspath(PROGRAM_FILE) == @__FILE__
    PvdWatcher.run_cli(ARGS)
end
