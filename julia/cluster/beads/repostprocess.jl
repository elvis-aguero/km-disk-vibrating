#!/usr/bin/env julia
#= Re-run the post-processing (gain_phase.jl) on already finished sweep directories, without
re-running any simulation. Needed after a change of the quality rules or of the columns.

Usage: julia julia/cluster/beads/repostprocess.jl <run dir> [<run dir> ...]
Each run dir is one that `run_bead_sweep.jl` wrote (it holds run_config.toml and one folder
per bead with cases/*.csv and summary.csv). The old gain_phase*.csv files are kept once as
*.v1.csv (never overwritten if a .v1 copy already exists). =#

if !@isdefined(FaradayDisk)
    include(joinpath(@__DIR__, "..", "..", "scripts", "_bootstrap.jl"))
end
using TOML
include(joinpath(@__DIR__, "run_bead_sweep.jl"))     # defines bead_params(...); does not run main

function repostprocess(outroot::AbstractString)
    cfg = TOML.parsefile(joinpath(outroot, "run_config.toml"))
    basis = String(cfg["basis"])
    scheme = Symbol(cfg["scheme"])
    tres = Int(cfg["temporalResolution"])
    sr = Float64(cfg["spatialResolution"])
    bd = Float64(cfg["bathDiameter"])
    gammas = Float64.(cfg["gammas"])
    freqs = Float64.(cfg["freqs_Hz"])
    domain = String(cfg["domain"])

    all_rows = Dict{String,Any}[]
    for key in String.(cfg["beads"])
        params = bead_params(key, basis, sr, bd, scheme, tres)
        beaddir = joinpath(outroot, key)
        old = joinpath(beaddir, "gain_phase.csv")
        isfile(old) && !isfile(old * ".v1") && cp(old, old * ".v1")
        append!(all_rows, postprocess_sweep(beaddir, params, gammas, freqs; bead = key, basis = basis,
                                             scheme = String(scheme), domain = domain))
    end
    old_all = joinpath(outroot, "gain_phase_all.csv")
    isfile(old_all) && !isfile(old_all * ".v1") && cp(old_all, old_all * ".v1")
    write_gain_phase(old_all, all_rows)

    counts = Dict{String,Int}()
    for r in all_rows
        counts[r["quality"]] = get(counts, r["quality"], 0) + 1
    end
    println(outroot, ": ", length(all_rows), " rows; quality counts ", sort(collect(counts)))
end

for d in ARGS
    try
        repostprocess(d)
    catch err
        @error "repostprocess failed for one directory" dir = d exception = (err, catch_backtrace())
    end
end
