#!/usr/bin/env julia
#= Dumps Julia arrays for the MATLAB parity tests (matlab/1_code/tests/DtnParityTest.m), one text file per array in the
directory given as ARGS[1]: whitespace-separated Float64 (shortest round-trip repr), one matrix row per line.
Stage 1: dtn_sl_nr50_D20_o4.txt, dtn_sl_nr50_D20_o6.txt (generate_dtn_single_layer(50, 20; order)).
Block comment, not a docstring: see _bootstrap.jl. Project choice as in parity_reference.jl. =#
import Pkg
let proj = get(ENV, "KM_JULIA_PROJECT", get(ENV, "GITHUB_ACTIONS", "") == "true" ? joinpath(@__DIR__, "..") : @__DIR__)
    Pkg.activate(proj)
    if !haskey(Pkg.project().dependencies, "FaradayDisk")
        Pkg.develop(path = joinpath(@__DIR__, ".."))
    end
    abspath(proj) == abspath(@__DIR__) && Pkg.instantiate()
end
using FaradayDisk
using Logging
global_logger(Logging.NullLogger())

outdir = isempty(ARGS) ? error("usage: parity_dump.jl <outdir>") : ARGS[1]
mkpath(outdir)

function dump_matrix(name, A)
    open(joinpath(outdir, name), "w") do io
        for i in 1:size(A, 1)
            println(io, join(A[i, :], " "))
        end
    end
end

for order in (4, 6)
    dump_matrix("dtn_sl_nr50_D20_o$(order).txt", generate_dtn_single_layer(50, 20.0; order = order))
end
println("done")
