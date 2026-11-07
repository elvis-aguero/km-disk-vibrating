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

# Stage 3 (graded grid). Bead fluid, 40 Hz, graded_grid_for defaults (ppw 32, q 1.04, r_in 2, ramp 20, smooth 0.1, order 2, wall_zone 2).
dump_vec(name, v) = dump_matrix(name, reshape(collect(Float64, v), :, 1))
pbead = SimulationParams(bathDensity = 1.148, bathSurfaceTension = 68.0, bathViscosity = 0.079, g = 981.0, diskRadius = 0.07938,
                         bathDiameter = 700.0, spatialResolution = 25.0)
for (tag, g) in ("bead40" => graded_grid_for(pbead, 40.0), "bead40_wz0" => graded_grid_for(pbead, 40.0; wall_zone = 0.0))
    dump_vec("grid_$(tag)_r.txt", g.r)
    dump_vec("grid_$(tag)_s.txt", g.s)
    dump_vec("grid_$(tag)_sp.txt", g.sp)
end
# smaller domain (R = 100 disk radii) for the Laplacian and the DtN matrix
psmall = with_overrides(pbead; bathDiameter = 200.0)
gs = graded_grid_for(psmall, 40.0)
dump_vec("grid_small_r.txt", gs.r)
dump_matrix("grid_small_laplacian.txt", Matrix(build_domain(gs).Laplacian))
dump_matrix("grid_small_dtn.txt", generate_dtn(gs))
println("done")
