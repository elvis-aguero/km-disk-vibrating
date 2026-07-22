#!/usr/bin/env julia
#= Runs the (gamma, frequency) parameter sweep — port of `matlab/1_code/sweeper.m` — using
the same fixed parameters and Cartesian grid (gamma in {0.05, 0.2, 0.5}, frequency
10:10:100 Hz, 30 cases). Output goes to `julia/data/results/sweeps/<timestamp>/`.

Usage: `julia -t auto julia/scripts/run_sweep.jl`

Only runs `main()` when executed directly — `default_sweep_spec()` is reused by
`run_matlab_baseline_check.jl` and `test/slow/test_matlab_baseline_validation.jl` via
`include`, so the exact sweep grid/parameters are defined in exactly one place.

Block comment, not a docstring: see _bootstrap.jl's header for why a bare top-level string
directly followed by an `if` here would fail to parse. =#

if !@isdefined(FaradayDisk)
    include("_bootstrap.jl")
end
using Dates  # Dates.now()/format below; not re-exported by `using FaradayDisk`

"""
    default_sweep_spec() -> SweepSpec

The same fixed parameters and Cartesian grid as `matlab/1_code/sweeper.m`'s `fixed`
struct and `sweep.gamma` / `sweep.bathFrequency` axes.
"""
function default_sweep_spec()
    fixed = SimulationParams(
        diskRadius = 0.2, diskMass = 0.0283,
        forceAmplitude = 0.0, forceFrequency = 90.0,
        phaseDifference = -90.0,
        bathDensity = 1.175, bathSurfaceTension = 66.5, bathViscosity = 0.18 / 1.175,
        bathDiameter = 100.0, spatialResolution = 50.0, temporalResolution = 60,
        simulationTime = 30 / 90,
        solverType = :auto, startStatic = true, earlyStop = true,
    )
    return SweepSpec(gammas = [0.05, 0.2, 0.5], bathFrequenciesHz = collect(10.0:10.0:100.0), fixed = fixed)
end

function main()
    ts = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
    outdir = joinpath(@__DIR__, "..", "data", "results", "sweeps", "sweep_$(ts)")

    global_logger(setup_logging(name = "sweep", jsonlines = true))

    results = run_sweep(default_sweep_spec(), outdir)
    println("Sweep complete. Results in: ", outdir)
    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
