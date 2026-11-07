#= One bath-driven run at nr = 2500 (sr 25, D 200) with the chosen solver; prints setup/step time and peak RSS, and
writes the CoM trace for the comparison. Usage: julia -t N julia/scripts/bench_schur.jl <lu|schur> <out.csv> =#
include(joinpath(@__DIR__, "_bootstrap.jl"))
using Printf
solver = Symbol(ARGS[1]); out = ARGS[2]
t_dtn = @elapsed dtn = FaradayDisk.generate_dtn(2500, 200.0)
p = SimulationParams(bathDiameter = 200.0, spatialResolution = 25.0, temporalResolution = 60, timeScheme = :bdf2,
                     forceAmplitude = 0.0, bathAmplitude = 0.5 * 981 / (2pi * 30)^2, bathFrequency = 30.0,
                     phaseDifference = -90.0, simulationTime = 5 / 30, startStatic = true, earlyStop = false,
                     solverType = solver, diskRadius = 0.07938)
t_run = @elapsed r = run_simulation(p; dtn = dtn)
steps = length(r.CoM_cm) - 1
@printf("solver=%s  dtn %.1f s  run %.1f s (%d steps)  peak RSS %.1f GB  hit_rate %.3f  status %s\n",
        solver, t_dtn, t_run, steps, Sys.maxrss() / 1e9, r.lu_cache_hit_rate, r.status)
open(out, "w") do io
    for (t, c) in zip(r.t_s, r.CoM_cm); @printf(io, "%.12e,%.17e\n", t, c); end
end
