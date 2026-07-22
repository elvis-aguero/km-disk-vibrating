#!/usr/bin/env julia
#= Tier-4 validation: runs the full 30-case sweep (the same grid as `default_sweep_spec()`
in `run_sweep.jl`/`sweeper.m`) with the *corrected* convergence check, computes
amplitude/phase RMSE against the digitized experimental data, and checks the result
against thresholds anchored to the known-good pre-regression MATLAB baseline
(`sweep_original` / `sweep_20260401_141625`: ~0.10 amplitude RMSE, ~5.7 deg phase RMSE)
rather than the buggy `sweep_20260721_211738` run audited earlier in this project
(~1.09 amplitude RMSE, ~285 deg phase RMSE).

This is deliberately NOT part of the default (fast) test tier — the full grid is a
genuinely multi-hour-scale computation even multithreaded (each case in the original
MATLAB sweeps took ~200-1500s). Run manually, or from
`test/slow/test_matlab_baseline_validation.jl` under `KMDISK_RUN_SLOW_TESTS=1`. This
script and that test share `run_matlab_baseline_check()` so there is exactly one
implementation of "run sweep, fit, compare, compute RMSE" to keep correct.

Usage: `julia -t auto julia/scripts/run_matlab_baseline_check.jl`

Block comment, not a docstring: see _bootstrap.jl's header for why a bare top-level string
directly followed by an `if` here would fail to parse. =#

if !@isdefined(FaradayDisk)
    include("_bootstrap.jl")
end
using Dates  # Dates.now()/format below; not re-exported by `using FaradayDisk`
using Logging  # global_logger() below; not re-exported by `using FaradayDisk`
include("run_sweep.jl")
include("run_validation_overlay.jl")

const AMPLITUDE_RMSE_THRESHOLD = 0.10 * 1.5   # 50% slack over the known-good MATLAB baseline (~0.10)
const PHASE_RMSE_DEG_THRESHOLD = 5.7 * 1.5     # 50% slack over the known-good MATLAB baseline (~5.7 deg)
const AMPLITUDE_RMSE_MUST_BEAT = 1.09          # the buggy sweep's amplitude RMSE — must not reproduce this
const PHASE_RMSE_DEG_MUST_BEAT = 100.0         # the buggy sweep's phase RMSE was ~285 deg — must not be anywhere near it

"""
    run_matlab_baseline_check(; outdir=nothing) -> NamedTuple

Runs the full sweep, computes the overlay RMSE against experiment, and returns
`(amp_rmse, phase_rmse_deg, n_amp, n_phase, outdir, pass_primary, pass_incident_guard)`.
Does not throw on threshold failure — the caller (script `main()` or the Tier-4 test)
decides what to do with the result.
"""
function run_matlab_baseline_check(; outdir::Union{Nothing,String} = nothing)
    if outdir === nothing
        ts = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
        outdir = joinpath(@__DIR__, "..", "data", "results", "sweeps", "baseline_check_$(ts)")
    end

    global_logger(setup_logging(name = "matlab_baseline_check", jsonlines = true))

    run_sweep(default_sweep_spec(), outdir)
    cases = compute_sim_overlay(outdir)
    isempty(cases) && error("baseline check produced no valid steady-state cases — cannot compute RMSE")

    stats = rmse_against_experiment(cases, joinpath(EXPERIMENTAL_DIR, "ampiezza_solo_misure_3.csv"),
                                     joinpath(EXPERIMENTAL_DIR, "fase_solo_misure_3.csv"))

    sweep_name = basename(rstrip(outdir, '/'))
    plots = plot_validation_overlay(cases, joinpath(EXPERIMENTAL_DIR, "ampiezza_solo_misure_3.csv"),
                                     joinpath(EXPERIMENTAL_DIR, "fase_solo_misure_3.csv"), outdir, sweep_name)

    pass_primary = stats.amp_rmse <= AMPLITUDE_RMSE_THRESHOLD && stats.phase_rmse_deg <= PHASE_RMSE_DEG_THRESHOLD
    pass_incident_guard = stats.amp_rmse < AMPLITUDE_RMSE_MUST_BEAT && stats.phase_rmse_deg < PHASE_RMSE_DEG_MUST_BEAT

    @info "baseline check complete" amp_rmse = stats.amp_rmse phase_rmse_deg = stats.phase_rmse_deg pass_primary = pass_primary pass_incident_guard = pass_incident_guard outdir = outdir amp_plot = plots.amp_png phase_plot = plots.phase_png

    return (amp_rmse = stats.amp_rmse, phase_rmse_deg = stats.phase_rmse_deg, n_amp = stats.n_amp,
            n_phase = stats.n_phase, outdir = outdir, pass_primary = pass_primary, pass_incident_guard = pass_incident_guard,
            amp_plot = plots.amp_png, phase_plot = plots.phase_png)
end

"""
    main()

CLI entry point. Takes one optional positional argument, the output directory (e.g. for a
cluster job to point results at a job-specific `output/<job_id>/sweep` path instead of the
default auto-timestamped `julia/data/results/sweeps/baseline_check_<ts>/`); with no
argument, behaves exactly as before.
"""
function main()
    outdir = isempty(ARGS) ? nothing : ARGS[1]
    result = run_matlab_baseline_check(; outdir = outdir)
    println("Amplitude RMSE: ", result.amp_rmse, " (threshold ", AMPLITUDE_RMSE_THRESHOLD, ")")
    println("Phase RMSE (deg): ", result.phase_rmse_deg, " (threshold ", PHASE_RMSE_DEG_THRESHOLD, ")")
    println("Passes primary threshold: ", result.pass_primary)
    println("Passes incident guard (not reproducing the ~1.09/~285deg regression): ", result.pass_incident_guard)
    if !result.pass_primary
        @error "baseline check FAILED primary threshold"
    end
    return result
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
