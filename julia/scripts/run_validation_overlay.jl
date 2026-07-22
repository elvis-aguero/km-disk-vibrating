#!/usr/bin/env julia
"""
Computes steady-state amplitude/phase for every case in a Julia sweep output directory
and compares against the digitized experimental measurements — the numerical core of
`overlay_validation.py`, ported to reuse the shared `fit_oscillation` (see
`src/convergence.jl`) instead of re-deriving the least-squares regression a third time.

Deliberately does NOT reproduce the PNG/PDF plotting from `overlay_validation.py` — doing
so would mean adding a plotting dependency (Plots.jl/Makie) to this port without being
able to install or run it to verify it actually works (see the network-access caveat in
`src/solver/gmres_solver.jl`). Instead, this writes a plain CSV table of
(gamma, freq, ampNorm, phaseDiff) per case plus the summary RMSE numbers, which is
sufficient for both human inspection and the automated MATLAB-baseline check in
`run_matlab_baseline_check.jl` / `test/slow/test_matlab_baseline_validation.jl` (which
`include`s this file rather than duplicating its logic). Adding plot output is a natural
follow-up once Julia is available to verify a plotting library end-to-end.

This file only *defines* functions when `include`d by another script/test; it only runs
`main()` when executed directly (`julia run_validation_overlay.jl <sweep_dir>`).
"""

if !@isdefined(FaradayDisk)
    include("_bootstrap.jl")
end
using Printf  # @printf below; not re-exported by `using FaradayDisk` even though FaradayDisk itself depends on it

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const EXPERIMENTAL_DIR = joinpath(REPO_ROOT, "matlab", "0_data", "external", "experimental_measurements")

"""
    read_csv_columns(path, wanted) -> Dict{String,Vector{Float64}}

Minimal, dependency-free CSV reader for the simple numeric CSVs this codebase produces
and consumes (comma-separated, one header row, no quoting) — reads only the named
`wanted` columns (by header lookup), ignoring any other (possibly non-numeric) columns
in the file, such as the experimental CSVs' `marker_face_color` column.
"""
function read_csv_columns(path::AbstractString, wanted::Vector{String})
    lines = readlines(path)
    isempty(lines) && error("empty CSV: $path")
    header = String.(split(lines[1], ','))
    idx = Dict(h => i for (i, h) in enumerate(header))
    for w in wanted
        haskey(idx, w) || error("column \"$w\" not found in $path (header: $header)")
    end
    data = Dict(w => Float64[] for w in wanted)
    for line in @view lines[2:end]
        isempty(strip(line)) && continue
        parts = split(line, ',')
        for w in wanted
            push!(data[w], parse(Float64, parts[idx[w]]))
        end
    end
    return data
end

struct OverlayCase
    gamma::Float64
    freq::Float64
    ampNorm::Float64
    phaseDiff::Float64
end

"""
    compute_sim_overlay(sweep_dir; g_cgs=981.0, phase_rad_bath=-pi/2) -> Vector{OverlayCase}

Direct port of `overlay_validation.py`'s per-case fit: detects forcing-driven vs
bath-driven from `summary.csv`, truncates each case's trace at the first sign of
finite-domain wave reflection (`eta_boundary_cm` crossing 10% of the reference
amplitude), requires at least 2 full periods of usable data, fits the last
`min(3, n_periods_available)` periods via the shared `fit_oscillation`, and normalizes
amplitude by `A_ref = gamma*g_cgs/omega^2`.
"""
function compute_sim_overlay(sweep_dir::AbstractString; g_cgs::Float64 = 981.0, phase_rad_bath::Float64 = -pi / 2)
    summary_path = joinpath(sweep_dir, "summary.csv")
    is_bath_driven = true
    if isfile(summary_path)
        s = read_csv_columns(summary_path, ["bathAmplitude_cm"])
        is_bath_driven = !all(iszero, s["bathAmplitude_cm"])
    end

    cases_dir = joinpath(sweep_dir, "cases")
    isdir(cases_dir) || error("no cases/ directory found under $sweep_dir")

    results = OverlayCase[]
    for fname in readdir(cases_dir)
        m = match(r"^gamma([\d.]+)_f([\d.]+)Hz\.csv$", fname)
        m === nothing && continue
        g_val = parse(Float64, m.captures[1])
        f_val = parse(Float64, m.captures[2])

        data = read_csv_columns(joinpath(cases_dir, fname), ["time_s", "CoM_cm", "eta_boundary_cm"])
        t = data["time_s"]
        CoM = data["CoM_cm"]
        eta_boundary = data["eta_boundary_cm"]
        length(t) < 10 && continue

        omega = 2 * pi * f_val
        A_ref = g_val * g_cgs / omega^2
        T_period = 1 / f_val

        wave_idx = findfirst(x -> abs(x) > 0.10 * A_ref, eta_boundary)
        if wave_idx !== nothing
            idx_cutoff = max(1, wave_idx - 1)
            t = t[1:idx_cutoff]
            CoM = CoM[1:idx_cutoff]
        end
        isempty(t) && continue

        total_time = t[end] - t[1]
        n_periods_available = floor(Int, total_time / T_period)
        n_periods_available < 2 && continue

        n_eval = min(3, n_periods_available)
        t_start_eval = t[end] - n_eval * T_period
        eval_idx = t .>= t_start_eval
        t_eval = t[eval_idx]
        CoM_eval = CoM[eval_idx]

        if is_bath_driven
            z_lab = CoM_eval .+ A_ref .* cos.(omega .* t_eval .+ phase_rad_bath)
            phase_ref = phase_rad_bath
        else
            z_lab = CoM_eval
            phase_ref = 0.0
        end

        fit = fit_oscillation(t_eval, z_lab, omega)
        dphi = mod(rad2deg(phase_ref - fit.phase), 360.0)
        push!(results, OverlayCase(g_val, f_val, fit.amplitude / A_ref, dphi))
    end
    return results
end

"""
    rmse_against_experiment(sim, exp_amp_path, exp_phase_path) -> (amp_rmse, phase_rmse_deg)

Matches each simulated case to the experimental point at the same `(gamma, frequency)`
and computes the RMSE of `ampNorm` and `phaseDiff` against the experimental `value`
columns.
"""
function rmse_against_experiment(sim::Vector{OverlayCase}, exp_amp_path::AbstractString, exp_phase_path::AbstractString)
    exp_amp = read_csv_columns(exp_amp_path, ["frequency_Hz", "value", "gamma"])
    exp_phase = read_csv_columns(exp_phase_path, ["frequency_Hz", "value", "gamma"])

    sq_amp = Float64[]
    sq_phase = Float64[]
    for s in sim
        for i in eachindex(exp_amp["gamma"])
            if isapprox(exp_amp["gamma"][i], s.gamma; atol = 1e-3) && isapprox(exp_amp["frequency_Hz"][i], s.freq; atol = 1e-6)
                push!(sq_amp, (exp_amp["value"][i] - s.ampNorm)^2)
            end
        end
        for i in eachindex(exp_phase["gamma"])
            if isapprox(exp_phase["gamma"][i], s.gamma; atol = 1e-3) && isapprox(exp_phase["frequency_Hz"][i], s.freq; atol = 1e-6)
                push!(sq_phase, (exp_phase["value"][i] - s.phaseDiff)^2)
            end
        end
    end

    amp_rmse = isempty(sq_amp) ? NaN : sqrt(sum(sq_amp) / length(sq_amp))
    phase_rmse = isempty(sq_phase) ? NaN : sqrt(sum(sq_phase) / length(sq_phase))
    return (amp_rmse = amp_rmse, phase_rmse_deg = phase_rmse, n_amp = length(sq_amp), n_phase = length(sq_phase))
end

function write_overlay_csv(path::AbstractString, cases::Vector{OverlayCase})
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "gamma,freq,ampNorm,phaseDiff")
        for c in sort(cases, by = c -> (c.gamma, c.freq))
            @printf(io, "%.4f,%g,%.6f,%.4f\n", c.gamma, c.freq, c.ampNorm, c.phaseDiff)
        end
    end
    return path
end

function main(args)
    length(args) >= 1 || error("usage: run_validation_overlay.jl <sweep_dir>")
    sweep_dir = args[1]
    cases = compute_sim_overlay(sweep_dir)
    isempty(cases) && error("no valid steady-state cases found in $sweep_dir")

    out_csv = joinpath(sweep_dir, "validation_overlay.csv")
    write_overlay_csv(out_csv, cases)

    stats = rmse_against_experiment(cases, joinpath(EXPERIMENTAL_DIR, "ampiezza_solo_misure_3.csv"),
                                     joinpath(EXPERIMENTAL_DIR, "fase_solo_misure_3.csv"))
    println("Validation overlay written to: ", out_csv)
    println("Amplitude RMSE: ", stats.amp_rmse, " (n=", stats.n_amp, ")")
    println("Phase RMSE (deg): ", stats.phase_rmse_deg, " (n=", stats.n_phase, ")")
    return stats
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
