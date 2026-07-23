#!/usr/bin/env julia
#= Computes steady-state amplitude/phase for every case in a Julia sweep output directory,
compares against the digitized experimental measurements, and plots both — the full
functional equivalent of `overlay_validation.py`, reusing the shared `fit_oscillation`
(see `src/convergence.jl`) instead of re-deriving the least-squares regression a third
time.

Plotting uses CairoMakie (a pure rasterizing/vector backend — no display/OpenGL context
needed, so it runs headlessly in CI; see `test/integration/test_plotting.jl`, which
exercises `plot_validation_overlay` directly and is the actual verification this code has
had, since it could not be run locally when written — see julia/README.md). For a truly
*live*, per-step view during a running simulation (MATLAB's per-step `drawnow`), see
`scripts/live_plot.jl` instead, which needs GLMakie and a real display and is NOT covered
by CI.

This file only *defines* functions when `include`d by another script/test; it only runs
`main()` when executed directly (`julia run_validation_overlay.jl <sweep_dir>`).

Block comment, not a docstring: see _bootstrap.jl's header for why a bare top-level string
directly followed by an `if` here would fail to parse. =#

if !@isdefined(FaradayDisk)
    include("_bootstrap.jl")
end
using Printf  # @printf below; not re-exported by `using FaradayDisk` even though FaradayDisk itself depends on it
using TOML    # TOML.parsefile below; same reason
using CairoMakie

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
bath-driven, truncates each case's trace at the first sign of finite-domain wave
reflection (`eta_boundary_cm` crossing 10% of the reference amplitude), requires at least
2 full periods of usable data, fits the last `min(3, n_periods_available)` periods via the
shared `fit_oscillation`, and normalizes amplitude by `A_ref = gamma*g_cgs/omega^2`.

This experimental study is bath-driven (the bath itself physically oscillates; the disk
carries no direct forcing) — confirmed against the real experimental setup, not merely
inferred from data. Forcing-mode detection prefers `sweep_dir/sweep_metadata.toml`'s
`is_bath_driven` flag (written by `run_sweep`, which knows the real physics parameters it
used) over inferring it from `summary.csv`'s `bathAmplitude_cm` column: that column holds
the same `gamma*g/omega^2` quantity regardless of convention (see `write_summary_csv`),
never literally zero either way, so an all-zero check on it can never actually distinguish
bath-driven from disk-forced — confirmed against real sweep data to always misfire
"bath-driven" regardless of which convention was actually used. The column-based
heuristic (matching `overlay_validation.py`) is kept only as a fallback for sweep output
predating `sweep_metadata.toml`.
"""
function compute_sim_overlay(sweep_dir::AbstractString; g_cgs::Float64 = 981.0, phase_rad_bath::Float64 = -pi / 2)
    metadata_path = joinpath(sweep_dir, "sweep_metadata.toml")
    is_bath_driven = true
    if isfile(metadata_path)
        is_bath_driven = Bool(TOML.parsefile(metadata_path)["is_bath_driven"])
    else
        summary_path = joinpath(sweep_dir, "summary.csv")
        if isfile(summary_path)
            s = read_csv_columns(summary_path, ["bathAmplitude_cm"])
            is_bath_driven = !all(iszero, s["bathAmplitude_cm"])
        end
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

# --- Plotting (port of overlay_validation.py's matplotlib figures to CairoMakie) ---

const _EXP_COLORS = Dict(0.05 => RGBf(0.6, 1.0, 0.6), 0.20 => RGBf(0.35, 0.75, 0.35), 0.50 => RGBf(0.1, 0.5, 0.1))
const _SIM_PALETTE = [:royalblue, :darkorange, :seagreen, :crimson, :purple, :saddlebrown, :teal, :magenta]

function _exp_color(g::Real)
    for (k, c) in _EXP_COLORS
        isapprox(k, g; atol = 1e-3) && return c
    end
    return :gray  # unrecognized gamma — still plots, just uncolor-coded
end

function _sim_color(u_gammas::Vector{Float64}, g::Real)
    idx = findfirst(x -> isapprox(x, g; atol = 1e-3), u_gammas)
    return _SIM_PALETTE[mod1(idx === nothing ? 1 : idx, length(_SIM_PALETTE))]
end

"""
    _plot_overlay_panel!(fig_pos, cases, exp_data, value_field, u_gammas, xlims_, ylims_,
                          xticks_, yticks_, title_, ylabel_)

Draws one panel (amplitude or phase): experimental error-bar scatter in the background,
per-gamma simulation `scatterlines!` on top, matching `overlay_validation.py`'s layering
and color scheme (light/medium/dark green for gamma 0.05/0.20/0.50).
"""
function _plot_overlay_panel!(fig_pos, cases::Vector{OverlayCase}, exp_data::Dict{String,Vector{Float64}},
                               value_field::Symbol, u_gammas::Vector{Float64}, xlims_::Tuple, ylims_::Tuple,
                               xticks_, yticks_, title_::AbstractString, ylabel_::AbstractString)
    ax = Axis(fig_pos; xlabel = "f (Hz)", ylabel = ylabel_, title = title_)
    xlims!(ax, xlims_...)
    ylims!(ax, ylims_...)
    ax.xticks = xticks_
    ax.yticks = yticks_

    for g in sort(unique(exp_data["gamma"]))
        idx = findall(x -> isapprox(x, g; atol = 1e-3), exp_data["gamma"])
        isempty(idx) && continue
        errorbars!(ax, exp_data["frequency_Hz"][idx], exp_data["value"][idx],
                   exp_data["error_lower"][idx], exp_data["error_upper"][idx];
                   color = :black, whiskerwidth = 6)
        scatter!(ax, exp_data["frequency_Hz"][idx], exp_data["value"][idx];
                 color = _exp_color(g), strokecolor = :black, strokewidth = 1, markersize = 12,
                 label = "Exp Γ=$(round(g, digits = 2))")
    end

    for g in u_gammas
        gcases = sort(filter(c -> isapprox(c.gamma, g; atol = 1e-3), cases), by = c -> c.freq)
        isempty(gcases) && continue
        freqs = [c.freq for c in gcases]
        values = [getfield(c, value_field) for c in gcases]
        scatterlines!(ax, freqs, values; color = _sim_color(u_gammas, g), linewidth = 2, markersize = 8,
                      label = "Sim Γ=$(round(g, digits = 2))")
    end

    axislegend(ax; position = :rt, framevisible = true, labelsize = 11)
    return ax
end

"""
    plot_validation_overlay(cases, exp_amp_path, exp_phase_path, out_dir, sweep_name) -> NamedTuple

Writes `val_amp_<sweep_name>.png/.pdf` and `val_phase_<sweep_name>.png/.pdf` under
`out_dir`, matching `overlay_validation.py`'s output naming and figure layout (axis
limits/ticks, title, legend). Returns the four output paths.
"""
function plot_validation_overlay(cases::Vector{OverlayCase}, exp_amp_path::AbstractString,
                                  exp_phase_path::AbstractString, out_dir::AbstractString, sweep_name::AbstractString)
    isempty(cases) && error("plot_validation_overlay: no cases to plot")
    mkpath(out_dir)
    u_gammas = sort(unique(c.gamma for c in cases))

    exp_amp = read_csv_columns(exp_amp_path, ["frequency_Hz", "value", "error_lower", "error_upper", "gamma"])
    exp_phase = read_csv_columns(exp_phase_path, ["frequency_Hz", "value", "error_lower", "error_upper", "gamma"])

    fig_amp = Figure()
    _plot_overlay_panel!(fig_amp[1, 1], cases, exp_amp, :ampNorm, u_gammas, (0, 105), (0, 1.2), 0:10:100, 0:0.2:1.2,
                         "Amplitude Validation Overlay ($sweep_name)", "A_disk / A_base")
    amp_png = joinpath(out_dir, "val_amp_$(sweep_name).png")
    amp_pdf = joinpath(out_dir, "val_amp_$(sweep_name).pdf")
    save(amp_png, fig_amp; px_per_unit = 3)
    save(amp_pdf, fig_amp)

    fig_phase = Figure()
    _plot_overlay_panel!(fig_phase[1, 1], cases, exp_phase, :phaseDiff, u_gammas, (0, 105), (0, 30), 0:10:100, 0:5:30,
                         "Phase Validation Overlay ($sweep_name)", "Phase Difference (deg)")
    phase_png = joinpath(out_dir, "val_phase_$(sweep_name).png")
    phase_pdf = joinpath(out_dir, "val_phase_$(sweep_name).pdf")
    save(phase_png, fig_phase; px_per_unit = 3)
    save(phase_pdf, fig_phase)

    return (amp_png = amp_png, amp_pdf = amp_pdf, phase_png = phase_png, phase_pdf = phase_pdf)
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

    sweep_name = basename(rstrip(sweep_dir, '/'))
    plots = plot_validation_overlay(cases, joinpath(EXPERIMENTAL_DIR, "ampiezza_solo_misure_3.csv"),
                                     joinpath(EXPERIMENTAL_DIR, "fase_solo_misure_3.csv"), sweep_dir, sweep_name)
    println("Amplitude plot: ", plots.amp_png, " / ", plots.amp_pdf)
    println("Phase plot: ", plots.phase_png, " / ", plots.phase_pdf)

    return (stats..., plots...)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
