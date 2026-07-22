#!/usr/bin/env julia
"""
Opt-in live visualization of a running simulation, via GLMakie.

STATUS / CAVEAT: this is the single least-verified file in the whole port. Every other
module and script here at least got real feedback from GitHub Actions CI once pushed; this
one cannot be, since GLMakie needs a real display/GPU context (unlike CairoMakie, which is
headless-safe and is exercised by `test/integration/test_plotting.jl`) and so is
deliberately excluded from `julia-ci.yml` entirely. It has also never been run locally even
once during development — no working Julia install was available while writing this port
(see the repo-level README section on how this port was authored). Treat every line below
as unverified until someone with local Julia and a display runs it and reports back.

Usage:

    using FaradayDisk
    include("scripts/live_plot.jl")

    dtn = ...  # e.g. via generate_dtn or load_dtn
    problem = build_problem(params, dtn)
    cb = make_live_plot_callback(problem)
    run_simulation(params; dtn = dtn, on_step = cb)

Never pass a callback built by this file as the `dtn_registry`/sweep path's `on_step` —
`run_sweep`/`run_sweep_case` run cases concurrently via `Threads.@threads`, and opening or
mutating a GLMakie window from more than one thread at once is not safe. This is exactly
why `on_step` defaults to `nothing` and is never set by the sweep orchestrator itself (see
`simulate.jl`'s docstring).
"""

if !@isdefined(FaradayDisk)
    include("_bootstrap.jl")
end
using GLMakie

"""
    make_live_plot_callback(problem; width_disk_radii=6, refresh_every=1) -> Function

Builds a closure `(state, problem, step_count) -> nothing` suitable for `run_simulation`'s
`on_step` keyword. Opens one GLMakie window immediately (via `display`) and, every
`refresh_every`-th call thereafter, updates it in place by mutating `Observable`s rather
than replotting from scratch: the bath surface profile `eta(r)`, the disk drawn as a short
horizontal segment at its current height spanning its own radius, and horizontal guide
lines at `+-A_bath_adim` giving a visual sense of the drive amplitude's scale. This mirrors
the debug plot embedded directly in `solve_motion.m`'s stepping loop (its `drawnow` block).

`A_bath_adim`/`A_forcing_adim` reuse fields already present on `Problem` — no new fields
were needed on the struct for this. `A_bath_adim = problem.bath_forcing_amplitude /
(problem.Fr * problem.bath_frequency^2)` is algebraically the nondimensional bath-drive
amplitude (inverting how `bath_forcing_amplitude = A*omega_bath^2/g` and `Fr = L/(g*T^2)`
were built from the physical `bathAmplitude` in `parameters.jl`), matching
`solve_motion.m`'s own `A_bath_adim = abs(bathAmplitude/L_unit)` debug-plot line. Likewise
`A_forcing_adim = problem.force_amplitude / (problem.force_frequency^2 + 1e-6)` matches
MATLAB's `force_adim/(freq_adim^2+1e-6)` (the `1e-6` guards the zero-forcing-amplitude
case, preserved as-is from the source rather than special-cased away).

`width_disk_radii` sets how many disk radii of the bath (from `r=0` outward) are shown —
the full domain (`problem.nr*problem.dr`, e.g. 50 disk radii at the default
`bathDiameter=100`) is far wider than the interesting near-disk region and would make the
disk itself and the near-field waves invisible at a normal window size. `refresh_every`
skips most redraws (default: every step) since a `display`-driven GLMakie redraw is
wall-clock-expensive relative to a single implicit solve of the step system.
"""
function make_live_plot_callback(problem::Problem; width_disk_radii::Real = 6, refresh_every::Integer = 1)
    refresh_every >= 1 || throw(ArgumentError("refresh_every must be >= 1, got $refresh_every"))

    nr = problem.nr
    dr = problem.dr
    r = collect(0:(nr - 1)) .* dr
    n_show = clamp(ceil(Int, width_disk_radii / dr), 2, nr)
    r_show = r[1:n_show]

    A_bath_adim = problem.bath_forcing_amplitude / (problem.Fr * problem.bath_frequency^2)
    A_forcing_adim = problem.force_amplitude / (problem.force_frequency^2 + 1e-6)
    H_limit_adim = 3 * max(A_bath_adim, A_forcing_adim, 1e-6)

    fig = Figure(size = (900, 600))
    ax = Axis(fig[1, 1]; xlabel = "r (disk radii)", ylabel = "eta, z (dimensionless)",
              title = "live simulation state (t = 0.000)")
    xlims!(ax, 0, r_show[end])
    ylims!(ax, -H_limit_adim, H_limit_adim)

    eta_obs = Observable(zeros(n_show))
    disk_y_obs = Observable([0.0, 0.0])

    lines!(ax, r_show, eta_obs; color = :royalblue, label = "bath surface eta")
    lines!(ax, [-1.0, 1.0], disk_y_obs; color = :crimson, linewidth = 6, label = "disk")
    hlines!(ax, [A_bath_adim, -A_bath_adim]; color = (:seagreen, 0.4), linestyle = :dash, label = "bath drive amplitude")
    axislegend(ax)

    display(fig)

    calls = 0
    return function live_plot_callback(state::SimulationState, problem::Problem, step_count::Int)
        calls += 1
        calls % refresh_every == 0 || return nothing

        eta_obs[] = state.bath_surface[1:n_show]
        disk_y_obs[] = [state.center_of_mass, state.center_of_mass]
        ax.title[] = "live simulation state (t = $(round(state.time, digits = 3)))"

        return nothing
    end
end
