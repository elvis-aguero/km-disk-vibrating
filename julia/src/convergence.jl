"""
Shared least-squares oscillation fit and the corrected early-stop convergence check.

## Why the MATLAB check was replaced, not patched

`matlab/1_code/simulation_code/solve_motion.m` (lines ~308-337) stopped a simulation once
the standard deviation of a rolling 1-period average of the raw center-of-mass trace,
divided by its mean, dropped below a tolerance. This was audited and found to be broken:
it measures whether the *smoothed baseline* has stopped drifting, which is a different
condition from "the oscillation has reached periodic steady state." Because the
static-equilibrium initial condition removes most of the transient drift up front, the
ratio's denominator (a fairly large, non-oscillating equilibrium offset) makes the ratio
trivially small almost immediately — verified directly against a real sweep run, where
all 30 (gamma, frequency) cases stopped at exactly the earliest allowed check point (3
forcing periods), before the oscillation had even completed one full cycle, degrading the
amplitude match to experimental data by roughly 10x and the phase match by roughly 50x.

The check below instead fits oscillation amplitude and phase (via the same least-squares
regression already used by the validation/overlay tooling) over two consecutive
multi-period windows, and only declares convergence when *both* are stable across those
windows — i.e., the state is actually periodic, not merely non-drifting. A hard
minimum-periods floor (independent of the tolerance check) prevents firing before the
system has had a chance to complete enough cycles to fit a meaningful window at all.
"""

Base.@kwdef struct ConvergenceOptions
    enabled::Bool = true
    minPeriods::Float64 = 10.0        # never declare convergence before this many periods
    windowPeriods::Float64 = 5.0      # length of each of the two compared windows, in periods
    checkEveryPeriods::Float64 = 1.0  # cadence of checks, in periods
    amplitudeRelTol::Float64 = 0.02   # relative amplitude change tolerated between windows
    phaseAbsTolDeg::Float64 = 2.0     # absolute phase change (degrees) tolerated between windows
end

struct ConvergenceResult
    converged::Bool
    periods_elapsed::Float64
    amplitude::Float64
    phase_deg::Float64
    amplitude_rel_change::Float64
    phase_change_deg::Float64
end

not_converged(periods_elapsed::Real) = ConvergenceResult(false, periods_elapsed, NaN, NaN, Inf, Inf)

"""
    fit_oscillation(t, z, ω) -> (amplitude, phase, offset)

Least-squares fit of `z(t) ≈ amplitude*cos(ω*t + phase) + offset` via the design matrix
`[cos(ωt) sin(ωt) 1]` (note the `+`: `X = [amplitude*cos(phase), -amplitude*sin(phase),
offset]`, so `phase = atan2(-X1,X0)` recovers the `+phase` convention below, not `-phase`
— confirmed against a synthetic-signal test after the first real test run caught the sign
flipped in this docstring and in `test/unit/test_fit.jl`). Ported from the fit already
used in `overlay_validation.py` (`M = [cos(ωt), sin(ωt), 1]`, `amplitude = hypot(X0,X1)`,
`phase = atan2(-X1,X0)`). The constant column absorbs any DC offset (e.g. the
static-equilibrium depth), so the recovered amplitude/phase are insensitive to that
offset — precisely the property the old rolling-mean check lacked.
"""
function fit_oscillation(t::AbstractVector{<:Real}, z::AbstractVector{<:Real}, ω::Real)
    length(t) == length(z) || throw(ArgumentError("t and z must have the same length"))
    length(t) >= 3 || throw(ArgumentError("fit_oscillation needs at least 3 samples, got $(length(t))"))
    M = [cos.(ω .* t) sin.(ω .* t) ones(length(t))]
    X = M \ collect(z)
    amplitude = hypot(X[1], X[2])
    phase = atan(-X[2], X[1])
    return (amplitude = amplitude, phase = phase, offset = X[3])
end

"""
    check_convergence(t, z, ω, opts) -> ConvergenceResult

Fits amplitude/phase over two consecutive, non-overlapping `opts.windowPeriods`-period
windows ending at `t[end]`, and declares convergence only when both are stable across
those windows within tolerance. Never returns `converged=true` before
`opts.minPeriods` periods have elapsed, regardless of tolerances.
"""
function check_convergence(t::AbstractVector{<:Real}, z::AbstractVector{<:Real}, ω::Real, opts::ConvergenceOptions)
    T = 2 * pi / ω
    periods_elapsed = t[end] / T

    if !opts.enabled || periods_elapsed < opts.minPeriods || periods_elapsed < 2 * opts.windowPeriods
        return not_converged(periods_elapsed)
    end

    t_end = t[end]
    in_window_b = t .>= (t_end - opts.windowPeriods * T)
    in_window_a = (t .>= (t_end - 2 * opts.windowPeriods * T)) .& .!in_window_b

    (count(in_window_a) >= 3 && count(in_window_b) >= 3) || return not_converged(periods_elapsed)

    fit_a = fit_oscillation(view(t, in_window_a), view(z, in_window_a), ω)
    fit_b = fit_oscillation(view(t, in_window_b), view(z, in_window_b), ω)

    amp_rel_change = abs(fit_b.amplitude - fit_a.amplitude) / max(fit_a.amplitude, eps())
    phase_change_deg = abs(mod(rad2deg(fit_b.phase - fit_a.phase) + 180, 360) - 180)

    converged = amp_rel_change < opts.amplitudeRelTol && phase_change_deg < opts.phaseAbsTolDeg

    return ConvergenceResult(converged, periods_elapsed, fit_b.amplitude, rad2deg(fit_b.phase),
                              amp_rel_change, phase_change_deg)
end
