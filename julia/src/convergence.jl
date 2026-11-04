"""
Least-squares oscillation fit and early-stop convergence check.

Convergence = amplitude and phase fitted over two consecutive multi-period windows agree within tolerance,
with a hard floor of `minPeriods`. (Replaces the MATLAB check, which tested drift of a smoothed baseline and
stopped every case after 3 periods.)
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

Least squares z ≈ amplitude cos(ωt + phase) + offset, using the basis [cos(ωt), sin(ωt), 1] = X;
amplitude = hypot(X1, X2), phase = atan(-X2, X1).
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

Fits amplitude and phase over two consecutive non-overlapping `opts.windowPeriods` windows ending at `t[end]`;
converged only if both change less than the tolerances. Never converged before `opts.minPeriods` periods.
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
