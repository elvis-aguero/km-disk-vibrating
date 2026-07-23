function result = check_convergence(t, z, omega, opts)
%CHECK_CONVERGENCE Amplitude/phase-stability early-stop check.
%   Mirrors julia/src/convergence.jl's check_convergence, which replaced (not
%   patched) the old rolling-mean-ratio check: that check measured whether the
%   *smoothed baseline* had stopped drifting, a different condition from "the
%   oscillation has reached periodic steady state." Because the
%   static-equilibrium initial condition removes most of the transient drift up
%   front, its ratio was trivially small almost immediately, regardless of
%   whether the oscillation itself had built up -- confirmed against real sweep
%   data, where weak-forcing cases were still drifting by tens of percent in
%   amplitude and tens of degrees in phase at the point the old check declared
%   convergence.
%
%   This check instead fits oscillation amplitude and phase (via
%   FIT_OSCILLATION, the same least-squares regression used by the
%   validation/overlay tooling) over two consecutive, non-overlapping
%   opts.windowPeriods-period windows ending at t(end), and only declares
%   convergence when *both* are stable across those windows -- i.e. the state
%   is actually periodic, not merely non-drifting. A hard minimum-periods floor
%   (independent of the tolerance check) prevents firing before the system has
%   had a chance to complete enough cycles to fit a meaningful window at all.
%
%   opts fields: enabled, minPeriods, windowPeriods, checkEveryPeriods,
%   amplitudeRelTol, phaseAbsTolDeg.
%
%   result fields: converged, periods_elapsed, amplitude, phase_deg,
%   amplitude_rel_change, phase_change_deg.

T = 2 * pi / omega;
periods_elapsed = t(end) / T;

result = struct( ...
    'converged', false, ...
    'periods_elapsed', periods_elapsed, ...
    'amplitude', NaN, ...
    'phase_deg', NaN, ...
    'amplitude_rel_change', Inf, ...
    'phase_change_deg', Inf);

if ~opts.enabled || periods_elapsed < opts.minPeriods || periods_elapsed < 2 * opts.windowPeriods
    return
end

t_end = t(end);
in_window_b = t >= (t_end - opts.windowPeriods * T);
in_window_a = (t >= (t_end - 2 * opts.windowPeriods * T)) & ~in_window_b;

if nnz(in_window_a) < 3 || nnz(in_window_b) < 3
    return
end

[amp_a, phase_a] = fit_oscillation(t(in_window_a), z(in_window_a), omega);
[amp_b, phase_b] = fit_oscillation(t(in_window_b), z(in_window_b), omega);

amp_rel_change = abs(amp_b - amp_a) / max(amp_a, eps);
phase_change_deg = abs(mod(rad2deg(phase_b - phase_a) + 180, 360) - 180);

result.converged = (amp_rel_change < opts.amplitudeRelTol) && (phase_change_deg < opts.phaseAbsTolDeg);
result.amplitude = amp_b;
result.phase_deg = rad2deg(phase_b);
result.amplitude_rel_change = amp_rel_change;
result.phase_change_deg = phase_change_deg;

end
