"""
Synthetic-signal tests for the corrected early-stop convergence check — these directly
encode the fix for the audited incident (MATLAB's rolling-mean check declared every one
of a 30-case sweep "converged" at exactly the earliest allowed check point, 3 forcing
periods, because it measured whether a smoothed baseline had stopped drifting rather
than whether the oscillation had reached periodic steady state).

`z(t) = z0 + A * envelope(t) * cos(omega*t - phi)`, where `envelope(t) = 1 -
exp(-t/tau)` models an oscillation whose amplitude is still *building up* toward its
steady-state value `A` (as opposed to MATLAB's failure mode, which was fooled by a
merely-decaying additive offset while the oscillation itself hadn't started) — this is a
strictly harder, more realistic version of the failure mode: a check that only looks at
`std(rolling mean)/|mean|` would be fooled by this just as badly, since the running mean
of a growing-envelope oscillation is also smooth.
"""
function synthetic_growing_oscillation(; A::Float64 = 1.0, phi::Float64 = 0.4, omega::Float64 = 2 * pi * 10.0,
                                         tau::Float64 = Inf, z0::Float64 = 0.3, t_end::Float64,
                                         steps_per_period::Int = 60)
    T = 2 * pi / omega
    dt = T / steps_per_period
    t = collect(0.0:dt:t_end)
    envelope = isinf(tau) ? ones(length(t)) : (1.0 .- exp.(-t ./ tau))
    z = z0 .+ A .* envelope .* cos.(omega .* t .- phi)
    return t, z
end

@testset "check_convergence" begin
    omega = 2 * pi * 10.0
    T = 2 * pi / omega
    opts = ConvergenceOptions(minPeriods = 10.0, windowPeriods = 5.0, amplitudeRelTol = 0.02, phaseAbsTolDeg = 2.0)

    @testset "does not false-converge on a slow-building oscillation (the audited failure mode)" begin
        # tau = 8 periods: at the earliest allowed check (10 periods), the envelope is
        # still visibly different between the two 5-period comparison windows.
        t, z = synthetic_growing_oscillation(; omega = omega, tau = 8 * T, t_end = 10 * T)
        result = check_convergence(t, z, omega, opts)
        @test !result.converged
        @test result.amplitude_rel_change > opts.amplitudeRelTol
    end

    @testset "converges later for a slower transient, not at a fixed period count regardless of dynamics" begin
        # This is the direct test of what the MATLAB bug got wrong: convergence timing
        # must depend on the actual signal, not fire at the same period count every time
        # regardless of how fast the underlying transient decays.
        # tau must be small enough that the transient is negligible even within the
        # FIRST comparison window (which starts at t=0 at the earliest allowed check
        # point) — 0.3*T left a non-negligible fraction of that 5-period window still
        # ramping up and made this assertion flaky; 0.01*T (same order used by the
        # "never fires before minPeriods" case below, which passed) does not.
        t_fast, z_fast = synthetic_growing_oscillation(; omega = omega, tau = 0.01 * T, t_end = 10 * T)
        fast_result = check_convergence(t_fast, z_fast, omega, opts)
        @test fast_result.converged  # transient long decayed by the floor

        t_slow, z_slow = synthetic_growing_oscillation(; omega = omega, tau = 3.0 * T, t_end = 10 * T)
        slow_result_at_floor = check_convergence(t_slow, z_slow, omega, opts)
        @test !slow_result_at_floor.converged  # same floor, but this transient hasn't decayed yet

        t_slow_later, z_slow_later = synthetic_growing_oscillation(; omega = omega, tau = 3.0 * T, t_end = 30 * T)
        slow_result_later = check_convergence(t_slow_later, z_slow_later, omega, opts)
        @test slow_result_later.converged  # ...but does converge once it actually has
    end

    @testset "never fires before minPeriods, even for a trivially-already-steady signal" begin
        t_short, z_short = synthetic_growing_oscillation(; omega = omega, tau = 0.01 * T, t_end = 8 * T)
        @test check_convergence(t_short, z_short, omega, opts).periods_elapsed < opts.minPeriods
        @test !check_convergence(t_short, z_short, omega, opts).converged

        t_long, z_long = synthetic_growing_oscillation(; omega = omega, tau = 0.01 * T, t_end = 11 * T)
        @test check_convergence(t_long, z_long, omega, opts).converged
    end

    @testset "converged flag matches directly recomputing the amplitude/phase predicate" begin
        t, z = synthetic_growing_oscillation(; omega = omega, tau = 1.5 * T, t_end = 20 * T)
        result = check_convergence(t, z, omega, opts)

        T_period = 2 * pi / omega
        t_end = t[end]
        win_b = t .>= (t_end - opts.windowPeriods * T_period)
        win_a = (t .>= (t_end - 2 * opts.windowPeriods * T_period)) .& .!win_b
        fit_a = fit_oscillation(t[win_a], z[win_a], omega)
        fit_b = fit_oscillation(t[win_b], z[win_b], omega)
        amp_rel = abs(fit_b.amplitude - fit_a.amplitude) / fit_a.amplitude
        phase_deg = abs(mod(rad2deg(fit_b.phase - fit_a.phase) + 180, 360) - 180)
        expected = amp_rel < opts.amplitudeRelTol && phase_deg < opts.phaseAbsTolDeg

        @test result.converged == expected
    end

    @testset "disabled options never converge" begin
        t, z = synthetic_growing_oscillation(; omega = omega, tau = 0.01 * T, t_end = 20 * T)
        disabled = ConvergenceOptions(enabled = false)
        @test !check_convergence(t, z, omega, disabled).converged
    end
end
