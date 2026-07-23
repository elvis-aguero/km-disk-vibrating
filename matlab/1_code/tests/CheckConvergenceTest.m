classdef CheckConvergenceTest < matlab.unittest.TestCase
    % Direct port of julia/test/unit/test_convergence_check.jl -- same synthetic
    % signals, same parameters, same expected outcomes (already verified passing
    % there; the underlying math, a plain M\z least-squares fit, is identical).
    %
    % z(t) = z0 + A*envelope(t)*cos(omega*t - phi), envelope(t) = 1 - exp(-t/tau),
    % models an oscillation whose amplitude is still *building up* toward its
    % steady-state value -- this is what the old rolling-mean-ratio check was
    % fooled by (it only measured whether a smoothed baseline had stopped
    % drifting, not whether the oscillation had reached periodic steady state).

    properties (Constant)
        omega = 2 * pi * 10.0
    end

    methods (Static)
        function [t, z] = growingOscillation(A, phi, tau, t_end, stepsPerPeriod)
            omega_ = CheckConvergenceTest.omega;
            T = 2 * pi / omega_;
            dt = T / stepsPerPeriod;
            t = (0:dt:t_end)';
            if isinf(tau)
                envelope = ones(size(t));
            else
                envelope = 1.0 - exp(-t / tau);
            end
            z0 = 0.3;
            z = z0 + A * envelope .* cos(omega_ * t - phi);
        end
    end

    methods (TestMethodSetup)
        function setup(testCase) %#ok<MANU>
        end
    end

    methods (Test)
        function doesNotFalseConvergeOnSlowBuildingOscillation(testCase)
            omega_ = testCase.omega;
            T = 2 * pi / omega_;
            opts = struct('enabled', true, 'minPeriods', 10.0, 'windowPeriods', 5.0, ...
                'checkEveryPeriods', 1.0, 'amplitudeRelTol', 0.02, 'phaseAbsTolDeg', 2.0);

            [t, z] = CheckConvergenceTest.growingOscillation(1.0, 0.4, 8 * T, 10 * T, 60);
            result = check_convergence(t, z, omega_, opts);
            testCase.verifyFalse(result.converged);
            testCase.verifyGreaterThan(result.amplitude_rel_change, opts.amplitudeRelTol);
        end

        function convergesLaterForSlowerTransientNotFixedPeriodCount(testCase)
            omega_ = testCase.omega;
            T = 2 * pi / omega_;
            opts = struct('enabled', true, 'minPeriods', 10.0, 'windowPeriods', 5.0, ...
                'checkEveryPeriods', 1.0, 'amplitudeRelTol', 0.02, 'phaseAbsTolDeg', 2.0);

            [t_fast, z_fast] = CheckConvergenceTest.growingOscillation(1.0, 0.4, 0.01 * T, 10 * T, 60);
            fast_result = check_convergence(t_fast, z_fast, omega_, opts);
            testCase.verifyTrue(fast_result.converged);

            [t_slow, z_slow] = CheckConvergenceTest.growingOscillation(1.0, 0.4, 3.0 * T, 10 * T, 60);
            slow_result_at_floor = check_convergence(t_slow, z_slow, omega_, opts);
            testCase.verifyFalse(slow_result_at_floor.converged);

            [t_slow_later, z_slow_later] = CheckConvergenceTest.growingOscillation(1.0, 0.4, 3.0 * T, 30 * T, 60);
            slow_result_later = check_convergence(t_slow_later, z_slow_later, omega_, opts);
            testCase.verifyTrue(slow_result_later.converged);
        end

        function neverFiresBeforeMinPeriods(testCase)
            omega_ = testCase.omega;
            T = 2 * pi / omega_;
            opts = struct('enabled', true, 'minPeriods', 10.0, 'windowPeriods', 5.0, ...
                'checkEveryPeriods', 1.0, 'amplitudeRelTol', 0.02, 'phaseAbsTolDeg', 2.0);

            [t_short, z_short] = CheckConvergenceTest.growingOscillation(1.0, 0.4, 0.01 * T, 8 * T, 60);
            result_short = check_convergence(t_short, z_short, omega_, opts);
            testCase.verifyLessThan(result_short.periods_elapsed, opts.minPeriods);
            testCase.verifyFalse(result_short.converged);

            [t_long, z_long] = CheckConvergenceTest.growingOscillation(1.0, 0.4, 0.01 * T, 11 * T, 60);
            testCase.verifyTrue(check_convergence(t_long, z_long, omega_, opts).converged);
        end

        function convergedFlagMatchesDirectPredicate(testCase)
            omega_ = testCase.omega;
            T = 2 * pi / omega_;
            opts = struct('enabled', true, 'minPeriods', 10.0, 'windowPeriods', 5.0, ...
                'checkEveryPeriods', 1.0, 'amplitudeRelTol', 0.02, 'phaseAbsTolDeg', 2.0);

            [t, z] = CheckConvergenceTest.growingOscillation(1.0, 0.4, 1.5 * T, 20 * T, 60);
            result = check_convergence(t, z, omega_, opts);

            t_end = t(end);
            in_b = t >= (t_end - opts.windowPeriods * T);
            in_a = (t >= (t_end - 2 * opts.windowPeriods * T)) & ~in_b;
            [amp_a, phase_a] = fit_oscillation(t(in_a), z(in_a), omega_);
            [amp_b, phase_b] = fit_oscillation(t(in_b), z(in_b), omega_);
            amp_rel = abs(amp_b - amp_a) / amp_a;
            phase_deg = abs(mod(rad2deg(phase_b - phase_a) + 180, 360) - 180);
            expected = (amp_rel < opts.amplitudeRelTol) && (phase_deg < opts.phaseAbsTolDeg);

            testCase.verifyEqual(result.converged, expected);
        end

        function disabledOptionsNeverConverge(testCase)
            omega_ = testCase.omega;
            T = 2 * pi / omega_;
            [t, z] = CheckConvergenceTest.growingOscillation(1.0, 0.4, 0.01 * T, 20 * T, 60);
            disabled = struct('enabled', false, 'minPeriods', 10.0, 'windowPeriods', 5.0, ...
                'checkEveryPeriods', 1.0, 'amplitudeRelTol', 0.02, 'phaseAbsTolDeg', 2.0);
            testCase.verifyFalse(check_convergence(t, z, omega_, disabled).converged);
        end
    end
end
