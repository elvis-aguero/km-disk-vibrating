classdef FitOscillationTest < matlab.unittest.TestCase
    % Mirrors julia/test/unit/test_fit.jl: fit_oscillation must recover the
    % known amplitude/phase/offset of a synthetic sinusoid exactly (up to
    % floating-point tolerance), including when there's a large DC offset.

    methods (Test)
        function recoversKnownSinusoid(testCase)
            omega = 2 * pi * 3.7;
            true_amp = 0.42;
            true_phase = deg2rad(35);
            true_offset = -0.0264; % mimics a static-equilibrium DC offset

            t = linspace(0, 10 / (omega / (2 * pi)), 2000)';
            z = true_amp * cos(omega * t + true_phase) + true_offset;

            [amp, phase, offset] = fit_oscillation(t, z, omega);

            testCase.verifyEqual(amp, true_amp, 'AbsTol', 1e-8);
            testCase.verifyEqual(phase, true_phase, 'AbsTol', 1e-8);
            testCase.verifyEqual(offset, true_offset, 'AbsTol', 1e-8);
        end

        function rejectsTooFewSamples(testCase)
            testCase.verifyError(@() fit_oscillation([0; 1], [0; 1], 1), ...
                'fit_oscillation:tooFewSamples');
        end
    end
end
