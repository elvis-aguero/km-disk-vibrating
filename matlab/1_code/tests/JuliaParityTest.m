classdef JuliaParityTest < matlab.unittest.TestCase
    % Cross-implementation parity check: shells out to julia/scripts/parity_reference.jl
    % to run the exact same bath-driven case through the Julia port's run_simulation, then
    % runs the identical case through solve_motion.m in THIS process, and compares the two
    % live. Both sides are computed fresh on every run -- deliberately not a
    % precomputed/hardcoded reference, which would silently go stale the moment either
    % implementation's physics changed without this test being re-derived.
    %
    % Uses the cheapest domain provably shared between the two implementations --
    % D5Quant20 (nr=50) -- whose DTN operator was cross-validated byte-for-byte between
    % Julia's native generator and this MATLAB cache earlier this session (~2e-15
    % relative difference), so any disagreement here points at the physics/stepping code,
    % not the DTN operator. Requires `julia` on PATH (see .github/workflows/matlab-ci.yml,
    % which sets up both runtimes in the same job) and julia/ instantiated.
    %
    % Tolerance is looser than floating-point-identical because the two implementations
    % use different BLAS/LAPACK libraries under the hood (MATLAB's vs OpenBLAS), so
    % per-step rounding differs slightly and accumulates over the run's ~384 sequential
    % direct solves; 0.1% relative is intended to catch a real algorithmic divergence, not
    % flag ordinary cross-platform floating-point noise.

    methods (Test)
        function matchesLiveJuliaRun(testCase)
            repoRoot = fullfile(fileparts(mfilename('fullpath')), '..', '..', '..');
            juliaScript = fullfile(repoRoot, 'julia', 'scripts', 'parity_reference.jl');
            testCase.assertTrue(isfile(juliaScript), sprintf('not found: %s', juliaScript));

            % MATLAB's system() inherits MATLAB's own LD_LIBRARY_PATH, which points at
            % MATLAB's bundled (often older) shared libraries. A JIT-compiled runtime
            % like Julia picking those up instead of the system libraries crashes at
            % startup (segfault, well before reaching user code) -- a known MATLAB
            % system()-with-external-binaries gotcha, not a bug in the Julia side.
            % `env -u` strips just those two variables for this subprocess only.
            [status, cmdout] = system(sprintf('env -u LD_LIBRARY_PATH -u LD_PRELOAD julia "%s" 2>/dev/null', juliaScript));
            testCase.assertEqual(status, 0, sprintf('parity_reference.jl failed (see stderr suppressed above): %s', cmdout));

            ref = sscanf(cmdout, '%f');
            testCase.assertEqual(numel(ref), 5, sprintf('unexpected parity_reference.jl output: "%s"', cmdout));
            ref_t_end = ref(1);
            ref_CoM_first = ref(2);
            ref_amplitude = ref(3);
            ref_phase = ref(4);
            ref_offset = ref(5);

            % Must match parity_reference.jl's SimulationParams exactly.
            freqHz = 30.0;
            gamma = 0.2;
            omega = 2 * pi * freqHz;
            A_bath = gamma * 981.0 / omega^2;
            nPeriods = 8;

            [t_s, CoM_cm, ~] = solve_motion( ...
                'diskRadius', 0.2, 'diskMass', 0.0283, ...
                'forceAmplitude', 0.0, 'forceFrequency', 90.0, ...
                'bathAmplitude', A_bath, 'bathFrequency', freqHz, ...
                'phaseDifference', -90.0, ...
                'bathDensity', 1.175, 'bathSurfaceTension', 66.5, 'bathViscosity', 0.18 / 1.175, ...
                'bathDiameter', 20.0, 'spatialResolution', 5.0, 'temporalResolution', 48, ...
                'simulationTime', nPeriods / freqHz, ...
                'debug_flag', false, 'solverType', 'lu', ...
                'startStatic', true, 'earlyStop', false);

            tol = 1e-3; % relative

            testCase.verifyEqual(t_s(end), ref_t_end, 'RelTol', tol);
            testCase.verifyEqual(CoM_cm(1), ref_CoM_first, 'RelTol', tol);

            T = 1 / freqHz;
            tail = t_s >= (t_s(end) - 3 * T);
            [amp, phase, offset] = fit_oscillation(t_s(tail), CoM_cm(tail), omega);

            testCase.verifyEqual(amp, ref_amplitude, 'RelTol', tol);
            testCase.verifyEqual(offset, ref_offset, 'RelTol', tol);
            % Phase is an angle, not a plain magnitude -- compare via wrapped difference
            % rather than RelTol (which would misbehave near phase == 0).
            phase_diff_deg = abs(mod(rad2deg(phase - ref_phase) + 180, 360) - 180);
            testCase.verifyLessThan(phase_diff_deg, 0.5);
        end
    end
end
