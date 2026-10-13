classdef JuliaParityTest < matlab.unittest.TestCase
    % Cross-implementation parity check: shells out to julia/scripts/parity_reference.jl,
    % which runs a set of cases through the Julia port's run_simulation, then runs the
    % identical cases through solve_motion.m in THIS process and compares the two live.
    % Both sides are computed fresh on every run -- deliberately not a precomputed
    % reference, which would silently go stale the moment either implementation changed.
    %
    % Cases (must match parity_reference.jl): bath-driven on the LU path, disk-forced,
    % bath-driven on the GMRES path, and bath-driven with early stopping on (same stop step).
    %
    % Uses the cheapest domain shared by both implementations, D5Quant20 (nr=50), so any
    % disagreement points at the physics/stepping code, not the DTN operator. Requires
    % `julia` on PATH (see .github/workflows/matlab-ci.yml).
    %
    % Tolerances: both sides run the same algorithm in double precision, so the direct (LU)
    % cases should agree to rounding accumulated over a few hundred solves -- far tighter
    % than a physical tolerance. The GMRES case is looser because each step is only solved
    % to gmresTolerance (and MATLAB's GMRES is unpreconditioned while Julia's is
    % preconditioned, so they stop at different iterates within that tolerance). Observed
    % differences are printed to the CI log.

    properties (Constant)
        RelTolDirect = 1e-8
        RelTolGmres = 1e-6
        PhaseTolDegDirect = 1e-6
        PhaseTolDegGmres = 1e-4
    end

    methods (Test)
        function matchesLiveJuliaRun(testCase)
            ref = JuliaParityTest.runJulia(testCase);

            freqHz = 30.0;
            omega = 2 * pi * freqHz;
            base = {'diskRadius', 0.2, 'diskMass', 0.0283, ...
                'forceAmplitude', 0.0, 'forceFrequency', 90.0, ...
                'bathAmplitude', 0.2 * 981.0 / omega^2, 'bathFrequency', freqHz, ...
                'phaseDifference', -90.0, ...
                'bathDensity', 1.175, 'bathSurfaceTension', 66.5, 'bathViscosity', 0.18 / 1.175, ...
                'bathDiameter', 20.0, 'spatialResolution', 5.0, 'temporalResolution', 48, ...
                'simulationTime', 8 / freqHz, 'debug_flag', false, 'solverType', 'lu', ...
                'startStatic', true, 'earlyStop', false};
            % Must match parity_reference.jl's cases exactly (later pairs override earlier).
            cases = {
                'bath_lu',        {}
                'disk_lu',        {'forceAmplitude', 0.1 * 0.0283 * 981.0, 'forceFrequency', freqHz, 'bathAmplitude', 0.0}
                'bath_gmres',     {'solverType', 'gmres', 'gmresTolerance', 1e-12}
                'bath_earlystop', {'earlyStop', true, 'simulationTime', 30 / freqHz}
            };

            for c = 1:size(cases, 1)
                name = cases{c, 1};
                testCase.assertTrue(isKey(ref, name), sprintf('parity_reference.jl printed no line for case %s', name));
                r = ref(name);

                args = JuliaParityTest.mergeArgs(base, cases{c, 2});
                [t_s, CoM_cm, ~] = solve_motion(args{:});

                tail = t_s >= (t_s(end) - 3 / freqHz);
                [amp, phase, offset] = fit_oscillation(t_s(tail), CoM_cm(tail), omega);

                if strcmp(name, 'bath_gmres')
                    relTol = testCase.RelTolGmres; phaseTol = testCase.PhaseTolDegGmres;
                else
                    relTol = testCase.RelTolDirect; phaseTol = testCase.PhaseTolDegDirect;
                end
                phase_diff_deg = abs(mod(rad2deg(phase - r(4)) + 180, 360) - 180);
                fprintf('[parity] %-15s t_end %.1e | CoM(1) %.1e | amplitude %.1e | offset %.1e (relative) | phase %.1e deg\n', ...
                    name, relDiff(t_s(end), r(1)), relDiff(CoM_cm(1), r(2)), relDiff(amp, r(3)), ...
                    relDiff(offset, r(5)), phase_diff_deg);

                % The stop step must match exactly (early stop decides on the same fits).
                testCase.verifyEqual(t_s(end), r(1), 'RelTol', 1e-12, sprintf('%s: t_end', name));
                testCase.verifyEqual(CoM_cm(1), r(2), 'RelTol', relTol, sprintf('%s: CoM(1)', name));
                testCase.verifyEqual(amp, r(3), 'RelTol', relTol, sprintf('%s: amplitude', name));
                testCase.verifyEqual(offset, r(5), 'RelTol', relTol, sprintf('%s: offset', name));
                testCase.verifyLessThan(phase_diff_deg, phaseTol, sprintf('%s: phase', name));
            end
        end
    end

    methods (Static, Access = private)
        function ref = runJulia(testCase)
            repoRoot = fullfile(fileparts(mfilename('fullpath')), '..', '..', '..');
            juliaScript = fullfile(repoRoot, 'julia', 'scripts', 'parity_reference.jl');
            testCase.assertTrue(isfile(juliaScript), sprintf('not found: %s', juliaScript));

            % MATLAB's system() inherits MATLAB's own LD_LIBRARY_PATH, which points at
            % MATLAB's bundled (often older) shared libraries; Julia picking those up crashes
            % at startup. `env -u` strips them for this subprocess only. stderr goes to a
            % temp file so a Julia-side error ends up in the test diagnostic.
            stderrFile = [tempname(), '.stderr.txt'];
            cleanupStderr = onCleanup(@() safeDelete(stderrFile));
            cmd = sprintf('env -u LD_LIBRARY_PATH -u LD_PRELOAD julia "%s" 2>"%s"', juliaScript, stderrFile);
            [status, cmdout] = system(cmd);
            stderrText = '';
            if isfile(stderrFile)
                stderrText = fileread(stderrFile);
            end
            testCase.assertEqual(status, 0, sprintf('parity_reference.jl failed (exit %d)\nstdout: %s\nstderr: %s', status, cmdout, stderrText));

            ref = containers.Map();
            lines = strsplit(strtrim(cmdout), newline);
            for k = 1:numel(lines)
                parts = strsplit(strtrim(lines{k}));
                vals = str2double(parts(2:end));
                testCase.assertEqual(numel(vals), 5, sprintf('unexpected parity_reference.jl line: "%s"', lines{k}));
                testCase.assertTrue(all(isfinite(vals)), sprintf('non-numeric parity_reference.jl line: "%s"', lines{k}));
                ref(parts{1}) = vals;
            end
        end

        function args = mergeArgs(base, overrides)
            args = base;
            for k = 1:2:numel(overrides)
                idx = find(strcmp(args(1:2:end), overrides{k}), 1);
                if isempty(idx)
                    args(end + 1:end + 2) = overrides(k:k + 1);
                else
                    args{2 * idx} = overrides{k + 1};
                end
            end
        end
    end
end

function d = relDiff(a, b)
    d = abs(a - b) / max(abs(b), eps);
end

function safeDelete(filePath)
    if isfile(filePath)
        delete(filePath);
    end
end
