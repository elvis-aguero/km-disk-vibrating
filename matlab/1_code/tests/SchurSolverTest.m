classdef SchurSolverTest < matlab.unittest.TestCase
    % Stage 2 MATLAB-internal tests (T3). schur must equal lu to 1e-9 relative (norm of the CoM trajectory difference
    % over the norm of the LU trajectory), for euler and bdf2, bath-driven and disk-forced, on D5Quant20 (nr = 50)
    % with the single-layer matrix. "auto" must resolve to schur on a tridiagonal Laplacian and fall back otherwise.
    % Tolerances fixed before any run; not to be loosened.

    properties (Constant)
        RelTol = 1e-9
    end

    methods (Test)
        function schurEqualsLuBathEuler(testCase)
            testCase.compareSolvers('euler', false);
        end

        function schurEqualsLuDiskEuler(testCase)
            testCase.compareSolvers('euler', true);
        end

        function schurEqualsLuBathBdf2(testCase)
            testCase.compareSolvers('bdf2', false);
        end

        function schurEqualsLuDiskBdf2(testCase)
            testCase.compareSolvers('bdf2', true);
        end

        function autoEqualsSchurEndToEnd(testCase)
            a = SchurSolverTest.run('euler', false, 'auto');
            s = SchurSolverTest.run('euler', false, 'schur');
            testCase.verifyEqual(a, s, 'AbsTol', 0, 'auto must run the schur path on the uniform grid');
        end

        function autoResolvesSchurOnTridiagonalAndFallsBack(testCase)
            [~, L] = domainMaker(20, 5);                       % uniform Laplacian, tridiagonal
            nr = size(L, 1);
            big = 1e12;
            testCase.verifyEqual(resolve_solver_type("auto", L, nr, 31, big), "schur");
            Lfull = L;  Lfull(1, 5) = 1e-3;                    % not tridiagonal any more
            testCase.verifyEqual(resolve_solver_type("auto", Lfull, nr, 31, big), "lu");
            testCase.verifyEqual(resolve_solver_type("auto", Lfull, nr, 31, 0), "gmres");
            testCase.verifyEqual(resolve_solver_type("auto", L, nr, 31, 0), "gmres");
            [~, info] = resolve_solver_type("auto", L, nr, 31, big);
            testCase.verifyLessThan(info.schurBytes, info.luBytes);
            testCase.verifyEqual(resolve_solver_type("lu", L, nr, 31, big), "lu");       % explicit passes through
            testCase.verifyEqual(resolve_solver_type("gmres", L, nr, 31, big), "gmres");
            testCase.verifyEqual(resolve_solver_type("schur", L, nr, 31, 0), "schur");
        end

        function luFactorCountMatchesJuliaDocs(testCase)
            % Julia lu_cache.jl docstring: 31 distinct of 60 steps at -90 deg; disk-forced: 1.
            testCase.verifyEqual(lu_factor_count(60, -pi / 2, true), 31);
            testCase.verifyEqual(lu_factor_count(60, -pi / 2, false), 1);
        end

        function schurRejectsNonTridiagonalLaplacian(testCase)
            PC = struct('nr', 20, 'contact_points', 6, 'reynolds', 1, 'froude', 1, 'weber', 1, 'obj_mass', 1, ...
                'laplacian', ones(20), 'DTN', eye(20), 'pressure_integral', ones(12, 1));
            testCase.verifyError(@() schur_system(PC, 0.1), 'schur_system:notTridiagonal');
        end
    end

    methods
        function compareSolvers(testCase, scheme, diskForced)
            ref = SchurSolverTest.run(scheme, diskForced, 'lu');
            sch = SchurSolverTest.run(scheme, diskForced, 'schur');
            testCase.assertEqual(size(ref), size(sch));
            d = norm(sch - ref) / norm(ref);
            fprintf('[schur-vs-lu] %s %s: relative difference %.2e\n', scheme, ...
                char(ternary(diskForced, 'disk', 'bath')), d);
            testCase.verifyLessThanOrEqual(d, testCase.RelTol);
        end
    end

    methods (Static)
        function CoM = run(scheme, diskForced, solver)
            freqHz = 30.0;  omega = 2 * pi * freqHz;
            if diskForced
                fA = 0.1 * 0.0283 * 981.0;  fF = freqHz;  bA = 0.0;
            else
                fA = 0.0;  fF = 90.0;  bA = 0.2 * 981.0 / omega^2;
            end
            args = {'diskRadius', 0.2, 'diskMass', 0.0283, ...
                'forceAmplitude', fA, 'forceFrequency', fF, ...
                'bathAmplitude', bA, 'bathFrequency', freqHz, 'phaseDifference', -90.0, ...
                'bathDensity', 1.175, 'bathSurfaceTension', 66.5, 'bathViscosity', 0.18 / 1.175, ...
                'bathDiameter', 20.0, 'spatialResolution', 5.0, 'temporalResolution', 48, ...
                'simulationTime', 8 / freqHz, 'debug_flag', false, 'startStatic', true, 'earlyStop', false, ...
                'dtnMethod', 'single_layer', 'timeScheme', scheme, 'solverType', solver};
            [~, CoM, ~] = solve_motion(args{:});
            CoM = CoM(:);
        end
    end
end

function v = ternary(c, a, b)
if c, v = a; else, v = b; end
end
