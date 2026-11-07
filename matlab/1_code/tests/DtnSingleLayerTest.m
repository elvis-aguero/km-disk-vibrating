classdef DtnSingleLayerTest < matlab.unittest.TestCase
    % Stage 1 DtN tests. T1a: MATLAB generate_dtn_single_layer equals Julia's (live dump by julia/scripts/parity_dump.jl),
    % max|diff| <= 1e-10 max|A|, orders 4 and 6. T1b: point-source solution with an algebraic tail
    % (port of julia/test/integration/test_dtn_pointsource.jl, uniform part): errors <= 1.5 x the Julia numbers.
    % Tolerances were fixed before any run and are not to be loosened.

    properties (Constant)
        ParityTol = 1e-10
    end

    methods (Test)
        function matchesJuliaSingleLayer(testCase)
            outdir = tempname();
            cleanup = onCleanup(@() rmdir_if(outdir)); %#ok<NASGU>
            run_julia_script(testCase, 'parity_dump.jl', outdir);
            for order = [4 6]
                J = readmatrix(fullfile(outdir, sprintf('dtn_sl_nr50_D20_o%d.txt', order)), ...
                    'FileType', 'text', 'Delimiter', ' ');
                M = generate_dtn_single_layer(50, 20, 'order', order);
                testCase.assertEqual(size(J), size(M));
                d = max(abs(M(:) - J(:))) / max(abs(J(:)));
                fprintf('[dtn-parity] order %d: max|diff|/max|A| = %.2e\n', order, d);
                testCase.verifyLessThanOrEqual(d, testCase.ParityTol, sprintf('order %d', order));
            end
        end

        function generateDtnDispatchesToSingleLayer(testCase)
            A1 = generate_dtn(20, 20);
            A2 = generate_dtn_single_layer(20, 20, 'order', 6);
            testCase.verifyEqual(A1, A2);
            A4 = generate_dtn(20, 20, 'order', 4);
            testCase.verifyEqual(A4, generate_dtn_single_layer(20, 20, 'order', 4));
        end

        function pointSourceUniformOrder6(testCase)
            A = generate_dtn_single_layer(400, 800, 'order', 6);
            r = (0:399).';
            e10 = DtnSingleLayerTest.globErr(A, r, 10);
            e25 = DtnSingleLayerTest.globErr(A, r, 25);
            fprintf('[dtn-pointsource] a=10: %.3e (Julia 4.62e-4), a=25: %.3e (Julia 1.77e-5)\n', e10, e25);
            testCase.verifyLessThanOrEqual(e10, 1.5 * 4.62e-4, 'a = 10');
            testCase.verifyLessThanOrEqual(e25, 1.5 * 1.77e-5, 'a = 25');
        end
    end

    methods (Static)
        function e = globErr(A, r, a)
            f = a ./ (r.^2 + a^2).^(3/2);
            Nf = (2*a^2 - r.^2) ./ (r.^2 + a^2).^(5/2);
            keep = r <= max(r) / 2;
            Af = A * f;
            e = max(abs(Af(keep) - Nf(keep))) / max(abs(Nf(keep)));
        end
    end
end

function rmdir_if(d)
if isfolder(d), rmdir(d, 's'); end
end
