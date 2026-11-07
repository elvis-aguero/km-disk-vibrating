classdef GradedGridTest < matlab.unittest.TestCase
    % Stage 3 tests (T4), graded radial grid. Criteria fixed before any run, never loosened:
    %  (a) radii, s, s', the Laplacian and the DtN matrix equal Julia's (live dump, parity_dump.jl): radii and s to
    %      1e-12 relative to max, the Laplacian to 1e-12 max|L|, the DtN matrix to 1e-10 max|A|;
    %  (b) graded DtN routine on uniform nodes equals the uniform routine to 1e-8 max|A|;
    %  (c) point-source test on the graded grid (bead fluid, 700 x 25, 40 Hz, defaults): errors <= 1.5 x the Julia
    %      numbers 1.96e-5, 2.17e-6, 3.82e-6 for a = 1, 3, 10 disk radii;
    %  (d) uniform domainMaker call unchanged; (e) the end-to-end graded run vs Julia is in JuliaParityTest (1e-8).

    methods (Test)
        function matchesJuliaGridOperators(testCase)
            outdir = tempname();
            cleanup = onCleanup(@() rmdir_if(outdir)); %#ok<NASGU>
            run_julia_script(testCase, 'parity_dump.jl', outdir);
            rd = @(f) readmatrix(fullfile(outdir, f), 'FileType', 'text', 'Delimiter', ' ');
            p = GradedGridTest.beadArgs(700);
            for tag = ["bead40", "bead40_wz0"]
                if tag == "bead40"
                    g = graded_grid_for(40, p{:});
                else
                    g = graded_grid_for(40, p{:}, 'wallZone', 0);
                end
                rJ = rd(sprintf('grid_%s_r.txt', tag));
                sJ = rd(sprintf('grid_%s_s.txt', tag));
                spJ = rd(sprintf('grid_%s_sp.txt', tag));
                testCase.assertEqual(numel(g.r), numel(rJ), sprintf('%s: node count', tag));
                dr = max(abs(g.r(:) - rJ(:))) / max(abs(rJ));
                ds = max(abs(g.s(:) - sJ(:))) / max(abs(sJ));
                dsp = max(abs(g.sp(:) - spJ(:))) / max(abs(sJ));
                fprintf('[grid-parity] %s nr=%d: r %.2e, s %.2e, sp %.2e (relative to max)\n', tag, g.nr, dr, ds, dsp);
                testCase.verifyLessThanOrEqual(dr, 1e-12, sprintf('%s r', tag));
                testCase.verifyLessThanOrEqual(ds, 1e-12, sprintf('%s s', tag));
                testCase.verifyLessThanOrEqual(dsp, 1e-12, sprintf('%s sp', tag));
            end
            ps = GradedGridTest.beadArgs(200);
            gs = graded_grid_for(40, ps{:});
            testCase.assertEqual(numel(gs.r), numel(rd('grid_small_r.txt')));
            testCase.verifyLessThanOrEqual(max(abs(gs.r(:) - rd('grid_small_r.txt'))) / max(gs.r), 1e-12, 'small r');
            [~, L] = domainMaker(gs);
            LJ = rd('grid_small_laplacian.txt');
            dL = max(abs(L(:) - LJ(:))) / max(abs(LJ(:)));
            fprintf('[grid-parity] small Laplacian: %.2e max|L|\n', dL);
            testCase.verifyLessThanOrEqual(dL, 1e-12, 'Laplacian');
            A = generate_dtn_grid(gs);
            AJ = rd('grid_small_dtn.txt');
            dA = max(abs(A(:) - AJ(:))) / max(abs(AJ(:)));
            fprintf('[grid-parity] small graded DtN: %.2e max|A|\n', dA);
            testCase.verifyLessThanOrEqual(dA, 1e-10, 'DtN');
        end

        function gradedDtnOnUniformNodesEqualsUniformRoutine(testCase)
            g = uniform_grid(5, 10);                       % nr = 50, D = 20
            Au = generate_dtn_single_layer(g.nr, 2 * g.R_req, 'order', 6);
            Ag = generate_dtn_grid(g, 'forceGeneral', true);
            d = max(abs(Au(:) - Ag(:))) / max(abs(Au(:)));
            fprintf('[grid-uniform] general vs uniform DtN: %.2e max|A|\n', d);
            testCase.verifyLessThanOrEqual(d, 1e-8);
            testCase.verifyEqual(generate_dtn_grid(g), Au, 'uniform grid must call the uniform routine unchanged');
        end

        function uniformGridDomainIsBitIdenticalToDomainMaker(testCase)
            [dr0, L0, I0] = domainMaker(20, 5);
            [dr1, L1, I1, nr, r] = domainMaker(uniform_grid(5, 10));
            testCase.verifyEqual(dr1, dr0, 'AbsTol', 0);
            testCase.verifyEqual(L1, L0, 'AbsTol', 0);
            testCase.verifyEqual(I1, I0, 'AbsTol', 0);
            testCase.verifyEqual(nr, 50);
            testCase.verifyEqual(r, dr0 * (0:49).', 'AbsTol', 0);
        end

        function gradedGridGuardsAndShape(testCase)
            g = graded_grid(5, 10, 'drMax', 0.4, 'q', 1.08);
            testCase.verifyFalse(g.uniform);
            testCase.verifyEqual(numel(g.r), g.nr + 4);
            testCase.verifyGreaterThanOrEqual(g.r(g.nr + 1), 10);        % the wall is the first node with r >= R
            testCase.verifyLessThan(g.r(g.nr), 10);
            testCase.verifyTrue(all(diff(g.r) > 0));
            testCase.verifyLessThanOrEqual(max(g.s), 0.4 * (1 + 1e-12)); % s <= dr_max
            testCase.verifyEqual(g.r(1:g.xin + 1), g.dr0 * (0:g.xin).', 'AbsTol', 0);   % identity region
            u = graded_grid(5, 10, 'drMax', 0.2);                        % drMax <= dr0: uniform fallback
            testCase.verifyTrue(u.uniform);
            testCase.verifyError(@() graded_grid(25, 100, 'drMax', 0.5, 'rIn', 1.0), 'graded_grid:rIn');
            testCase.verifyError(@() graded_grid(5, 10, 'drMax', 0.4, 'laplacianOrder', 3), 'graded_grid:order');
            testCase.verifyError(@() graded_grid(5, 10), 'graded_grid:drMax');
        end

        function pointSourceOnGradedGrid(testCase)
            p = GradedGridTest.beadArgs(700);
            g = graded_grid_for(40, p{:});
            A = generate_dtn_grid(g);
            r = g.r(1:g.nr);
            ref = [1 1.96e-5; 3 2.17e-6; 10 3.82e-6];
            for k = 1:3
                e = DtnSingleLayerTest.globErr(A, r, ref(k, 1));
                fprintf('[grid-pointsource] a=%g: %.3e (Julia %.3e), nr=%d\n', ref(k, 1), e, ref(k, 2), g.nr);
                testCase.verifyLessThanOrEqual(e, 1.5 * ref(k, 2), sprintf('a = %g', ref(k, 1)));
            end
        end
    end

    methods (Static)
        function c = beadArgs(bathDiameter)
            % bead fluid and 25 nodes per disk radius, as test_dtn_pointsource.jl
            c = {'bathDensity', 1.148, 'bathSurfaceTension', 68.0, 'bathViscosity', 0.079, 'g', 981.0, ...
                'diskRadius', 0.07938, 'bathDiameter', bathDiameter, 'spatialResolution', 25};
        end
    end
end

function rmdir_if(d)
if isfolder(d), rmdir(d, 's'); end
end
