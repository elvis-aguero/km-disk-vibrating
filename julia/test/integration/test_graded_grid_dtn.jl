# Graded grid: DtN matrix against the exact operator, consistency with the uniform routine, and a small end-to-end case.
# Design and criteria (fixed before the first run): docs/graded_grid.md, sections 4 and 5.

isdefined(Main, :_gg_exact_packet) || include(joinpath(@__DIR__, "..", "..", "scripts", "graded_grid_checks.jl"))

@testset "graded grid DtN and end-to-end" begin
    @testset "uniform grid: generate_dtn(grid) is generate_dtn(nr, D) bit for bit" begin
        ug = uniform_grid(1.0, 40.0)                        # nr = 40, D = 80, dr = 1
        A = generate_dtn(40, 80.0)
        @test generate_dtn(ug) == A
        @test generate_dtn_single_layer(ug; order = 4) == generate_dtn_single_layer(40, 80.0; order = 4)
        @test generate_dtn(ug; method = :polar) == generate_dtn_polar(40, 80.0)
        @test_throws ArgumentError generate_dtn(graded_grid(25.0, 20.0; dr_max = 0.5); method = :polar)
    end

    @testset "the non-uniform routine reproduces the uniform one on uniform nodes" begin
        ug = uniform_grid(1.0, 40.0)
        A = generate_dtn(40, 80.0)
        Ag = generate_dtn(ug; force_general = true)
        @test maximum(abs, Ag .- A) <= 1e-8 * maximum(abs, A)
        ug2 = uniform_grid(2.0, 20.0)                       # dr = 0.5: physical units, no rescale
        @test maximum(abs, generate_dtn(ug2; force_general = true) .- generate_dtn(ug2)) <= 1e-8 * maximum(abs, generate_dtn(ug2))
    end

    @testset "C1: graded DtN vs exact wave packet, compared with a uniform grid of spacing dr_max" begin
        sr, R, q, dmax, rin, ppw = 4.0, 40.0, 1.05, 1.0, 2.25, 10.0   # r_in >= 1 + 5 dr0 = 2.25
        g = graded_grid(sr, R; q = q, dr_max = dmax, r_in = rin)
        ug = uniform_grid(1 / dmax, R)
        lam = ppw * dmax
        k0 = 2pi / lam
        s = R / 4
        eg = eu = NaN
        for (grid, tag) in ((g, :graded), (ug, :uniform))
            A = generate_dtn(grid)
            @test size(A) == (grid.nr, grid.nr) && all(isfinite, A)
            rs = grid.r[1:grid.nr]
            f = @. exp(-rs^2 / (2 * s^2)) * _gg_besselj0(k0 * rs)
            e = _gg_local_error(A, f, _gg_exact_packet(s, k0, rs), rs, lam)
            tag === :graded ? (eg = e) : (eu = e)
        end
        @info "C1 (unit-test size): local error" graded = eg uniform = eu ratio = eg / eu nr_graded = g.nr nr_uniform = ug.nr
        @test eg <= 1.5 * eu                                # criterion C1; Python prototype of the same algorithm: 9.98e-3 vs 7.82e-3
        @test eg <= 2e-2                                    # absolute sanity (prototype 1.0e-2)
        @test size(generate_dtn(g; order = 4)) == (g.nr, g.nr)
    end

    @testset "small end-to-end case: graded vs uniform grid (smoke bound; C3 of docs/graded_grid.md is the real criterion)" begin
        # bead fluid and an Al 1/16" sphere on a 120-radius-wide bath, bath-driven at 40 Hz (viscous: waves die before the wall)
        a = 0.07938
        f = 40.0
        p = SimulationParams(bathDensity = 1.148, bathSurfaceTension = 68.0, bathViscosity = 0.079, g = 981.0,
                             diskRadius = a, diskMass = 4 / 3 * pi * a^3 * 2.79, forceAmplitude = 0.0, forceFrequency = 90.0,
                             phaseDifference = -90.0, bathAmplitude = 0.5 * 981 / (2pi * f)^2, bathFrequency = f,
                             bathDiameter = 120.0, spatialResolution = 10.0, temporalResolution = 48,
                             simulationTime = 12 / f, startStatic = true, earlyStop = false, solverType = :lu, timeScheme = :bdf2)
        function gain_phase(params, dtn)
            res = run_simulation(params; dtn = dtn)
            @test res.status == :ok
            z = lab_frame_CoM(res)
            keep = res.t_s .>= res.t_s[end] - 3 / f
            amp, ph, _ = fit_oscillation(res.t_s[keep], z[keep], 2pi * f)
            return amp / (0.5 * 981 / (2pi * f)^2), rad2deg(ph)
        end
        ug = uniform_grid(10.0, 60.0)
        gg = graded_grid_for(p, f; ppw = 16.0, q = 1.05, laplacian_order = 4)
        @test gg.nr < ug.nr ÷ 2
        t0 = time()
        gu, pu = gain_phase(p, generate_dtn(ug))
        t_u = time() - t0
        t0 = time()
        gg_gain, gg_ph = gain_phase(with_overrides(p; grid = gg), generate_dtn(gg))
        t_g = time() - t0
        dphi = mod(gg_ph - pu + 180, 360) - 180
        @info "end-to-end 40 Hz" nr_uniform = ug.nr nr_graded = gg.nr gain_uniform = gu gain_graded = gg_gain gain_rel_diff = gg_gain / gu - 1 phase_diff_deg = dphi time_uniform_s = t_u time_graded_s = t_g
        @test isapprox(gg_gain, gu; rtol = 0.05)            # smoke bound, fixed before the first run
        @test abs(dphi) <= 5.0                              # degrees, same
    end
end
