# Radial grids (src/grid.jl, build_domain(grid)); design and criteria: docs/graded_grid.md.

using LinearAlgebra
isdefined(Main, :_gg_exact_packet) || include(joinpath(@__DIR__, "..", "..", "scripts", "graded_grid_checks.jl"))

@testset "graded grid" begin
    @testset "uniform grid is today's grid, bit for bit" begin
        sr, D = 2.0, 40.0
        ug = uniform_grid(sr, D / 2)
        old = build_domain(D, sr)
        new = build_domain(ug)
        @test ug.uniform
        @test ug.nr == old.nr
        @test new.dr === old.dr
        @test new.Laplacian isa Tridiagonal
        @test new.Laplacian == old.Laplacian
        @test new.IntMat == old.IntMat
        @test new.r == [old.dr * i for i in 0:(old.nr - 1)]
        @test wall_radius(ug) ≈ D / 2

        p = SimulationParams(bathDiameter = 20.0, spatialResolution = 5.0)
        dtn = zeros(build_domain(20.0, 5.0).nr, build_domain(20.0, 5.0).nr)
        pr0 = build_problem(p, dtn)
        pr1 = build_problem(with_overrides(p; grid = uniform_grid(5.0, 10.0)), dtn)
        @test pr1.laplacian == pr0.laplacian
        @test pr1.dr === pr0.dr
        @test pr1.nr == pr0.nr
        @test pr1.contact_weights == pr0.contact_weights
        @test pr1.r == pr0.r
        @test_throws ArgumentError build_problem(with_overrides(p; grid = uniform_grid(6.0, 10.0)), dtn)
        # a default (grid = nothing) run keeps its configuration digest unchanged by the new field: no `grid` line is hashed
        @test config_digest(p) == config_digest(with_overrides(p; grid = nothing))
        @test config_digest(p) != config_digest(with_overrides(p; grid = uniform_grid(5.0, 10.0)))
    end

    @testset "map: identity near the disk, growth and cap, wall" begin
        sr, R, q, dmax, rin = 25.0, 40.0, 1.04, 0.5, 2.0
        g = graded_grid(sr, R; q = q, dr_max = dmax, r_in = rin)
        dr0 = 1 / sr
        @test !g.uniform
        @test g.dr0 == dr0
        @test length(g.r) == g.nr + 4 && length(g.s) == g.nr + 3 && length(g.sp) == g.nr + 3
        @test g.r[1] == 0.0
        @test all(g.r[i + 1] == dr0 * i for i in 0:g.xin)                         # exactly today's nodes up to r_in
        @test g.xin * dr0 >= rin
        steps = diff(g.r)
        @test all(steps .> 0)
        @test all(isapprox.(steps[1:g.xin], dr0; rtol = 1e-12))
        @test maximum(steps) <= dmax * (1 + 1e-12)                                  # cap
        @test steps[end] > 0.95 * dmax                                              # and the cap is reached
        @test all(steps[2:end] ./ steps[1:(end - 1)] .<= q + 1e-9)                  # growth per node
        @test wall_radius(g) >= R && g.r[g.nr] < R
        @test g.nr < 0.4 * ceil(Int, R * sr)                                        # far fewer nodes than the uniform grid
        # s = dr/dxi and s' = ds/dxi agree with finite differences of the nodes
        for i in (g.xin + 5, g.xin + 30, g.nr ÷ 2, g.nr - 3)
            @test isapprox((g.r[i + 2] - g.r[i]) / 2, g.s[i + 1]; rtol = 2e-3)
            @test isapprox((g.s[i + 2] - g.s[i]) / 2, g.sp[i + 1]; rtol = 2e-2, atol = 5e-4 * g.s[i + 1])   # measured in the prototype: 1e-4 s at the start of the ramp
        end
        # degenerate requests give the uniform grid; bad ones throw
        @test graded_grid(sr, R; dr_max = dr0).uniform
        @test graded_grid(sr, R; dr_max = dmax, q = 1.0).uniform
        @test_throws ArgumentError graded_grid(sr, R; dr_max = dmax, r_in = 1.0)
        @test_throws ArgumentError graded_grid(sr, R; dr_max = dmax, laplacian_order = 3)
        let lam = capillary_wavelength_cm(50.0, 1.0, 72.0, 981.0), k = 2pi / lam
            @test 981.0 * k + 72.0 * k^3 ≈ (2pi * 50.0)^2 rtol = 1e-9
        end
        @test boundary_node_count(uniform_grid(sr, R)) == ceil(Int, 2 * sr)
        @test boundary_node_count(g) >= 3
    end

    @testset "Laplacian: contact region unchanged, constants, tridiagonality" begin
        g = graded_grid(25.0, 40.0; q = 1.04, dr_max = 0.5, r_in = 2.0)
        L = build_domain(g).Laplacian
        @test L isa Tridiagonal
        old = build_domain(100.0, 25.0).Laplacian                    # dr = 0.04 as dr0: same rows inside r_in
        for i in 1:g.xin
            @test L[i, i] == old[i, i]
            @test L[i, i + 1] == old[i, i + 1]
            i > 1 && @test L[i, i - 1] == old[i, i - 1]
        end
        Lones = L * ones(g.nr)
        @test all(abs.(Lones[1:(g.nr - 1)]) .<= 1e-8 * maximum(abs, L))      # the last row feels the wall
        for order in (4, 6)
            g2 = graded_grid(25.0, 40.0; q = 1.04, dr_max = 0.5, r_in = 2.0, laplacian_order = order)
            L2 = build_domain(g2).Laplacian
            @test L2 isa Matrix{Float64}
            @test all(isapprox.((L2 * ones(g2.nr))[1:(g2.nr - 4)], 0.0; atol = 1e-8 * maximum(abs, L2)))
        end
    end

    @testset "Laplacian order under refinement (C2 of docs/graded_grid.md)" begin
        errs = Dict(form => [_gg_laplacian_error(p, form) for p in (1, 2, 4)] for form in (:cons2, :fd4, :fd6))
        ord(form) = log2(errs[form][2] / errs[form][3])            # between p = 2 and p = 4
        @test ord(:cons2) >= 1.9                                   # measured in a Python prototype: 1.99
        @test ord(:fd4) >= 3.5                                     # prototype: 3.88
        @test errs[:fd6][3] < errs[:fd4][3]                        # 6th order is more accurate than 4th (order itself not asserted)
        @info "Laplacian errors, p = 1, 2, 4" errs orders = (cons2 = ord(:cons2), fd4 = ord(:fd4), fd6 = ord(:fd6))
    end

    @testset "eta_boundary at a fixed physical radius" begin
        sr = 25.0
        ug = uniform_grid(sr, 20.0)
        eta = [sin(0.3 * i + 0.1 * j) for i in 1:ug.nr, j in 1:5]
        @test compute_eta_boundary(eta, sr; grid = ug) == compute_eta_boundary(eta, sr)         # uniform: bit-identical
        g = graded_grid(sr, 40.0; q = 1.04, dr_max = 1.5, r_in = 2.0)
        j0, w, radii = boundary_stencil(g)
        @test length(radii) == ceil(Int, 2 * sr)
        @test radii ≈ [wall_radius(g) - k / sr for k in 1:ceil(Int, 2 * sr)]
        @test all(sum(w; dims = 2) .≈ 1.0)
        # a cubic in r is reproduced exactly by the stencil: eta_boundary = max |p(x_k)|
        p(r) = 1.0 - 0.3r + 0.02r^2 - 0.0007r^3
        r = g.r[1:g.nr]
        eta = [c * p(ri) for ri in r, c in (1.0, -2.0, 0.5)]
        got = compute_eta_boundary(eta, sr; grid = g)
        @test got ≈ [abs(c) * maximum(abs.(p.(radii))) for c in (1.0, -2.0, 0.5)] rtol = 1e-9
    end
end
