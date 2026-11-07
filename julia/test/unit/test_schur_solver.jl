using LinearAlgebra, Random

@testset "schur_solver" begin
    p = SimulationParams(bathDiameter = 20.0, spatialResolution = 5.0, temporalResolution = 30, timeScheme = :bdf2)
    dtn = generate_dtn(50, 20.0)
    problem = build_problem(p, dtn)
    nr = problem.nr
    n = 2nr + 2

    # Residual of x in the ORIGINAL (2nr+2) system. A correct reduction must satisfy the original equations as well as
    # a direct LU solve does; the bound is tied to the LU residual, not to a fixed loose number.
    relres(M, x, b) = norm(M * x - b) / norm(b)

    @testset "reduced solve = full solve, several g and right-hand sides (implicit coefficient dt and 2dt/3)" begin
        for dtc in (problem.dt, 2 * problem.dt / 3), g in (0.0, 0.37, 1.0, 1.63, -1.2)
            sys = build_schur_system(problem, dtc)
            F = lu(sys.S0 .+ g .* sys.G1)
            M = assemble_system_matrix(problem, g; dt = dtc)
            Mf = lu(M)
            for seed in 1:3
                b = randn(MersenneTwister(100seed + round(Int, 10g)), n)
                x_full = Mf \ b
                x_red = schur_solve!(sys, F, g, b, zeros(n))
                @test relres(M, x_red, b) <= max(10 * relres(M, x_full, b), 1e-14)
                @test norm(x_red - x_full) / norm(x_full) <= 1e-9
            end
        end
    end

    @testset "S(g) is affine in g: S(g) = S0 + g G1 equals the reduction built at that g" begin
        sys = build_schur_system(problem, problem.dt)
        # an independent g: solve with S0 + 0.8 G1 and compare to the full system at g = 0.8
        b = randn(MersenneTwister(7), n)
        M = assemble_system_matrix(problem, 0.8)
        x_red = schur_solve!(sys, lu(sys.S0 .+ 0.8 .* sys.G1), 0.8, b, zeros(n))
        @test relres(M, x_red, b) <= max(10 * relres(M, lu(M) \ b, b), 1e-14)
    end

    @testset "implicit-Euler start-up step of the Schur solver = dense start-up step" begin
        solver, resolved = build_step_solver(problem, :schur, 1e-12; ram_budget_bytes = typemax(Int), scheme = :bdf2)
        @test resolved === :schur
        dense, _ = build_step_solver(problem, :lu, 1e-12; ram_budget_bytes = typemax(Int), scheme = :bdf2)
        b = randn(MersenneTwister(3), n)
        xs = solve_startup!(solver, problem, 0.9, b, zeros(n))
        xd = solve_startup!(dense, problem, 0.9, b, zeros(n))
        @test norm(xs - xd) / norm(xd) <= 1e-9
    end

    @testset "a banded (non-tridiagonal) Laplacian is refused with a clear error" begin
        pg = SimulationParams(bathDiameter = 40.0, spatialResolution = 5.0, temporalResolution = 20)
        gg = graded_grid_for(pg, 40.0; ppw = 12.0, q = 1.05, laplacian_order = 4)
        prob4 = build_problem(with_overrides(pg; grid = gg), generate_dtn(gg))
        @test !(prob4.laplacian isa Tridiagonal)
        @test_throws ArgumentError build_schur_system(prob4, prob4.dt)
    end

    @testset "memory estimate: about one quarter of the full-system LU cache" begin
        r = estimated_schur_bytes(31, 8750) / estimated_lu_cache_bytes(31, 2 * 8750 + 2)
        @test 0.24 < r < 0.30
    end
end
