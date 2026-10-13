@testset "solver internals" begin
    fixed = SimulationParams(bathDiameter = 20.0, spatialResolution = 5.0, temporalResolution = 20,
                              forceAmplitude = 0.1 * 0.0283 * 981, forceFrequency = 90.0, bathAmplitude = 0.0,
                              bathFrequency = 90.0, simulationTime = 3 / 90, startStatic = true,
                              earlyStop = false)
    dtn = generate_dtn(50, 20.0)

    @testset "LU cache is actually reused (hit rate) rather than silently refactorizing every step" begin
        # Not a correctness invariant but a real, easy-to-miss performance regression: the
        # whole point of the LU-caching architecture is that only the first cycle's steps
        # are cache misses. 3 periods at 20 steps/period = 60 steps; the first 20 are
        # necessarily misses, so hit rate should be comfortably above 0.5.
        p = with_overrides(fixed; solverType = :lu)
        result = run_simulation(p; dtn = dtn)
        @test result.resolved_solver == :lu
        @test !isnan(result.lu_cache_hit_rate)
        @test result.lu_cache_hit_rate > 0.5
    end

    @testset "LU-cached and GMRES solvers agree on the same problem" begin
        # Two independent solution methods for the same linear system must agree —
        # catches bugs specific to either path (e.g. a caching/indexing bug in the LU
        # path, since GMRES rebuilds the matrix fresh every step with no caching at all).
        p_lu = with_overrides(fixed; solverType = :lu)
        p_gmres = with_overrides(fixed; solverType = :gmres, gmresTolerance = 1e-12)
        r_lu = run_simulation(p_lu; dtn = dtn)
        r_gmres = run_simulation(p_gmres; dtn = dtn)

        @test r_lu.resolved_solver == :lu
        @test r_gmres.resolved_solver == :gmres
        @test length(r_lu.CoM_cm) == length(r_gmres.CoM_cm)
        @test isapprox(r_lu.CoM_cm, r_gmres.CoM_cm; rtol = 1e-6, atol = 1e-9)
    end

    # Bath-driven fixture: g_prefactor sweeps 1 ± gamma over each cycle, so (unlike the
    # disk-forced `fixed` above, where g ≡ 1) these exercise both the LU-cache dedup and
    # the GMRES preconditioner on genuinely different per-step matrices.
    bath = with_overrides(fixed; forceAmplitude = 0.0, bathAmplitude = 0.5 * 981 / (2pi * 30)^2,
                          bathFrequency = 30.0, phaseDifference = -90.0, simulationTime = 3 / 30)

    function run_steps(p, solver_type, nsteps)
        problem = build_problem(p, dtn)
        solver, _ = build_step_solver(problem, solver_type, 1e-12; ram_budget_bytes = typemax(Int))
        state = init_state(problem, p)
        buffers = StepBuffers(2 * problem.nr + 2)
        for _ in 1:nsteps
            advance_one_step!(state, problem, solver, buffers)
        end
        return solver
    end

    @testset "LU cache shares one factorization between phase-mirrored steps" begin
        # g = 1 - gamma*cos(theta_k): at -90 deg phase, steps mirrored about the forcing
        # phase give the same matrix, so a 20-step cycle needs 20/2 + 1 = 11 distinct
        # factorizations. Exact-float matching found more (mirrored values differ by a
        # few eps); this pins the tolerant match.
        solver = run_steps(bath, :lu, 20)
        @test n_factorizations(solver.cache) == 11
        # ...and never merges genuinely different matrices: at a phase with no mirror
        # symmetry on the step grid, every one of the 20 steps is distinct.
        solver_37 = run_steps(with_overrides(bath; phaseDifference = 37.0), :lu, 20)
        @test n_factorizations(solver_37.cache) == 20
    end

    @testset "LU-preconditioned GMRES converges in a handful of iterations" begin
        # A performance invariant, not a correctness one: without the preconditioner GMRES
        # still converges (same answer) but needs ~n iterations per step.
        solver = run_steps(bath, :gmres, 40)
        @test 0 < mean_iterations(solver.gmres) <= 10
    end

    @testset "LU-cached and GMRES solvers agree on a bath-driven problem" begin
        r_lu = run_simulation(with_overrides(bath; solverType = :lu); dtn = dtn)
        r_gmres = run_simulation(with_overrides(bath; solverType = :gmres, gmresTolerance = 1e-12); dtn = dtn)
        @test r_lu.status == :ok && r_gmres.status == :ok
        @test length(r_lu.CoM_cm) == length(r_gmres.CoM_cm)
        @test isapprox(r_lu.CoM_cm, r_gmres.CoM_cm; rtol = 1e-8)
    end
end
