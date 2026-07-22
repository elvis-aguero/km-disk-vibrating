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
end
