using LinearAlgebra

@testset "schur solver end to end" begin
    fixed = SimulationParams(bathDiameter = 20.0, spatialResolution = 5.0, temporalResolution = 20,
                              forceAmplitude = 0.1 * 0.0283 * 981, forceFrequency = 90.0, bathAmplitude = 0.0,
                              bathFrequency = 90.0, simulationTime = 3 / 90, startStatic = true, earlyStop = false)
    bath = with_overrides(fixed; forceAmplitude = 0.0, bathAmplitude = 0.5 * 981 / (2pi * 30)^2,
                          bathFrequency = 30.0, phaseDifference = -90.0, simulationTime = 4 / 30)
    dtn = generate_dtn(50, 20.0)

    # Same trajectory as the dense LU path: the reduction is exact, so the only differences are rounding.
    for (name, base) in (("disk-forced", fixed), ("bath-driven", bath)), scheme in (:euler, :bdf2)
        @testset "$name, $scheme: :schur trajectory = :lu trajectory" begin
            r_lu = run_simulation(with_overrides(base; solverType = :lu, timeScheme = scheme); dtn = dtn)
            r_sc = run_simulation(with_overrides(base; solverType = :schur, timeScheme = scheme); dtn = dtn)
            @test r_lu.resolved_solver == :lu
            @test r_sc.resolved_solver == :schur
            @test length(r_lu.CoM_cm) == length(r_sc.CoM_cm)
            scale = maximum(abs, r_lu.CoM_cm)
            @test maximum(abs, r_sc.CoM_cm .- r_lu.CoM_cm) / scale <= 1e-9
            @test norm(r_sc.eta_history_cm - r_lu.eta_history_cm) / norm(r_lu.eta_history_cm) <= 1e-9
            @test r_sc.lu_cache_hit_rate > 0.5
        end
    end
end
