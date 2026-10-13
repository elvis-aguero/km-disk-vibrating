# Top-level, guarded include: see test_plotting.jl.
if !@isdefined(OverlayCase)
    include(joinpath(@__DIR__, "..", "..", "scripts", "run_validation_overlay.jl"))
end

@testset "experiment comparison reproduces the committed reference run" begin
    # Re-scores the committed nr=2500 reference sweep (output/4193275_*) against the lab data
    # with the current overlay/fit code — no simulation, seconds. Catches regressions in
    # the fitting, lab-frame correction or matching logic, which otherwise only the
    # hours-long cluster run would reveal. Expected values are what that run logged.
    sweep_dir = joinpath(REPO_ROOT, "output", "4193275_20260722_214351", "sweep")
    sim = compute_sim_overlay(sweep_dir)
    @test length(sim) == 30
    stats = rmse_against_experiment(sim, joinpath(EXPERIMENTAL_DIR, "ampiezza_solo_misure_3.csv"),
                                    joinpath(EXPERIMENTAL_DIR, "fase_solo_misure_3.csv"))
    @test isapprox(stats.amp_rmse, 0.08482509815240952; rtol = 1e-9)
    @test isapprox(stats.phase_rmse_deg, 6.001592386470637; rtol = 1e-9)
    @test stats.n_amp > 0 && stats.n_phase > 0
end
