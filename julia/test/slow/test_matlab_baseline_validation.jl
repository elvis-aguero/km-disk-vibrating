@testset "MATLAB-baseline validation (full 30-case sweep, slow)" begin
    # Requires the D50Quant100 DTN cache to be registered (run
    # scripts/migrate_dtn_caches.jl first if julia/data/dtn_cache/manifest.toml doesn't
    # have it yet). Gated behind KMDISK_RUN_SLOW_TESTS: this is a genuinely multi-hour
    # computation even multithreaded (30 cases, each historically ~200-1500s in MATLAB).
    registry = default_registry()
    if !haskey(registry.entries, (50.0, 100.0))
        @test_skip "D50Quant100 DTN cache not registered — run scripts/migrate_dtn_caches.jl first"
    else
        include(joinpath(@__DIR__, "..", "..", "scripts", "run_matlab_baseline_check.jl"))
        result = run_matlab_baseline_check()

        @info "MATLAB baseline validation result" amp_rmse = result.amp_rmse phase_rmse_deg = result.phase_rmse_deg

        @test result.pass_incident_guard  # must not reproduce the ~1.09 amp / ~285deg phase regression
        @test result.pass_primary          # within 50% slack of the known-good ~0.10 amp / ~5.7deg baseline
    end
end
