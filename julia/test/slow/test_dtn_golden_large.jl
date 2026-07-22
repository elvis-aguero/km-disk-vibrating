using MAT

@testset "DTN generator matches legacy MATLAB caches (nr=2500 domains, slow)" begin
    # Gated behind KMDISK_RUN_SLOW_TESTS: generating nr=2500 DTN matrices took MATLAB's
    # parfor-based generator on the order of hours historically for the largest domains;
    # even with Threads.jl this must never run on every push.
    repo_root = normpath(joinpath(@__DIR__, "..", "..", ".."))

    domains = [
        (spatialResolution = 50.0, bathDiameter = 100.0, nr = 2500,
         mat_path = joinpath(repo_root, "matlab", "1_code", "D50Quant100", "DTNnew345nr2500D100refp10.mat")),
        (spatialResolution = 25.0, bathDiameter = 200.0, nr = 2500,
         mat_path = joinpath(repo_root, "matlab", "1_code", "D25Quant200", "DTNnew345nr2500D25refp10.mat")),
    ]

    for d in domains
        @testset "nr=$(d.nr), bathDiameter=$(d.bathDiameter)" begin
            if !isfile(d.mat_path)
                @test_skip "legacy DTN cache not found at $(d.mat_path)"
                continue
            end
            DTN_legacy = Matrix{Float64}(MAT.matread(d.mat_path)["DTNnew345"])
            t0 = time()
            DTN_native = generate_dtn(d.nr, d.bathDiameter)
            elapsed = time() - t0
            rel_diff = maximum(abs, DTN_native .- DTN_legacy) / maximum(abs, DTN_legacy)
            @info "large DTN golden comparison" nr = d.nr bathDiameter = d.bathDiameter elapsed_s = elapsed rel_diff = rel_diff
            @test rel_diff < 1e-6
        end
    end
end
