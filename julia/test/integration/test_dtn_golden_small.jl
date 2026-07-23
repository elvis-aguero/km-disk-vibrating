using MAT

@testset "DTN generator matches the MATLAB cache (nr=50)" begin
    # This is the load-bearing test for the highest-risk file in this port
    # (src/dtn/dtn_generator.jl — see its module docstring for the full story).
    #
    # History: this cache (matlab/1_code/D5Quant20/) was originally generated at the
    # literal argument D=5, not its folder-implied bathDiameter=20 — a real, if previously
    # benign (this domain is never used by any production/experimental sweep — confirmed by
    # audit), data bug. `scripts/regenerate_d5quant20_matlab_cache.jl` fixed it by
    # regenerating natively at the correct D=20 and removing solve_motion.m's
    # "Machine-specific patch for D5Quant20" special case. See migrate_dtn_caches.jl's
    # header for the full history.
    #
    # Still checks both D=20 and D=5 (rather than hardcoding D=20) so a future regression in
    # either direction is caught and reported with the actual numbers, not just "test
    # failed". Tolerance is tight (1e-6) since machine precision is what's actually expected
    # once the domain is right.
    repo_root = normpath(joinpath(@__DIR__, "..", "..", ".."))
    mat_path = joinpath(repo_root, "matlab", "1_code", "D5Quant20", "DTNnew345nr50D20refp10.mat")

    if !isfile(mat_path)
        @test_skip "MATLAB DTN cache not found at $mat_path (expected if the repo layout changed)"
    else
        DTN_legacy = Matrix{Float64}(MAT.matread(mat_path)["DTNnew345"])

        results = NamedTuple[]
        for D_candidate in (20.0, 5.0)
            DTN_native = generate_dtn(50, D_candidate)
            max_abs_diff = maximum(abs, DTN_native .- DTN_legacy)
            rel_diff = max_abs_diff / maximum(abs, DTN_legacy)
            @info "DTN golden comparison" D_candidate = D_candidate max_abs_diff = max_abs_diff rel_diff = rel_diff
            push!(results, (D = D_candidate, rel_diff = rel_diff))
        end

        best = results[argmin([r.rel_diff for r in results])]
        @info "best-matching D for the D5Quant20 cache" D = best.D rel_diff = best.rel_diff

        @test best.D == 20.0
        @test size(generate_dtn(50, best.D)) == size(DTN_legacy)
        @test best.rel_diff < 1e-6
    end
end
