using MAT

@testset "DTN generator matches the legacy MATLAB cache (nr=50)" begin
    # This is the load-bearing test for the highest-risk file in this port
    # (src/dtn/dtn_generator.jl — see its module docstring for the full story). It used to
    # show a "2.15% unexplained residual" at D=5 (and 75% at D=20, the value the cache's
    # folder name `D5Quant20` implies) — that gap turned out to be a real bug (`far_part!`
    # was missing a `/refp` division, silently zeroing almost every far-off-diagonal DTN
    # entry), not a property of the legacy cache. With that fixed, D=5 matches to machine
    # precision (~2e-15 relative difference) and D=20 still disagrees by 75% — i.e. this
    # cache genuinely was generated at D=5, consistent with `solve_motion.m`'s own
    # "Machine-specific patch for D5Quant20" (see scripts/migrate_dtn_caches.jl's header for
    # the full explanation of why D=5 is correct despite the D5Quant20 folder name implying
    # bathDiameter=20).
    #
    # Still checks both D=20 and D=5 (rather than hardcoding D=5) so a future regression in
    # either direction is caught and reported with the actual numbers, not just "test
    # failed". The tolerance is now tight (1e-6, not the old 5% headroom) since machine
    # precision is what's actually expected once the domain is right.
    repo_root = normpath(joinpath(@__DIR__, "..", "..", ".."))
    mat_path = joinpath(repo_root, "matlab", "1_code", "D5Quant20", "DTNnew345nr50D5refp10.mat")

    if !isfile(mat_path)
        @test_skip "legacy DTN cache not found at $mat_path (expected if the repo layout changed)"
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
        @info "best-matching D for the legacy D5Quant20 cache" D = best.D rel_diff = best.rel_diff

        @test best.D == 5.0
        @test size(generate_dtn(50, best.D)) == size(DTN_legacy)
        @test best.rel_diff < 1e-6
    end
end
