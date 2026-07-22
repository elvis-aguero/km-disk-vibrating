using MAT

@testset "DTN generator matches the legacy MATLAB cache (nr=50, D=20)" begin
    # This is the load-bearing test for the highest-risk file in this port
    # (src/dtn/dtn_generator.jl — see its module docstring): it was transcribed from
    # DTNVectorized.m without the ability to run either MATLAB or Julia to check it.
    repo_root = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
    mat_path = joinpath(repo_root, "matlab", "1_code", "D5Quant20", "DTNnew345nr50D5refp10.mat")

    if !isfile(mat_path)
        @test_skip "legacy DTN cache not found at $mat_path (expected if the repo layout changed)"
    else
        DTN_legacy = Matrix{Float64}(MAT.matread(mat_path)["DTNnew345"])
        DTN_native = generate_dtn(50, 20.0)

        @test size(DTN_native) == size(DTN_legacy)
        max_abs_diff = maximum(abs, DTN_native .- DTN_legacy)
        rel_diff = max_abs_diff / maximum(abs, DTN_legacy)
        @info "DTN golden comparison (nr=50, D=20)" max_abs_diff = max_abs_diff rel_diff = rel_diff

        # Loose-ish tolerance: batched (MATLAB) vs. unbatched (Julia) angular quadrature
        # summation order differs, so exact bitwise equality isn't expected — but a real
        # transcription bug should produce an error many orders of magnitude above this.
        @test rel_diff < 1e-6
    end
end
