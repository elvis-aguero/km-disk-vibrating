using MAT

@testset "DTN caches" begin
    repo_root = normpath(joinpath(@__DIR__, "..", "..", ".."))
    legacy(folder, file) = Matrix{Float64}(MAT.matread(joinpath(repo_root, "matlab", "1_code", folder, file))["DTNnew345"])

    @testset "dtn_row is exactly a row of generate_dtn" begin
        M = generate_dtn(50, 20.0)
        @test all(dtn_row(50, 20.0, k) == M[k, :] for k in 1:50)
    end

    @testset "check_dtn rejects incomplete or non-finite matrices" begin
        M = generate_dtn(50, 20.0)
        @test check_dtn(M) === M
        M_partial = copy(M); M_partial[30:end, :] .= 0
        @test_throws ErrorException check_dtn(M_partial)
        M_nan = copy(M); M_nan[5, 7] = NaN
        @test_throws ErrorException check_dtn(M_nan)
    end

    @testset "production cache D50Quant100 (nr=2500) matches the generator on spot rows" begin
        # Full regeneration takes hours; rows are independent, so check two: row 1 pins the
        # domain size (every row scales as 1/dr), row nr/2 exercises the full far-field
        # quadrature (~40 s).
        M = legacy("D50Quant100", "DTNnew345nr2500D100refp10.mat")
        @test check_dtn(M) === M
        for k in (1, 1250)
            ref = dtn_row(2500, 100.0, k)
            @test maximum(abs, ref - M[k, :]) / maximum(abs, ref) < 1e-10
        end
    end

    @testset "legacy D25Quant200 cache is incomplete and is rejected" begin
        # Rows 564:2500 are all zero in the committed file. If someone regenerates it, this
        # test should be updated to spot-check it like D50Quant100 instead.
        @test_throws ErrorException check_dtn(legacy("D25Quant200", "DTNnew345nr2500D25refp10.mat"))
    end
end
