@testset "domain" begin
    D = 20.0
    quant = 5.0
    dom = build_domain(D, quant)

    expected_nr = ceil(Int, D * quant / 2)
    @test dom.nr == expected_nr
    @test dom.dr ≈ D / (2 * dom.nr)
    @test size(dom.Laplacian) == (dom.nr, dom.nr)

    @testset "Laplacian of a constant is zero away from the truncated outer boundary" begin
        # (1/r) d/dr(r d/dr) of a constant is zero everywhere it's well-defined. The
        # discretization preserves this at the origin (row 1's mirror-boundary stencil)
        # and every interior row; only the truncated outer edge (row nr, which has no
        # far-field neighbor to complete the stencil) is expected to deviate.
        ones_vec = ones(dom.nr)
        Lones = dom.Laplacian * ones_vec
        @test all(isapprox.(Lones[1:(dom.nr - 1)], 0.0; atol = 1e-9))
    end

    @testset "row 1 mirror-boundary stencil" begin
        @test dom.Laplacian[1, 1] ≈ -4.0 / dom.dr^2
        @test dom.Laplacian[1, 2] ≈ 4.0 / dom.dr^2
        dom.nr > 2 && @test dom.Laplacian[1, 3] == 0.0
    end

    @testset "interior row stencil matches the closed-form radial Laplacian weights" begin
        i = dom.nr ÷ 2  # an interior row, away from both boundaries
        expected_sub = (1.0 / dom.dr^2) * (1.0 - 1.0 / (2.0 * (i - 1)))
        expected_diag = -2.0 / dom.dr^2
        expected_super = (1.0 / dom.dr^2) * (1.0 + 1.0 / (2.0 * (i - 1)))
        @test dom.Laplacian[i, i - 1] ≈ expected_sub
        @test dom.Laplacian[i, i] ≈ expected_diag
        @test dom.Laplacian[i, i + 1] ≈ expected_super
    end

    @testset "pressure-integral matrix" begin
        nlmax = Int(quant) + 1
        @test size(dom.IntMat) == (2 * nlmax, 2 * nlmax)
        @test dom.IntMat[1, 1] ≈ pi * dom.dr^2 / 12
        @test dom.IntMat[2, 1] ≈ pi * dom.dr^2 / 3
        # Upper-triangular part (jj > ii) is untouched (stays zero) by construction.
        @test dom.IntMat[2, 3] == 0.0
    end
end
