# The single-layer DtN construction (src/dtn/dtn_single_layer.jl, derivation in docs/dtn_single_layer.md).
# Reference for accuracy: the exact DtN of a Gaussian, from its Hankel transform (computed here by direct quadrature),
# the same kind of test as the paper's Appendix B.1. Reference for consistency: the paper's scheme (generate_dtn_polar).

# J0 by the periodic trapezoid rule of (1/π) ∫_0^π cos(x sin t) dt (spectrally accurate for N > |x|)
function _slt_besselj0(x)
    N = max(64, ceil(Int, 2 * abs(x)) + 64)
    return sum(cos(x * sin(pi * (i + 0.5) / N)) for i in 0:(N - 1)) / N
end

# (N f)(r) for f = exp(-r^2/(2 s^2)), unit node spacing: s^2 ∫_0^∞ k^2 exp(-k^2 s^2/2) J0(k r) dk
function _slt_exact_gaussian_dtn(s, rs)
    dk = 0.4 / maximum(rs)
    ks = 0:dk:(12 / s)
    return [s^2 * dk * sum((ik == 1 ? 0.5 : 1.0) * k^2 * exp(-k^2 * s^2 / 2) * _slt_besselj0(k * r) for (ik, k) in enumerate(ks)) for r in rs]
end

@testset "single-layer DtN construction" begin
    nr, s = 100, 12.0                      # envelope 12 nodes wide, negligible at the last node
    rs = collect(0.0:(nr - 1))
    f = @. exp(-rs^2 / (2 * s^2))
    exact = _slt_exact_gaussian_dtn(s, rs)
    scale = maximum(abs, exact)

    @testset "accurate on a smooth function (both finite-difference orders)" begin
        for order in (4, 6)
            A = generate_dtn_single_layer(nr, 2.0 * nr; order = order)     # D = 2 nr: node spacing 1
            @test size(A) == (nr, nr)
            @test check_dtn(A) === A
            @test maximum(abs, A * f .- exact) / scale < 1e-4              # measured: 2e-5 (docs/dtn_single_layer.md, section 5)
        end
    end

    @testset "agrees with the paper's scheme on a smooth function" begin
        Apolar = generate_dtn_polar(nr, 2.0 * nr)
        A = generate_dtn(nr, 2.0 * nr)
        @test maximum(abs, A * f .- Apolar * f) / scale < 1e-4
    end

    @testset "scales as 1/dr, like the paper's matrix" begin
        A1 = generate_dtn_single_layer(40, 80.0)
        A2 = generate_dtn_single_layer(40, 160.0)
        @test maximum(abs, A1 .- 2 .* A2) <= 1e-12 * maximum(abs, A1)
    end

    @testset "method dispatch" begin
        @test generate_dtn(40, 80.0) == generate_dtn_single_layer(40, 80.0)
        @test generate_dtn(40, 80.0; method = :polar) == generate_dtn_polar(40, 80.0)
        @test_throws ArgumentError generate_dtn(40, 80.0; method = :nope)
        @test_throws ArgumentError generate_dtn_single_layer(40, 80.0; order = 5)
    end

    @testset "registry keeps both methods per domain" begin
        path = tempname() * ".toml"
        m = DTNManifest()
        register!(m, DTNCacheEntry(5, 20, 50, 10, "a.jld2", "polar test"))
        register!(m, DTNCacheEntry(5, 20, 50, 10, "b.jld2", "single-layer test", :single_layer))
        save_manifest(path, m)
        m2 = load_manifest(path)
        @test length(m2.entries) == 2
        @test basename(resolve_dtn(m2, 5, 20; method = :polar).path) == "a.jld2"
        @test basename(resolve_dtn(m2, 5, 20; method = :single_layer).path) == "b.jld2"
        @test_throws ErrorException resolve_dtn(m2, 5, 99; method = :single_layer)
        # a manifest entry without a `method` key (written before the field existed) is a polar entry
        write(path, "[[entry]]\nspatialResolution = 5\nbathDiameter = 20\nnr = 50\nrefp = 10\npath = \"c.jld2\"\nsource_file = \"\"\n")
        @test only(values(load_manifest(path).entries)).method === :polar
        rm(path)
    end
end
