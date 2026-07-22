@testset "fit_oscillation" begin
    @testset "recovers known amplitude/phase from a clean sinusoid" begin
        for (A, phi, f) in ((1.0, 0.3, 5.0), (2.5, -1.1, 12.0), (0.02, 2.9, 30.0))
            omega = 2 * pi * f
            t = collect(0.0:0.001:1.0)
            z = A .* cos.(omega .* t .- phi) .+ 0.7  # offset column should absorb the +0.7
            fit = fit_oscillation(t, z, omega)
            @test isapprox(fit.amplitude, A; rtol = 1e-6)
            @test isapprox(mod(fit.phase - phi + pi, 2pi) - pi, 0.0; atol = 1e-6)
            @test isapprox(fit.offset, 0.7; atol = 1e-6)
        end
    end

    @testset "robust to modest noise" begin
        A, phi, f = 1.0, 0.5, 8.0
        omega = 2 * pi * f
        t = collect(0.0:0.0005:2.0)
        noise = [0.01 * sin(97 * ti) for ti in t]  # deterministic pseudo-noise, no RNG dependency
        z = A .* cos.(omega .* t .- phi) .+ noise
        fit = fit_oscillation(t, z, omega)
        @test isapprox(fit.amplitude, A; atol = 0.02)
        @test isapprox(mod(fit.phase - phi + pi, 2pi) - pi, 0.0; atol = 0.02)
    end

    @testset "input validation" begin
        @test_throws ArgumentError fit_oscillation([1.0, 2.0], [1.0], 1.0)
        @test_throws ArgumentError fit_oscillation([1.0, 2.0], [1.0, 2.0], 1.0)
    end
end
