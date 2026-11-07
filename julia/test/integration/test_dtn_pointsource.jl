# DtN matrix against an exact solution with an ALGEBRAIC tail and nonzero mean (not an oscillating packet): a point source at height a above the
# surface, phi = (a - z)/(r^2 + (a - z)^2)^(3/2), f(r) = a/(r^2 + a^2)^(3/2), (N f)(r) = (2a^2 - r^2)/(r^2 + a^2)^(5/2).
# Error = max over nodes r <= R/2 of |A f - exact| / max|exact|. Bounds are 1.5 x the errors measured in job 8558665 (not looser).
@testset "DtN against a point-source solution (algebraic tail)" begin
    f(r, a) = a / (r^2 + a^2)^(3 / 2)
    Nf(r, a) = (2a^2 - r^2) / (r^2 + a^2)^(5 / 2)
    function glob_err(A, r, a)
        keep = r .<= maximum(r) / 2
        return maximum(abs.((A*f.(r, a))[keep] .- Nf.(r, a)[keep])) / maximum(abs.(Nf.(r, a)[keep]))
    end
    @testset "uniform grid, order 6, nr = 400" begin
        A = generate_dtn_single_layer(400, 800; order = 6)
        r = collect(0.0:399.0)
        @test glob_err(A, r, 10.0) <= 1.5 * 4.62e-4
        @test glob_err(A, r, 25.0) <= 1.5 * 1.77e-5
    end
    @testset "graded grid (defaults: ppw 32, wall zone 2), radii in disk radii" begin
        p = SimulationParams(bathDensity = 1.148, bathSurfaceTension = 68.0, bathViscosity = 0.079, g = 981.0, diskRadius = 0.07938,
                             bathDiameter = 700.0, spatialResolution = 25.0)
        g = graded_grid_for(p, 40.0)
        A = generate_dtn(g); r = collect(g.r[1:g.nr])
        @test glob_err(A, r, 1.0) <= 1.5 * 1.96e-5
        @test glob_err(A, r, 3.0) <= 1.5 * 2.17e-6
        @test glob_err(A, r, 10.0) <= 1.5 * 3.82e-6
    end
end
