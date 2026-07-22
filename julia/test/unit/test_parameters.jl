@testset "parameters" begin
    p = SimulationParams(diskRadius = 0.2, diskMass = 0.0283, forceFrequency = 90.0,
                          bathDensity = 1.0, bathSurfaceTension = 72.20, bathViscosity = 0.978e-2,
                          g = 981.0, bathDiameter = 20.0, spatialResolution = 5.0)

    dom = build_domain(p.bathDiameter, p.spatialResolution)
    dtn = zeros(dom.nr, dom.nr)
    problem = build_problem(p, dtn)

    T_unit = 1 / (p.forceFrequency * 2 * pi)
    L_unit = p.diskRadius

    @test problem.units.length ≈ L_unit
    @test problem.units.time ≈ T_unit
    @test problem.units.mass ≈ p.bathDensity * L_unit^3
    @test problem.units.velocity ≈ L_unit / T_unit
    @test problem.Fr ≈ L_unit / (p.g * T_unit^2)
    @test problem.We ≈ p.bathDensity * L_unit^3 / (p.bathSurfaceTension * T_unit^2)
    @test problem.Re ≈ L_unit^2 / (p.bathViscosity * T_unit)

    # force_frequency (freq_adim in the MATLAB source) is always exactly 1.0 by
    # construction: T_unit is defined as 1/forceFrequency, so the disk's own forcing
    # frequency is always the reference clock. This looks redundant but is intentional
    # (see parameters.jl docstring) — pin it so a future "simplification" doesn't
    # accidentally change what effective_w means when bath forcing is absent.
    @test problem.force_frequency ≈ 1.0

    @testset "zero bath forcing uses disk frequency as effective_w" begin
        p0 = with_overrides(p; bathAmplitude = 0.0)
        problem0 = build_problem(p0, dtn)
        @test problem0.bath_forcing_amplitude == 0.0
        @test problem0.effective_w ≈ problem0.force_frequency
    end

    @testset "with_overrides" begin
        p2 = with_overrides(p; forceAmplitude = 5.0, bathFrequency = 42.0)
        @test p2.forceAmplitude == 5.0
        @test p2.bathFrequency == 42.0
        @test p2.diskRadius == p.diskRadius
        @test p2.diskMass == p.diskMass
        @test_throws ArgumentError with_overrides(p; not_a_real_field = 1.0)
    end
end
