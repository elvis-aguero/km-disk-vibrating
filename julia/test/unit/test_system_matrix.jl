@testset "system_matrix" begin
    p = SimulationParams(bathDiameter = 20.0, spatialResolution = 5.0, temporalResolution = 30)
    dom = build_domain(p.bathDiameter, p.spatialResolution)
    # Deterministic, non-trivial (not just zeros) stand-in DTN matrix — these tests exercise
    # the assembly/patch plumbing, not DTN correctness (that's test_dtn_golden_small.jl's job).
    dtn = [0.001 * (i + j) for i in 1:dom.nr, j in 1:dom.nr]
    problem = build_problem(p, dtn)

    @testset "materialize! matches a full rebuild for several g_prefactor values" begin
        tmpl = build_template(problem)
        n = 2 * problem.nr + 2
        Mat_patched = Matrix{Float64}(undef, n, n)
        for g in (0.0, 0.37, -1.2, 5.0)
            Mat_full = assemble_system_matrix(problem, g)
            materialize!(Mat_patched, tmpl, g)
            @test Mat_patched ≈ Mat_full
        end
    end

    @testset "pressure-under-disk block is exactly dt*I" begin
        Mat0 = assemble_system_matrix(problem, 0.0)
        nr, cPoints = problem.nr, problem.cPoints
        for k in 1:cPoints
            @test Mat0[nr + k, (2 * nr - cPoints) + k] ≈ problem.dt
        end
    end

    @testset "force-balance and kinematic rows have the expected fixed entries" begin
        Mat0 = assemble_system_matrix(problem, 0.0)
        nr = problem.nr
        @test Mat0[2 * nr + 1, 2 * nr + 1] ≈ 1.0
        @test Mat0[2 * nr + 1, 2 * nr + 2] ≈ problem.surface_force_constant * problem.dt / problem.dr
        @test Mat0[2 * nr + 2, 2 * nr + 1] ≈ -problem.dt
        @test Mat0[2 * nr + 2, 2 * nr + 2] ≈ 1.0
    end

    @testset "assemble_rhs! layout" begin
        nr = problem.nr
        rhs = zeros(2 * nr + 2)
        bath_surface = collect(1.0:nr)
        bath_potential = collect((nr + 1.0):(2 * nr))
        assemble_rhs!(rhs, bath_surface, bath_potential, 3.0, 7.0, problem.dt, problem.Fr, 0.5, 0.1)
        @test rhs[1:nr] == bath_surface
        @test rhs[(nr + 1):(2 * nr)] == bath_potential
        @test rhs[2 * nr + 1] ≈ 3.0 - problem.dt / problem.Fr * 0.5 - 0.1
        @test rhs[2 * nr + 2] ≈ 7.0
    end

    @testset "repatch! is bit-identical to materialize!" begin
        tmpl = build_template(problem)
        M_full = similar(tmpl.Mat0)
        M_patched = similar(tmpl.Mat0)
        materialize!(M_patched, tmpl, 1.0)
        for g in (0.5, 1.37, 1.0, 0.91, 1.5)
            materialize!(M_full, tmpl, g)
            repatch!(M_patched, tmpl, g)
            @test M_patched == M_full
        end
    end
end
