#= Physics checks against limits known independently of this code: Archimedes buoyancy,
the quasi-static (slow shaking) limit, and convergence under refinement. Sweep fluid,
small domains. =#
@testset "physics limits" begin
    fluid = (bathDensity = 1.175, bathSurfaceTension = 66.5, bathViscosity = 0.18 / 1.175)
    labfit(r, f; periods = 3) = begin
        tail = r.t_s .>= r.t_s[end] - periods / f
        fit = fit_oscillation(r.t_s[tail], lab_frame_CoM(r)[tail], 2π * f)
        (ratio = fit.amplitude / r.params.bathAmplitude,
         lag = mod(rad2deg(deg2rad(r.params.phaseDifference) - fit.phase) + 180, 360) - 180)
    end

    @testset "no surface tension: the disk floats at the Archimedes depth" begin
        # σ → 0 leaves buoyancy alone to carry the weight: depth h = m / (ρ π R²). The solver
        # takes the contact pressure at the midpoint of the last two grid points (a deliberate
        # ghost-point assumption, see domain.jl), so a uniform pressure integrates to
        # π R² (1 + 1/(2q) + 1/(12q²)) and the computed depth is h divided by that factor —
        # checked to 1e-6 here; the difference must vanish as q grows.
        errs = Float64[]
        for q in (5, 10, 20)
            p = SimulationParams(; fluid..., bathSurfaceTension = 1e-8, bathDiameter = 20.0, spatialResolution = q)
            problem = build_problem(p, generate_dtn(10q, 20.0))
            z = FaradayDisk.solve_static_equilibrium(problem).z * problem.units.length
            h = p.diskMass / (p.bathDensity * π * p.diskRadius^2)
            @test isapprox(-z * (1 + 1 / (2q) + 1 / (12q^2)), h; rtol = 1e-6)
            push!(errs, abs(z + h) / h)
        end
        @test errs[1] > errs[2] > errs[3]
    end

    @testset "slow shaking: the disk rides the bath (quasi-static limit)" begin
        # Far below every wave frequency of the domain, inertia is negligible: in the lab
        # frame the disk moves with the bath (amplitude ratio → 1, lag → 0), and the leftover
        # deviation is inertial, so it scales as f² (×4 per doubling of f).
        dtn = generate_dtn(150, 30.0)
        dev = Dict{Float64,Float64}()
        for f in (0.5, 1.0, 2.0)
            p = SimulationParams(; fluid..., bathFrequency = f, bathAmplitude = 0.5 * 981 / (2π * f)^2,
                                 simulationTime = 20 / f, earlyStop = false, bathDiameter = 30.0, spatialResolution = 10.0)
            r = run_simulation(p; dtn = dtn)
            @test r.status == :ok
            fit = labfit(r, f; periods = 1)
            dev[f] = fit.ratio - 1
            f == 0.5 && @test abs(fit.ratio - 1) < 1e-3 && abs(fit.lag) < 0.1
            if f == 0.5
                # Relative to the bath, the disk follows the static equilibrium for the
                # instantaneous effective gravity g(1 - Γ cos(ωt + φ)). Agreement is ~4% of
                # that swing, not exact: the viscous surface terms make the dynamic model's
                # rest state differ slightly from the (inviscid) static solve.
                tail = r.t_s .>= r.t_s[end] - 1 / f
                z_qs = map(r.t_s[tail]) do t
                    geff = p.g * (1 - 0.5 * cos(2π * f * t + deg2rad(p.phaseDifference)))
                    pq = build_problem(with_overrides(p; g = geff), dtn)
                    FaradayDisk.solve_static_equilibrium(pq).z * pq.units.length
                end
                swing = maximum(z_qs) - minimum(z_qs)
                @test maximum(abs, r.CoM_cm[tail] - z_qs) < 0.06 * swing
            end
        end
        @test 3 < dev[2.0] / dev[1.0] < 5
        @test 3 < dev[1.0] / dev[0.5] < 5
    end

    @testset "refinement: results converge as the grid and time step shrink" begin
        f, Γ = 30.0, 0.2
        run(q, spc, scheme) = labfit(run_simulation(
            SimulationParams(; fluid..., bathFrequency = f, bathAmplitude = Γ * 981 / (2π * f)^2,
                             temporalResolution = spc, simulationTime = 25 / f, earlyStop = false,
                             bathDiameter = 20.0, spatialResolution = q, timeScheme = scheme);
            dtn = generate_dtn(10q, 20.0)), f)
        spcs = (60, 120, 240, 480)

        # Implicit Euler: first order — successive changes shrink, about 2x per halving of dt.
        e = [run(5, spc, :euler) for spc in spcs]
        de_amp = abs.(diff([x.ratio for x in e]))
        de_lag = abs.(diff([x.lag for x in e]))
        @test issorted(de_amp; rev = true) && issorted(de_lag; rev = true)

        # BDF2: second order — phase changes shrink ≥3x per halving (4x asymptotically), and
        # the amplitude is converged to <0.1% already at 60 steps/period.
        b = [run(5, spc, :bdf2) for spc in spcs]
        db_lag = abs.(diff([x.lag for x in b]))
        @test all(db_lag[1:end-1] ./ db_lag[2:end] .> 3)
        @test abs(b[1].ratio - b[end].ratio) < 1e-3 * b[end].ratio

        # Both converge to the same answer: Richardson-extrapolating Euler's first-order error
        # (2·E(dt/2) − E(dt)) lands on the BDF2 result.
        @test isapprox(2 * e[end].ratio - e[end-1].ratio, b[end].ratio; rtol = 5e-3)
        @test isapprox(2 * e[end].lag - e[end-1].lag, b[end].lag; atol = 0.1)

        # Grid: at a fixed time step, amplitude is already resolved to <1% at q = 5, and
        # the phase changes shrink with refinement.
        s = [run(q, 60, :bdf2) for q in (5, 10, 20)]
        @test maximum(x.ratio for x in s) - minimum(x.ratio for x in s) < 0.01
        @test abs(s[3].lag - s[2].lag) < abs(s[2].lag - s[1].lag)
    end
end
