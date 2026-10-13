@testset "lab_frame_CoM" begin
    dtn = generate_dtn(50, 20.0)
    f = 30.0
    p = SimulationParams(bathDiameter = 20.0, spatialResolution = 5.0, temporalResolution = 60,
                         bathAmplitude = 0.3 * 981 / (2pi * f)^2, bathFrequency = f,
                         simulationTime = 2 / f, earlyStop = false)
    r = run_simulation(p; dtn = dtn)

    @testset "no bath motion -> identical to the bath-frame trace" begin
        r0 = run_simulation(with_overrides(p; bathAmplitude = 0.0); dtn = dtn)
        @test lab_frame_CoM(r0) == r0.CoM_cm
    end

    @testset "bath displacement it adds is the one the solver's effective gravity implies" begin
        # g_eff = g*g_prefactor = g + (bath acceleration), so the second time-derivative of
        # (lab - bath-frame) must equal g*(g_prefactor - 1) at each step-end time. A sign or
        # phase error in lab_frame_CoM breaks this; a bathAmplitude=0 check never would.
        disp = lab_frame_CoM(r) .- r.CoM_cm
        dt = r.t_s[2] - r.t_s[1]
        problem = build_problem(p, dtn)
        for k in 2:(length(r.t_s) - 1)
            accel_fd = (disp[k + 1] - 2disp[k] + disp[k - 1]) / dt^2
            # t_s[k] is the end of global step k-2 (t_s[1] is t=0, before any step).
            accel_solver = p.g * (g_prefactor(problem, step_phase(k - 2, problem.stepsPerCycle), problem.dt) - 1)
            @test isapprox(accel_fd, accel_solver; rtol = 1e-2, atol = 1e-2 * p.bathAmplitude * (2pi * f)^2)
        end
    end
end
