"""
Physically-grounded integration tests, run on a small, fast domain (nr=50) so the whole
group runs in seconds, not minutes. These test properties expected from the physics
itself, not frozen "the answer was N last time" golden numbers — this suite was written
without the ability to run Julia to generate golden reference output (see the repo-level
notes on how this port was authored), so correctness has to be established this way for
anything beyond the DTN-generator comparison (`test_dtn_golden_small.jl`) and the
synthetic-signal convergence-check tests (`unit/test_convergence_check.jl`).

A shared DTN matrix (bathDiameter=20, spatialResolution=5 -> nr=50) is generated once and
reused across every case in this file to avoid recomputing it repeatedly.
"""

const _PG_DTN = generate_dtn(50, 20.0)

_mean(x) = sum(x) / length(x)

"""
    naive_dft_peak_frequency(t, x) -> Float64

Direct (O(n^2), no FFTW dependency) discrete Fourier transform: returns the positive
frequency (Hz, excluding DC) with maximum power in a mean-subtracted signal `x` sampled
at times `t` (assumed uniformly spaced). Used to check that the simulated center-of-mass
trace actually oscillates at the driving frequency, without needing to add and verify an
FFT library in an environment where Julia itself could not be run to test it.
"""
function naive_dft_peak_frequency(t::AbstractVector{Float64}, x::AbstractVector{Float64})
    n = length(t)
    dt = t[2] - t[1]
    fs = 1 / dt
    xd = x .- _mean(x)
    n_bins = n ÷ 2
    power = zeros(n_bins)
    for k in 1:n_bins
        f = k * fs / n
        omega = 2 * pi * f
        c = 0.0
        s = 0.0
        @inbounds for j in 1:n
            c += xd[j] * cos(omega * t[j])
            s += xd[j] * sin(omega * t[j])
        end
        power[k] = c^2 + s^2
    end
    best_k = argmax(power)
    return best_k * fs / n
end

@testset "physically grounded" begin

    @testset "1: CoM oscillates at the forcing frequency (peak of a direct DFT)" begin
        freqHz = 20.0
        p = SimulationParams(bathDiameter = 20.0, spatialResolution = 5.0, temporalResolution = 48,
                              forceAmplitude = 0.1 * 0.0283 * 981, forceFrequency = freqHz, bathAmplitude = 0.0,
                              bathFrequency = freqHz, simulationTime = 15 / freqHz, startStatic = true,
                              earlyStop = false, solverType = :lu)
        result = run_simulation(p; dtn = _PG_DTN)
        @test result.status == :ok

        peak_hz = naive_dft_peak_frequency(result.t_s, result.CoM_cm)
        @test isapprox(peak_hz, freqHz; rtol = 0.1)
    end

    @testset "2: CoM amplitude is not orders of magnitude larger than the bath wave amplitude" begin
        freqHz = 20.0
        gamma_like = 0.1
        p = SimulationParams(bathDiameter = 20.0, spatialResolution = 5.0, temporalResolution = 48,
                              forceAmplitude = gamma_like * 0.0283 * 981, forceFrequency = freqHz, bathAmplitude = 0.0,
                              bathFrequency = freqHz, simulationTime = 15 / freqHz, startStatic = true,
                              earlyStop = false, solverType = :lu)
        result = run_simulation(p; dtn = _PG_DTN)
        @test result.status == :ok

        T = 1 / freqHz
        tail = result.t_s .>= (result.t_s[end] - 3 * T)
        com_fit = fit_oscillation(result.t_s[tail], result.CoM_cm[tail], 2 * pi * freqHz)
        wave_amp = maximum(abs, @view result.eta_history_cm[:, tail])

        @test com_fit.amplitude > 0
        @test com_fit.amplitude < 10 * wave_amp
    end

    @testset "3: static equilibrium is a fixed point when nothing forces the system" begin
        p = SimulationParams(bathDiameter = 20.0, spatialResolution = 5.0, temporalResolution = 40,
                              forceAmplitude = 0.0, bathAmplitude = 0.0, startStatic = true,
                              earlyStop = false, simulationTime = 5 / 90, solverType = :lu)
        result = run_simulation(p; dtn = _PG_DTN)
        @test result.status == :ok

        z_eq = result.CoM_cm[1]
        @test all(x -> isapprox(x, z_eq; atol = 1e-6), result.CoM_cm)
    end

    @testset "4: disk relaxes toward equilibrium when started away from it (damping, not anti-damping)" begin
        p_flat = SimulationParams(bathDiameter = 20.0, spatialResolution = 5.0, temporalResolution = 40,
                                   forceAmplitude = 0.0, bathAmplitude = 0.0, startStatic = false,
                                   earlyStop = false, simulationTime = 8 / 90, solverType = :lu)
        result_flat = run_simulation(p_flat; dtn = _PG_DTN)
        @test result_flat.status == :ok

        p_eq = with_overrides(p_flat; startStatic = true)
        result_eq = run_simulation(p_eq; dtn = _PG_DTN)
        z_eq = result_eq.CoM_cm[1]

        dist = abs.(result_flat.CoM_cm .- z_eq)
        n = length(dist)
        first_third = _mean(@view dist[1:(n ÷ 3)])
        last_third = _mean(@view dist[(2 * n ÷ 3 + 1):end])
        @test last_third < first_third
    end

    @testset "5: linear response — doubling forcing amplitude ~doubles steady-state CoM amplitude" begin
        freqHz = 25.0
        small_F = 0.02 * 0.0283 * 981
        p1 = SimulationParams(bathDiameter = 20.0, spatialResolution = 5.0, temporalResolution = 48,
                               forceAmplitude = small_F, forceFrequency = freqHz, bathAmplitude = 0.0,
                               bathFrequency = freqHz, simulationTime = 12 / freqHz, startStatic = true,
                               earlyStop = false, solverType = :lu)
        p2 = with_overrides(p1; forceAmplitude = 2 * small_F)

        r1 = run_simulation(p1; dtn = _PG_DTN)
        r2 = run_simulation(p2; dtn = _PG_DTN)
        @test r1.status == :ok && r2.status == :ok

        T = 1 / freqHz
        tail1 = r1.t_s .>= (r1.t_s[end] - 3 * T)
        tail2 = r2.t_s .>= (r2.t_s[end] - 3 * T)
        a1 = fit_oscillation(r1.t_s[tail1], r1.CoM_cm[tail1], 2 * pi * freqHz).amplitude
        a2 = fit_oscillation(r2.t_s[tail2], r2.CoM_cm[tail2], 2 * pi * freqHz).amplitude

        @test isapprox(a2 / a1, 2.0; rtol = 0.15)
    end

    @testset "6: CoM stays bounded (no blow-up) for a typical, stable forcing configuration" begin
        freqHz = 40.0
        p = SimulationParams(bathDiameter = 20.0, spatialResolution = 5.0, temporalResolution = 48,
                              forceAmplitude = 0.2 * 0.0283 * 981, forceFrequency = freqHz, bathAmplitude = 0.0,
                              bathFrequency = freqHz, simulationTime = 15 / freqHz, startStatic = true,
                              earlyStop = false, solverType = :lu)
        result = run_simulation(p; dtn = _PG_DTN)
        @test result.status == :ok
        @test maximum(abs, result.CoM_cm) < 1.0  # cm; disk radius is 0.2 cm, so this is a generous bound
    end

    @testset "7: the divergence guard actually triggers on a genuinely blown-up state" begin
        # Rather than hunting for physical parameters empirically known to cause blow-up
        # (which cannot be verified without a working Julia environment — see the
        # repo-level notes on how this port was authored), this exercises the real
        # advance_one_step! guard directly against a state already far past the
        # threshold, as if a previous (unshown) step had gone unstable.
        p = SimulationParams(bathDiameter = 20.0, spatialResolution = 5.0, temporalResolution = 40,
                              forceAmplitude = 0.0, bathAmplitude = 0.0, startStatic = true, solverType = :lu)
        problem = build_problem(p, _PG_DTN)
        state = init_state(problem, p)
        state.center_of_mass = 1e8
        state.center_of_mass_velocity = 1e8

        solver, _ = build_step_solver(problem, p.solverType, p.gmresTolerance)
        buffers = StepBuffers(2 * problem.nr + 2)
        outcome = advance_one_step!(state, problem, solver, buffers)
        @test outcome == :diverged
    end

    @testset "8: fast-test configurations don't trigger boundary-reflection contamination" begin
        # Meta-check on the fixtures used throughout this file: several tests above share
        # this small domain/short duration, and would be silently invalid if the outer
        # boundary were already reflecting waves back within that window. Uses the same
        # 10%-of-A_ref threshold as the validation/overlay tooling.
        freqHz = 30.0
        gamma_like = 0.1
        p = SimulationParams(bathDiameter = 20.0, spatialResolution = 5.0, temporalResolution = 48,
                              forceAmplitude = gamma_like * 0.0283 * 981, forceFrequency = freqHz, bathAmplitude = 0.0,
                              bathFrequency = freqHz, simulationTime = 15 / freqHz, startStatic = true,
                              earlyStop = false, solverType = :lu)
        result = run_simulation(p; dtn = _PG_DTN)
        @test result.status == :ok

        eta_boundary = compute_eta_boundary(result.eta_history_cm, p.spatialResolution)
        omega = 2 * pi * freqHz
        A_ref = gamma_like * 981.0 / omega^2
        @test maximum(eta_boundary) < 0.10 * A_ref
    end
end
