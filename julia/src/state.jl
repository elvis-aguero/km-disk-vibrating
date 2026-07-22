"""
Mutable per-step simulation state and the single-step advance function.

Deliberately separate from `Problem` (immutable, never mutates once a run starts) — see
`parameters.jl`'s docstring for why MATLAB's `PROBLEM_CONSTANTS`/`current_conditions`
split doesn't draw this line as cleanly.
"""
mutable struct SimulationState
    time::Float64
    step_counter::Int
    center_of_mass::Float64
    center_of_mass_velocity::Float64
    bath_surface::Vector{Float64}     # length nr
    bath_potential::Vector{Float64}   # length nr
    pressure::Vector{Float64}         # length cPoints
end

"""
    init_state(problem, params) -> SimulationState

Builds the initial condition: static equilibrium (see `solver/static_equilibrium.jl`) if
`params.startStatic`, else flat/zero — matching `solve_motion.m`'s `startStatic` branch.
The initial center-of-mass velocity is always `problem.v_bath_0` (the offset that makes
the *lab-frame* initial velocity zero when the bath itself is oscillating), regardless of
`startStatic`.
"""
function init_state(problem::Problem, params::SimulationParams)
    if params.startStatic
        eq = solve_static_equilibrium(problem)
        eta_initial, pressure_initial, z_initial = eq.eta, eq.pressure, eq.z
    else
        eta_initial = zeros(problem.nr)
        pressure_initial = zeros(problem.cPoints)
        z_initial = 0.0
    end
    phi_initial = zeros(problem.nr)
    return SimulationState(0.0, 0, z_initial, problem.v_bath_0, eta_initial, phi_initial, pressure_initial)
end

"""
Preallocated scratch buffers for `advance_one_step!`, reused every step so the
steady-state hot loop performs no heap allocation (`advance_one_step.m`'s MATLAB
equivalent reallocates `indep` via array concatenation every single step).
"""
struct StepBuffers
    rhs::Vector{Float64}
    sol::Vector{Float64}
end

StepBuffers(n::Int) = StepBuffers(zeros(n), zeros(n))

"""
    advance_one_step!(state, problem, solver, buffers) -> :ok | :diverged

Port of `advance_one_step.m`: assembles the RHS, solves the (LU-cached or GMRES) linear
system for this step, and unpacks the solution back into `state` in place. Returns
`:diverged` (instead of MATLAB's `error(...)`) if the solution contains a non-finite
value or `|CoM| > 10` (dimensionless) — the caller (`simulate.jl`) decides how to react,
so an expected/anticipated failure mode doesn't need exception-based control flow.
"""
function advance_one_step!(state::SimulationState, problem::Problem, solver::AbstractStepSolver, buffers::StepBuffers)
    dt = problem.dt
    phase = step_phase(state.step_counter, problem.stepsPerCycle)
    state.step_counter += 1

    g_pf = g_prefactor(problem, phase, dt)
    f_term = force_term(problem, phase, dt)

    assemble_rhs!(buffers.rhs, state.bath_surface, state.bath_potential,
                  state.center_of_mass_velocity, state.center_of_mass, dt, problem.Fr, g_pf, f_term)

    solve_step!(solver, phase, g_pf, buffers.rhs, buffers.sol)

    nr = problem.nr
    cPoints = problem.cPoints
    sol = buffers.sol
    z = sol[2 * nr + 2]

    @views state.bath_surface[1:cPoints] .= z
    @views state.bath_surface[(cPoints + 1):nr] .= sol[1:(nr - cPoints)]
    @views state.bath_potential .= sol[(nr - cPoints + 1):(2 * nr - cPoints)]
    @views state.pressure .= sol[(2 * nr - cPoints + 1):(2 * nr)]
    state.center_of_mass_velocity = sol[2 * nr + 1]
    state.center_of_mass = z
    state.time += dt

    diverged = abs(z) > 10 || !isfinite(z) || !isfinite(state.center_of_mass_velocity) || any(!isfinite, sol)
    return diverged ? :diverged : :ok
end
