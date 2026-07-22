"""
    SimulationResult

Everything a single `run_simulation` call produces, in physical units (matching the
`nargout>0` return values of `solve_motion.m`: `t_s`, `CoM_cm`, `eta_history_cm`), plus
metadata (`resolved_solver`, the final `convergence` state, wall-clock `elapsed_s`, and
the LU-cache hit rate — `NaN` if GMRES was used — for the cache-reuse invariant check).
"""
struct SimulationResult
    params::SimulationParams
    resolved_solver::Symbol
    status::Symbol                  # :ok | :diverged
    t_s::Vector{Float64}
    CoM_cm::Vector{Float64}
    eta_history_cm::Matrix{Float64} # nr x nsteps
    convergence::ConvergenceResult
    elapsed_s::Float64
    lu_cache_hit_rate::Float64
end

"""
    run_simulation(params; dtn_registry=default_registry(), dtn=nothing, ram_budget_bytes=nothing,
                   on_step=nothing) -> SimulationResult

Runs one simulation to completion (or until early-stop convergence, or divergence).
Entry point equivalent to `solve_motion.m`. `ram_budget_bytes` defaults to
`Sys.free_memory()`; a sweep orchestrator running many cases concurrently should pass a
per-case share of that budget instead (see `sweep.jl`).

If `dtn` is provided directly (an already-loaded matrix), it is used as-is and
`dtn_registry` is never consulted — this lets a sweep over a fixed `(spatialResolution,
bathDiameter)` load the (potentially tens-of-MB) DTN matrix once and share it read-only
across every case, instead of every case re-reading it from disk.

`on_step`, if given, is called as `on_step(state, problem, step_count)` after every
completed step (including the one that triggers divergence, but not called for the
initial pre-loop state). This is the *only* hook for per-step visualization — MATLAB's
per-step `drawnow` debug plot has no direct equivalent baked into the loop itself, since
that would force a plotting dependency onto every user of this function. Instead, an
opt-in live-plotting callback lives entirely outside the core package in
`scripts/live_plot.jl` (`make_live_plot_callback`) and is passed in through this
parameter; `run_sweep`/`run_sweep_case` never set it, since a multithreaded sweep must
never try to open plot windows from multiple threads at once.

The internal time-stepping loop, and the convergence check applied to it, both operate
entirely in the solver's dimensionless units (matching `problem.effective_w`) — the
result is only converted to physical units (seconds, cm) once, at the end, for output.
"""
function run_simulation(params::SimulationParams; dtn_registry::Union{Nothing,DTNManifest} = nothing,
                         dtn::Union{Nothing,AbstractMatrix{<:Real}} = nothing,
                         ram_budget_bytes::Union{Nothing,Integer} = nothing,
                         on_step::Union{Nothing,Function} = nothing)
    t_wall_start = time()

    if dtn === nothing
        # Registry lookup (and the TOML file I/O it entails) is only attempted when no
        # DTN matrix was supplied directly — callers that already have one in hand (e.g.
        # a sweep sharing one matrix across cases, or a test using a freshly-generated
        # small domain) never touch the manifest at all.
        registry = dtn_registry === nothing ? default_registry() : dtn_registry
        entry = resolve_dtn(registry, params.spatialResolution, params.bathDiameter)
        dtn = load_dtn(entry)
    end
    problem = build_problem(params, dtn)
    state = init_state(problem, params)

    @info "starting simulation" diskRadius = params.diskRadius diskMass = params.diskMass forceAmplitude = params.forceAmplitude forceFrequency = params.forceFrequency bathAmplitude = params.bathAmplitude bathFrequency = params.bathFrequency nr = problem.nr solverType = params.solverType startStatic = params.startStatic

    if params.startStatic
        @info "static equilibrium solved" z_eq = state.center_of_mass z_eq_cm = state.center_of_mass * problem.units.length
    end

    solver, resolved_type = build_step_solver(problem, params.solverType, params.gmresTolerance;
                                                ram_budget_bytes = ram_budget_bytes)
    n = 2 * problem.nr + 2
    @info "solver path selected" requested = params.solverType resolved = resolved_type n = n stepsPerCycle = problem.stepsPerCycle
    buffers = StepBuffers(n)

    steps_estimate = max(1, ceil(Int, params.simulationTime / (problem.dt * problem.units.time)))
    time_hist = Float64[]
    com_hist = Float64[]
    eta_hist = Vector{Vector{Float64}}()
    sizehint!(time_hist, steps_estimate + 1)
    sizehint!(com_hist, steps_estimate + 1)
    sizehint!(eta_hist, steps_estimate + 1)

    push!(time_hist, state.time)
    push!(com_hist, state.center_of_mass)
    push!(eta_hist, copy(state.bath_surface))

    convergence_result = not_converged(0.0)
    status = :ok
    check_every_steps = max(1, round(Int, params.convergence.checkEveryPeriods * problem.stepsPerCycle))

    step_count = 0
    while state.time * problem.units.time < params.simulationTime
        outcome = advance_one_step!(state, problem, solver, buffers)
        step_count += 1
        push!(time_hist, state.time)
        push!(com_hist, state.center_of_mass)
        push!(eta_hist, copy(state.bath_surface))

        on_step === nothing || on_step(state, problem, step_count)

        if outcome == :diverged
            status = :diverged
            @error "numerical divergence detected" step = step_count com = state.center_of_mass
            break
        end

        if params.earlyStop && step_count % check_every_steps == 0
            convergence_result = check_convergence(time_hist, com_hist, problem.effective_w, params.convergence)
            if convergence_result.converged
                @info "converged" step = step_count periods = convergence_result.periods_elapsed amplitude = convergence_result.amplitude phase_deg = convergence_result.phase_deg amp_rel_change = convergence_result.amplitude_rel_change phase_change_deg = convergence_result.phase_change_deg
                break
            end
        end
    end

    if status == :ok && params.simulationTime > 0 && step_count == 0
        @warn "simulation ran zero steps" simulationTime = params.simulationTime dt_seconds = problem.dt * problem.units.time
    elseif status == :ok && params.earlyStop && !convergence_result.converged
        @warn "reached simulationTime without meeting the convergence criterion" periods_elapsed = convergence_result.periods_elapsed
    end

    T_unit = problem.units.time
    L_unit = problem.units.length
    t_s = time_hist .* T_unit
    CoM_cm = com_hist .* L_unit
    eta_history_cm = isempty(eta_hist) ? zeros(problem.nr, 0) : reduce(hcat, eta_hist) .* L_unit

    hit_rate_value = solver isa LUStepSolver ? hit_rate(solver.cache) : NaN

    elapsed_s = time() - t_wall_start
    @info "finished simulation" status = status steps = step_count elapsed_s = elapsed_s

    return SimulationResult(params, resolved_type, status, t_s, CoM_cm, eta_history_cm,
                             convergence_result, elapsed_s, hit_rate_value)
end
