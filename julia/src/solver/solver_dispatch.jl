"""
Ties `system_matrix.jl`, `lu_cache.jl`, and `gmres_solver.jl` together: resolves
`solverType=:auto/:lu/:gmres` into a concrete step solver, and dispatches a single step's
solve through whichever one was chosen.

## The physics-step / cache-index split, made structurally safe

MATLAB's `advance_one_step.m` uses one integer (`current_step`, a monotonically
increasing counter) for two conceptually different purposes: (1) the *time-phase
argument* fed into the oscillating-gravity and disk-forcing cosines, which must be
continuous across cycle boundaries for the physics to be correct, and (2) the *cache
lookup key* for the LU library, which must instead *wrap* every `stepsPerCycle` steps
(`mod(current_step, stepsPerCycle)+1`), since the system matrix itself is exactly
periodic with that period. The git history of the MATLAB port shows this exact
conflation was once introduced as a bug and had to be fixed (commit e89f9bd, "fixed
boundary discontinuity by using global step for physics and cycle index for caching").

`StepPhase` makes the two uses distinct fields with different names, so a future edit
that tries to pass `cycle_index` into a physics function (or `global_step` into the cache
lookup) is a call-site type/name mismatch rather than a silent copy-paste bug.
"""
struct StepPhase
    global_step::Int   # monotonic; feeds g_prefactor / force_term (physical continuity)
    cycle_index::Int   # mod(global_step, stepsPerCycle) + 1; LU-cache lookup key ONLY
end

"""
    step_phase(global_step, stepsPerCycle) -> StepPhase

Constructs a `StepPhase` from a global step counter and the cycle length, computing
`cycle_index` from them. Deliberately a plain function, not a second `StepPhase(::Int,
::Int)` outer constructor — that would have the exact same argument signature as the
struct's own default (auto-generated) positional constructor `StepPhase(global_step,
cycle_index)` and silently overwrite it.
"""
step_phase(global_step::Int, stepsPerCycle::Int) =
    StepPhase(global_step, mod(global_step, stepsPerCycle) + 1)

"""
    g_prefactor(problem, phase, dt) -> Float64

The time-varying effective-gravity term `1 - gamma*cos(bath_w*(global_step+1)*dt +
phase)`. Uses `global_step+1` (the step-END time), matching `advance_one_step.m`'s
implicit-Euler-style evaluation of the forcing at the new time level, not the old one.
"""
g_prefactor(problem::Problem, phase::StepPhase, dt::Float64) =
    1.0 - problem.bath_forcing_amplitude *
    cos(problem.bath_frequency * (phase.global_step + 1) * dt + problem.phase_difference)

"""
    force_term(problem, phase, dt) -> Float64

The disk-forcing contribution to the RHS, `dt*F*cos(w*(global_step+1)*dt)`.
"""
force_term(problem::Problem, phase::StepPhase, dt::Float64) =
    dt * problem.force_amplitude * cos(problem.force_frequency * (phase.global_step + 1) * dt)

# --- solverType resolution ---

"""
    estimated_lu_cache_bytes(stepsPerCycle, n) -> Int

Estimated peak memory for a full LU cache of `stepsPerCycle` factorizations of an `n x n`
matrix. Julia's `LinearAlgebra.lu` stores one `n x n` factors matrix plus an `O(n)` pivot
vector, roughly a 3x reduction versus MATLAB's separate dense `[L,U,P]` triples (MATLAB:
`stepsPerCycle * 3 * n^2 * 8` bytes) — so `:auto` will favor `:lu` at larger `nr`/
`stepsPerCycle` in this port than it would in MATLAB for the same machine. That is
intentional, not a discrepancy to "fix" back to MATLAB's more conservative threshold.
"""
function estimated_lu_cache_bytes(stepsPerCycle::Integer, n::Integer)
    return stepsPerCycle * (n^2 * sizeof(Float64) + n * sizeof(Int))
end

# Applied to whatever available_memory_bytes() auto-detects (never to an explicit
# KMDISK_RAM_BUDGET_BYTES override, which the caller has presumably already sized
# deliberately) — the LU cache is the dominant but not the only consumer of a step's
# resident memory (the DTN matrix, system-matrix template, GC overhead, and Julia's own
# runtime all add up too), and a real cluster run got OOM-killed sizing :auto's decision
# right up against the detected total with no headroom at all.
const MEMORY_SAFETY_FACTOR = 0.7

"""
    available_memory_bytes() -> Int

The RAM budget `:auto` solver-type resolution should size itself against. In order:

1. `\$KMDISK_RAM_BUDGET_BYTES`, if set — an explicit override, used as-is (no safety
   margin applied; size it however you intend).
2. SLURM's `\$SLURM_MEM_PER_NODE` (megabytes; set automatically when a job requests
   `--mem`, e.g. from `julia/cluster/run_baseline_sweep.sh`), times `MEMORY_SAFETY_FACTOR`.
3. `Sys.free_memory()`, times `MEMORY_SAFETY_FACTOR` — correct on an unshared, non-SLURM
   machine, but **not cgroup-aware**: under SLURM (or any container/cgroup), it reports
   the whole node's free memory via libuv, not this job's actual `--mem` allocation. On a
   shared cluster node this can wildly overstate what a job can actually use — exactly
   what caused a real OOM kill: `:auto` sized a 16-thread sweep's LU caches against the
   apparent free memory of the *whole node*, not the job's real 32GB cgroup limit, and
   every thread happily tried to fit its own multi-GB cache within that illusory budget.
"""
function available_memory_bytes()
    override = get(ENV, "KMDISK_RAM_BUDGET_BYTES", "")
    isempty(override) || return parse(Int, override)

    slurm_mem_mb = get(ENV, "SLURM_MEM_PER_NODE", "")
    detected = isempty(slurm_mem_mb) ? Sys.free_memory() : parse(Int, slurm_mem_mb) * 1024 * 1024
    return floor(Int, detected * MEMORY_SAFETY_FACTOR)
end

"""
    resolve_solver_type(requested, stepsPerCycle, n, ram_budget_bytes) -> Symbol

`requested` is `:auto`, `:lu`, or `:gmres`. An explicit `:lu`/`:gmres` passes through
unchanged (an explicit override by the caller); `:auto` picks `:lu` if the estimated cache
fits in `ram_budget_bytes`, else `:gmres`.
"""
function resolve_solver_type(requested::Symbol, stepsPerCycle::Integer, n::Integer, ram_budget_bytes::Integer)
    requested === :auto || return requested
    return estimated_lu_cache_bytes(stepsPerCycle, n) <= ram_budget_bytes ? :lu : :gmres
end

# --- step solvers ---

abstract type AbstractStepSolver end

struct LUStepSolver <: AbstractStepSolver
    template::SystemMatrixTemplate
    cache::LUCache
    scratch::Matrix{Float64}
end

struct GMRESStepSolver <: AbstractStepSolver
    template::SystemMatrixTemplate
    gmres::GMRESSolver
    scratch::Matrix{Float64}
end

"""
    build_step_solver(problem, requested_type, gmres_tol; ram_budget_bytes=nothing) -> (solver, resolved_type)

Builds the system-matrix template once and constructs whichever concrete step solver
`resolve_solver_type` chooses. `ram_budget_bytes` defaults to `available_memory_bytes()`;
a sweep orchestrator running many cases concurrently should instead pass a per-case share
of that budget (see `sweep.jl`), so that concurrently-running cases don't each
independently assume they have the whole budget to themselves.
"""
function build_step_solver(problem::Problem, requested_type::Symbol, gmres_tol::Float64;
                            ram_budget_bytes::Union{Nothing,Integer} = nothing)
    n = 2 * problem.nr + 2
    budget = ram_budget_bytes === nothing ? available_memory_bytes() : ram_budget_bytes
    resolved = resolve_solver_type(requested_type, problem.stepsPerCycle, n, budget)
    tmpl = build_template(problem)
    scratch = Matrix{Float64}(undef, n, n)

    if resolved === :lu
        return LUStepSolver(tmpl, LUCache(problem.stepsPerCycle), scratch), resolved
    elseif resolved === :gmres
        return GMRESStepSolver(tmpl, build_gmres_solver(n, gmres_tol), scratch), resolved
    else
        error("resolve_solver_type returned unexpected symbol :$resolved (requested :$requested_type)")
    end
end

"""
    solve_step!(solver, phase, g_pf, rhs, sol) -> sol

Solves the linear system for this step's `g_prefactor` (`g_pf`) and RHS, writing the
solution into `sol` (and returning it). Dispatches on the concrete solver type.
"""
function solve_step!(solver::LUStepSolver, phase::StepPhase, g_pf::Float64,
                      rhs::AbstractVector{Float64}, sol::AbstractVector{Float64})
    F = get_or_factorize!(solver.cache, solver.template, solver.scratch, g_pf, phase.cycle_index)
    ldiv!(sol, F, rhs)
    return sol
end

function solve_step!(solver::GMRESStepSolver, phase::StepPhase, g_pf::Float64,
                      rhs::AbstractVector{Float64}, sol::AbstractVector{Float64})
    materialize!(solver.scratch, solver.template, g_pf)
    x = solve_gmres!(solver.gmres, solver.scratch, rhs)
    copyto!(sol, x)
    return sol
end
