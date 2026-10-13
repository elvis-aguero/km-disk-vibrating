"""
GMRES solver path, used when `solverType=:gmres` (explicit) or when `solverType=:auto`
decides the LU cache would not fit in available RAM.

Preconditioned with ONE LU factorization of the system matrix at `g_prefactor = 1` (the
centre of the `1 ± gamma` range the matrix sweeps over a cycle), applied as a right
preconditioner on every iteration. Each step still solves that step's own matrix
`M(g)`; the preconditioner only changes how fast GMRES gets there, not the answer.
Because `M(g) = M(1) + (g-1)*(dt/Fr)*P`, the preconditioned operator is a small
perturbation of the identity: measured 1-5 iterations per step at nr = 100-300 for
gamma up to 0.5, versus ~n iterations (one per unknown) unpreconditioned. Memory is one
`n x n` factor instead of the LU cache's one per cycle step.

GMRES is not restarted (`memory` is only an initial basis size; Krylov.jl grows it if
needed) and `itmax = n`, so convergence is still guaranteed in exact arithmetic even if
the preconditioner were poor — matching MATLAB's `gmres(Mat, indep, [], tol, n, ...)`.
"""
mutable struct GMRESSolver
    workspace::Krylov.GmresWorkspace
    precond::LU{Float64,Matrix{Float64},Vector{Int}}
    tol::Float64
    x0::Vector{Float64}
    last_converged::Bool
    last_niter::Int
    total_niter::Int
    nsolves::Int
end

const GMRES_INITIAL_MEMORY = 20

"""
    build_gmres_solver(precond, tol) -> GMRESSolver

`precond` is the LU factorization used as the right preconditioner (see the file
docstring); the system size `n` is taken from it.
"""
function build_gmres_solver(precond::LU{Float64,Matrix{Float64},Vector{Int}}, tol::Float64)
    n = size(precond, 1)
    # S must be the storage (array) type, not the element type — GmresWorkspace
    # allocates internal buffers via `S(undef, n)`, which only an array type supports.
    workspace = Krylov.GmresWorkspace(n, n, Vector{Float64}; memory = min(n, GMRES_INITIAL_MEMORY))
    return GMRESSolver(workspace, precond, tol, zeros(n), true, 0, 0, 0)
end

"""
    solve_gmres!(solver, Mat, rhs) -> Vector{Float64}

Solves `Mat * x = rhs`, warm-started from the previous call's solution (zero vector on
the first call). Returns the solution vector (a view into the solver's internal
workspace — copy it if you need it to outlive the next call). Logs a warning (does not
throw) if GMRES fails to converge within `tol`/`itmax`, matching MATLAB's non-fatal
`warning('advance_one_step:gmresNoConverge', ...)`.
"""
function solve_gmres!(solver::GMRESSolver, Mat::AbstractMatrix{Float64}, rhs::AbstractVector{Float64})
    Krylov.gmres!(solver.workspace, Mat, rhs, solver.x0;
                  N = solver.precond, ldiv = true,
                  atol = solver.tol, rtol = solver.tol, itmax = length(rhs))
    stats = solver.workspace.stats
    solver.last_converged = stats.solved
    solver.last_niter = stats.niter
    solver.total_niter += stats.niter
    solver.nsolves += 1
    if !stats.solved
        @warn "GMRES did not converge" tol = solver.tol itmax = length(rhs) niter = stats.niter
    end
    copyto!(solver.x0, solver.workspace.x)
    return solver.workspace.x
end

"""
    mean_iterations(solver) -> Float64

Average GMRES iterations per solve so far (`NaN` before the first solve).
"""
mean_iterations(solver::GMRESSolver) = solver.nsolves == 0 ? NaN : solver.total_niter / solver.nsolves
