"""
GMRES fallback solver path, used when `solverType=:gmres` (explicit) or when
`solverType=:auto` decides the LU cache would not fit in available RAM.

Port of the MATLAB call `gmres(Mat, indep, [], gmresTolerance, size(Mat,1), [], [],
gmres_x0)`: empty restart parameter and `maxit = size(Mat,1)` together mean *full,
non-restarted* GMRES — this must be preserved (Krylov.jl's `memory` keyword set to the
full system size reproduces this), not silently swapped for a restarted variant. The
matrix is rebuilt fresh every step (no caching, unlike the `:lu` path) via
`materialize!`, and the solve is warm-started from the previous step's solution.

CAVEAT: this file was written without the ability to run Julia or install Krylov.jl in
the environment it was authored in (network access to the Julia package registry was
blocked), so the exact in-place `gmres!`/workspace keyword names below are a best-effort
match to the Krylov.jl API and have NOT been executed. If `Pkg.instantiate()` or the test
suite reports a `MethodError`/`UndefKeywordError` originating here, this is the first
place to check against whatever Krylov.jl version actually resolves — the fix is local to
this file (`GMRESSolver`/`solve_gmres!`); nothing else in the package depends on Krylov's
exact call signature.
"""
mutable struct GMRESSolver
    workspace::Krylov.GmresWorkspace
    tol::Float64
    x0::Vector{Float64}
    last_converged::Bool
end

"""
    build_gmres_solver(n, tol) -> GMRESSolver

`n` is the system size (`2*nr+2`). `memory = n` gives full (non-restarted) GMRES,
matching MATLAB's `gmres(...,[],tol,n,...)`.
"""
function build_gmres_solver(n::Int, tol::Float64)
    workspace = Krylov.GmresWorkspace(n, n, Float64; memory = n)
    return GMRESSolver(workspace, tol, zeros(n), true)
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
                  atol = solver.tol, rtol = solver.tol, itmax = length(rhs))
    stats = solver.workspace.stats
    solver.last_converged = stats.solved
    if !stats.solved
        @warn "GMRES did not converge" tol = solver.tol itmax = length(rhs) niter = stats.niter
    end
    copyto!(solver.x0, solver.workspace.x)
    return solver.workspace.x
end
