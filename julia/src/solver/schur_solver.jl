"""
Schur-complement step solver (`solverType = :schur`): the same linear system as the LU path, reduced from
`(2nr+2) x (2nr+2)` to `(nr+1) x (nr+1)` unknowns before it is factorized.

Why it is exact. In `system_matrix.jl` the bath potential `phi` is determined row by row by the lower block
`B phi = rhs - K eta - dt [p; 0]`, with `B = I - 2 dt Δ/Re` and `K = -(dt/We) Δ + (dt g/Fr) I`, both TRIDIAGONAL
(`Δ` is the radial Laplacian). Only the upper block carries the dense DtN matrix `N`: `A eta - dt N phi = r_top`
(`A = B`). Solving the lower block for `phi` with the tridiagonal `B` and substituting gives, for the unknowns
`y = [eta_out; z; p]` (surface height outside the disk, disk height, pressure under the disk; `eta` under the disk
equals `z`),

    S(g) y = [r_top + dt Y r_bot ; r_force + r_kin/dt],     Y = N B^-1,
    S(g)   = S0 + g G1,        (nr+1) x (nr+1),

where the velocity `v = (z - r_kin)/dt` is eliminated with the kinematic row, and `g` (the oscillating effective
gravity) enters only through `K`, so `S` is affine in `g`. `Y` costs `O(nr^2)` (one tridiagonal solve per row of
`N`), not `O(nr^3)`. The matrices cached per distinct `g` are `(nr+1)^2`: one quarter of the memory and about one
eighth of the factorization flops of the full system. `phi`, `v` follow by back substitution. Algebra checked
against the full solve at machine precision (`julia/test/unit/test_schur_solver.jl`).

Memory: `S0` and `G1` (2 matrices) plus one LU per distinct `g` (31 for the standard bath-driven cycle) plus a
transient `Y` while building; see `estimated_schur_bytes`.
"""

"""
    estimated_schur_bytes(n_factors, nr) -> Int

Peak memory of the Schur path: `n_factors` factorizations of `(nr+1)^2`, the two stored matrices `S0`, `G1`, and the
transient `Y` (+ its transpose) while they are built.
"""
function estimated_schur_bytes(n_factors::Integer, nr::Integer)
    m = nr + 1
    return n_factors * (m^2 * sizeof(Float64) + m * sizeof(Int)) + 2 * m^2 * sizeof(Float64) + 2 * nr^2 * sizeof(Float64)
end

"""
Reduced system for one implicit coefficient `dt`: `S(g) = S0 + g G1` and the tridiagonal pieces for the back
substitution. See the file header.
"""
struct SchurSystem
    nr::Int
    cPoints::Int
    dt::Float64
    Fr::Float64
    Bfact::Factorization{Float64}                 # lu(B), B = I - 2 dt Δ/Re (tridiagonal)
    K0::Tridiagonal{Float64,Vector{Float64}}      # -(dt/We) Δ  (the g-independent part of K)
    S0::Matrix{Float64}
    G1::Matrix{Float64}
    DTN::Matrix{Float64}                          # reference to problem.DTN (not copied)
    w_over_Ma::Vector{Float64}
end

"""
    build_schur_system(problem, dt) -> SchurSystem

`dt` is the implicit coefficient (`implicit_dt(problem, scheme)`). Cost `O(nr^2)`.
"""
function build_schur_system(problem::Problem, dt::Float64)
    nr, c = problem.nr, problem.cPoints
    no = nr - c
    Δ = problem.laplacian
    Δ isa Tridiagonal || throw(ArgumentError(
        "solverType = :schur needs a tridiagonal Laplacian (uniform grid, or graded grid with laplacian_order = 2); " *
        "got $(typeof(Δ)). Use solverType = :lu, or laplacian_order = 2."))
    Re, Fr, We, Ma = problem.Re, problem.Fr, problem.We, problem.obj_mass
    N = problem.DTN
    m = nr + 1

    B = Tridiagonal(Vector{Float64}(-(2dt / Re) .* Δ.dl), 1.0 .- (2dt / Re) .* Δ.d, Vector{Float64}(-(2dt / Re) .* Δ.du))
    Bfact = lu(B)

    # Y = N B^-1  <=>  Y' = B'^-1 N'  (tridiagonal solves, O(nr^2))
    Bt = Tridiagonal(copy(B.du), copy(B.d), copy(B.dl))
    Yt = Matrix{Float64}(undef, nr, nr)
    transpose!(Yt, N)
    ldiv!(lu(Bt), Yt)
    Y = Matrix{Float64}(undef, nr, nr)
    transpose!(Y, Yt)
    Yt = nothing

    # Hη(g) = A + dt Y K(g) = [A - (dt^2/We) Y Δ] + g (dt^2/Fr) Y,  columns j of the full η vector
    d, dlow, dup = Δ.d, Δ.dl, Δ.du     # Δ[j,j], Δ[j+1,j], Δ[j-1,j] = du[j-1]
    S0 = zeros(m, m)
    G1 = zeros(m, m)
    cz = (dt^2 / We)
    cg = (dt^2 / Fr)
    function eta_column!(col::Int, j::Int, scale::Float64)   # adds column j of Hη0 / G1 into column `col`
        # A[:, j] = e_j - (2dt/Re) Δ[:, j]  = B[:, j]; nonzeros rows j-1, j, j+1
        @inbounds begin
            j > 1 && (S0[j - 1, col] += scale * B.du[j - 1])
            S0[j, col] += scale * B.d[j]
            j < nr && (S0[j + 1, col] += scale * B.dl[j])
            # - (dt^2/We) (Y Δ)[:, j] = - (dt^2/We) (Y[:, j-1] Δ[j-1, j] + Y[:, j] Δ[j, j] + Y[:, j+1] Δ[j+1, j])
            for i in 1:nr
                acc = Y[i, j] * d[j]
                j > 1 && (acc += Y[i, j - 1] * dup[j - 1])
                j < nr && (acc += Y[i, j + 1] * dlow[j])
                S0[i, col] -= scale * cz * acc
                G1[i, col] += scale * cg * Y[i, j]
            end
        end
    end
    for j in (c + 1):nr            # eta_out columns
        eta_column!(j - c, j, 1.0)
    end
    for q in 1:c                   # the z column = sum of the first c eta columns
        eta_column!(no + 1, q, 1.0)
    end
    @inbounds for p in 1:c, i in 1:nr   # pressure columns: dt^2 Y[:, 1:c]
        S0[i, no + 1 + p] = dt^2 * Y[i, p]
    end
    S0[m, no + 1] = 1 / dt             # (z - r_kin)/dt - dt w.p/Ma = r_force
    @inbounds for p in 1:c
        S0[m, no + 1 + p] = -dt * problem.contact_weights[p] / Ma
    end

    K0 = Tridiagonal(Vector{Float64}(-(dt / We) .* Δ.dl), -(dt / We) .* Δ.d, Vector{Float64}(-(dt / We) .* Δ.du))
    return SchurSystem(nr, c, dt, Fr, Bfact, K0, S0, G1, N, problem.contact_weights ./ Ma)
end

"Solves one step with the reduced system and its factorization `F` of `S(g)`; writes the full-layout solution."
function schur_solve!(sys::SchurSystem, F::Factorization{Float64}, g::Float64, rhs::AbstractVector{Float64},
                      sol::AbstractVector{Float64})
    nr, c, dt = sys.nr, sys.cPoints, sys.dt
    no = nr - c
    u = Vector{Float64}(undef, nr)
    @views u .= rhs[(nr + 1):(2nr)]
    ldiv!(sys.Bfact, u)                                   # u = B^-1 r_bot
    b = Vector{Float64}(undef, nr + 1)
    mul!(view(b, 1:nr), sys.DTN, u, dt, 0.0)              # dt N u
    @views b[1:nr] .+= rhs[1:nr]
    r_f, r_k = rhs[2nr + 1], rhs[2nr + 2]
    b[nr + 1] = r_f + r_k / dt
    ldiv!(F, b)                                           # y = [eta_out; z; p]
    z = b[no + 1]
    q = Vector{Float64}(undef, nr)                        # eta = [z 1_c; eta_out]
    @views q[1:c] .= z
    @views q[(c + 1):nr] .= b[1:no]
    t = sys.K0 * q
    @. t += (dt * g / sys.Fr) * q
    @views t[1:c] .+= dt .* b[(no + 2):(nr + 1)]
    ldiv!(sys.Bfact, t)                                   # B^-1 (K eta + dt [p;0])
    @views sol[1:no] .= b[1:no]
    @views sol[(no + 1):(no + nr)] .= u .- t              # phi
    @views sol[(2nr - c + 1):(2nr)] .= b[(no + 2):(nr + 1)]
    sol[2nr + 1] = (z - r_k) / dt
    sol[2nr + 2] = z
    return sol
end

"Per-run solver for `solverType = :schur`: reduced system plus a cache of the (nr+1)^2 factorizations."
struct SchurStepSolver <: AbstractStepSolver
    sys::SchurSystem
    cache::LUCache
    scheme::Symbol
end

function get_or_factorize_schur!(cache::LUCache, sys::SchurSystem, g::Float64, cycle_index::Int)
    1 <= cycle_index <= length(cache) || throw(BoundsError(cache.factors, cycle_index))
    f = cache.factors[cycle_index]
    if f === nothing
        f = _find_factorization(cache, g)
        if f === nothing
            S = sys.S0 .+ g .* sys.G1          # fresh matrix, factorized in place
            f = lu!(S)
            push!(cache.by_value, g => f)
        end
        cache.factors[cycle_index] = f
        cache.misses += 1
    else
        cache.hits += 1
    end
    return f
end

function solve_step!(solver::SchurStepSolver, phase::StepPhase, g_pf::Float64,
                      rhs::AbstractVector{Float64}, sol::AbstractVector{Float64})
    F = get_or_factorize_schur!(solver.cache, solver.sys, g_pf, phase.cycle_index)
    return schur_solve!(solver.sys, F, g_pf, rhs, sol)
end

"""
    solve_startup!(solver, problem, g_pf, rhs, sol)

The first BDF2 step is one implicit Euler step with its own matrix (coefficient `dt`, not `2dt/3`). The dense
solvers factorize the full system once; the Schur solver builds the reduced Euler system, solves, and discards it.
"""
function solve_startup!(::AbstractStepSolver, problem::Problem, g_pf::Float64, rhs::AbstractVector{Float64},
                        sol::AbstractVector{Float64})
    ldiv!(sol, lu!(assemble_system_matrix(problem, g_pf)), rhs)
    return sol
end

function solve_startup!(::SchurStepSolver, problem::Problem, g_pf::Float64, rhs::AbstractVector{Float64},
                        sol::AbstractVector{Float64})
    sys = build_schur_system(problem, problem.dt)
    F = lu!(sys.S0 .+ g_pf .* sys.G1)
    return schur_solve!(sys, F, g_pf, rhs, sol)
end

"""
    schur_applicable(params) -> Bool

Whether a run with these parameters uses the Schur solver: `solverType = :schur`, or `:auto` on a grid whose Laplacian is
tridiagonal (the uniform grid, or a graded grid with `laplacian_order = 2`). Used by the sweep to size its memory budget.
"""
schur_applicable(params::SimulationParams) =
    params.solverType === :schur ||
    (params.solverType === :auto && (params.grid === nothing || params.grid.laplacian_order == 2))
