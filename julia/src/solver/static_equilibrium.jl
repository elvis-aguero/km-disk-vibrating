"""
Port of the static-equilibrium initial condition (`solve_motion.m` lines ~159-202,
`startStatic` option): solves a coupled linear system for the equilibrium bath profile
outside the disk, the pressure distribution under the disk, and the disk's own
equilibrium depth `z`, so a simulation can start at rest at the physically correct depth
instead of `z=0` (which produced a large spurious transient in earlier versions of the
MATLAB code).

Unknown vector layout (length `nr+1`): `[eta_out (nr-cPoints); pressure_under_disk
(cPoints); z (1)]` — the same "eta under the disk collapses to a single unknown z"
reduction used by the main system matrix in `system_matrix.jl`.

This is a small, one-time dense solve (`(nr+1) x (nr+1)`), not a per-step hot path, so a
direct, literal translation is appropriate here (unlike `system_matrix.jl`'s
performance-sensitive per-step assembly).
"""
function solve_static_equilibrium(problem::Problem)
    nr = problem.nr
    cPoints = problem.cPoints
    Fr = problem.Fr
    We = problem.We
    Delta = problem.laplacian
    dr = problem.dr
    SF = problem.surface_force_constant
    Ma = problem.obj_mass
    pIntegral = problem.pressureIntegral

    # A_static(i,j) = (eye(nr)/Fr - laplacian/We)(i,j), evaluated lazily (Delta is
    # tridiagonal, so this is O(1) per entry) rather than materializing a dense nr x nr
    # matrix just to slice a handful of rows out of it.
    Aij(i, j) = (i == j ? 1.0 / Fr : 0.0) - Delta[i, j] / We

    n = nr + 1
    M = zeros(n, n)
    rhs = zeros(n)

    row = 1
    # 1. Equations for r > R_disk (pure bath equilibrium, no pressure coupling)
    for i in (cPoints + 1):nr
        for k in 1:(nr - cPoints)
            M[row, k] = Aij(i, cPoints + k)
        end
        s = 0.0
        for j in 1:cPoints
            s += Aij(i, j)
        end
        M[row, n] = s
        row += 1
    end

    # 2. Equations for r <= R_disk (solve for the pressure that keeps eta flat at z)
    for i in 1:cPoints
        for k in 1:(nr - cPoints)
            M[row, k] = Aij(i, cPoints + k)
        end
        M[row, nr - cPoints + i] = 1.0
        s = 0.0
        for j in 1:cPoints
            s += Aij(i, j)
        end
        M[row, n] = s
        row += 1
    end

    # 3. Disk force balance: capillary edge force + integrated pressure force = weight
    M[row, 1] = SF / dr
    @views M[row, (nr - cPoints + 1):nr] .= pIntegral ./ Ma
    M[row, n] = -SF / dr
    rhs[row] = 1.0 / Fr
    row == n || error("internal error: static equilibrium row count mismatch ($row != $n)")

    sol = M \ rhs

    eta_out = sol[1:(nr - cPoints)]
    pressure_initial = sol[(nr - cPoints + 1):nr]
    z_initial = sol[n]
    eta_initial = vcat(fill(z_initial, cPoints), eta_out)

    return (eta = eta_initial, pressure = pressure_initial, z = z_initial)
end
