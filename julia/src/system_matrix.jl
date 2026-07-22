"""
Port of `matlab/1_code/simulation_code/buildSystemMatrix.m`: assembles the dense
`(2nr+2) x (2nr+2)` implicit time-step operator.

Unknown vector layout (matches `advance_one_step.m`'s unpacking of `sol` exactly):
  [1 : nr-cPoints]              -> bath surface height outside the disk (eta_out)
  [nr-cPoints+1 : 2nr-cPoints]  -> surface velocity potential (phi), all nr nodes
  [2nr-cPoints+1 : 2nr]         -> pressure under the disk (cPoints nodes)
  [2nr+1]                        -> center-of-mass velocity
  [2nr+2]                        -> center-of-mass position (also the eta value under
                                     the disk, since the rigid disk forces eta to be flat
                                     and equal to its own depth there)

Only one physical quantity depends on the time-varying `g_prefactor` (the oscillating
effective-gravity term): the diagonal of `eye(nr)/Fr * g_prefactor` in the eta-evolution
block. Tracing that diagonal through the column-reduction below shows it touches exactly
`nr` entries of the assembled matrix (`cPoints` entries in the last column, `nr-cPoints`
entries in the eta_out block) — everything else is invariant across steps. `materialize!`
exploits this: it patches those `nr` entries onto a matrix built once with
`g_prefactor=0`, instead of rebuilding the whole `O(nr^2)` dense matrix from scratch every
step (which is what MATLAB's GMRES path does on every single step).
"""

"""
    assemble_system_matrix!(Mat, problem, g_prefactor)

Fills `Mat` (must be `(2nr+2) x (2nr+2)`) with the system matrix for the given
`g_prefactor`, from scratch. This is the reference implementation — `materialize!` is a
faster equivalent, verified against this function in the test suite.
"""
function assemble_system_matrix!(Mat::AbstractMatrix{Float64}, problem::Problem, g_prefactor::Float64)
    nr = problem.nr
    cPoints = problem.cPoints
    dt = problem.dt
    Re = problem.Re
    Fr = problem.Fr
    We = problem.We
    Delta = problem.laplacian
    DTN = problem.DTN
    pIntegral = problem.pressureIntegral
    SF = problem.surface_force_constant
    dr = problem.dr
    Ma = problem.obj_mass

    n = 2 * nr + 2
    size(Mat) == (n, n) || throw(DimensionMismatch("Mat must be $(n)x$(n), got $(size(Mat))"))
    fill!(Mat, 0.0)

    # --- eta_out columns: Sist columns (cPoints+1):nr -> Mat columns 1:(nr-cPoints) ---
    for j in 1:(nr - cPoints)
        sist_col = cPoints + j
        for i in 1:nr
            # top block: eye(nr) - dt*2*Delta/Re
            val = -dt * 2 * Delta[i, sist_col] / Re
            i == sist_col && (val += 1.0)
            Mat[i, j] = val
        end
        for i in 1:nr
            # bottom-left block: dt*(eye(nr)/Fr*g_prefactor - Delta/We)
            val = -dt * Delta[i, sist_col] / We
            i == sist_col && (val += dt / Fr * g_prefactor)
            Mat[nr + i, j] = val
        end
    end

    # --- phi columns: Sist columns (nr+1):2nr -> Mat columns (nr-cPoints+1):(2nr-cPoints) ---
    for j in 1:nr
        mat_col = (nr - cPoints) + j
        for i in 1:nr
            Mat[i, mat_col] = -dt * DTN[i, j]                         # top-right: -dt*DTN
        end
        for i in 1:nr
            val = -dt * 2 * Delta[i, j] / Re                          # bottom-right: eye(nr) - dt*2*Delta/Re
            i == j && (val += 1.0)
            Mat[nr + i, mat_col] = val
        end
    end

    # --- pressure-under-disk columns: Mat columns (2nr-cPoints+1):2nr = dt*I(cPoints), rows nr+1:nr+cPoints ---
    for p in 1:cPoints
        Mat[nr + p, (2 * nr - cPoints) + p] = dt
    end

    # --- last column (CoM position z) = row-sum of Sist(:, 1:cPoints) ---
    for i in 1:nr
        s_top = 0.0
        for q in 1:cPoints
            val = -dt * 2 * Delta[i, q] / Re
            i == q && (val += 1.0)
            s_top += val
        end
        Mat[i, 2 * nr + 2] = s_top

        s_bot = 0.0
        for q in 1:cPoints
            val = -dt * Delta[i, q] / We
            i == q && (val += dt / Fr * g_prefactor)
            s_bot += val
        end
        Mat[nr + i, 2 * nr + 2] = s_bot
    end

    # --- disk force balance row (2nr+1) ---
    Mat[2 * nr + 1, 1] = -SF * dt / dr
    for p in 1:cPoints
        Mat[2 * nr + 1, (2 * nr - cPoints) + p] = -dt * pIntegral[p] / Ma
    end
    Mat[2 * nr + 1, 2 * nr + 1] = 1.0
    Mat[2 * nr + 1, 2 * nr + 2] = SF * dt / dr

    # --- kinematic update row (2nr+2): z_{n+1} = z_n + dt*v_{n+1} ---
    Mat[2 * nr + 2, 2 * nr + 1] = -dt
    Mat[2 * nr + 2, 2 * nr + 2] = 1.0

    return Mat
end

function assemble_system_matrix(problem::Problem, g_prefactor::Float64)
    n = 2 * problem.nr + 2
    Mat = Matrix{Float64}(undef, n, n)
    return assemble_system_matrix!(Mat, problem, g_prefactor)
end

struct SystemMatrixTemplate
    Mat0::Matrix{Float64}
    patch_rows::Vector{Int}
    patch_cols::Vector{Int}
    patch_coeff::Float64
end

"""
    build_template(problem) -> SystemMatrixTemplate

Builds the `g_prefactor=0` base matrix once, plus the list of `(row,col)` entries that
scale linearly with `g_prefactor` (all with the same coefficient `dt/Fr`, so a single
scalar suffices rather than a per-entry coefficient vector).
"""
function build_template(problem::Problem)
    nr = problem.nr
    cPoints = problem.cPoints
    Mat0 = assemble_system_matrix(problem, 0.0)

    patch_rows = Int[]
    patch_cols = Int[]
    sizehint!(patch_rows, nr)
    sizehint!(patch_cols, nr)
    for p in 1:cPoints
        push!(patch_rows, nr + p)
        push!(patch_cols, 2 * nr + 2)
    end
    for j in 1:(nr - cPoints)
        push!(patch_rows, nr + cPoints + j)
        push!(patch_cols, j)
    end

    return SystemMatrixTemplate(Mat0, patch_rows, patch_cols, problem.dt / problem.Fr)
end

"""
    materialize!(Mat, tmpl, g_prefactor)

Fills `Mat` with `tmpl.Mat0` plus the `g_prefactor`-scaled patch — an `O(nr)` update
instead of the `O(nr^2)` full rebuild `assemble_system_matrix!` performs. Verified
equivalent to `assemble_system_matrix!` in the test suite.
"""
function materialize!(Mat::Matrix{Float64}, tmpl::SystemMatrixTemplate, g_prefactor::Float64)
    copyto!(Mat, tmpl.Mat0)
    @inbounds for k in eachindex(tmpl.patch_rows)
        r = tmpl.patch_rows[k]
        c = tmpl.patch_cols[k]
        Mat[r, c] += tmpl.patch_coeff * g_prefactor
    end
    return Mat
end

"""
    assemble_rhs!(rhs, bath_surface, bath_potential, center_of_mass_velocity, center_of_mass,
                  dt, Fr, g_prefactor, force_term)

Port of `indep = [b; CoM_vel - dt/Fr*g_prefactor - force_term; CoM]` from
`advance_one_step.m`, where `b = [bath_surface; bath_potential]`.
"""
function assemble_rhs!(rhs::Vector{Float64}, bath_surface::AbstractVector{Float64},
                        bath_potential::AbstractVector{Float64}, center_of_mass_velocity::Float64,
                        center_of_mass::Float64, dt::Float64, Fr::Float64, g_prefactor::Float64,
                        force_term::Float64)
    nr = length(bath_surface)
    length(bath_potential) == nr || throw(DimensionMismatch("bath_surface and bath_potential must have equal length"))
    length(rhs) == 2 * nr + 2 || throw(DimensionMismatch("rhs must have length 2*nr+2"))

    @views rhs[1:nr] .= bath_surface
    @views rhs[(nr + 1):(2 * nr)] .= bath_potential
    rhs[2 * nr + 1] = center_of_mass_velocity - dt / Fr * g_prefactor - force_term
    rhs[2 * nr + 2] = center_of_mass
    return rhs
end
