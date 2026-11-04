"""
    build_domain(bathDiameter, spatialResolution) -> (dr, Laplacian, IntMat)

Port of `matlab/1_code/simulation_code/domainMaker.m`.

Builds the radial finite-difference discretization used by the rest of the solver:
- `dr`: mesh spacing (in units of disk radii).
- `Laplacian`: the discrete radial Laplacian `(1/r) d/dr(r d/dr)`, as a `Tridiagonal`.
  MATLAB's `domainMaker.m` builds this densely from several intermediate `nr x nr`
  matrices (`subdiag`, `Dr`, `DrOverR`, `Derivrr`), but the whole construction reduces
  algebraically to a tridiagonal operator (verified by hand against the MATLAB source),
  so it is built directly here as one, without ever allocating an `nr x nr` dense matrix.
- `IntMat`: the pressure-integral matrix (only ever indexed at a single row elsewhere
  in the codebase — the full matrix is returned so callers can slice whichever row they
  need, matching `pressureIntegral(cPoints, 1:cPoints)` in the MATLAB source).

Row 1 (r = 0) uses a mirror/Neumann boundary stencil (`-4/dr^2, 4/dr^2`), the standard
finite-difference treatment of the removable singularity in `(1/r) d/dr(r d/dr)` at the
axis of a cylindrically symmetric domain — this is intentional physics, not an artifact,
and is preserved exactly from the MATLAB source.
"""
function build_domain(bathDiameter::Real, spatialResolution::Real)
    D = float(bathDiameter)
    quant = float(spatialResolution)

    nr = ceil(Int, D * quant / 2)
    dr = D / (2 * nr)

    # --- Laplacian (tridiagonal) ---
    d  = fill(-2.0 / dr^2, nr)
    du = zeros(nr - 1)   # du[i] = entry (i, i+1)
    dl = zeros(nr - 1)   # dl[i] = entry (i+1, i)

    # Row 1: mirror-boundary stencil at r = 0.
    d[1]  = -4.0 / dr^2
    du[1] = 4.0 / dr^2

    # Interior rows i = 2:(nr-1) contribute both an upper and lower neighbor entry;
    # row nr only ever contributes a lower neighbor entry (no column nr+1 exists).
    for i in 2:(nr - 1)
        du[i] = (1.0 / dr^2) * (1.0 + 1.0 / (2.0 * (i - 1)))
    end
    for i in 2:nr
        dl[i - 1] = (1.0 / dr^2) * (1.0 - 1.0 / (2.0 * (i - 1)))
    end

    Laplacian = Tridiagonal(dl, d, du)

    # --- Pressure-integral matrix (legacy; no longer used for the disk force) ---
    # Exact integral of a piecewise-linear pressure that decays to zero half a cell beyond
    # the edge node (from the original kinematic-match solver, where the contact area can
    # change). The disk force now uses `contact_weights` instead: the edge pressure here
    # carries the meniscus pull, which these weights would count only ~3/4 of.
    nlmax = Int(quant) + 1
    n_int = 2 * nlmax
    IntMat = zeros(n_int, n_int)
    IntMat[1, 1] = 1.0 / 12.0
    for ii in 2:n_int
        IntMat[ii, 1] = 1.0 / 3.0
        for jj in 2:(ii - 1)
            IntMat[ii, jj] = 2.0 * jj - 2.0
        end
        IntMat[ii, ii] = 1.5 * ii - 21.0 / 12.0
    end
    IntMat .*= pi * dr^2

    return (dr = dr, Laplacian = Laplacian, IntMat = IntMat, nr = nr)
end

"""
    contact_weights(cPoints, dr) -> Vector{Float64}

Area weights for the pressure under the disk, as the fluid applies it: the ring
`2π r_j dr` at contact node `j` (`π (dr/2)²` at the centre), in units of the disk radius
squared. These are the weights under which the discrete Laplacian is symmetric, so the
disk receives exactly the load the water feels. Summed, the pressure rows then give the
buoyancy plus the meniscus pull on the rim exactly — the edge node's pressure carries that
pull, because its Laplacian row reaches the first free node — so the disk force uses these
weights and no separate contact-line term (see `system_matrix.jl`).
"""
function contact_weights(cPoints::Integer, dr::Real)
    w = [2pi * (j - 1) * dr * dr for j in 1:cPoints]
    w[1] = pi * (dr / 2)^2
    return w
end
