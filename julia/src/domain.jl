# Pressure-integral matrix (legacy; no longer used for the disk force)
# Exact integral of a piecewise-linear pressure that decays to zero half a cell beyond
# the edge node (from the original kinematic-match solver, where the contact area can
# change). The disk force now uses `contact_weights` instead: the edge pressure here
# carries the meniscus pull, which these weights would count only ~3/4 of.
function _pressure_integral_matrix(quant::Float64, dr::Float64)
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
    return IntMat
end

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

    IntMat = _pressure_integral_matrix(quant, dr)

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

# Conservative second-order Laplacian on a mapped grid (TRIDIAGONAL):
#   Delta g_i = (1/(r_i s_i)) [ c_{i+1/2} (g_{i+1} - g_i) - c_{i-1/2} (g_i - g_{i-1}) ],  c = r/s at the half nodes,
# s = dr/dxi. Rows inside the identity region (i < xin) are the rows of `build_domain(D, sr)`, bit for bit.
function _laplacian_conservative(g::RadialGrid)
    nr, dr0, xin = g.nr, g.dr0, g.xin
    d = fill(-2.0 / dr0^2, nr)
    du = zeros(nr - 1)
    dl = zeros(nr - 1)
    d[1] = -4.0 / dr0^2
    du[1] = 4.0 / dr0^2
    c = [_grid_r_half(g, j) / _grid_spacing_at(g, j + 0.5)[1] for j in 0:(nr - 1)]     # c[j+1] = (r/s) at xi = j + 1/2
    for i in 1:(nr - 1)                                                                  # node index (0-based)
        if i < xin
            i <= nr - 2 && (du[i + 1] = (1.0 / dr0^2) * (1.0 + 1.0 / (2.0 * i)))
            dl[i] = (1.0 / dr0^2) * (1.0 - 1.0 / (2.0 * i))
        else
            f = 1.0 / (g.r[i + 1] * g.s[i + 1])
            d[i + 1] = -f * (c[i + 1] + c[i])
            dl[i] = f * c[i]
            i <= nr - 2 && (du[i + 1] = f * c[i + 1])
        end
    end
    return Tridiagonal(dl, d, du)
end

# Finite-difference Laplacian of order 4 or 6 on a mapped grid, banded but returned dense:
#   Delta g = g_xixi / s^2 + (1/(r s) - s'/s^3) g_xi ,  stencils for g_xixi, g_xi at spacing 1.
# Axis row: 2 g_xixi / s^2 with even extension g_{-k} = g_k. Wall: g_j = 0 for j >= nr.
function _laplacian_banded(g::RadialGrid, order::Int)
    haskey(_FD_STENCIL, order) || throw(ArgumentError("order must be 4 or 6, got $order"))
    c2, c1 = _FD_STENCIL[order]
    h = length(c2) ÷ 2
    nr = g.nr
    L = zeros(nr, nr)
    for k in -h:h
        col = abs(k)
        col < nr && (L[1, col + 1] += 2 * c2[k + h + 1] / g.s[1]^2)
    end
    for i in 1:(nr - 1)
        s = g.s[i + 1]
        A2 = 1 / s^2
        A1 = 1 / (g.r[i + 1] * s) - g.sp[i + 1] / s^3
        for k in -h:h
            col = abs(i + k)
            col < nr && (L[i + 1, col + 1] += A2 * c2[k + h + 1] + A1 * c1[k + h + 1])
        end
    end
    return L
end

"""
    build_domain(grid::RadialGrid) -> (dr, Laplacian, IntMat, nr, r)

Domain operators for a radial grid. `dr` is the spacing `dr0` of the disk region, `r` the node radii (disk radii).
A uniform grid with `laplacian_order = 2` returns exactly `build_domain(D, sr)` (same `Tridiagonal`, bit for bit).
Otherwise `Laplacian` is the conservative second-order `Tridiagonal` (`laplacian_order = 2`) or a dense banded matrix
(order 4 or 6). Formulas: `docs/graded_grid.md`.
"""
function build_domain(grid::RadialGrid)
    if grid.uniform && grid.laplacian_order == 2
        dom = build_domain(2 * grid.R_req, grid.sr)
        dom.nr == grid.nr || error("internal error: uniform grid nr=$(grid.nr) differs from build_domain nr=$(dom.nr)")
        return (dr = dom.dr, Laplacian = dom.Laplacian, IntMat = dom.IntMat, nr = dom.nr, r = grid.r[1:dom.nr])
    end
    Lap = grid.laplacian_order == 2 ? _laplacian_conservative(grid) : _laplacian_banded(grid, grid.laplacian_order)
    IntMat = _pressure_integral_matrix(grid.sr, grid.dr0)
    return (dr = grid.dr0, Laplacian = Lap, IntMat = IntMat, nr = grid.nr, r = grid.r[1:grid.nr])
end
