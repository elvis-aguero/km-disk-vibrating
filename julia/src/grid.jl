#=
Radial grids: the uniform grid used so far, and a graded grid (fine near the disk, coarse far away).

Design, formulas, pass/fail criteria: `docs/graded_grid.md`. In short, node `i` has the index coordinate `xi = i` (spacing 1) and
the radius `r(xi)` with `s = r_xi`:

  * `r <= r_in`: identity map, `s = dr0` (the disk, the contact nodes and `contact_weights` are untouched).
  * outside: `s = dr0 exp(l)`, the log-spacing `l` grows with slope `ln q` (growth `q` per node, after a smooth ramp of `ramp`
    nodes) and is capped smoothly at `ln(dr_max/dr0)`. So `s <= dr_max` and `s_{i+1}/s_i <= q`.
  * `r(xi)` is the integral of `s` (10-point Gauss-Legendre on each unit interval).

Lengths are in disk radii. Nodes `0 .. nr-1` are unknowns; node `nr` is the wall (`phi = 0`). Three more nodes `nr+1 .. nr+3`
exist because the DtN stencil rows reach them.
=#

const _GRID_GL = let n = 10
    E = eigen(SymTridiagonal(zeros(n), [k / sqrt(4k^2 - 1) for k in 1:(n - 1)]))
    (E.values, 2 .* E.vectors[1, :] .^ 2)
end

# centred finite-difference stencils (offsets -h:h) for g_xixi and g_xi at spacing 1
const _FD_STENCIL = Dict(
    4 => (c2 = [-1.0, 16.0, -30.0, 16.0, -1.0] ./ 12, c1 = [1.0, -8.0, 0.0, 8.0, -1.0] ./ 12),
    6 => (c2 = [2.0, -27.0, 270.0, -490.0, 270.0, -27.0, 2.0] ./ 180, c1 = [-1.0, 9.0, -45.0, 0.0, 45.0, -9.0, 1.0] ./ 60),
)

"""
    RadialGrid

Radial nodes of a simulation. Build it with [`uniform_grid`](@ref) or [`graded_grid`](@ref).

Fields: `nr` unknown nodes; `r[i+1] = r_i` for `i = 0 .. nr+3` (disk radii; `r[nr+1]` is the wall); `s[i+1] = dr/dxi`
and `sp[i+1] = d²r/dxi²` at the nodes `0 .. nr+2`; `dr0` the spacing in the disk/contact region; `sr = 1/dr0` nodes per
disk radius (`spatialResolution`); `R_req` the requested domain radius; `uniform`; the map parameters `q`, `dr_max`, `r_in`,
`xin`, `ramp`, `smooth`; `laplacian_order` (2 = conservative and tridiagonal, 4 or 6 = finite difference, banded, stored dense).
"""
struct RadialGrid
    nr::Int
    r::Vector{Float64}
    s::Vector{Float64}
    sp::Vector{Float64}
    dr0::Float64
    sr::Float64
    R_req::Float64
    uniform::Bool
    q::Float64
    dr_max::Float64
    r_in::Float64
    xin::Int
    ramp::Float64
    smooth::Float64
    laplacian_order::Int
end

function Base.show(io::IO, g::RadialGrid)
    if g.uniform
        print(io, "RadialGrid(uniform, nr=", g.nr, ", dr=", g.dr0, ")")
    else
        print(io, "RadialGrid(graded, nr=", g.nr, ", dr0=", g.dr0, ", q=", g.q, ", dr_max=", g.dr_max, ", r_in=", g.r_in,
              ", ramp=", g.ramp, ", smooth=", g.smooth, ", R=", g.r[g.nr + 1], ", laplacian_order=", g.laplacian_order, ")")
    end
end

"Radius of the wall node (the first node with `phi = 0`), in disk radii."
wall_radius(g::RadialGrid) = g.r[g.nr + 1]

"""
    uniform_grid(a_over_dr0, R) -> RadialGrid

Today's grid: `nr = ceil(a_over_dr0 * R)` nodes, spacing `dr = R/nr` (disk radii). `R = bathDiameter/2`,
`a_over_dr0 = spatialResolution`. `build_domain(uniform_grid(sr, D/2))` and `generate_dtn(uniform_grid(sr, D/2))` are
bit-identical to `build_domain(D, sr)` and `generate_dtn(nr, D)`.
"""
function uniform_grid(a_over_dr0::Real, R::Real)
    sr = float(a_over_dr0)
    Rf = float(R)
    D = 2 * Rf
    nr = ceil(Int, D * sr / 2)
    dr = D / (2 * nr)          # same expressions as build_domain
    r = [dr * i for i in 0:(nr + 3)]
    return RadialGrid(nr, r, fill(dr, nr + 3), zeros(nr + 3), dr, sr, Rf, true, 1.0, dr, 0.0, 0, 0.0, 0.0, 2)
end

# smooth step from 0 (x <= 0) to 1 (x >= 1), C^3, and its integral
@inline function _grid_rho(x::Float64)
    x <= 0 && return 0.0
    x >= 1 && return 1.0
    return x^4 * (35 - 84x + 70x^2 - 20x^3)
end

@inline function _grid_Pbar(x::Float64)
    x <= 0 && return 0.0
    x >= 1 && return 0.5 + (x - 1)
    return 7x^5 - 14x^6 + 10x^7 - 2.5x^8
end

# log-spacing l and its slope l' at t = xi - xi_in
function _grid_log_spacing(t::Float64, lnq::Float64, m::Float64, lmax::Float64, beta::Float64)
    t <= 0 && return (0.0, 0.0)
    x = t / m
    lg = lnq * m * _grid_Pbar(x)
    z = beta * (lg - lmax)
    softplus = max(z, 0.0) + log1p(exp(-abs(z)))
    return (lg - softplus / beta, lnq * _grid_rho(x) / (1 + exp(z)))
end

# s = r_xi and s' = r_xixi at the (real) index xi, from the map parameters
function _grid_spacing(xi::Float64, dr0::Float64, xin::Int, lnq::Float64, m::Float64, lmax::Float64, beta::Float64)
    l, lp = _grid_log_spacing(xi - xin, lnq, m, lmax, beta)
    s = dr0 * exp(l)
    return (s, s * lp)
end

# integral of s over [a, b] (Gauss-Legendre, 10 points)
function _grid_integrate_s(a::Float64, b::Float64, dr0, xin, lnq, m, lmax, beta)
    xg, wg = _GRID_GL
    h = (b - a) / 2
    c = (a + b) / 2
    acc = 0.0
    for k in eachindex(xg)
        acc += wg[k] * _grid_spacing(c + h * xg[k], dr0, xin, lnq, m, lmax, beta)[1]
    end
    return h * acc
end

_grid_params(g::RadialGrid) = (g.dr0, g.xin, log(g.q), g.ramp, log(g.dr_max / g.dr0), 1 / g.smooth)

"""
    graded_grid(a_over_dr0, R; q = 1.04, dr_max, r_in = 2.0, ramp = 20, smooth = 0.1, laplacian_order = 2) -> RadialGrid

Graded grid for the domain radius `R` (disk radii), with `a_over_dr0` nodes per disk radius in the disk region. Identity map
(`dr0 = 1/a_over_dr0`) for `r <= r_in` (needs `r_in >= 1 + 5 dr0`); then spacing growth `q` per node, capped at `dr_max`
(disk radii). `ramp` is the length (nodes) of the smooth start of the growth, `smooth` the width (in `ln` spacing) of the cap.
The wall is the first node with `r >= R`. If `dr_max <= dr0`, `q <= 1` or `R <= r_in`, the result is `uniform_grid`.
`laplacian_order` is 2 (conservative, tridiagonal), 4 or 6 (finite difference, banded; not tridiagonal).
"""
function graded_grid(a_over_dr0::Real, R::Real; q::Real = 1.04, dr_max::Real, r_in::Real = 2.0, ramp::Real = 20.0,
                     smooth::Real = 0.1, laplacian_order::Int = 2)
    laplacian_order in (2, 4, 6) || throw(ArgumentError("laplacian_order must be 2, 4 or 6, got $laplacian_order"))
    sr, Rf, qf, dmax, rin, rampf, smf = float(a_over_dr0), float(R), float(q), float(dr_max), float(r_in), float(ramp), float(smooth)
    dr0 = 1 / sr
    if dmax <= dr0 * (1 + 1e-9) || qf <= 1 || Rf <= rin
        g = uniform_grid(sr, Rf)
        return RadialGrid(g.nr, g.r, g.s, g.sp, g.dr0, g.sr, g.R_req, true, 1.0, g.dr0, 0.0, 0, 0.0, 0.0, laplacian_order)
    end
    rin >= 1 + 5 * dr0 || throw(ArgumentError("r_in must be at least 1 + 5 dr0 = $(1 + 5 * dr0) (contact nodes stay in the identity region), got $rin"))
    rampf > 0 && smf > 0 || throw(ArgumentError("ramp and smooth must be positive"))
    xin = ceil(Int, rin / dr0)
    lnq, lmax, beta = log(qf), log(dmax / dr0), 1 / smf
    r = Float64[0.0]
    i = 0
    while r[end] < Rf
        i += 1
        i > 20_000_000 && error("graded_grid: more than 2e7 nodes (check R and dr_max)")
        push!(r, i <= xin ? dr0 * i : r[end] + _grid_integrate_s(i - 1.0, Float64(i), dr0, xin, lnq, rampf, lmax, beta))
    end
    nr = length(r) - 1
    for _ in 1:3
        i += 1
        push!(r, i <= xin ? dr0 * i : r[end] + _grid_integrate_s(i - 1.0, Float64(i), dr0, xin, lnq, rampf, lmax, beta))
    end
    s = Vector{Float64}(undef, nr + 3)
    sp = Vector{Float64}(undef, nr + 3)
    for k in 0:(nr + 2)
        s[k + 1], sp[k + 1] = _grid_spacing(Float64(k), dr0, xin, lnq, rampf, lmax, beta)
    end
    return RadialGrid(nr, r, s, sp, dr0, sr, Rf, false, qf, dmax, rin, xin, rampf, smf, laplacian_order)
end

# s and s' at a real index (works for both kinds of grid)
function _grid_spacing_at(g::RadialGrid, xi::Float64)
    g.uniform && return (g.dr0, 0.0)
    dr0, xin, lnq, m, lmax, beta = _grid_params(g)
    return _grid_spacing(xi, dr0, xin, lnq, m, lmax, beta)
end

# radius at the half node xi = j + 1/2
function _grid_r_half(g::RadialGrid, j::Int)
    (g.uniform || j < g.xin) && return g.dr0 * (j + 0.5)
    dr0, xin, lnq, m, lmax, beta = _grid_params(g)
    return g.r[j + 1] + _grid_integrate_s(Float64(j), j + 0.5, dr0, xin, lnq, m, lmax, beta)
end

"""
    capillary_wavelength_cm(freqHz, density, tension, g) -> Float64

Wavelength (cm) of the deep-water gravity-capillary wave of frequency `freqHz`: `omega^2 = g k + (sigma/rho) k^3`.
CGS units (g/cm^3, dyne/cm, cm/s^2).
"""
function capillary_wavelength_cm(freqHz::Real, density::Real, tension::Real, g::Real)
    w2 = (2pi * freqHz)^2
    s = tension / density
    f(k) = g * k + s * k^3 - w2
    lo, hi = 0.0, 1.0
    while f(hi) < 0
        hi *= 2
    end
    for _ in 1:200
        mid = (lo + hi) / 2
        f(mid) < 0 ? (lo = mid) : (hi = mid)
    end
    return 2pi / ((lo + hi) / 2)
end

"""
    boundary_node_count(grid; width = 2.0, minimum_nodes = 3) -> Int

Number of outermost nodes used for the wall-echo diagnostic `eta_boundary`: the nodes within `width` disk radii of the wall
(`2 * spatialResolution` nodes on the uniform grid, as before), but at least `minimum_nodes`.
"""
function boundary_node_count(g::RadialGrid; width::Real = 2.0, minimum_nodes::Int = 3)
    g.uniform && return min(g.nr, ceil(Int, 2 * g.sr))
    rw = g.r[g.nr + 1]
    n = count(i -> g.r[i + 1] >= rw - width, 0:(g.nr - 1))
    return min(g.nr, max(minimum_nodes, n))
end

"""
    boundary_stencil(grid) -> (j0, w, radii)

Cubic-interpolation stencil of the wall-echo diagnostic `eta_boundary` on a graded grid. The uniform grid samples `eta` at the
nodes `wall - k dr`, `k = 1 .. ceil(2 sr)` (the outermost two disk radii). Here the same physical radii `x_k = r_wall - k dr0`
are used, whatever the local spacing: `eta(x_k) = sum_m w[k, m] eta[j0[k] + m]`, `m = 1..4` (0-based node `j0 + m - 1`), the
Lagrange cubic through the 4 nodes `j0 .. j0+3` around `x_k` (the same piecewise cubic as the DtN: nodes `q-1 .. q+2` of the panel
containing `x_k`, shifted to stay inside `0 .. nr-1`; in the last panel it extrapolates from nodes `nr-4 .. nr-1`).
"""
function boundary_stencil(g::RadialGrid)
    nb = ceil(Int, 2 * g.sr)
    rw = g.r[g.nr + 1]
    radii = [rw - k * g.dr0 for k in 1:nb]
    j0 = zeros(Int, nb)
    w = zeros(nb, 4)
    for (k, x) in enumerate(radii)
        q = searchsortedlast(view(g.r, 1:g.nr), x) - 1          # 0-based panel with r_q <= x
        j = clamp(q - 1, 0, g.nr - 4)
        j0[k] = j
        for m in 1:4
            num = 1.0
            den = 1.0
            for n in 1:4
                n == m && continue
                num *= x - g.r[j + n]
                den *= g.r[j + m] - g.r[j + n]
            end
            w[k, m] = num / den
        end
    end
    return j0, w, radii
end
