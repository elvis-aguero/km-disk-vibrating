"""
Single-layer construction of the dense Dirichlet-to-Neumann (DtN) matrix, O(nr^2) kernel evaluations.

Derivation, error tables and checks: `docs/dtn_single_layer.md`. In short:

  * `N = -Δ_H S`, with `S` the single-layer operator `(Sψ)(r) = (1/2π) ∫ ψ(s)/|r - s| dA(s)` (paper Eq. A 3 gives
    `S N = I`; the identity itself follows from Laplace's equation).
  * For axisymmetric ψ the angle integral is exact: `(Sψ)(r) = (2/π) ∫ r'/(r + r') K(m) ψ(r') dr'`,
    `m = 4 r r'/(r + r')^2`, `K` the complete elliptic integral of the first kind (arithmetic-geometric mean).
  * φ is interpolated by piecewise four-point Lagrange cubics (even extension about r = 0, zero beyond the last node);
    panel integrals use 24-point Gauss-Legendre with the graded substitution `r' = r_i ± h u^4` on the two panels that
    touch the logarithmic singularity at `r' = r_i`.
  * `Δ_H g = g'' + g'/r` is applied to the nodal values of `S φ` by centred finite differences of order 4 or 6.

Units inside this file: node spacing 1. `generate_dtn_single_layer` rescales by `1/dr`, `dr = D/(2 nr)`, to match
`generate_dtn_polar`. Rows are independent, so the work is threaded over rows (small per-row allocations only).
"""

# arithmetic-geometric mean; returns K(m) for b = sqrt(1 - m), 0 < b <= 1
function _ellipK_b(b::Float64)
    a = 1.0
    for _ in 1:40
        a, b = (a + b) / 2, sqrt(a * b)
        abs(a - b) <= 1e-15 * a && break
    end
    return pi / (2 * a)
end

# Gauss-Legendre nodes and weights on [-1, 1] (Golub-Welsch)
function _gauss_legendre(n::Int)
    E = eigen(SymTridiagonal(zeros(n), [k / sqrt(4k^2 - 1) for k in 1:(n - 1)]))
    return E.values, 2 .* E.vectors[1, :] .^ 2
end

const _SL_GL = _gauss_legendre(24)

_sl_kernel(r::Float64, rp::Float64) = (2 / pi) * rp / (r + rp) * _ellipK_b(abs(r - rp) / (r + rp))

# I[n+1] = ∫_a^{a+1} kernel(r, r') (r' - a)^n dr', n = 0..N, `sing` = the panel end with the log singularity
function _sl_panel_moments!(I::Vector{Float64}, r::Float64, a::Float64, sing::Symbol, N::Int)
    fill!(I, 0.0)
    xs, ws = _SL_GL
    for k in eachindex(xs)
        u = (xs[k] + 1) / 2
        s, ds = sing === :left ? (u^4, 4u^3) : sing === :right ? (1 - (1 - u)^4, 4 * (1 - u)^3) : (u, 1.0)
        wt = (ws[k] / 2) * ds * _sl_kernel(r, a + s)
        p = 1.0
        for n in 0:N
            I[n + 1] += wt * p
            p *= s
        end
    end
    return I
end

# coefficients (powers t^0..t^3) of the Lagrange cubic for node tp among t = -1, 0, 1, 2
function _sl_lagrange_coeffs(tp::Int)
    c = [1.0]
    den = 1.0
    for tj in (-1, 0, 1, 2)
        tj == tp && continue
        c = [(i <= length(c) ? -tj * c[i] : 0.0) + (i > 1 ? c[i - 1] : 0.0) for i in 1:(length(c) + 1)]
        den *= (tp - tj)
    end
    return c ./ den
end

# Matrix of S on nodal values (rows r = 0..nr+2, columns node r = 0..nr-1), node spacing 1.
function _sl_single_layer_matrix(nr::Int)
    LC = [_sl_lagrange_coeffs(tp) for tp in (-1, 0, 1, 2)]
    G = zeros(nr + 3, nr)
    Threads.@threads :dynamic for m in 0:(nr + 2)
        I = zeros(4)
        for q in 0:(nr - 2)
            sing = m == q ? :left : (m == q + 1 ? :right : :none)
            _sl_panel_moments!(I, Float64(m), Float64(q), sing, 3)
            for (p, tp) in enumerate((-1, 0, 1, 2))
                col = abs(q + tp) + 1          # even extension about r = 0
                col > nr && continue           # φ = 0 beyond the last node
                lc = LC[p]
                G[m + 1, col] += lc[1] * I[1] + lc[2] * I[2] + lc[3] * I[3] + lc[4] * I[4]
            end
        end
    end
    return G
end

# stencils for Δ_H g = g'' + g'/r at node r = m >= 1: offsets -h:h
const _SL_STENCIL = Dict(
    4 => (c2 = [-1.0, 16.0, -30.0, 16.0, -1.0] ./ 12, c1 = [1.0, -8.0, 0.0, 8.0, -1.0] ./ 12),
    6 => (c2 = [2.0, -27.0, 270.0, -490.0, 270.0, -27.0, 2.0] ./ 180, c1 = [-1.0, 9.0, -45.0, 0.0, 45.0, -9.0, 1.0] ./ 60),
)

"""
    generate_dtn_single_layer(nr, D; order=6) -> Matrix{Float64}

DtN matrix for `nr` nodes and bath diameter `D` (in disk radii), `dr = D/(2 nr)`, same operator and units as
[`generate_dtn_polar`](@ref). `order` is the finite-difference order (4 or 6) of the horizontal Laplacian.
"""
function generate_dtn_single_layer(nr::Int, D::Real; order::Int = 6)
    nr >= 8 || error("generate_dtn_single_layer requires nr >= 8, got nr=$nr")
    haskey(_SL_STENCIL, order) || throw(ArgumentError("order must be 4 or 6, got $order"))
    G = _sl_single_layer_matrix(nr)
    st = _SL_STENCIL[order]
    h = length(st.c2) ÷ 2
    A = zeros(nr, nr)
    Threads.@threads :dynamic for m in 0:(nr - 1)
        if m == 0                                           # Δ_H g(0) = 2 g''(0), even symmetry
            for j in 1:nr
                A[1, j] = order == 4 ? -(-G[3, j] + 16 * G[2, j] - 15 * G[1, j]) / 3 :
                                       -(4 * G[4, j] - 54 * G[3, j] + 540 * G[2, j] - 490 * G[1, j]) / 90
            end
        else
            for j in 1:nr
                acc = 0.0
                for k in -h:h
                    acc += (st.c2[k + h + 1] + st.c1[k + h + 1] / m) * G[abs(m + k) + 1, j]
                end
                A[m + 1, j] = -acc
            end
        end
    end
    return A ./ (D / (2nr))
end

"""
    generate_dtn(nr, D; method=:single_layer, order=6, refp=10) -> Matrix{Float64}

Dense DtN matrix of the domain (`nr` radial nodes, bath diameter `D` disk radii). `method = :single_layer` (default)
is [`generate_dtn_single_layer`](@ref), O(nr^2); `method = :polar` is the paper's scheme,
[`generate_dtn_polar`](@ref), O(nr^3) with a large constant (`refp` applies to it only).
"""
function generate_dtn(nr::Int, D::Real; method::Symbol = :single_layer, order::Int = 6, refp::Int = DTN_DEFAULT_REFP)
    method === :single_layer && return generate_dtn_single_layer(nr, D; order = order)
    method === :polar && return generate_dtn_polar(nr, D; refp = refp)
    throw(ArgumentError("method must be :single_layer or :polar, got :$method"))
end

# ---------------------------------------------------------------------------------------------------------------------
# Non-uniform nodes (graded grid, `docs/graded_grid.md`). Same construction in physical units (disk radii): panels
# [r_q, r_{q+1}] of width h_q, cubic Lagrange interpolation through r_{q-1}..r_{q+2} written in t = (r' - r_q)/h_q, graded
# substitution on the two singular panels, N = -Δ_H G with the mapped stencil Δ g = g_ξξ/s² + (1/(r s) - s'/s³) g_ξ.
# ---------------------------------------------------------------------------------------------------------------------

# I[n+1] = ∫_{panel} kernel(r, r') t^n h dt, n = 0..3, `t = (r' - a)/h`. For the singular panels the distance |r' - r| comes from
# the substitution (h u^4 or h (1-u)^4), not from subtracting two nearby radii.
function _sl_panel_moments_grid!(I::Vector{Float64}, r::Float64, a::Float64, h::Float64, sing::Symbol)
    fill!(I, 0.0)
    xs, ws = _SL_GL
    for k in eachindex(xs)
        u = (xs[k] + 1) / 2
        if sing === :left                       # r = a
            t, dt, d = u^4, 4u^3, h * u^4
        elseif sing === :right                  # r = a + h
            v = (1 - u)^4
            t, dt, d = 1 - v, 4 * (1 - u)^3, -h * v
        else
            t, dt = u, 1.0
            d = a + h * t - r
        end
        rp = a + h * t
        kv = (2 / pi) * rp / (r + rp) * _ellipK_b(abs(d) / (r + rp))
        wt = (ws[k] / 2) * dt * kv * h
        p = 1.0
        for n in 1:4
            I[n] += wt * p
            p *= t
        end
    end
    return I
end

# Matrix of S on nodal values (rows: nodes 0..nr+2, columns: nodes 0..nr-1), physical units.
function _sl_single_layer_matrix_grid(g::RadialGrid)
    nr, r = g.nr, g.r
    node(j) = j < 0 ? -r[-j + 1] : r[j + 1]            # odd extension about the axis
    LC = zeros(4, 4, nr - 1)                           # LC[n, p, q+1]: coefficient of t^(n-1) of the cubic of node offset p-2
    hq = zeros(nr - 1)
    for q in 0:(nr - 2)
        a = r[q + 1]
        h = r[q + 2] - a
        hq[q + 1] = h
        x = ((node(q - 1) - a) / h, 0.0, (r[q + 2] - a) / h, (r[q + 3] - a) / h)
        for p in 1:4
            c = [1.0]
            den = 1.0
            for j in 1:4
                j == p && continue
                xj = x[j]
                c = [(i <= length(c) ? -xj * c[i] : 0.0) + (i > 1 ? c[i - 1] : 0.0) for i in 1:(length(c) + 1)]
                den *= (x[p] - xj)
            end
            LC[:, p, q + 1] .= c ./ den
        end
    end
    G = zeros(nr + 3, nr)
    Threads.@threads :dynamic for m in 0:(nr + 2)
        I = zeros(4)
        rm = r[m + 1]
        for q in 0:(nr - 2)
            sing = m == q ? :left : (m == q + 1 ? :right : :none)
            _sl_panel_moments_grid!(I, rm, r[q + 1], hq[q + 1], sing)
            for (p, tp) in enumerate((-1, 0, 1, 2))
                col = abs(q + tp) + 1                  # even extension about r = 0
                col > nr && continue                   # φ = 0 beyond the last node
                G[m + 1, col] += LC[1, p, q + 1] * I[1] + LC[2, p, q + 1] * I[2] + LC[3, p, q + 1] * I[3] + LC[4, p, q + 1] * I[4]
            end
        end
    end
    return G
end

"""
    generate_dtn_single_layer(grid::RadialGrid; order=6, force_general=false) -> Matrix{Float64}

DtN matrix for a [`RadialGrid`](@ref), physical units (disk radii), `order` the finite-difference order (4 or 6) of the
mapped Laplacian applied to `S φ`. For a uniform grid this calls `generate_dtn_single_layer(nr, 2R; order)` unchanged
(bit-identical); `force_general = true` runs the non-uniform routine on it instead (used by the tests).
"""
function generate_dtn_single_layer(grid::RadialGrid; order::Int = 6, force_general::Bool = false)
    haskey(_SL_STENCIL, order) || throw(ArgumentError("order must be 4 or 6, got $order"))
    if grid.uniform && !force_general
        return generate_dtn_single_layer(grid.nr, 2 * grid.R_req; order = order)
    end
    nr = grid.nr
    nr >= 8 || error("generate_dtn_single_layer requires nr >= 8, got nr=$nr")
    G = _sl_single_layer_matrix_grid(grid)
    st = _SL_STENCIL[order]
    h = length(st.c2) ÷ 2
    A = zeros(nr, nr)
    Threads.@threads :dynamic for m in 0:(nr - 1)
        s = grid.s[m + 1]
        if m == 0                                           # Δ_H g(0) = 2 g_ξξ(0)/s², even symmetry
            for j in 1:nr
                acc = 0.0
                for k in -h:h
                    acc += st.c2[k + h + 1] * G[abs(k) + 1, j]
                end
                A[1, j] = -2 * acc / s^2
            end
        else
            A2 = 1 / s^2
            A1 = 1 / (grid.r[m + 1] * s) - grid.sp[m + 1] / s^3
            for j in 1:nr
                acc = 0.0
                for k in -h:h
                    acc += (A2 * st.c2[k + h + 1] + A1 * st.c1[k + h + 1]) * G[abs(m + k) + 1, j]
                end
                A[m + 1, j] = -acc
            end
        end
    end
    return A
end

"""
    generate_dtn(grid::RadialGrid; method=:single_layer, order=6, force_general=false) -> Matrix{Float64}

DtN matrix of a radial grid (in process, no cache). `method = :polar` needs a uniform grid.
"""
function generate_dtn(grid::RadialGrid; method::Symbol = :single_layer, order::Int = 6, force_general::Bool = false)
    method === :single_layer && return generate_dtn_single_layer(grid; order = order, force_general = force_general)
    if method === :polar
        grid.uniform || throw(ArgumentError("method :polar needs a uniform grid"))
        return generate_dtn_polar(grid.nr, 2 * grid.R_req)
    end
    throw(ArgumentError("method must be :single_layer or :polar, got :$method"))
end
