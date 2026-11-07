# Shared numerical checks for the graded grid (docs/graded_grid.md): used by the tests and by julia/scripts/graded_pilot.jl.
# Plain definitions (no module): they use `graded_grid`, `build_domain` of FaradayDisk, which the includer has loaded.

# J0 by the periodic trapezoid rule of (1/π) ∫_0^π cos(x sin t) dt
function _gg_besselj0(x)
    N = max(64, ceil(Int, 2 * abs(x)) + 64)
    return sum(cos(x * sin(pi * (i + 0.5) / N)) for i in 0:(N - 1)) / N
end

# exp(-x) I0(x), periodic midpoint rule of (1/π) ∫_0^π exp(x (cos t - 1)) dt
function _gg_scaled_i0(x)
    N = max(2000, ceil(Int, 40 * sqrt(x)) + 200)
    return sum(exp(x * (cos(pi * (i + 0.5) / N) - 1)) for i in 0:(N - 1)) / N
end

# exact (N f)(r) for f = exp(-r²/2s²) J0(k0 r): ∫_0^∞ k² F(k) J0(k r) dk, F = s² exp(-s²(k-k0)²/2) exp(-x) I0(x), x = s² k k0
function _gg_exact_packet(s, k0, rs)
    dk = 0.2 / maximum(rs)
    ks = collect(0.0:dk:(k0 + 14 / s))
    F = [s^2 * exp(-s^2 * (k - k0)^2 / 2) * _gg_scaled_i0(s^2 * k * k0) for k in ks]
    return [dk * sum(ks[i]^2 * F[i] * _gg_besselj0(ks[i] * r) for i in eachindex(ks)) for r in rs]
end

# E_loc: max_i |(A f)_i - exact_i| / max_{|r_j - r_i| <= λ} |exact_j|, over nodes whose local scale is >= 1e-3 max|exact|
function _gg_local_error(A, f, exact, rs, lam)
    glob = maximum(abs, exact)
    worst = 0.0
    for i in eachindex(rs)
        sc = maximum(abs(exact[j]) for j in eachindex(rs) if abs(rs[j] - rs[i]) <= lam)
        sc >= 1e-3 * glob || continue
        worst = max(worst, abs(sum(A[i, :] .* f) - exact[i]) / sc)
    end
    return worst
end

# Laplacian error on g = exp(-r^2/2w^2) for the refinement level p: the map is the same continuous function for every p
# (sr -> p sr, q -> q^(1/p), dr_max -> dr_max/p, ramp -> p ramp), only the sampling gets finer.
function _gg_laplacian_error(p::Int, form::Symbol)
    w = 5.0
    order = form === :cons2 ? 2 : (form === :fd4 ? 4 : 6)
    g = graded_grid(5 * p, 50.0; q = 1.15^(1 / p), dr_max = 1.0 / p, r_in = 2.0, ramp = 10.0 * p, smooth = 0.1,
                    laplacian_order = order)
    L = build_domain(g).Laplacian
    r = g.r[1:g.nr]
    gf = @. exp(-r^2 / (2 * w^2))
    ex = @. (r^2 / w^4 - 2 / w^2) * gf
    e = L * gf .- ex
    return maximum(abs, e[1:(g.nr - 6)]) / maximum(abs, ex)       # the last rows feel the wall stencil
end
