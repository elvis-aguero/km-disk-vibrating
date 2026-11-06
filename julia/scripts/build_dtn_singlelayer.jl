#= Build a DtN cache with the single-layer method (cubic interpolation, finite-difference order 4 or 6) and write it
in the production format (physical units: A_prod = A_unit / dr, dr = D / (2 nr)). TEST ONLY.
Usage: julia julia/scripts/build_dtn_singlelayer.jl <spatialResolution> <bathDiameter> <ord> <out.jld2> =#

include(joinpath(@__DIR__, "_bootstrap.jl"))
using LinearAlgebra, Printf, JLD2

function agm(a, b)
    for _ in 1:40
        a, b = (a + b) / 2, sqrt(a * b)
        abs(a - b) <= 1e-15 * a && break
    end
    return a
end
ellipK_b(b) = pi / (2 * agm(1.0, b))            # K(m) with b = sqrt(1 - m)

function gauss_legendre(n)
    E = eigen(SymTridiagonal(zeros(n), [k / sqrt(4k^2 - 1) for k in 1:(n - 1)]))
    return E.values, 2 .* E.vectors[1, :] .^ 2
end

kern(r, rp) = (2 / pi) * rp / (r + rp) * ellipK_b(abs(r - rp) / (r + rp))

# moments I_n = Int kern (r' - a)^n dr', n = 0..N, over [a, a+1]; `sing` = end with the log singularity
function panel_moments(r, a, xg, wg, sing, N)
    I = zeros(N + 1)
    for (x, w) in zip(xg, wg)
        u = (x + 1) / 2
        s, ds = sing === :left ? (u^4, 4u^3) : sing === :right ? (1 - (1 - u)^4, 4 * (1 - u)^3) : (u, 1.0)
        kv = kern(r, a + s)
        wt = (w / 2) * ds * kv
        for n in 0:N
            I[n + 1] += wt * s^n
        end
    end
    return I
end

# coefficients (powers t^0..t^3) of the Lagrange polynomial for node tp among t = -1, 0, 1, 2
function lagrange_coeffs(tp)
    c = [1.0]
    den = 1.0
    for tj in (-1, 0, 1, 2)
        tj == tp && continue
        c = [(i <= length(c) ? -tj * c[i] : 0.0) + (i > 1 ? c[i - 1] : 0.0) for i in 1:(length(c) + 1)]
        den *= (tp - tj)
    end
    return c ./ den
end

# rows m = 0..nr+1 (r = m), columns j = 1..nr (node r = j - 1)
function single_layer(nr, basis; ng = 24)
    xg, wg = gauss_legendre(ng)
    N = basis === :cubic ? 3 : 1
    G = zeros(nr + 3, nr)
    LC = [lagrange_coeffs(tp) for tp in (-1, 0, 1, 2)]
    for q in 0:(nr - 2), m in 0:(nr + 2)
        sing = m == q ? :left : (m == q + 1 ? :right : :none)
        I = panel_moments(Float64(m), Float64(q), xg, wg, sing, N)
        if basis === :linear
            G[m + 1, q + 1] += I[1] - I[2]
            G[m + 1, q + 2] += I[2]
        else
            for (p, tp) in enumerate((-1, 0, 1, 2))
                col = abs(q + tp) + 1                 # even extension about r = 0
                col > nr && continue                  # f = 0 beyond the last node
                G[m + 1, col] += sum(LC[p][n + 1] * I[n + 1] for n in 0:3)
            end
        end
    end
    return G
end

# radial Laplacian, node m, finite-difference order 4 or 6; `gr(m)` returns g at r = m (a row of G or a number)
function lap(gr, m, ord = 4)
    if ord == 4
        m == 0 && return (-gr(2) .+ 16 .* gr(1) .- 15 .* gr(0)) ./ 3
        g2 = (-gr(m + 2) .+ 16 .* gr(m + 1) .- 30 .* gr(m) .+ 16 .* gr(m - 1) .- gr(m - 2)) ./ 12
        g1 = (gr(m - 2) .- 8 .* gr(m - 1) .+ 8 .* gr(m + 1) .- gr(m + 2)) ./ 12
        return g2 .+ g1 ./ m
    end
    m == 0 && return (4 .* gr(3) .- 54 .* gr(2) .+ 540 .* gr(1) .- 490 .* gr(0)) ./ 90
    g2 = (2 .* gr(m - 3) .- 27 .* gr(m - 2) .+ 270 .* gr(m - 1) .- 490 .* gr(m) .+ 270 .* gr(m + 1) .- 27 .* gr(m + 2) .+ 2 .* gr(m + 3)) ./ 180
    g1 = (-gr(m - 3) .+ 9 .* gr(m - 2) .- 45 .* gr(m - 1) .+ 45 .* gr(m + 1) .- 9 .* gr(m + 2) .+ gr(m + 3)) ./ 60
    return g2 .+ g1 ./ m
end

function dtn_from_single_layer(nr, basis; ord = 4)
    G = single_layer(nr, basis)
    A = zeros(nr, nr)
    for m in 0:(nr - 1)
        A[m + 1, :] .= -lap(k -> view(G, abs(k) + 1, :), m, ord)
    end
    return A, G
end


sr, D, ord, out = parse(Float64, ARGS[1]), parse(Float64, ARGS[2]), parse(Int, ARGS[3]), ARGS[4]
nr = ceil(Int, sr * D / 2)
t = @elapsed ((A, _) = dtn_from_single_layer(nr, :cubic; ord = ord))
A ./= (D / (2nr))
FaradayDisk.check_dtn(A; source = "single-layer")
JLD2.jldsave(out; DTN = A)
@printf("built nr=%d ord=%d in %.1f s -> %s\n", nr, ord, t, out)
