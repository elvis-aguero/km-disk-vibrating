#= Error of a DtN matrix as a function of WAVELENGTH, against the exact operator.
Test functions: wave packets f(r) = exp(-r^2/(2 s^2)) J0(k0 r), the shape of the paper's own test (Appendix B.1).
Exact: the Hankel transform of f is F(k) = s^2 exp(-s^2 (k-k0)^2/2) * exp(-x) I0(x), x = s^2 k k0, and
(N f)(r) = Int_0^inf k^2 F(k) J0(k r) dk. Units: node spacing 1 (generate_dtn(nr, 2nr)); wavelength in nodes.
Usage: julia julia/scripts/pilot_dtn_spectrum.jl <nr> =#

include(joinpath(@__DIR__, "_bootstrap.jl"))
using LinearAlgebra, Printf

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

# J0 by the periodic trapezoid rule of (1/pi) Int_0^pi cos(x sin t) dt (spectrally accurate for N > x)
function besselj0(x)
    N = max(64, ceil(Int, 2 * abs(x)) + 64)
    s = 0.0
    for i in 0:(N - 1)
        s += cos(x * sin(pi * (i + 0.5) / N))
    end
    return s / N
end


function scaled_I0(x)                      # exp(-x) I0(x), periodic midpoint rule of (1/pi) Int_0^pi exp(x (cos t - 1)) dt
    N = max(2000, ceil(Int, 40 * sqrt(x)) + 200)
    s = 0.0
    for i in 0:(N - 1)
        s += exp(x * (cos(pi * (i + 0.5) / N) - 1))
    end
    return s / N
end

function exact_packet(s, k0, rs)
    dk = 0.4 / maximum(rs)
    ks = collect(0:dk:(k0 + 12 / s))
    Fk = [s^2 * exp(-s^2 * (k - k0)^2 / 2) * scaled_I0(s^2 * k * k0) for k in ks]
    out = zeros(length(rs))
    for (ir, r) in enumerate(rs)
        acc = 0.0
        for (ik, k) in enumerate(ks)
            acc += (ik == 1 ? 0.5 : 1.0) * k^2 * Fk[ik] * besselj0(k * r)
        end
        out[ir] = acc * dk
    end
    return out
end

function main(nr)
    s = nr / 8
    rs = collect(0.0:(nr - 1))
    t_ref = @elapsed Aref = FaradayDisk.generate_dtn_polar(nr, 2.0 * nr)
    (A4, _) = dtn_from_single_layer(nr, :cubic; ord = 4)
    (A6, _) = dtn_from_single_layer(nr, :cubic; ord = 6)
    # ring-area weights (unit spacing): W A is symmetric if A is self-adjoint in the L2 inner product of the plane
    W = [m == 0 ? pi * 0.25 : 2pi * m for m in 0:(nr - 1)]
    symW(A) = (A .+ (A' .* W') ./ W) ./ 2          # (A + W^-1 A' W) / 2
    asym(A) = (M = W .* A; norm(M - M') / norm(M))
    Aref_s, A6_s = symW(Aref), symW(A6)
    @printf("nr = %d, envelope s = %g nodes, reference build %.0f s\n", nr, s, t_ref)
    @printf("asymmetry |WA-(WA)'|/|WA|: ref %.2e   new(6) %.2e   ref sym %.2e   new(6) sym %.2e\n", asym(Aref), asym(A6), asym(Aref_s), asym(A6_s))
    @printf("%-12s %-10s | %-10s %-10s | %-10s %-10s\n", "wavelength", "k*dr", "ref", "ref sym", "new(6)", "new(6) sym")
    for lam in (2.5, 3.0, 4.0, 6.0, 8.0, 12.0, 16.0, 24.0, 32.0, 48.0, 64.0, 96.0, 150.0)
        k0 = 2pi / lam
        f = @. exp(-rs^2 / (2 * s^2)) * besselj0(k0 * rs)
        ex = exact_packet(s, k0, rs)
        nrm = maximum(abs, ex)
        e(A) = maximum(abs, A * f .- ex) / nrm
        @printf("%-12g %-10.3f | %-10.2e %-10.2e | %-10.2e %-10.2e\n", lam, k0, e(Aref), e(Aref_s), e(A6), e(A6_s))
        flush(stdout)
    end
end

main(isempty(ARGS) ? 200 : parse(Int, ARGS[1]))
