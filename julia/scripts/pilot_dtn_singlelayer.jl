#= PILOT (option D): build the DtN matrix by a different, much cheaper method and compare it with
the library generator (`generate_dtn`, the reference).

The generator's kernel is the free-space, infinite-depth operator |k| (symbol of the 2D
Dirichlet-to-Neumann map), restricted to an axisymmetric function that is zero outside r < R.
It is NOT Toeplitz in the radial index, but |k| = k^2 * (1/|k|), so
    DtN f = -Laplacian( S f ),   (S f)(x) = (1/2pi) Int f(x') / |x - x'| d2x'.
S is a weakly (log) singular integral. Axisymmetric reduction (K = complete elliptic integral):
    (S f)(r) = (2/pi) Int_0^R f(r') r' / (r + r') K(m) dr',   m = 4 r r' / (r + r')^2.
Steps: (1) G = nodal matrix of S for piecewise-linear f (panel Gauss-Legendre, graded map at
the log singularity, K by AGM); (2) A = -L G with the radial Laplacian by 4th-order finite
differences (even extension at r = 0). Cost O(nr^2) kernel evaluations, no angular quadrature.

Units: node spacing dr = 1 (so generate_dtn(nr, D = 2 nr)).
Test: Gaussian f = exp(-r^2 / (2 s^2)), exact DtN by direct k-integration
    (DtN f)(r) = s^2 Int_0^inf k^2 exp(-k^2 s^2 / 2) J0(k r) dk.
The truncation convention at r = R does not matter (the Gaussian is negligible there).
Variants: piecewise-LINEAR or piecewise-CUBIC (4-point Lagrange, even extension at r=0, f=0 beyond the last
node) interpolation of f; Richardson extrapolation between nr and nr/2 (even nodes only).
Diagnostic: error of S f alone (G f vs exact), and of the finite differences alone (applied to the exact S f).
Usage: julia julia/scripts/pilot_dtn_singlelayer.jl [nr1 nr2 ...] =#

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

# exact Hankel integrals for the Gaussian: pow = 2 -> DtN f, pow = 0 -> S f
function exact_gaussian(s, rs, pow)
    dk = 0.4 / maximum(rs)
    ks = 0:dk:(12 / s)
    out = zeros(length(rs))
    for (ir, r) in enumerate(rs)
        acc = 0.0
        for (ik, k) in enumerate(ks)
            acc += (ik == 1 ? 0.5 : 1.0) * k^pow * exp(-k^2 * s^2 / 2) * besselj0(k * r)
        end
        out[ir] = s^2 * acc * dk
    end
    return out
end

relmax(a, b, nrm) = maximum(abs, a .- b) / nrm

# matrix-level checks: ring-weighted symmetry (W A should be symmetric positive semi-definite), regional difference to the reference
function matrix_report(nr, Aref, Ac4, Ac6)
    w = [m == 0 ? pi * 0.25 : 2pi * m for m in 0:(nr - 1)]
    println("\nC) matrix level, nr=$nr (units dr=1)")
    for (name, A) in (("ref", Aref), ("cub4", Ac4), ("cub6", Ac6))
        M = w .* A
        S = (M .+ M') ./ 2
        @printf("  %-5s asym |M-M'|/|M| = %.2e   min eig of sym(M) / max = %.2e\n", name, norm(M - M') / norm(M), minimum(eigvals(Symmetric(S))) / maximum(eigvals(Symmetric(S))))
    end
    for (name, A) in (("cub4", Ac4), ("cub6", Ac6))
        d(rows) = norm((A .- Aref)[rows, :]) / norm(Aref[rows, :])
        @printf("  %-5s rel diff to ref: all %.2e | rows 0-4 %.2e | middle %.2e | last 5 rows %.2e\n", name, d(1:nr), d(1:5), d(div(nr, 2):(div(nr, 2) + 5)), d((nr - 4):nr))
    end
end

function main(nrs)
    println("A) error vs EXACT DtN (relative to max|exact|)")
    @printf("%-4s %-3s | %-9s | %-9s %-9s %-9s | %-9s %-9s | %-8s %-8s\n", "nr", "s", "ref", "lin tot", "lin Gpart", "FD part", "cub tot", "cub Gpart", "t_ref", "t_new")
    rich = String[]
    for nr in nrs
        t_ref = nr <= 200 ? (@elapsed Aref = FaradayDisk.generate_dtn(nr, 2.0 * nr)) : NaN
        t_new = @elapsed ((Al, Gl) = dtn_from_single_layer(nr, :linear))
        (Ac, Gc) = dtn_from_single_layer(nr, :cubic)
        (Ac6, _) = dtn_from_single_layer(nr, :cubic; ord = 6)
        nr <= 200 && matrix_report(nr, Aref, Ac, Ac6)
        nc = nr ÷ 2
        (Alc, _) = dtn_from_single_layer(nc, :linear)
        (Acc, _) = dtn_from_single_layer(nc, :cubic)
        for s in (6.0, 12.0)
            8s <= nr - 1 || continue
            rs = collect(0.0:(nr + 2))
            f = @. exp(-rs[1:nr]^2 / (2 * s^2))
            ex = exact_gaussian(s, rs[1:nr], 2)
            gx = exact_gaussian(s, rs, 0)
            nrm = maximum(abs, ex)
            fd = [-lap(k -> gx[abs(k) + 1], m, 4) for m in 0:(nr - 1)]
            fd6 = [-lap(k -> gx[abs(k) + 1], m, 6) for m in 0:(nr - 1)]
            e_ref = nr <= 200 ? relmax(Aref * f, ex, nrm) : NaN
            e_l = relmax(Al * f, ex, nrm)
            e_lG = relmax(Gl[1:nr, :] * f, gx[1:nr], maximum(abs, gx))
            e_fd = relmax(fd, ex, nrm)
            e_c = relmax(Ac * f, ex, nrm)
            e_cG = relmax(Gc[1:nr, :] * f, gx[1:nr], maximum(abs, gx))
            @printf("%-4d %-3g | %-9.2e | %-9.2e %-9.2e %-9.2e | %-9.2e %-9.2e | %-8.1f %-8.2f\n", nr, s, e_ref, e_l, e_lG, e_fd, e_c, e_cG, t_ref, t_new)
            push!(rich, @sprintf("C0 nr=%d s=%g: cubic + 6th-order FD: total %.2e  (FD6 part alone %.2e)", nr, s, relmax(Ac6 * f, ex, nrm), relmax(fd6, ex, nrm)))
            # Richardson on even nodes: fine (spacing 1) vs coarse (spacing 2, matrix/2)
            ev = 1:2:(2nc - 1)
            fc = f[ev]
            for (name, A, Ac_, pth) in (("lin", Al, Alc, 2), ("cub", Ac, Acc, 4))
                af = (A*f)[ev]
                ac = (Ac_ * fc) ./ 2
                ef = relmax(af, ex[ev], nrm)
                ecs = relmax(ac, ex[ev], nrm)
                R = (2^pth .* af .- ac) ./ (2^pth - 1)
                push!(rich, @sprintf("nr=%d s=%g %s: err fine %.2e  err coarse %.2e  p_obs %.2f  Richardson(p=%d) %.2e", nr, s, name, ef, ecs, log2(ecs / ef), pth, relmax(R, ex[ev], nrm)))
            end
        end
        flush(stdout)
    end
    println("\nB) Richardson between nr and nr/2 (even nodes only)")
    foreach(println, rich)
end

main(isempty(ARGS) ? [100, 200] : parse.(Int, ARGS))
