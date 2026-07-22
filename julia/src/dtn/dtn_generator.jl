"""
Native Julia port of `matlab/1_code/simulation_code/DTNVectorized.m`, which computes the
dense Dirichlet-to-Neumann (DTN) operator matrix for a given domain (`nr` radial nodes,
diameter `D` disk radii) via a singular boundary-integral quadrature.

This is a line-for-line-faithful translation of the MATLAB algorithm's index arithmetic
and polynomial coefficients (not a re-derivation), with two structural simplifications
that do not change the result: (1) MATLAB's `num_batches=20` chunking of the angular
quadrature (`l_vals`) existed only to bound per-*process* memory inside `parfor` workers
(each a separate OS process) — Julia's `Threads.@threads` runs in one shared address
space, so the angular quadrature is evaluated in one pass instead of 20 batches
(confirmed mathematically identical, not just "should be", now that this is verified
against a real MATLAB-generated reference — see STATUS below); (2) MATLAB's
`accumarray(idx, val, [nr,1])` scatter-add is replaced by `accumulate_masked!`, a small
helper with the same semantics.

Correctness is checked by `test/integration/test_dtn_golden_small.jl` against the legacy
`D5Quant20` MATLAB reference cache — treat that test as load-bearing, not incidental.

STATUS: fully verified. This port was originally written blind (no Julia install
available), and an initial comparison against the legacy `D5Quant20` cache
(`DTNnew345nr50D5refp10.mat`) showed a 75% relative discrepancy at `generate_dtn(50,
20.0)` — the `D` value the cache's *folder name* (`D5Quant20`) implies (spatialResolution
5, bathDiameter 20). A first (still-blind) investigation guessed the cache might instead
have been generated at `D=5` (matching the `.mat` filename literally) and found that
dropped the gap to ~2.15%, but left that residual "unexplained." Once Julia actually
became available in this environment, direct numerical comparison found the real bug:
`far_part!` computed its outer-product radial position without dividing `i_far` by
`refp` (unlike the structurally identical "near" part code, which did), inflating
`radn`/`idxs` roughly `refp`-fold and pushing nearly every far-field contribution outside
the valid `1:nr` range where it was silently masked out — i.e. every off-diagonal DTN
entry beyond a small near-diagonal band was coming out as exact `0.0` regardless of `D`.
Fixing that one missing `/refp` makes `generate_dtn(50, 5.0)` match the legacy
`D5Quant20` cache to ~2e-15 relative difference (machine precision) — not an
approximation, and not a coincidence. `generate_dtn(50, 20.0)` still disagrees (as it
should: this genuinely is a `D=5` matrix), consistent with `solve_motion.m`'s own
"Machine-specific patch for D5Quant20", which loads this exact `.mat` file by its literal
name instead of computing the filename from `bathDiameter` the way every other domain
does — i.e. the real MATLAB pipeline already knowingly uses a `D=5`-generated matrix for
the nominally-bathDiameter=20 domain. That is a genuine, pre-existing MATLAB-side
data/usage inconsistency (not introduced by this port, and not something to silently
"fix" by swapping in a properly-regenerated bathDiameter=20 matrix, which would change
simulation behavior relative to what the real pipeline has always produced) — see
`scripts/migrate_dtn_caches.jl` for how the registry key and the validation `D` are kept
separate to reflect this. `D25Quant200`'s cache (`DTNnew345nr2500D25refp10.mat`) shows the
identical filename-vs-folder mismatch pattern (D25 vs implied bathDiameter=200) but its
`nr=2500` is too expensive to regenerate natively here to confirm either way — flagged as
an open question in `migrate_dtn_caches.jl`, not assumed resolved by analogy.
"""

const DTN_DEFAULT_REFP = 10

@inline function _term0(i::Float64)
    (-i^2 / 2 - i - 1 / 3) * log((i + 1) / i) +
    (i^2 / 6 + i / 2 + 1 / 3) * (1 - i / (i + 1)) -
    (i + 0.5) / 6 + i / 2 + 1 / 2
end

@inline function _term1(i::Float64)
    (3 * i^2 / 2 + 2 * i - 1 / 2) * log((i + 1) / i) +
    (-i^2 / 2 - i + 1 / 2 + 1 / i) * (1 - i / (i + 1)) +
    (i + 0.5) / 2 - 3 * i / 2 - 1
end

@inline function _term2(i::Float64)
    (-3 * i^2 / 2 - i + 1) * log((i + 1) / i) +
    (i^2 / 2 + i / 2 - 1) * (1 - i / (i + 1)) -
    (i + 0.5) / 2 + 3 * i / 2 + 1 / 2
end

@inline function _term3(i::Float64)
    (i^2 / 2 - 1 / 6) * log((i + 1) / i) +
    (-i^2 + 1) / 6 * (1 - i / (i + 1)) +
    (i + 0.5) / 6 - i / 2
end

# 2*(1/(i-1/2) - 1/(i+1/2)) elementwise, matching MATLAB's `Kern*` vectors.
_near_far_kernel(i_range::AbstractVector{<:Real}) = @. 2.0 * (1.0 / (i_range - 0.5) - 1.0 / (i_range + 0.5))

"""
    accumulate_masked!(Line, vals, mask, idx, nr)

Scatter-add `vals[j]` into `Line[Int(idx[j])]` for every `j` with `mask[j]` true,
mirroring MATLAB's `accumarray(idx(mask), vals(mask), [nr,1])'`. `vals`, `mask`, `idx`
must be same-shaped arrays (a `length(i_far) x length(l_vals)` matrix in every call site
below). Bounds are checked defensively — the masks at every call site are constructed so
`idx` always lands in `1:nr`, but a silently-dropped or silently-wrapped out-of-range
write would be a far worse failure mode than an explicit bounds error.
"""
function accumulate_masked!(Line::Vector{Float64}, vals::AbstractArray, mask::AbstractArray{Bool}, idx::AbstractArray, nr::Int)
    @inbounds for j in eachindex(mask)
        if mask[j]
            ix = Int(idx[j])
            1 <= ix <= nr || error("DTN quadrature index $ix out of bounds 1:$nr (this indicates a transcription bug, not bad input data)")
            Line[ix] += vals[j]
        end
    end
    return Line
end

"""
    far_part!(Line, rn_k, i_far, Kern_far, cos_l, sin_l, nr, refp)

The "far" (piecewise-cubic-interpolant) contribution to a DTN row, shared verbatim across
rows 2, 3, and every row k>=4 in the MATLAB source (`LineUpdate2`, `LineUpdate3`,
`LineUpdatek` all use the identical polynomial coefficients — only the row-dependent
`rn_k`, `i_far`, `Kern_far` differ). Adds in place to `Line` (length `nr`).
"""
function far_part!(Line::Vector{Float64}, rn_k::Float64, i_far::Vector{Float64}, Kern_far::Vector{Float64},
                    cos_l::Vector{Float64}, sin_l::Vector{Float64}, nr::Int, refp::Int)
    # Outer products: (length(i_far)) x (length(l_vals)). MATLAB: `(rn(k) + i_vals*cos_l/refp)`
    # / `(i_vals*sin_l/refp)` — the `/refp` here is not optional bookkeeping, it converts
    # `i_far` (an index into the refp-fold-refined radial sub-grid) back into physical
    # radial-node units before adding to `rn_k`; omitting it inflates `re`/`im` (and hence
    # `radn`/`idxs`) by a factor of `refp`, pushing nearly every far-field contribution
    # outside the valid `1:nr` index range where it gets masked out as out-of-domain —
    # which is exactly what made every far-off-diagonal DTN entry silently come out as
    # exact 0.0 (caught by comparing against the legacy MATLAB DTN cache with real Julia
    # execution, after this port had only ever been checked by static re-reading before).
    re = rn_k .+ i_far .* cos_l' ./ refp
    im = i_far .* sin_l' ./ refp
    radn = abs.(sqrt.(re .^ 2 .+ im .^ 2))
    idxs = floor.(radn)
    w1 = clamp.(radn .- idxs, 0.0, 1.0)
    w2 = w1 .^ 2
    w3 = w1 .^ 3

    cond1 = idxs .< 0.5
    accumulate_masked!(Line, @.(-(3 * w3 / 4 - 7 * w2 / 4 + 1) * Kern_far), cond1, idxs .+ 1, nr)
    accumulate_masked!(Line, @.(-(-w3 + 2 * w2) * Kern_far), cond1, idxs .+ 2, nr)
    accumulate_masked!(Line, @.(-(w3 - w2) / 4 * Kern_far), cond1, idxs .+ 3, nr)

    cond2 = (idxs .>= 0.5) .& (idxs .< nr)
    accumulate_masked!(Line, @.(-(-w3 / 6 + w2 / 2 - w1 / 3) * Kern_far), cond2, idxs, nr)
    accumulate_masked!(Line, @.(-(w3 / 2 - w2 - w1 / 2 + 1) * Kern_far), cond2, idxs .+ 1, nr)

    cond3 = (idxs .< nr - 1) .& (idxs .>= 0.5)
    cond4 = (idxs .< nr - 2) .& (idxs .>= 0.5)
    accumulate_masked!(Line, @.(-(-w3 / 2 + w2 / 2 + w1) * Kern_far), cond3, idxs .+ 2, nr)
    accumulate_masked!(Line, @.(-(w3 - w1) / 6 * Kern_far), cond4, idxs .+ 3, nr)
    return Line
end

"""
    generate_dtn(nr, D; refp=10) -> Matrix{Float64}

Computes the `nr x nr` dense Dirichlet-to-Neumann operator for a domain of diameter `D`
(in disk radii), matching `DTNVectorized.m`'s output (see module-level warning above
about verification status). Rows `4:nr` are computed in parallel via
`Threads.@threads :dynamic` — deliberately dynamic scheduling, not static, because the
per-row cost grows with `k` (the far-field quadrature range widens with row index), so a
naive equal-sized static split would load-balance poorly across threads.
"""
function generate_dtn(nr::Int, D::Real; refp::Int = DTN_DEFAULT_REFP)
    nr >= 4 || error("generate_dtn requires nr >= 4, got nr=$nr")
    dr = D / (2 * nr)
    drp = dr / refp

    numer = ceil(Int, pi * D / drp)
    isodd(numer) && (numer += 1)
    dtheta = 2 * pi / numer

    DTN = zeros(nr, nr)

    rn(k) = Float64(k - 1)  # rn = 0:nr+1 in MATLAB (1-indexed array); rn(k) == k-1.

    # ---------- Row 1 (r = 0): near-origin singular integral, closed form ----------
    for i in 2:nr
        DTN[1, i] -= _term0(Float64(i))
    end
    for i in 2:(nr - 1)
        DTN[1, i + 1] -= _term1(Float64(i))
    end
    for i in 2:(nr - 2)
        DTN[1, i + 2] -= _term2(Float64(i))
    end
    for i in 2:(nr - 3)
        DTN[1, i + 3] -= _term3(Float64(i))
    end
    DTN[1, 1] += 0.5
    DTN[1, :] ./= dr
    DTN[1, 1] += 209 / (54 * dr)
    DTN[1, 2] += -29 / (6 * dr)
    DTN[1, 3] += 7 / (6 * dr)
    DTN[1, 4] += -11 / (54 * dr)

    # ---------- Shared angular quadrature nodes ----------
    l_vals = collect((dtheta / 2):dtheta:(pi - dtheta / 4))
    cos_l = cos.(l_vals)
    sin_l = sin.(l_vals)

    i_near = collect(1.0:(2 * refp))
    Kern_near = _near_far_kernel(i_near)

    diag_correction = 2 / (4 * dr + drp)

    # ---------- Row 2 ----------
    let k = 2
        rn2 = rn(k)
        re = rn2 .+ i_near .* cos_l' ./ refp
        im = i_near .* sin_l' ./ refp
        radn = abs.(sqrt.(re .^ 2 .+ im .^ 2))
        x1 = i_near .* cos_l' ./ refp
        posr = radn .- rn2

        DTN[2, 2] -= sum(@.((-6 - 3 * posr + 9 * posr^2 + 3 * posr^3 - 3 * posr^4) * posr + 6 * x1) / 72 .* Kern_near)
        DTN[2, 1] -= sum(@.((-88 + 48 * posr + 62 * posr^2 - 12 * posr^3 - 10 * posr^4) * posr + 88 * x1) / 72 .* Kern_near)
        DTN[2, 2] -= sum(@.((72 - 90 * posr - 90 * posr^2 + 18 * posr^3 + 18 * posr^4) * posr - 72 * x1) / 72 .* Kern_near)
        if k < nr - 0.5
            DTN[2, 3] -= sum(@.((24 + 48 * posr + 18 * posr^2 - 12 * posr^3 - 6 * posr^4) * posr - 24 * x1) / 72 .* Kern_near)
        end
        if k < nr - 1.5
            DTN[2, 4] -= sum(@.((-2 - 3 * posr + posr^2 + 3 * posr^3 + posr^4) * posr + 2 * x1) / 72 .* Kern_near)
        end

        i_far = collect((2.0 * refp + 1):((rn2 + nr) * refp))
        Kern_far = _near_far_kernel(i_far)
        LineUpdate2 = zeros(nr)
        far_part!(LineUpdate2, rn2, i_far, Kern_far, cos_l, sin_l, nr, refp)

        DTN[2, :] .= (dtheta / (2 * pi * drp)) .* (DTN[2, :] .+ LineUpdate2)
        DTN[2, 2] += diag_correction
    end

    # ---------- Row 3 ----------
    let k = 3
        rn3 = rn(k)
        re = rn3 .+ i_near .* cos_l' ./ refp
        im = i_near .* sin_l' ./ refp
        radn = abs.(sqrt.(re .^ 2 .+ im .^ 2))
        x1 = i_near .* cos_l' ./ refp
        posr = radn .- rn3

        DTN[3, 1] -= sum(@.((124 - 12 * posr - 149 * posr^2 + 12 * posr^3 + 25 * posr^4) * posr - 124 * x1) / 288 .* Kern_near)
        DTN[3, 2] -= sum(@.((-384 + 192 * posr + 288 * posr^2 - 48 * posr^3 - 48 * posr^4) * posr + 384 * x1) / 288 .* Kern_near)
        DTN[3, 3] += sum(@.((-4 + 10 * posr + 5 * posr^2 - 2 * posr^3 - posr^4) * posr + 4 * x1) / 8 .* Kern_near)
        if k < nr - 0.5
            DTN[3, 4] -= sum(@.((128 + 192 * posr + 32 * posr^2 - 48 * posr^3 - 16 * posr^4) * posr - 128 * x1) / 288 .* Kern_near)
        end
        if k < nr - 1.5
            DTN[3, 5] -= sum(@.((-12 - 12 * posr + 9 * posr^2 + 12 * posr^3 + 3 * posr^4) * posr + 12 * x1) / 288 .* Kern_near)
        end

        i_far = collect((2.0 * refp + 1):((rn3 + nr) * refp))
        Kern_far = _near_far_kernel(i_far)
        LineUpdate3 = zeros(nr)
        far_part!(LineUpdate3, rn3, i_far, Kern_far, cos_l, sin_l, nr, refp)

        DTN[3, :] .= (dtheta / (2 * pi * drp)) .* (DTN[3, :] .+ LineUpdate3)
        DTN[3, 3] += diag_correction
    end

    # ---------- Rows 4:nr (parallel) ----------
    progress = Threads.Atomic{Int}(0)
    progress_every = max(1, cld(nr, 20))
    n_rows_total = nr - 3
    Threads.@threads :dynamic for k in 4:nr
        rnk = rn(k)
        re = rnk .+ i_near .* cos_l' ./ refp
        im = i_near .* sin_l' ./ refp
        radn = abs.(sqrt.(re .^ 2 .+ im .^ 2))
        x1 = i_near .* cos_l' ./ refp
        posr = radn .- rnk

        Line = zeros(nr)
        Line[k - 2] -= sum(@.((2 - posr - 2 * posr^2 + posr^3) * posr - 2 * x1) / 24 .* Kern_near)
        Line[k - 1] -= sum(@.(4 * (-4 + 4 * posr + posr^2 - posr^3) * posr + 16 * x1) / 24 .* Kern_near)
        Line[k] -= sum(@.((posr^2 - 5) * posr^2 / 4) .* Kern_near)
        if k < nr - 0.5
            Line[k + 1] -= sum(@.(4 * (4 + 4 * posr - posr^2 - posr^3) * posr - 16 * x1) / 24 .* Kern_near)
            if k < nr - 1.5
                Line[k + 2] -= sum(@.((-2 - posr + 2 * posr^2 + posr^3) * posr + 2 * x1) / 24 .* Kern_near)
            end
        end

        i_far = collect((2.0 * refp + 1):((rnk + nr) * refp))
        Kern_far = _near_far_kernel(i_far)
        far_part!(Line, rnk, i_far, Kern_far, cos_l, sin_l, nr, refp)

        Line .*= dtheta / (2 * pi * drp)
        @views DTN[k, :] .+= Line
        DTN[k, k] += diag_correction

        n = Threads.atomic_add!(progress, 1) + 1
        if n % progress_every == 0 || n == n_rows_total
            @info "DTN generation progress" rows_done = n rows_total = n_rows_total nr = nr
        end
    end

    return DTN
end
