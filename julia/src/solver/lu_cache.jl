"""
LU-factorization cache for the `:lu` solver path.

Port of the caching behaviour actually shipped in `advance_one_step.m` (the *lazy-sync*
path — a comment block in `solve_motion.m` describes an async `parfeval`-based
pre-computation of every cycle's factorization ahead of time, but no code in that file
actually fires those futures, so lazy-on-first-use is the real current behavior, not a
simplification of it). Since `dt` is chosen so the forcing is exactly periodic over
`stepsPerCycle` steps, the system matrix only ever takes `stepsPerCycle` distinct values
over the whole run — factorize each one once, on first encounter, and reuse it forever
after.

Julia's `LinearAlgebra.lu` returns a single factors object (`LU`) rather than MATLAB's
separate dense `L`, `U`, `P` matrices, so the same cache is roughly 3x more memory
efficient per cached factorization than MATLAB's `[L,U,P]` cell-array storage — this
matters for `solverType=:auto`'s RAM-budget decision, see `solver_dispatch.jl`.

`hits`/`misses` are cheap counters left in production code (not test-only
instrumentation) specifically so a test can assert the cache is actually being reused
(`hits / total_steps` should be close to 1 after the first cycle) — silently falling back
to refactorizing every step would be a real, easy-to-miss performance regression with no
correctness symptom, so it needs its own invariant check rather than relying on
correctness tests to catch it.
"""
mutable struct LUCache
    factors::Vector{Union{Nothing,LU}}
    hits::Int
    misses::Int
end

LUCache(stepsPerCycle::Integer) = LUCache(Vector{Union{Nothing,LU}}(nothing, stepsPerCycle), 0, 0)

Base.length(cache::LUCache) = length(cache.factors)

"""
    get_or_factorize!(cache, tmpl, Mat, g_prefactor, cycle_index) -> LU

Returns the cached factorization for `cycle_index`, computing (and caching) it via
`materialize!(Mat, tmpl, g_prefactor)` followed by `lu(Mat)` on first use. `Mat` is reused
scratch space — it is only actually read on a cache miss.
"""
function get_or_factorize!(cache::LUCache, tmpl::SystemMatrixTemplate, Mat::Matrix{Float64},
                            g_prefactor::Float64, cycle_index::Int)
    1 <= cycle_index <= length(cache) || throw(BoundsError(cache.factors, cycle_index))
    f = cache.factors[cycle_index]
    if f === nothing
        materialize!(Mat, tmpl, g_prefactor)
        f = lu(Mat)
        cache.factors[cycle_index] = f
        cache.misses += 1
    else
        cache.hits += 1
    end
    return f
end

"""
    hit_rate(cache) -> Float64

Fraction of `get_or_factorize!` calls that were cache hits. `NaN` if never called.
"""
function hit_rate(cache::LUCache)
    total = cache.hits + cache.misses
    return total == 0 ? NaN : cache.hits / total
end
