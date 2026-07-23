#!/usr/bin/env julia
#= One-off migration: convert the existing MATLAB-generated DTN `.mat` caches into
Julia-native JLD2 files and register them in `julia/data/dtn_cache/manifest.toml`.

Run from the repo root: `julia julia/scripts/migrate_dtn_caches.jl`

Each legacy cache is registered under `(spatialResolution, bathDiameter)` taken from its
*folder name* (`D<spatialResolution>Quant<bathDiameter>`) — this is the key `solve_motion.m`
actually resolves DTN caches by, and what the rest of the physics pipeline (`domainMaker`,
called as `domainMaker(bathDiameter, spatialResolution)`) uses for this domain regardless
of what DTN matrix ends up loaded.

That's a distinct question from what `D` value each `.mat` file's *contents* were actually
generated with — `DTNVectorized.m`'s last line, `save(['DTNnew345nr',num2str(nr),'D',
num2str(D),...])`, saves the literal `D` argument it was called with, so a filename is a
mechanically exact record of that, not a guess. For `D50Quant100`
(`DTNnew345nr2500D100refp10.mat`) that literal value (100) matches the folder's implied
bathDiameter, so the cache should be exactly what the domain needs.

`D5Quant20` (originally `DTNnew345nr50D5refp10.mat`, literal D=5) did NOT match its folder's
implied bathDiameter=20 — confirmed by regenerating natively at D=5 and matching the legacy
cache to ~2e-15 relative difference (machine precision), while D=20 disagreed by 75%.
`solve_motion.m` had a "Machine-specific patch for D5Quant20" that loaded that mismatched
file by its literal name regardless, meaning the DTN operator and the bulk-domain operators
(Laplacian, pressure integral — both built from the nominal bathDiameter=20) were computed
on physically different domain sizes (domain radius = bathDiameter/2, so D=5 vs D=20 is a
4x difference, not a rounding-level one) for the same simulation. A real usage audit
confirmed D5Quant20 is never used by any production/experimental sweep in either MATLAB or
Julia (both use bathDiameter=100/spatialResolution=50) — it exists solely as a cheap
test/parity fixture — so this was fixed directly: `scripts/regenerate_d5quant20_matlab_cache.jl`
regenerated the cache natively at the correct D=20 (`DTNnew345nr50D20refp10.mat`, matching
the standard naming convention) and `solve_motion.m`'s special case was removed. The
`validation_D` field below is now just `bathDiameter` for this entry, same as every other
self-consistent cache.

`D25Quant200` (`DTNnew345nr2500D25refp10.mat`, literal D=25) does NOT match its folder's
implied bathDiameter=200 either, and is NOT independently verified the way `D5Quant20` was
— its `nr=2500` is too expensive to regenerate natively in a reasonable time (the per-row
quadrature cost scales with the far-field range, which itself scales with `nr`), so whether
it has the same kind of mismatch is currently unknown — flagged, not assumed either way.
`D50Quant100` is self-consistent by filename but likewise has never actually been
cross-validated against the native generator (same cost problem).

If the *validation* comparison for a small (`nr<=200`) domain disagrees beyond a tight
(1e-8) tolerance, that is now treated as a genuine, unexplained bug to surface loudly —
not silently patched over by substituting the native matrix as "probably fine" (that
overly-charitable fallback is exactly what masked the real `far_part!` `/refp` bug this
port shipped with for a while; see `src/dtn/dtn_generator.jl`'s history).

Block comment, not a docstring: see _bootstrap.jl's header for why a bare top-level string
directly followed by `include(...)` here would fail to parse. =#

include("_bootstrap.jl")
using MAT
using JLD2

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const CACHE_DIR = joinpath(@__DIR__, "..", "data", "dtn_cache")
const MANIFEST_PATH = joinpath(CACHE_DIR, "manifest.toml")

struct LegacySource
    spatialResolution::Float64
    bathDiameter::Float64       # registry key component; what solve_motion.m/domainMaker use
    nr::Int
    mat_path::String
    jld2_name::String
    validation_D::Float64       # the D this cache's contents were *actually* generated with
end

const LEGACY_SOURCES = [
    LegacySource(50, 100, 2500, joinpath(REPO_ROOT, "matlab", "1_code", "D50Quant100", "DTNnew345nr2500D100refp10.mat"), "D50_Quant100_nr2500.jld2", 100),
    LegacySource(25, 200, 2500, joinpath(REPO_ROOT, "matlab", "1_code", "D25Quant200", "DTNnew345nr2500D25refp10.mat"), "D25_Quant200_nr2500.jld2", 25),
    LegacySource(5, 20, 50, joinpath(REPO_ROOT, "matlab", "1_code", "D5Quant20", "DTNnew345nr50D20refp10.mat"), "D5_Quant20_nr50.jld2", 20),
]

function main()
    mkpath(CACHE_DIR)
    manifest = isfile(MANIFEST_PATH) ? load_manifest(MANIFEST_PATH) : DTNManifest()

    for src in LEGACY_SOURCES
        if !isfile(src.mat_path)
            @warn "legacy DTN cache not found, skipping" path = src.mat_path
            continue
        end

        @info "migrating DTN cache" spatialResolution = src.spatialResolution bathDiameter = src.bathDiameter nr = src.nr from = src.mat_path
        DTN = MAT.matread(src.mat_path)["DTNnew345"]
        size(DTN) == (src.nr, src.nr) || error("expected $(src.nr)x$(src.nr), got $(size(DTN)) in $(src.mat_path)")

        out_path = joinpath(CACHE_DIR, src.jld2_name)
        source_note = src.mat_path
        DTN_to_store = Matrix{Float64}(DTN)

        if src.nr <= 200
            @info "cross-validating against native Julia generator" nr = src.nr validation_D = src.validation_D
            DTN_native = generate_dtn(src.nr, src.validation_D)
            max_abs_diff = maximum(abs, DTN_native .- DTN_to_store)
            rel_diff = max_abs_diff / max(maximum(abs, DTN_to_store), eps())
            if rel_diff < 1e-8
                @info "native generator matches migrated cache" max_abs_diff = max_abs_diff rel_diff = rel_diff
            else
                @error "native generator DISAGREES with migrated legacy cache at its own validation_D — this is now treated as a real bug, not a presumed mislabeling" max_abs_diff = max_abs_diff rel_diff = rel_diff validation_D = src.validation_D
                error("DTN cache validation failed for $(src.mat_path) (rel_diff=$rel_diff at D=$(src.validation_D)) — investigate before trusting this cache")
            end
        else
            @warn "large domain (nr>200) — DTN cache NOT independently verified against the native generator (too expensive to regenerate here); registering the legacy matrix as-is" spatialResolution = src.spatialResolution bathDiameter = src.bathDiameter nr = src.nr
        end

        JLD2.jldsave(out_path; DTN = DTN_to_store)
        entry = DTNCacheEntry(src.spatialResolution, src.bathDiameter, src.nr, 10, out_path, source_note)
        register!(manifest, entry)
        @info "registered DTN cache entry" spatialResolution = src.spatialResolution bathDiameter = src.bathDiameter path = out_path
    end

    save_manifest(MANIFEST_PATH, manifest)
    @info "migration complete" manifest = MANIFEST_PATH n_entries = length(manifest.entries)
end

main()
