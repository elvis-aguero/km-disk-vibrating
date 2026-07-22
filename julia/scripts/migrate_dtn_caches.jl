#!/usr/bin/env julia
#= One-off migration: convert the existing MATLAB-generated DTN `.mat` caches into
Julia-native JLD2 files and register them in `julia/data/dtn_cache/manifest.toml`.

Run from the repo root: `julia julia/scripts/migrate_dtn_caches.jl`

For each legacy cache, `(spatialResolution, bathDiameter)` is taken from the *folder
name* (`D<spatialResolution>Quant<bathDiameter>`), which is how `solve_motion.m` actually
used the domain at runtime — NOT from the value embedded in the `.mat` filename itself,
which is ambiguous/wrong for two of the three legacy folders (`D25Quant200` and
`D5Quant20` both encode `spatialResolution` in the filename where the normal
`DTNnew345nr<nr>D<bathDiameter>refp10.mat` convention expects `bathDiameter`).

For the small `D5Quant20` domain (nr=50, fast to regenerate), this script also generates
the DTN matrix natively and compares it against the migrated cache. Because that legacy
cache's filename is ambiguous, there is a real possibility the *wrong* `D` argument was
originally passed into the MATLAB generator itself (not just that the output file was
mislabeled) — `DTNVectorized(nr, D)` uses `D` to compute the mesh spacing `dr = D/(2nr)`,
so a wrong `D` would mean a numerically wrong matrix, not just a misnamed one. If the
native and migrated matrices disagree, the manifest is pointed at the freshly generated
(native) matrix instead of the migrated one, and the discrepancy is recorded in its
`source_file` provenance field so it's auditable later.

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
    bathDiameter::Float64
    nr::Int
    mat_path::String
    jld2_name::String
end

const LEGACY_SOURCES = [
    LegacySource(50, 100, 2500, joinpath(REPO_ROOT, "matlab", "1_code", "D50Quant100", "DTNnew345nr2500D100refp10.mat"), "D50_Quant100_nr2500.jld2"),
    LegacySource(25, 200, 2500, joinpath(REPO_ROOT, "matlab", "1_code", "D25Quant200", "DTNnew345nr2500D25refp10.mat"), "D25_Quant200_nr2500.jld2"),
    LegacySource(5, 20, 50, joinpath(REPO_ROOT, "matlab", "1_code", "D5Quant20", "DTNnew345nr50D5refp10.mat"), "D5_Quant20_nr50.jld2"),
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
            @info "cross-validating against native Julia generator (small domain, cheap)" nr = src.nr bathDiameter = src.bathDiameter
            DTN_native = generate_dtn(src.nr, src.bathDiameter)
            max_abs_diff = maximum(abs, DTN_native .- DTN_to_store)
            rel_diff = max_abs_diff / max(maximum(abs, DTN_to_store), eps())
            if rel_diff < 1e-6
                @info "native generator matches migrated cache" max_abs_diff = max_abs_diff rel_diff = rel_diff
            else
                @warn "native generator DISAGREES with migrated legacy cache — legacy .mat may have been generated with the wrong D argument, not just a mislabeled filename; registering the NATIVE matrix as canonical" max_abs_diff = max_abs_diff rel_diff = rel_diff
                DTN_to_store = DTN_native
                source_note = "generated natively (disagreed with legacy $(src.mat_path): rel_diff=$(rel_diff))"
            end
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
