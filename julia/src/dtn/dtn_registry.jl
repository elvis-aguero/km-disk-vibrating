"""
Explicit registry for precomputed Dirichlet-to-Neumann (DTN) operator caches.

Replaces MATLAB's filename-formatting lookup (`sprintf('DTNnew345nr%dD%drefp10.mat', nr,
bathDiameter)`), which is fragile — two of the three existing legacy caches
(`D25Quant200`, `D5Quant20`) have filenames that encode `spatialResolution` rather than
`bathDiameter`, silently breaking the naive formula and requiring a hardcoded special
case in `solve_motion.m`. Here, resolving a domain to a cache file is a plain dictionary
lookup keyed by the physically meaningful `(spatialResolution, bathDiameter)` pair,
loaded from a manifest file — adding a new domain means adding a manifest entry, not
hoping a string-formatting convention holds.
"""

struct DTNCacheEntry
    spatialResolution::Float64
    bathDiameter::Float64
    nr::Int
    refp::Int
    path::String
    source_file::String
    method::Symbol      # :polar (the paper's scheme, generate_dtn_polar) or :single_layer (docs/dtn_single_layer.md)
end

# Entries written before the `method` field existed, and callers that build polar caches, omit it.
DTNCacheEntry(sr, bd, nr, refp, path, source_file) = DTNCacheEntry(sr, bd, nr, refp, path, source_file, :polar)

# One entry per (spatialResolution, bathDiameter, method): both constructions of a domain can be registered side by side.
struct DTNManifest
    entries::Dict{Tuple{Float64,Float64,Symbol},DTNCacheEntry}
end

DTNManifest() = DTNManifest(Dict{Tuple{Float64,Float64,Symbol},DTNCacheEntry}())

"""
    default_dtn_method() -> Symbol

The construction a domain is resolved with: `:single_layer`, unless the environment variable `KMDISK_DTN_METHOD`
is set to `polar` (the paper's scheme, the matrices used before the single-layer construction became the default).
"""
function default_dtn_method()
    m = Symbol(get(ENV, "KMDISK_DTN_METHOD", "single_layer"))
    m in (:single_layer, :polar) || error("KMDISK_DTN_METHOD must be single_layer or polar, got $m")
    return m
end

"""
    load_manifest(path) -> DTNManifest

Reads a `manifest.toml` of the form:

    [[entry]]
    spatialResolution = 50
    bathDiameter = 100
    nr = 2500
    refp = 10
    path = "D50_Quant100_nr2500.jld2"
    source_file = "matlab/1_code/D50Quant100/DTNnew345nr2500D100refp10.mat"

`path` is resolved relative to the manifest file's own directory.

Exercised by the cluster workflow (`migrate_dtn_caches.jl` writes the manifest, the sweep
reads it); not by the fast CI test tier, which passes DTN matrices directly.
"""
function load_manifest(path::AbstractString)
    manifest = DTNManifest()
    isfile(path) || return manifest
    doc = TOML.parsefile(path)
    dir = dirname(abspath(path))
    for e in get(doc, "entry", [])
        entry = DTNCacheEntry(
            float(e["spatialResolution"]),
            float(e["bathDiameter"]),
            Int(e["nr"]),
            Int(get(e, "refp", 10)),
            joinpath(dir, e["path"]),
            get(e, "source_file", ""),
            Symbol(get(e, "method", "polar")),
        )
        manifest.entries[(entry.spatialResolution, entry.bathDiameter, entry.method)] = entry
    end
    return manifest
end

function save_manifest(path::AbstractString, manifest::DTNManifest)
    dir = dirname(abspath(path))
    doc = Dict(
        "entry" => [
            Dict(
                "spatialResolution" => e.spatialResolution,
                "bathDiameter" => e.bathDiameter,
                "nr" => e.nr,
                "refp" => e.refp,
                "path" => relpath(abspath(e.path), dir),
                "source_file" => e.source_file,
                "method" => String(e.method),
            ) for e in values(manifest.entries)
        ],
    )
    open(path, "w") do io
        TOML.print(io, doc)
    end
    return nothing
end

"""
    default_manifest_path() -> String

Path to the manifest shipped with the package, `julia/data/dtn_cache/manifest.toml`,
resolved relative to this source file so it works regardless of the caller's working
directory.
"""
default_manifest_path() = joinpath(@__DIR__, "..", "..", "data", "dtn_cache", "manifest.toml")

default_registry() = load_manifest(default_manifest_path())

"""
    resolve_dtn(registry, spatialResolution, bathDiameter; method=default_dtn_method()) -> DTNCacheEntry

Look up the cache entry for a domain and construction method. Throws a clear, actionable error (naming the
script to run) rather than silently guessing a filename or switching method, if no entry is registered.
"""
function resolve_dtn(registry::DTNManifest, spatialResolution::Real, bathDiameter::Real; method::Symbol = default_dtn_method())
    key = (float(spatialResolution), float(bathDiameter), method)
    entry = get(registry.entries, key, nothing)
    if entry === nothing
        others = sort!([String(k[3]) for k in keys(registry.entries) if k[1] == key[1] && k[2] == key[2]])
        error(
            "No DTN cache (method $method) registered for spatialResolution=$spatialResolution, " *
            "bathDiameter=$bathDiameter" * (isempty(others) ? "" : " (registered for this domain: $(join(others, ", ")); " *
            "set KMDISK_DTN_METHOD to use one of those)") *
            ". Run `julia scripts/generate_dtn.jl $spatialResolution $bathDiameter` " *
            "(or `scripts/migrate_dtn_caches.jl` if a MATLAB .mat cache already exists " *
            "for it), which adds an entry to $(default_manifest_path()).",
        )
    end
    return entry
end

"""
    load_dtn(entry) -> Matrix{Float64}

Loads the DTN matrix referenced by a registry entry from its JLD2 file, and rejects it if
it is visibly incomplete (`check_dtn`).
"""
load_dtn(entry::DTNCacheEntry) = check_dtn(JLD2.load(entry.path, "DTN"); source = entry.path)

"""
    check_dtn(DTN; source="") -> DTN

Throws if the matrix is not finite or has an all-zero row. Every row of a valid DTN
operator has a nonzero diagonal (`diag_correction` in `dtn_generator.jl`), so an all-zero
row means the cache was only partly generated (as a former `D25Quant200` cache was:
rows 564:2500 were zero).
"""
function check_dtn(DTN::AbstractMatrix; source::AbstractString = "")
    all(isfinite, DTN) || error("DTN matrix $source contains non-finite values")
    for k in axes(DTN, 1)
        if all(iszero, view(DTN, k, :))
            nz = count(j -> all(iszero, view(DTN, j, :)), axes(DTN, 1))
            error("DTN matrix $source is incomplete: row $k (and $nz rows in total) is all zeros. " *
                  "Regenerate it (scripts/generate_dtn.jl) rather than using it.")
        end
    end
    return DTN
end

"""
    register!(manifest, entry)

Insert or replace a manifest entry, keyed by `(spatialResolution, bathDiameter, method)`.
"""
function register!(manifest::DTNManifest, entry::DTNCacheEntry)
    manifest.entries[(entry.spatialResolution, entry.bathDiameter, entry.method)] = entry
    return manifest
end
