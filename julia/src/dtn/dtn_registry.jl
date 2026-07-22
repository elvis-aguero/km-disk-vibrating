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
end

struct DTNManifest
    entries::Dict{Tuple{Float64,Float64},DTNCacheEntry}
end

DTNManifest() = DTNManifest(Dict{Tuple{Float64,Float64},DTNCacheEntry}())

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

CAVEAT: `load_manifest`/`save_manifest` use `TOML.parsefile`/`TOML.print`, written from
memory without the ability to confirm `TOML.print` specifically is available/stable in
the resolved Julia 1.10 TOML stdlib (see the repo-level notes on how this port was
authored). Not exercised by the default (fast) test tier, which passes DTN matrices
directly and never touches the manifest file.
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
        )
        manifest.entries[(entry.spatialResolution, entry.bathDiameter)] = entry
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
    resolve_dtn(registry, spatialResolution, bathDiameter) -> DTNCacheEntry

Look up the cache entry for a domain. Throws a clear, actionable error (naming the
script to run) rather than silently guessing a filename, if no entry is registered.
"""
function resolve_dtn(registry::DTNManifest, spatialResolution::Real, bathDiameter::Real)
    key = (float(spatialResolution), float(bathDiameter))
    entry = get(registry.entries, key, nothing)
    if entry === nothing
        error(
            "No DTN cache registered for spatialResolution=$spatialResolution, " *
            "bathDiameter=$bathDiameter. Run `scripts/generate_dtn.jl` for this domain " *
            "(or `scripts/migrate_dtn_caches.jl` if a MATLAB .mat cache already exists " *
            "for it), then add an entry to $(default_manifest_path()).",
        )
    end
    return entry
end

"""
    load_dtn(entry) -> Matrix{Float64}

Loads the DTN matrix referenced by a registry entry from its JLD2 file.
"""
load_dtn(entry::DTNCacheEntry) = JLD2.load(entry.path, "DTN")

"""
    register!(manifest, entry)

Insert or replace a manifest entry, keyed by `(spatialResolution, bathDiameter)`.
"""
function register!(manifest::DTNManifest, entry::DTNCacheEntry)
    manifest.entries[(entry.spatialResolution, entry.bathDiameter)] = entry
    return manifest
end
