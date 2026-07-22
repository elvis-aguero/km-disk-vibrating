#!/usr/bin/env julia
"""
Generate a DTN cache natively in Julia for a given domain, and register it.

Usage: `julia -t auto julia/scripts/generate_dtn.jl <spatialResolution> <bathDiameter>`

Equivalent to MATLAB's `generate_dtn.m` (which calls `DTNVectorized.m`), but threaded via
`Threads.@threads :dynamic` instead of `parfor` — run with `-t auto` (or `-t N`) to get
any parallelism at all; this script is single-threaded otherwise.

WARNING: the largest domains in this repo (nr=2500) took MATLAB's `parfor`-based
generator on the order of hours; some legacy variants (`ParRadDTNStops.m`) even had
ad hoc checkpoint/restart logic specifically because a single run might not finish in one
sitting. Expect this to be slow for large `nr` regardless of the language.
"""

include("_bootstrap.jl")
using JLD2

function main(args)
    length(args) == 2 || error("usage: generate_dtn.jl <spatialResolution> <bathDiameter>")
    spatialResolution = parse(Float64, args[1])
    bathDiameter = parse(Float64, args[2])
    nr = ceil(Int, spatialResolution * bathDiameter / 2)

    @info "generating DTN" spatialResolution = spatialResolution bathDiameter = bathDiameter nr = nr threads = Threads.nthreads()
    t0 = time()
    DTN = generate_dtn(nr, bathDiameter)
    elapsed = time() - t0
    @info "DTN generation complete" nr = nr elapsed_s = elapsed

    cache_dir = joinpath(@__DIR__, "..", "data", "dtn_cache")
    mkpath(cache_dir)
    jld2_name = "D$(Int(spatialResolution))_Quant$(Int(bathDiameter))_nr$(nr).jld2"
    out_path = joinpath(cache_dir, jld2_name)
    JLD2.jldsave(out_path; DTN = DTN)

    manifest_path = joinpath(cache_dir, "manifest.toml")
    manifest = isfile(manifest_path) ? load_manifest(manifest_path) : DTNManifest()
    register!(manifest, DTNCacheEntry(spatialResolution, bathDiameter, nr, 10, out_path,
                                        "generated natively via scripts/generate_dtn.jl"))
    save_manifest(manifest_path, manifest)

    @info "registered DTN cache" path = out_path manifest = manifest_path
end

main(ARGS)
