#!/usr/bin/env julia
#= Regenerate one DTN cache with write -> verify -> rename, then register it.
Rows are computed by SEPARATE WORKER PROCESSES (Distributed), not by threads.

Usage: `DTN_WORKERS=24 julia julia/cluster/beads/regenerate_dtn.jl <spatialResolution> <bathDiameter>`
(run julia with its default 1 thread: the parallelism is the worker processes.)
Optional environment:
  DTN_WORKERS    number of worker processes (default: 4)
  DTN_CACHE_DIR  write here instead of julia/data/dtn_cache. TEST MODE: the manifest is not
                 touched and, if the default registry already has a cache for the same
                 (spatialResolution, bathDiameter), the new matrix must match it.

Why processes, not threads (measured, job 8262774): `generate_dtn` with `-t 32` ran at
exactly 1.00 busy core while its memory grew. Each row allocates ~10 large temporaries (the
"near" part) and then runs a long allocation-free loop (the "far" part); Julia's GC needs
every thread at a safepoint, so one thread's allocation waits for the others' long loops.
Separate processes each have their own GC, so they do not wait for each other. Rows are
independent (`dtn_row`), and the row cost grows with the row index, so rows are handed out
from the last to the first.

Why not `scripts/generate_dtn.jl` as is: it writes the final file name directly and registers
it at once. A job killed during the write would leave a bad file at the final name -- the
failure that produced the old, incomplete D25Quant200 cache (rows 564:2500 all zero). This
wrapper:
  1. refuses to overwrite an existing cache;
  2. checks the matrix in memory (`check_dtn`: finite, no all-zero row);
  3. writes `<name>.jld2.partial`, reloads it, and requires an exact round trip;
  4. recomputes a few rows (including 563/564, where the old cache broke) and compares;
  5. only then renames the file and adds the manifest entry.
Any failure stops the script with a non-zero exit code and leaves no cache and no manifest
entry. =#

using Distributed

include(joinpath(@__DIR__, "..", "..", "scripts", "_bootstrap.jl"))   # activates env, `using FaradayDisk` on the master
using JLD2

const SCRIPTS_PROJECT = abspath(joinpath(@__DIR__, "..", "..", "scripts"))

# Workers are started at TOP LEVEL, not inside a function: `@everywhere using ...` expands to
# a `:toplevel` expression, which Julia rejects inside a function body (this was the error
# of the first self-test, job 8264010). The job function is defined with `@everywhere` AFTER
# the workers exist: workers added later do not receive earlier `@everywhere` definitions.
const NWORKERS = parse(Int, get(ENV, "DTN_WORKERS", "4"))
const T_WORKERS = time()
addprocs(NWORKERS; exeflags = ["--project=$SCRIPTS_PROJECT", "-t1", "--startup-file=no"])
@everywhere using FaradayDisk
@everywhere begin
    function _dtn_row_job(nr::Int, bd::Float64, k::Int, progress)
        row = FaradayDisk.dtn_row(nr, bd, k)
        put!(progress, k)
        return row
    end
end
@info "workers running" n = nworkers() start_up_s = round(time() - T_WORKERS; digits = 1)
flush(stderr)

function compute_rows(nr::Int, bd::Float64, rows::AbstractVector{Int}; log_progress::Bool = true)
    progress = RemoteChannel(() -> Channel{Int}(length(rows) + 16))
    logger = @async begin
        done = 0
        every = max(1, cld(length(rows), 20))
        t0 = time()
        while done < length(rows)
            take!(progress)
            done += 1
            if log_progress && (done % every == 0 || done == length(rows))
                @info "DTN generation progress" rows_done = done rows_total = length(rows) elapsed_s = round(time() - t0; digits = 1)
                flush(stderr)
            end
        end
    end
    results = pmap(k -> _dtn_row_job(nr, bd, k, progress), rows; retry_delays = zeros(2))
    wait(logger)
    return results
end

function main(args)
    length(args) == 2 || error("usage: regenerate_dtn.jl <spatialResolution> <bathDiameter>")
    sr = parse(Float64, args[1])
    bd = parse(Float64, args[2])
    nr = ceil(Int, sr * bd / 2)
    nworkers_req = NWORKERS

    test_mode = haskey(ENV, "DTN_CACHE_DIR")
    cache_dir = abspath(get(ENV, "DTN_CACHE_DIR", joinpath(@__DIR__, "..", "..", "data", "dtn_cache")))
    mkpath(cache_dir)
    name = "D$(Int(sr))_Quant$(Int(bd))_nr$(nr).jld2"
    final = joinpath(cache_dir, name)
    partial = final * ".partial"
    isfile(final) && error("refusing to overwrite existing cache $final")
    rm(partial; force = true)

    @info "generating DTN" spatialResolution = sr bathDiameter = bd nr = nr workers = nworkers_req test_mode = test_mode
    t0 = time()
    flush(stderr)

    order = collect(nr:-1:1)                      # heaviest rows first
    row_data = compute_rows(nr, bd, order)
    DTN = zeros(nr, nr)
    for (i, k) in enumerate(order)
        DTN[k, :] .= row_data[i]
    end
    @info "generation done" elapsed_s = round(time() - t0; digits = 1)

    size(DTN) == (nr, nr) || error("unexpected DTN size $(size(DTN))")
    FaradayDisk.check_dtn(DTN; source = "in-memory matrix")

    JLD2.jldsave(partial; DTN = DTN)
    back = JLD2.load(partial, "DTN")
    back == DTN || error("round trip through $partial changed the matrix")
    FaradayDisk.check_dtn(back; source = partial)

    check_rows = unique(filter(k -> 1 <= k <= nr, [4, 563, 564, 1500, nr]))
    again = compute_rows(nr, bd, check_rows; log_progress = false)
    rel = [maximum(abs, back[k, :] .- again[i]) / maximum(abs, again[i]) for (i, k) in enumerate(check_rows)]
    @info "independent row check" rows = check_rows max_rel_diff = rel
    all(<(1e-12), rel) || error("independent row check failed: $rel")

    if test_mode
        # Reference: the library's own serial generator (checked against the corrected MATLAB
        # reference by test/integration/test_dtn_golden_small.jl). NOT the registered cache of
        # the same domain: the local D5_Quant20 cache is stale (migrated from the old D=5
        # matrix before commit 5391c98), which made the first test fail with a 0.75 difference.
        nr <= 200 || error("test mode compares with the serial generator and is meant for small domains (nr <= 200), got nr=$nr")
        ref = FaradayDisk.generate_dtn_polar(nr, bd)
        relerr = maximum(abs, back .- ref) / maximum(abs, ref)
        @info "test mode: comparison with the serial generator" relerr = relerr
        flush(stderr)
        relerr < 1e-12 || error("test mode: parallel matrix differs from the serial generator (relative difference $relerr)")
        mv(partial, final)
        @info "TEST MODE OK; manifest not touched" path = final total_elapsed_s = round(time() - t0; digits = 1)
        return nothing
    end

    mv(partial, final)

    sha = try
        readchomp(`git -C $(@__DIR__) rev-parse --short HEAD`)
    catch
        "unknown"
    end
    manifest_path = joinpath(cache_dir, "manifest.toml")
    manifest = isfile(manifest_path) ? load_manifest(manifest_path) : DTNManifest()
    register!(manifest, DTNCacheEntry(sr, bd, nr, 10, final,
        "generated natively by julia/cluster/beads/regenerate_dtn.jl (git $sha, slurm job $(get(ENV, "SLURM_JOB_ID", "none")))"))
    save_manifest(manifest_path, manifest)

    @info "DTN cache verified and registered" path = final manifest = manifest_path total_elapsed_s = round(time() - t0; digits = 1)
end

main(ARGS)
