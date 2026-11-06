#!/usr/bin/env julia
#= Bead campaign driver: for each requested bead, run a (gamma x frequency) bath-driven sweep
with the library's `run_sweep`, then post-process every case into gain and phase lag in the
lab and shaker frames (`gain_phase.jl`).

Usage (key=value arguments, order free):
  julia -t N julia/cluster/beads/run_bead_sweep.jl \
      outdir=<dir> gammas=0.1,0.5 [beads=ptfe,al32,al16,nylon] [basis=R|rc] [scheme=bdf2|euler] \
      [freqs=grid|<comma list in Hz>] [domain=25,200] [tres=60] [stage=<label>] [require_finite=true]

  gammas     required: there is deliberately no default.
  basis      R  = diskRadius is the actual bead radius; rc = the contact-line radius.
  domain     "<spatialResolution>,<bathDiameter>": resolves the DTN cache in the manifest.
  tres       steps per bath period (reference runs used 60).
  require_finite=true  makes the script fail if any case is error/missing/nonfinite (smoke tests).

Failure policy: a failing case never stops the campaign. `run_sweep_case` catches errors per
case; a failing bead never stops the other beads; post-processing catches errors per case.
Only a missing or corrupt DTN cache stops the job, and it stops it before any compute.

Output under `outdir`:
  run_config.toml          every parameter, git commit, host, DTN file
  <bead>/cases/*.csv       raw traces (time_s, CoM_cm, eta_boundary_cm)
  <bead>/summary.csv       library summary (status ok/error/skipped)
  <bead>/gain_phase.csv    gain/phase per case with quality flags
  gain_phase_all.csv       all beads together
Re-running with the same `outdir` resumes: cases whose CSV exists are skipped. =#

if !@isdefined(FaradayDisk)
    include(joinpath(@__DIR__, "..", "..", "scripts", "_bootstrap.jl"))
end
using Dates
using Logging
using TOML

include(joinpath(@__DIR__, "bead_params.jl"))
include(joinpath(@__DIR__, "gain_phase.jl"))

function parse_kv(args)
    kv = Dict{String,String}()
    for a in args
        i = findfirst('=', a)
        i === nothing && error("argument '$a' is not key=value")
        kv[a[1:(i - 1)]] = a[(i + 1):end]
    end
    return kv
end

parse_floats(s::AbstractString) = [parse(Float64, strip(x)) for x in split(s, ',') if !isempty(strip(x))]

function bead_params(key::AbstractString, basis::AbstractString, sr::Float64, bd::Float64, scheme::Symbol, tres::Int)
    b = BEADS[key]
    a = basis == "R" ? b.R : b.rc
    return SimulationParams(
        diskRadius = a, diskMass = bead_mass_g(b),
        forceAmplitude = 0.0, forceFrequency = 90.0,      # forceFrequency only sets the time unit (see sweep.jl)
        phaseDifference = -90.0,
        bathDensity = BEAD_FLUID.density, bathSurfaceTension = BEAD_FLUID.surface_tension,
        bathViscosity = BEAD_FLUID.viscosity, g = BEAD_FLUID.g,
        bathDiameter = bd, spatialResolution = sr, temporalResolution = tres,
        simulationTime = 30 / 90,                          # replaced per case by 30/f in run_sweep_case
        solverType = :auto, startStatic = true, earlyStop = true,
        timeScheme = scheme,
    )
end

git_info() = try
    sha = readchomp(`git -C $(@__DIR__) rev-parse --short HEAD`)
    dirty = !isempty(readchomp(`git -C $(@__DIR__) status --porcelain`))
    (sha = sha, dirty = dirty)
catch
    (sha = "unknown", dirty = false)
end

function main(args)
    kv = parse_kv(args)
    haskey(kv, "outdir") || error("outdir=<dir> is required")
    haskey(kv, "gammas") || error("gammas=... is required (no default: choose it deliberately)")
    outroot = kv["outdir"]
    beads = String.(split(get(kv, "beads", "ptfe,al32,al16,nylon"), ','))
    for b in beads
        haskey(BEADS, b) || error("unknown bead '$b'; known: $(sort(collect(keys(BEADS))))")
    end
    basis = get(kv, "basis", "R")
    basis in ("R", "rc") || error("basis must be R or rc")
    scheme = Symbol(get(kv, "scheme", "bdf2"))
    scheme in (:bdf2, :euler) || error("scheme must be bdf2 or euler")
    gammas = parse_floats(kv["gammas"])
    freqs = (!haskey(kv, "freqs") || kv["freqs"] == "grid") ? DEFAULT_FREQS : parse_floats(kv["freqs"])
    sr, bd = parse_floats(get(kv, "domain", "25,200"))
    tres = parse(Int, get(kv, "tres", "60"))
    stage = get(kv, "stage", "unnamed")
    require_finite = get(kv, "require_finite", "false") == "true"
    domain = "$(Int(sr))x$(Int(bd))"

    mkpath(outroot)
    global_logger(setup_logging(name = "bead_sweep", jsonlines = true))

    if any(>(1.0), gammas)
        @warn "gamma > 1 requested: effective gravity changes sign during the cycle; the model (pinned contact line, flat interface) may not be valid. Rows are flagged gamma_gt_1." gammas = gammas
    end

    # Fail early, before any compute, if the DTN cache is missing or incomplete.
    # KMDISK_DTN_FILE=<file.jld2 with key "DTN"> replaces the registered cache (TEST ONLY, e.g. an alternative DtN builder).
    dtn_override = get(ENV, "KMDISK_DTN_FILE", "")
    entry = isempty(dtn_override) ? resolve_dtn(default_registry(), sr, bd) :
            DTNCacheEntry(sr, bd, ceil(Int, sr * bd / 2), 10, abspath(dtn_override), "override via KMDISK_DTN_FILE")
    dtn = load_dtn(entry)      # load_dtn runs check_dtn
    @info "DTN cache loaded" path = entry.path nr = size(dtn, 1)

    git = git_info()
    config = Dict{String,Any}(
        "stage" => stage, "basis" => basis, "scheme" => String(scheme), "temporalResolution" => tres,
        "domain" => domain, "spatialResolution" => sr, "bathDiameter" => bd, "dtn_file" => entry.path,
        "gammas" => gammas, "freqs_Hz" => freqs, "beads" => beads,
        "fluid" => Dict{String,Any}("density_g_cm3" => BEAD_FLUID.density, "surface_tension_dyne_cm" => BEAD_FLUID.surface_tension,
                                    "viscosity_St" => BEAD_FLUID.viscosity, "g_cm_s2" => BEAD_FLUID.g),
        "git_commit" => git.sha, "git_dirty" => git.dirty, "julia_version" => string(VERSION),
        "threads" => Threads.nthreads(), "host" => gethostname(),
        "slurm_job_id" => get(ENV, "SLURM_JOB_ID", "none"), "started" => string(Dates.now()),
        "bead" => Dict{String,Any}(k => Dict{String,Any}(
            "name" => BEADS[k].name, "R_cm" => BEADS[k].R, "rc_cm" => BEADS[k].rc, "rho_s_g_cm3" => BEADS[k].rho_s,
            "mass_g" => bead_mass_g(BEADS[k]), "diskRadius_cm" => (basis == "R" ? BEADS[k].R : BEADS[k].rc)) for k in beads),
    )
    open(joinpath(outroot, "run_config.toml"), "w") do io
        TOML.print(io, config)
    end

    all_rows = Dict{String,Any}[]
    for key in beads
        params = bead_params(key, basis, sr, bd, scheme, tres)
        beaddir = joinpath(outroot, key)
        @info "bead sweep starting" bead = key basis = basis diskRadius_cm = params.diskRadius diskMass_g = params.diskMass scheme = scheme n_cases = length(gammas) * length(freqs)
        try
            run_sweep(SweepSpec(gammas = gammas, bathFrequenciesHz = freqs, fixed = params), beaddir; dtn = dtn)
        catch err
            @error "sweep aborted for this bead; continuing with the next bead" bead = key exception = (err, catch_backtrace())
        end
        rows = postprocess_sweep(beaddir, params, gammas, freqs; bead = key, basis = basis, scheme = String(scheme), domain = domain)
        append!(all_rows, rows)
        counts = Dict{String,Int}()
        for r in rows
            counts[r["quality"]] = get(counts, r["quality"], 0) + 1
        end
        @info "bead finished" bead = key quality_counts = counts
    end

    write_gain_phase(joinpath(outroot, "gain_phase_all.csv"), all_rows)
    total = Dict{String,Int}()
    for r in all_rows
        total[r["quality"]] = get(total, r["quality"], 0) + 1
    end
    @info "campaign stage finished" stage = stage outdir = outroot n_rows = length(all_rows) quality_counts = total
    println("Quality counts: ", total)

    if require_finite
        # Pipeline-integrity check only: a case may legitimately have no fit (wave_cut_short,
        # too_short), but it must not be an error, missing, nonfinite or a post-processing failure.
        bad = filter(r -> r["quality"] in ("error", "missing", "nonfinite", "postproc_error"), all_rows)
        isempty(bad) || error("require_finite=true but $(length(bad)) of $(length(all_rows)) cases failed: " *
                              join(unique(r["quality"] for r in bad), ", "))
    end
    return nothing
end

# Only run when executed directly, so that repostprocess.jl can `include` this file.
# (`... == @__FILE__ && main(ARGS)` is a parse error: a macro call swallows what follows it.)
if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
