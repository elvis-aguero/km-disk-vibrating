"""
Threads.jl-based parameter sweep orchestration — port of `sweeper.m`'s Cartesian-grid
`parfor` loop. Each thread writes only its own index into a preallocated results vector,
mirroring MATLAB's `summaryRows{ii} = {...}` pattern, so there is no shared-mutable-state
race to reason about.
"""

Base.@kwdef struct SweepSpec
    gammas::Vector{Float64}
    bathFrequenciesHz::Vector{Float64}
    fixed::SimulationParams
end

struct SweepCaseResult
    gamma::Float64
    freqHz::Float64
    forcingAmplitude::Float64
    CoM_max_cm::Float64
    CoM_rms_cm::Float64
    elapsed_s::Float64
    status::Symbol   # :ok | :error | :skipped
end

"""
    case_id(gamma, freqHz) -> String

Deterministic case identifier, used consistently by the CSV writer, the resume/skip
check, and any downstream reader — a single source of truth instead of the
format-string-duplication risk of MATLAB's `sprintf` being repeated (with potentially
inconsistent precision) across `sweeper.m` and `overlay_validation.py`.
"""
case_id(gamma::Real, freqHz::Real) = "gamma$(round(gamma, digits = 4))_f$(round(freqHz, digits = 3))Hz"

case_csv_path(outdir::AbstractString, gamma::Real, freqHz::Real) =
    joinpath(outdir, "cases", case_id(gamma, freqHz) * ".csv")

"""
    run_sweep(spec, outdir) -> Vector{SweepCaseResult}

Runs every `(gamma, frequency)` case in `spec`'s Cartesian grid, in parallel via
`Threads.@threads :dynamic`, writing one CSV per case under `outdir/cases/` plus a
`summary.csv`. Resumes a previous partial sweep by skipping any case whose CSV already
exists (same semantics as `sweeper.m`).

The DTN matrix for `(spec.fixed.spatialResolution, spec.fixed.bathDiameter)` is resolved
and loaded exactly once, before the threaded loop starts, and shared read-only across
every case (a sweep holds the domain fixed across its whole grid) — avoiding redundant
per-thread loads of what can be a tens-of-MB matrix. `available_memory_bytes()` (not raw
`Sys.free_memory()` — see its docstring for why that matters under SLURM) is divided by
`Threads.nthreads()` up front so each case's `solverType=:auto` decision accounts for
everything else running concurrently, rather than each case assuming it alone owns the
whole budget (a real robustness gap `getAvailableRAM()`-per-case in MATLAB does not have
to contend with in quite the same way, since MATLAB's `parfor` workers are separate
processes).
"""
function run_sweep(spec::SweepSpec, outdir::AbstractString; dtn_registry::Union{Nothing,DTNManifest} = nothing,
                    dtn::Union{Nothing,AbstractMatrix{<:Real}} = nothing)
    mkpath(joinpath(outdir, "cases"))

    cases = vec(collect(Iterators.product(spec.gammas, spec.bathFrequenciesHz)))
    n = length(cases)
    results = Vector{SweepCaseResult}(undef, n)

    if dtn === nothing
        registry = dtn_registry === nothing ? default_registry() : dtn_registry
        entry = resolve_dtn(registry, spec.fixed.spatialResolution, spec.fixed.bathDiameter)
        dtn = load_dtn(entry)
    end

    nthreads = max(1, Threads.nthreads())
    ram_per_case = max(1, available_memory_bytes() ÷ nthreads)

    @info "sweep starting" n_cases = n gammas = spec.gammas freqs_hz = spec.bathFrequenciesHz threads = nthreads ram_per_case_bytes = ram_per_case outdir = outdir

    Threads.@threads :dynamic for i in 1:n
        gamma, freqHz = cases[i]
        results[i] = run_sweep_case(spec, gamma, freqHz, dtn, outdir, ram_per_case, i, n)
    end

    write_summary_csv(joinpath(outdir, "summary.csv"), results)

    n_ok = count(r -> r.status == :ok, results)
    n_err = count(r -> r.status == :error, results)
    n_skip = count(r -> r.status == :skipped, results)
    @info "sweep finished" ok = n_ok errors = n_err skipped = n_skip summary = joinpath(outdir, "summary.csv")

    return results
end

function run_sweep_case(spec::SweepSpec, gamma::Float64, freqHz::Float64, dtn::AbstractMatrix{Float64},
                         outdir::AbstractString, ram_budget_bytes::Integer, idx::Int, n::Int)
    csv_path = case_csv_path(outdir, gamma, freqHz)
    forcing_amplitude = (gamma * spec.fixed.g) / (2 * pi * freqHz)^2

    if isfile(csv_path)
        @info "case skipped (csv exists)" case = idx of = n gamma = gamma freqHz = freqHz
        return SweepCaseResult(gamma, freqHz, forcing_amplitude, NaN, NaN, 0.0, :skipped)
    end

    t0 = time()
    @info "case starting" case = idx of = n gamma = gamma freqHz = freqHz thread = Threads.threadid()
    try
        F = gamma * spec.fixed.diskMass * spec.fixed.g
        case_params = with_overrides(spec.fixed;
            forceAmplitude = F, forceFrequency = freqHz,
            bathAmplitude = 0.0, bathFrequency = freqHz,
            simulationTime = 30 / freqHz)

        result = run_simulation(case_params; dtn = dtn, ram_budget_bytes = ram_budget_bytes)

        eta_boundary_cm = compute_eta_boundary(result.eta_history_cm, spec.fixed.spatialResolution)
        write_case_csv(csv_path, result.t_s, result.CoM_cm, eta_boundary_cm)

        elapsed = time() - t0
        com_max = isempty(result.CoM_cm) ? NaN : maximum(abs, result.CoM_cm)
        com_rms = isempty(result.CoM_cm) ? NaN : sqrt(sum(abs2, result.CoM_cm) / length(result.CoM_cm))
        status = result.status == :ok ? :ok : :error

        @info "case finished" case = idx of = n gamma = gamma freqHz = freqHz status = status elapsed_s = elapsed CoM_max_cm = com_max CoM_rms_cm = com_rms
        return SweepCaseResult(gamma, freqHz, forcing_amplitude, com_max, com_rms, elapsed, status)
    catch err
        elapsed = time() - t0
        @error "case errored" case = idx of = n gamma = gamma freqHz = freqHz elapsed_s = elapsed exception = (err, catch_backtrace())
        return SweepCaseResult(gamma, freqHz, forcing_amplitude, NaN, NaN, elapsed, :error)
    end
end
