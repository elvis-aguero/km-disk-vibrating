"""
Result I/O: sweep CSVs (schema kept identical to `sweeper.m`'s output so existing or
ported validation/overlay tooling needs no format changes), and a flatter,
content-addressed layout for single-run results — replacing MATLAB's 4-level
parameter-string nested folders (`rho..gcm3-sigma.../diskRadius.../forceAmplitude.../
bathAmplitude...`), whose float-formatting inconsistencies (`%.2g` vs `%g`) can silently
produce different folder names for the same physical configuration.
"""

"""
    compute_eta_boundary(eta_history_cm, spatialResolution) -> Vector{Float64}

Port of `sweeper.m`'s post-processing step: for each recorded timestep (column), the max
`|eta|` over the outermost `ceil(2*spatialResolution)` radial nodes (rows) — used to
detect when a finite-domain wave reflection has contaminated the "steady state."
"""
function compute_eta_boundary(eta_history_cm::AbstractMatrix{<:Real}, spatialResolution::Real)
    nr, nsteps = size(eta_history_cm)
    boundary_nodes = min(nr, ceil(Int, 2 * spatialResolution))
    result = Vector{Float64}(undef, nsteps)
    @inbounds for j in 1:nsteps
        result[j] = maximum(abs, @view eta_history_cm[(nr - boundary_nodes + 1):nr, j])
    end
    return result
end

"""
    write_case_csv(path, t_s, CoM_cm, eta_boundary_cm)

Writes a sweep-case CSV with columns `time_s,CoM_cm,eta_boundary_cm`, matching
`sweeper.m`'s output format and column names exactly.
"""
function write_case_csv(path::AbstractString, t_s::AbstractVector{<:Real}, CoM_cm::AbstractVector{<:Real},
                         eta_boundary_cm::AbstractVector{<:Real})
    n = length(t_s)
    (length(CoM_cm) == n && length(eta_boundary_cm) == n) ||
        throw(DimensionMismatch("t_s, CoM_cm, eta_boundary_cm must have equal length"))
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "time_s,CoM_cm,eta_boundary_cm")
        for i in 1:n
            @printf(io, "%.6e,%.6e,%.6e\n", t_s[i], CoM_cm[i], eta_boundary_cm[i])
        end
    end
    return path
end

"""
    write_timeseries_csv(path, t_s, CoM_cm)

Two-column flat export (`time_s,CoM_cm`) for a single-run result, for quick external
plotting without needing to load the JLD2 file.
"""
function write_timeseries_csv(path::AbstractString, t_s::AbstractVector{<:Real}, CoM_cm::AbstractVector{<:Real})
    length(t_s) == length(CoM_cm) || throw(DimensionMismatch("t_s and CoM_cm must have equal length"))
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "time_s,CoM_cm")
        for i in eachindex(t_s)
            @printf(io, "%.6e,%.6e\n", t_s[i], CoM_cm[i])
        end
    end
    return path
end

"""
    write_summary_csv(path, results)

Writes the sweep summary CSV, columns `gamma,bathFrequency_Hz,bathAmplitude_cm,
CoM_max_cm,CoM_rms_cm,elapsed_s,status` — matching `sweeper.m` exactly. This sweep is
bath-driven (the actual physical setup of this experimental study): `bathAmplitude_cm` is
the real physical bath oscillation amplitude `gamma*g/omega^2` at each case's frequency,
not a derived stand-in for anything else.
"""
function write_summary_csv(path::AbstractString, results::AbstractVector{SweepCaseResult})
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "gamma,bathFrequency_Hz,bathAmplitude_cm,CoM_max_cm,CoM_rms_cm,elapsed_s,status")
        for r in results
            @printf(io, "%.4f,%g,%.6e,%.4e,%.4e,%.1f,%s\n",
                    r.gamma, r.freqHz, r.bathAmplitude, r.CoM_max_cm, r.CoM_rms_cm, r.elapsed_s, String(r.status))
        end
    end
    return path
end

"""
    write_sweep_metadata(path; is_bath_driven)

Records the real forcing convention `run_sweep`/`run_sweep_case` actually used for this
sweep, so `compute_sim_overlay` can key its lab-frame correction off the true physics
parameters instead of inferring it from `summary.csv`'s `bathAmplitude_cm` column. That
column is unreliable for this purpose: whether the sweep is bath-driven (physical bath
oscillation) or disk-forced (direct disk forcing), the column holds the same
`gamma*g/omega^2` quantity (see `write_summary_csv`'s docstring) — bath amplitude in the
former, a merely-derived display value in the latter — so it is never literally zero
either way, and an all-zero check on it can never actually tell the two apart.
"""
function write_sweep_metadata(path::AbstractString; is_bath_driven::Bool)
    mkpath(dirname(path))
    open(path, "w") do io
        TOML.print(io, Dict{String,Any}("is_bath_driven" => is_bath_driven))
    end
    return path
end

_toml_safe(x::AbstractFloat) = isfinite(x) ? x : string(x)
_toml_safe(x) = x

"""
    config_digest(params) -> String

16-hex-char SHA-256 digest of a canonical (sorted-field-name) serialization of every
`SimulationParams` field. Used to build short, reproducible run identifiers — the full
parameter set is always stored alongside in `metadata.toml`, so the digest never has to
be *decoded* to recover the configuration, only compared for equality.
"""
function config_digest(params::SimulationParams)
    io = IOBuffer()
    for f in sort(collect(fieldnames(SimulationParams)); by = string)
        v = getfield(params, f)
        entry = v isa ConvergenceOptions ?
                join(("$(cf)=$(getfield(v, cf))" for cf in fieldnames(ConvergenceOptions)), ";") :
                string(v)
        println(io, string(f), "=", entry)
    end
    return bytes2hex(sha256(take!(io)))[1:16]
end

"""
    run_id(params) -> String

`<timestamp>_<8-char config digest>`, e.g. `20260722_140512_a1b2c3d4`.
"""
function run_id(params::SimulationParams)
    ts = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
    return string(ts, "_", first(config_digest(params), 8))
end

function _params_to_dict(params::SimulationParams)
    d = Dict{String,Any}()
    for f in fieldnames(SimulationParams)
        v = getfield(params, f)
        d[string(f)] = if v isa ConvergenceOptions
            Dict{String,Any}(string(cf) => _toml_safe(getfield(v, cf)) for cf in fieldnames(ConvergenceOptions))
        elseif v isa Symbol
            string(v)
        else
            _toml_safe(v)
        end
    end
    return d
end

"""
    save_result(result, dir)

Saves a `SimulationResult` under `dir`: `result.jld2` (the full result, Julia-native),
`metadata.toml` (every parameter, derived solver/convergence info, and status — full
traceability without needing to decode a directory name), and `com_timeseries.csv` (a
quick flat export).
"""
function save_result(result::SimulationResult, dir::AbstractString)
    mkpath(dir)
    JLD2.jldsave(joinpath(dir, "result.jld2");
                 params = result.params, resolved_solver = String(result.resolved_solver),
                 status = String(result.status), t_s = result.t_s, CoM_cm = result.CoM_cm,
                 eta_history_cm = result.eta_history_cm, convergence = result.convergence,
                 elapsed_s = result.elapsed_s, lu_cache_hit_rate = result.lu_cache_hit_rate)

    write_timeseries_csv(joinpath(dir, "com_timeseries.csv"), result.t_s, result.CoM_cm)

    meta = Dict{String,Any}(
        "status" => String(result.status),
        "resolved_solver" => String(result.resolved_solver),
        "elapsed_s" => result.elapsed_s,
        "lu_cache_hit_rate" => _toml_safe(result.lu_cache_hit_rate),
        "convergence" => Dict{String,Any}(
            "converged" => result.convergence.converged,
            "periods_elapsed" => _toml_safe(result.convergence.periods_elapsed),
            "amplitude" => _toml_safe(result.convergence.amplitude),
            "phase_deg" => _toml_safe(result.convergence.phase_deg),
            "amplitude_rel_change" => _toml_safe(result.convergence.amplitude_rel_change),
            "phase_change_deg" => _toml_safe(result.convergence.phase_change_deg),
        ),
        "params" => _params_to_dict(result.params),
    )
    open(joinpath(dir, "metadata.toml"), "w") do io
        TOML.print(io, meta)
    end
    return dir
end
