#= Gain and phase lag per (bead, gamma, frequency), lab and shaker frame, plus quality flags.
(Same conventions as scripts/run_validation_overlay.jl.)

  bath displacement  b(t) = A_ref cos(ωt + φ_b),  A_ref = Γ g / ω²,  φ_b = deg2rad(phaseDifference) = -90°
  shaker frame       z_s(t) = saved trace `CoM_cm` (bead height relative to the bath)
  lab frame          z_l(t) = z_s(t) + b(t)
  fit                z ≈ a cos(ωt) + b sin(ωt) + c, least squares over the last min(3, available) periods (`fit_oscillation`)
  gain               fitted amplitude / A_ref
  phase lag          mod(φ_b - fitted phase, 360°): lag behind the bath displacement; same reference in both frames
  wave cut           trace cut at the first time |eta_boundary_cm| > 0.10 A_ref

quality: good | unconverged (full time, criterion unmet) | wave_cut | wave_cut_short (<2 periods left) | too_short (<2 periods)
         | error (threw or diverged) | nonfinite | missing (no CSV) | postproc_error
         | subharmonic (|f/2| > 0.10 |f| over the last <= 8 periods; gain and phase not reliable; overrides good/unconverged/wave_cut)
echo_pred: a-priori wall-echo estimate (see `echo_prediction`); echo_ok = echo_pred < 0.01. =#

using Printf

const GP_COLUMNS = ["bead", "basis", "scheme", "domain", "gamma", "freq_Hz", "A_ref_cm", "diskRadius_cm",
                    "sim_status", "quality", "periods_total", "periods_used", "wave_cut", "echo_pred",
                    "echo_ok", "converged", "amp_rel_change", "phase_change_deg",
                    "gain_lab", "phase_lag_lab_deg", "gain_shaker", "phase_lag_shaker_deg",
                    "max_rel_motion_over_a", "subharmonic_ratio", "gamma_gt_1"]

_fmt(x::AbstractString) = String(x)
_fmt(x::Bool) = x ? "true" : "false"
_fmt(x::Integer) = string(x)
_fmt(x::Real) = isfinite(x) ? @sprintf("%.8g", x) : "NaN"

function blank_row(; bead, basis, scheme, domain, gamma, freq, params::SimulationParams, sim_status, quality)
    ω = 2π * freq
    row = Dict{String,Any}()
    for c in GP_COLUMNS
        row[c] = NaN
    end
    row["bead"] = bead
    row["basis"] = basis
    row["scheme"] = scheme
    row["domain"] = domain
    row["gamma"] = gamma
    row["freq_Hz"] = freq
    row["A_ref_cm"] = gamma * params.g / ω^2
    row["diskRadius_cm"] = params.diskRadius
    row["sim_status"] = sim_status
    row["quality"] = quality
    row["wave_cut"] = false
    row["echo_pred"] = echo_prediction(freq, params)
    row["echo_ok"] = row["echo_pred"] < 0.01
    row["converged"] = false
    row["gamma_gt_1"] = gamma > 1
    return row
end

"""
echo_prediction(freq, params) -> Float64

A-priori (not measured) wall-echo amplitude factor exp(-2 R_domain / L_e); R_domain = bathDiameter/2 * diskRadius,
L_e = c_g / (2 ν k²): viscous decay length of a deep-water capillary-gravity wave at the drive frequency.
Geometric spreading ignored.
"""
function echo_prediction(freq::Real, params::SimulationParams)
    ω = 2π * freq
    s = params.bathSurfaceTension / params.bathDensity
    g = params.g
    k = ω^2 / g                       # starts to the right of the root: Newton converges monotonically
    for _ in 1:100
        k -= (g * k + s * k^3 - ω^2) / (g + 3 * s * k^2)
    end
    cg = (g + 3 * s * k^2) / (2ω)
    Le = cg / (2 * params.bathViscosity * k^2)
    R_domain = 0.5 * params.bathDiameter * params.diskRadius
    return exp(-2 * R_domain / Le)
end

function read_case_trace(path::AbstractString)
    lines = readlines(path)
    length(lines) >= 2 || error("empty case CSV: $path")
    t = Float64[]
    z = Float64[]
    eb = Float64[]
    for l in @view lines[2:end]
        isempty(strip(l)) && continue
        p = split(l, ',')
        push!(t, parse(Float64, p[1]))
        push!(z, parse(Float64, p[2]))
        push!(eb, parse(Float64, p[3]))
    end
    return t, z, eb
end

"Map (gamma, freq) -> status string from a sweep's `summary.csv`; empty if the file is absent."
function read_summary_status(path::AbstractString)
    out = Dict{Tuple{Float64,Float64},String}()
    isfile(path) || return out
    lines = readlines(path)
    for l in @view lines[2:end]
        isempty(strip(l)) && continue
        p = split(l, ',')
        length(p) >= 7 || continue
        out[(round(parse(Float64, p[1]); digits = 4), round(parse(Float64, p[2]); digits = 3))] = String(strip(p[7]))
    end
    return out
end

"""
subharmonic_ratio(t, z, ω) -> Float64

Least squares z ≈ a cos(ωt/2) + b sin(ωt/2) + c cos(ωt) + d sin(ωt) + e; returns |f/2 amplitude| / |f amplitude|.
"""
function subharmonic_ratio(t::AbstractVector{<:Real}, z::AbstractVector{<:Real}, ω::Real)
    M = [cos.(ω .* t ./ 2) sin.(ω .* t ./ 2) cos.(ω .* t) sin.(ω .* t) ones(length(t))]
    X = M \ collect(z)
    return hypot(X[1], X[2]) / max(hypot(X[3], X[4]), eps())
end

function gain_phase_row(params::SimulationParams, gamma::Real, freq::Real, csv_path::AbstractString,
                         sim_status::AbstractString; bead, basis, scheme, domain)
    row = blank_row(; bead, basis, scheme, domain, gamma, freq, params, sim_status, quality = "missing")
    isfile(csv_path) || return row
    sim_status == "error" && (row["quality"] = "error"; return row)

    t, z, eb = read_case_trace(csv_path)
    if !(all(isfinite, t) && all(isfinite, z) && all(isfinite, eb))
        row["quality"] = "nonfinite"
        return row
    end

    ω = 2π * freq
    T = 1 / freq
    A_ref = row["A_ref_cm"]
    φb = deg2rad(params.phaseDifference)
    row["periods_total"] = (t[end] - t[1]) / T
    if length(t) < 10
        row["quality"] = "too_short"
        return row
    end

    wave_idx = findfirst(x -> abs(x) > 0.10 * A_ref, eb)
    row["wave_cut"] = wave_idx !== nothing
    cut = wave_idx === nothing ? length(t) : max(1, wave_idx - 1)
    tt = t[1:cut]
    zz = z[1:cut]
    n_avail = floor(Int, (tt[end] - tt[1]) / T)
    row["periods_used"] = (tt[end] - tt[1]) / T
    if n_avail < 2
        row["quality"] = row["wave_cut"] ? "wave_cut_short" : "too_short"
        return row
    end

    n_eval = min(3, n_avail)
    keep = tt .>= (tt[end] - n_eval * T)
    te = tt[keep]
    ze = zz[keep]
    z_lab = ze .+ A_ref .* cos.(ω .* te .+ φb)
    fit_lab = fit_oscillation(te, z_lab, ω)
    fit_shk = fit_oscillation(te, ze, ω)

# The CSV stores time to 7 digits, so t[end]/T can read 9.999995 and fail the 10-period check; 1e-3 periods of slack (0.06 step).
    o = params.convergence
    conv_opts = ConvergenceOptions(enabled = o.enabled, minPeriods = o.minPeriods - 1e-3,
                                   windowPeriods = o.windowPeriods - 1e-3, checkEveryPeriods = o.checkEveryPeriods,
                                   amplitudeRelTol = o.amplitudeRelTol, phaseAbsTolDeg = o.phaseAbsTolDeg)
    conv = check_convergence(tt, zz, ω, conv_opts)
    row["converged"] = conv.converged
    row["amp_rel_change"] = conv.amplitude_rel_change
    row["phase_change_deg"] = conv.phase_change_deg

    row["gain_lab"] = fit_lab.amplitude / A_ref
    row["phase_lag_lab_deg"] = mod(rad2deg(φb - fit_lab.phase), 360.0)
    row["gain_shaker"] = fit_shk.amplitude / A_ref
    row["phase_lag_shaker_deg"] = mod(rad2deg(φb - fit_shk.phase), 360.0)
    row["max_rel_motion_over_a"] = (maximum(ze) - minimum(ze)) / 2 / params.diskRadius

    row["quality"] = row["wave_cut"] ? "wave_cut" : (conv.converged ? "good" : "unconverged")

    n_sub = min(8, n_avail)
    if n_sub >= 4        # at least two cycles of f/2 are needed to separate it from f
        keep_sub = tt .>= (tt[end] - n_sub * T)
        sub = subharmonic_ratio(tt[keep_sub], zz[keep_sub], ω)
        row["subharmonic_ratio"] = sub
        isfinite(sub) && sub > 0.10 && (row["quality"] = "subharmonic")
    end
    return row
end

function write_gain_phase(path::AbstractString, rows::AbstractVector)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, join(GP_COLUMNS, ","))
        for r in rows
            println(io, join((_fmt(r[c]) for c in GP_COLUMNS), ","))
        end
    end
    return path
end

"""
postprocess_sweep(outdir, params, gammas, freqs; bead, basis, scheme, domain) -> Vector{Dict}

One row per (gamma, freq), grid order. Missing CSV -> `missing`; a failure in one case never stops the others.
"""
function postprocess_sweep(outdir::AbstractString, params::SimulationParams, gammas, freqs; bead, basis, scheme, domain)
    status = read_summary_status(joinpath(outdir, "summary.csv"))
    rows = Dict{String,Any}[]
    for γ in gammas, f in freqs
        st = get(status, (round(γ; digits = 4), round(f; digits = 3)), "missing")
        row = try
            gain_phase_row(params, γ, f, case_csv_path(outdir, γ, f), st; bead, basis, scheme, domain)
        catch err
            @error "post-processing failed for one case" bead = bead gamma = γ freqHz = f exception = (err, catch_backtrace())
            blank_row(; bead, basis, scheme, domain, gamma = γ, freq = f, params, sim_status = st, quality = "postproc_error")
        end
        push!(rows, row)
    end
    write_gain_phase(joinpath(outdir, "gain_phase.csv"), rows)
    return rows
end
