#!/usr/bin/env julia
#= Graded-grid pilot: C1, C2, C3 and C6 of docs/graded_grid.md (criteria fixed there, before any run).

Usage: julia -t 8 --project=julia/scripts julia/scripts/graded_pilot.jl <outdir> [parts=ABCD] [budget_min=95]
  A  C1  DtN error on wave packets, graded vs uniform at equal wavelength / dr_max
  B  C2  Laplacian order under refinement
  C  C3  28 bead cases, graded vs the uniform rows of bead_gain_phase_bdf2.csv, for several (Laplacian order, ppw, q)
  D  C6  one uniform nr = 2500 case (25 x 200) timed in the same job, for the cost comparison
Prints tables; C3 rows are also written to <outdir>/c3_rows.csv (appended after every case). New configurations are not started
after `budget_min` minutes. =#

if !@isdefined(FaradayDisk)
    include("_bootstrap.jl")
end
using Printf, Dates, Logging, LinearAlgebra
include("graded_grid_checks.jl")
include(joinpath(@__DIR__, "..", "cluster", "beads", "run_bead_sweep.jl"))    # bead_params, postprocess_sweep (main() does not run)

const REF_CSV = "/oscar/data/dharri15/eaguerov/Github/km-disk-vibrating/output/beads/campaign_20261104_174348/deliverables/bead_gain_phase_bdf2.csv"
const BEAD_CSV_NAME = Dict("al32" => "Al_1/32in", "al16" => "Al_1/16in", "nylon" => "Nylon_1/16in", "ptfe" => "PTFE")

# (domain, bead, basis, freq Hz): all have quality = good in REF_CSV (checked below)
const C3_CASES = vcat(
    [("25x700", b, k, f) for (b, k, fs) in (("al32", "rc", (20.0, 30.0, 45.0, 60.0, 90.0)), ("al32", "R", (15.0, 25.0, 50.0)),
                                            ("al16", "rc", (25.0, 40.0)), ("al16", "R", (20.0,)), ("nylon", "rc", (15.0, 35.0, 52.0)),
                                            ("nylon", "R", (25.0,)), ("ptfe", "rc", (15.0, 30.0, 45.0)), ("ptfe", "R", (20.0,))) for f in fs],
    [("25x200", b, k, f) for (b, k, fs) in (("al16", "R", (40.0, 70.0, 90.0)), ("nylon", "R", (30.0, 60.0)), ("ptfe", "R", (45.0, 80.0)),
                                            ("ptfe", "rc", (60.0, 90.0)), ("al16", "rc", (55.0,)), ("nylon", "rc", (65.0,))) for f in fs],
)

# (laplacian_order, ppw, q)
const C3_CONFIGS = [(4, 12.0, 1.04), (4, 16.0, 1.04), (2, 16.0, 1.04), (2, 24.0, 1.04), (2, 32.0, 1.04), (4, 8.0, 1.04),
                    (4, 12.0, 1.08), (4, 12.0, 1.02), (2, 12.0, 1.04)]

function read_ref()
    lines = readlines(REF_CSV)
    hdr = split(lines[1], ',')
    ref = Dict{Tuple{String,String,String,Float64},Tuple{Float64,Float64}}()
    for ln in lines[2:end]
        isempty(strip(ln)) && continue
        v = split(ln, ',')
        row = Dict(hdr[i] => v[i] for i in eachindex(hdr))
        (row["role"] == "production" && row["quality"] == "good") || continue
        ref[(row["domain"], row["bead"], row["radius_basis"], parse(Float64, row["freq_Hz"]))] =
            (parse(Float64, row["gain_lab"]), parse(Float64, row["phase_lag_lab_deg"]))
    end
    return ref
end

wrap180(x) = mod(x + 180, 360) - 180

# ------------------------------------------------------------------------------------------------ A: C1
function part_C1()
    println("\n=== C1: DtN on wave packets exp(-r^2/2s^2) J0(k0 r), graded vs uniform spacing dr_max (E_loc, E_glob as in docs/graded_grid.md) ===")
    sr, R, dmax, rin = 25.0, 120.0, 0.5, 2.0
    s = R / 4
    @printf("%-5s %-6s | %-6s %-9s %-9s | %-6s %-9s %-9s | %-6s %-5s\n", "ppw", "q", "nr_g", "Eloc_g", "Eglob_g", "nr_u", "Eloc_u", "Eglob_u", "ratio", "C1")
    for ppw in (8.0, 12.0, 16.0, 24.0)
        lam = ppw * dmax
        k0 = 2pi / lam
        ug = uniform_grid(1 / dmax, R)
        Au = generate_dtn(ug)
        rsu = ug.r[1:ug.nr]
        fu = @. exp(-rsu^2 / (2 * s^2)) * _gg_besselj0(k0 * rsu)
        exu = _gg_exact_packet(s, k0, rsu)
        eu = _gg_local_error(Au, fu, exu, rsu, lam)
        gu = maximum(abs, Au * fu .- exu) / maximum(abs, exu)
        for q in (1.03, 1.05, 1.08)
            g = graded_grid(sr, R; q = q, dr_max = dmax, r_in = rin)
            A = generate_dtn(g)
            rs = g.r[1:g.nr]
            f = @. exp(-rs^2 / (2 * s^2)) * _gg_besselj0(k0 * rs)
            ex = _gg_exact_packet(s, k0, rs)
            eg = _gg_local_error(A, f, ex, rs, lam)
            gg = maximum(abs, A * f .- ex) / maximum(abs, ex)
            @printf("%-5g %-6g | %-6d %-9.2e %-9.2e | %-6d %-9.2e %-9.2e | %-6.2f %-5s\n", ppw, q, g.nr, eg, gg, ug.nr, eu, gu, eg / eu,
                    eg <= 1.5 * eu ? "PASS" : "FAIL")
            flush(stdout)
        end
    end
end

# ------------------------------------------------------------------------------------------------ B: C2
function part_C2()
    println("\n=== C2: Laplacian error under refinement p (g = exp(-r^2/50); sr = 5p, q = 1.15^(1/p), dr_max = 1/p, ramp = 10p) ===")
    forms = (:cons2, :fd4, :fd6)
    errs = Dict(f => [_gg_laplacian_error(p, f) for p in (1, 2, 4, 8)] for f in forms)
    @printf("%-6s %-12s %-12s %-12s %-12s | orders 1->2, 2->4, 4->8\n", "form", "p=1", "p=2", "p=4", "p=8")
    for f in forms
        e = errs[f]
        @printf("%-6s %-12.3e %-12.3e %-12.3e %-12.3e | %.2f %.2f %.2f\n", f, e[1], e[2], e[3], e[4], log2(e[1] / e[2]), log2(e[2] / e[3]), log2(e[3] / e[4]))
    end
    o2 = log2(errs[:cons2][2] / errs[:cons2][3])
    o4 = log2(errs[:fd4][2] / errs[:fd4][3])
    @printf("C2 conservative (tridiagonal) order %.2f (need >= 1.9): %s\n", o2, o2 >= 1.9 ? "PASS" : "FAIL")
    @printf("C2 4th-order form order %.2f (need >= 3.5): %s\n", o4, o4 >= 3.5 ? "PASS" : "FAIL")
end

# ------------------------------------------------------------------------------------------------ C: C3
function run_graded_case(dom, key, basis, f, cfg, outroot; bd = dom == "25x200" ? 200.0 : 700.0)
    order, ppw, q = cfg
    params = bead_params(key, basis, 25.0, bd, :bdf2, 60)
    t0 = time()
    grid = graded_grid_for(params, f; ppw = ppw, q = q, laplacian_order = order)
    t_grid = time() - t0
    params = with_overrides(params; grid = grid)
    t0 = time()
    dtn = generate_dtn(grid)
    t_dtn = time() - t0
    dir = joinpath(outroot, "o$(order)_ppw$(Int(ppw))_q$(q)", "$(key)_$(basis)_$(dom)_bd$(Int(bd))_f$(Int(f))")
    t0 = time()
    run_sweep(SweepSpec(gammas = [0.5], bathFrequenciesHz = [f], fixed = params), dir; dtn = dtn)
    t_run = time() - t0
    row = postprocess_sweep(dir, params, [0.5], [f]; bead = key, basis = basis, scheme = "bdf2", domain = dom)[1]
    n_fact = lu_factor_count(params.temporalResolution, deg2rad(params.phaseDifference), true)
    return (nr = grid.nr, wall_cm = wall_radius(grid) * params.diskRadius, t_grid = t_grid, t_dtn = t_dtn, t_run = t_run,
            lu_bytes = estimated_lu_cache_bytes(n_fact, 2 * grid.nr + 2), quality = String(row["quality"]),
            gain = Float64(row["gain_lab"]), phase = Float64(row["phase_lag_lab_deg"]), wave_cut = row["wave_cut"], echo = Float64(row["echo_pred"]))
end

function part_C3(outdir, t_start, budget_min)
    println("\n=== C3: graded vs uniform rows (criterion per row: |d gain| <= 0.3 %, |d phase| <= 0.3 deg; all 28 rows must pass) ===")
    ref = read_ref()
    cases = [c for c in C3_CASES if haskey(ref, (c[1], BEAD_CSV_NAME[c[2]], c[3], c[4]))]
    missing_ref = [c for c in C3_CASES if !haskey(ref, (c[1], BEAD_CSV_NAME[c[2]], c[3], c[4]))]
    isempty(missing_ref) || println("REFERENCE ROWS MISSING OR NOT good (skipped): ", missing_ref)
    println("cases: ", length(cases), " of ", length(C3_CASES))
    csv = joinpath(outdir, "c3_rows.csv")
    open(csv, "w") do io
        println(io, "order,ppw,q,domain,bead,basis,freq_Hz,nr,wall_cm,t_grid_s,t_dtn_s,t_run_s,lu_bytes,quality,gain,gain_ref,dgain_pct,phase,phase_ref,dphase_deg,pass")
    end
    summary = Tuple[]
    for cfg in C3_CONFIGS
        if (time() - t_start) / 60 > budget_min
            println("BUDGET REACHED: configuration ", cfg, " and later ones not started")
            break
        end
        order, ppw, q = cfg
        println("\n--- config: laplacian_order=$order ppw=$ppw q=$q ---")
        @printf("%-7s %-6s %-5s %-4s | %-5s %-8s %-8s | %-8s %-9s %-9s | %-9s %-9s %-8s %-5s\n", "domain", "bead", "basis", "f", "nr", "t_dtn_s", "t_run_s",
                "quality", "gain", "dgain_%", "phase", "dphase", "", "pass")
        npass = 0
        nrs = Int[]
        dg = Float64[]
        dp = Float64[]
        for (dom, key, basis, f) in cases
            gref, pref = ref[(dom, BEAD_CSV_NAME[key], basis, f)]
            try
                r = run_graded_case(dom, key, basis, f, cfg, joinpath(outdir, "runs"))
                dgp = 100 * (r.gain / gref - 1)
                dph = wrap180(r.phase - pref)
                ok = r.quality == "good" && isfinite(dgp) && isfinite(dph) && abs(dgp) <= 0.3 && abs(dph) <= 0.3
                npass += ok
                push!(nrs, r.nr)
                isfinite(dgp) && push!(dg, abs(dgp))
                isfinite(dph) && push!(dp, abs(dph))
                @printf("%-7s %-6s %-5s %-4g | %-5d %-8.2f %-8.1f | %-8s %-9.5f %-+9.3f | %-9.4f %-+9.3f %-8s %-5s\n", dom, key, basis, f, r.nr, r.t_dtn, r.t_run,
                        r.quality, r.gain, dgp, r.phase, dph, "", ok ? "PASS" : "FAIL")
                open(csv, "a") do io
                    @printf(io, "%d,%g,%g,%s,%s,%s,%g,%d,%.6g,%.4g,%.4g,%.4g,%d,%s,%.8g,%.8g,%.5g,%.8g,%.8g,%.5g,%s\n", order, ppw, q, dom, key, basis, f, r.nr, r.wall_cm,
                            r.t_grid, r.t_dtn, r.t_run, r.lu_bytes, r.quality, r.gain, gref, dgp, r.phase, pref, dph, ok)
                end
            catch err
                @printf("%-7s %-6s %-5s %-4g | ERROR %s\n", dom, key, basis, f, sprint(showerror, err))
            end
            flush(stdout)
        end
        push!(summary, (cfg, npass, length(cases), isempty(nrs) ? 0 : minimum(nrs), isempty(nrs) ? 0 : maximum(nrs),
                        isempty(dg) ? NaN : maximum(dg), isempty(dp) ? NaN : maximum(dp)))
    end
    println("\n=== C3 summary (config = laplacian_order, ppw, q) ===")
    @printf("%-22s %-10s %-8s %-8s %-12s %-12s %-8s\n", "config", "rows pass", "nr min", "nr max", "max|dgain|%", "max|dphase|", "C3")
    for (cfg, npass, n, nmin, nmax, mg, mp) in summary
        @printf("%-22s %-10s %-8d %-8d %-12.3f %-12.3f %-8s\n", string(cfg), "$npass/$n", nmin, nmax, mg, mp, npass == n ? "PASS" : "FAIL")
    end
end

# ------------------------------------------------------------------------------------------------ D: C6
function part_C6(outdir)
    println("\n=== C6: cost of one uniform nr = 2500 case (al16 R, 25 x 200, 40 Hz), measured here ===")
    println("peak RSS of this process before the uniform case (graded work only): ", round(Sys.maxrss() / 2^30; digits = 2), " GiB")
    params = bead_params("al16", "R", 25.0, 200.0, :bdf2, 60)
    t0 = time()
    dtn = generate_dtn(2500, 200.0)
    t_dtn = time() - t0
    n_fact = lu_factor_count(params.temporalResolution, deg2rad(params.phaseDifference), true)
    t0 = time()
    run_sweep(SweepSpec(gammas = [0.5], bathFrequenciesHz = [40.0], fixed = params), joinpath(outdir, "runs", "uniform_2500"); dtn = dtn)
    t_run = time() - t0
    @printf("uniform nr=2500: DtN build %.1f s, run_sweep %.1f s, LU cache estimate %.2f GiB (%d factorizations of n=5002), peak RSS now %.2f GiB\n",
            t_dtn, t_run, estimated_lu_cache_bytes(n_fact, 5002) / 2^30, n_fact, Sys.maxrss() / 2^30)
    @printf("(for scale) uniform nr=8750: LU cache estimate %.1f GiB; measured earlier (08_resim_small.sh): 9 to 19 min per case, peak 83 GiB RSS, DtN build 222 s\n",
            estimated_lu_cache_bytes(n_fact, 2 * 8750 + 2) / 2^30)
    row = postprocess_sweep(joinpath(outdir, "runs", "uniform_2500"), params, [0.5], [40.0]; bead = "al16", basis = "R", scheme = "bdf2", domain = "25x200")[1]
    @printf("uniform result: quality %s gain_lab %.6f phase_lab %.4f deg\n", row["quality"], row["gain_lab"], row["phase_lag_lab_deg"])
end

# ------------------------------------------------------------------------------------------------ E: C3b (5 Hz rows)
const C3B_CONFIGS = [(2, 24.0, 1.04), (2, 32.0, 1.04)]

function part_C3b(outdir, t_start, budget_min)
    println("\n=== C3b: every quality=good 5 Hz row of the 25x700 domain; same criterion as C3 (|d gain| <= 0.3 %, |d phase| <= 0.3 deg) ===")
    ref = read_ref()
    inv = Dict(v => k for (k, v) in BEAD_CSV_NAME)
    cases = [(d, inv[b], k, f) for (d, b, k, f) in sort(collect(keys(ref))) if d == "25x700" && f == 5.0]
    println("rows: ", length(cases), " ", cases)
    for cfg in C3B_CONFIGS
        (time() - t_start) / 60 > budget_min && (println("BUDGET REACHED: C3b config ", cfg, " not started"); break)
        println("\n--- config (order, ppw, q) = ", cfg, " ---")
        npass = 0
        for (dom, key, basis, f) in cases
            gref, pref = ref[(dom, BEAD_CSV_NAME[key], basis, f)]
            try
                r = run_graded_case(dom, key, basis, f, cfg, joinpath(outdir, "runs_c3b"))
                dgp = 100 * (r.gain / gref - 1)
                dph = wrap180(r.phase - pref)
                ok = r.quality == "good" && isfinite(dgp) && isfinite(dph) && abs(dgp) <= 0.3 && abs(dph) <= 0.3
                npass += ok
                @printf("%-5s %-3s nr=%-5d t_dtn=%-6.1f t_run=%-7.1f %-8s gain %.5f (ref %.5f, %+.3f %%) phase %.3f (ref %.3f, %+.3f deg) %s\n", key, basis, r.nr, r.t_dtn, r.t_run,
                        r.quality, r.gain, gref, dgp, r.phase, pref, dph, ok ? "PASS" : "FAIL")
            catch err
                println(key, " ", basis, " ERROR ", sprint(showerror, err))
            end
            flush(stdout)
        end
        println("C3b config ", cfg, ": ", npass, "/", length(cases), " rows pass")
    end
end

# ------------------------------------------------------------------------------------------------ F: C7
function part_C7(outdir, t_start, budget_min)
    println("\n=== C7 (informational, no pass/fail): al16 R, 5 Hz, domain convergence of the graded grid (order 2, ppw 32, q 1.04) ===")
    cfg = (2, 32.0, 1.04)
    for bd in (700.0, 2000.0, 6000.0, 16000.0)
        (time() - t_start) / 60 > budget_min && (println("BUDGET REACHED: R = ", bd / 2, " not started"); break)
        try
            r = run_graded_case("25x$(Int(bd))", "al16", "R", 5.0, cfg, joinpath(outdir, "runs_c7"); bd = bd)
            @printf("R = %-6g a (bathDiameter %-6g): nr=%-5d wall=%.1f cm  gain %.5f phase %.4f quality %s wave_cut %s  t_dtn %.1f s  t_run %.1f s  LU est %.2f GiB\n",
                    bd / 2, bd, r.nr, r.wall_cm, r.gain, r.phase, r.quality, r.wave_cut, r.t_dtn, r.t_run, r.lu_bytes / 2^30)
        catch err
            println("R = ", bd / 2, " ERROR ", sprint(showerror, err))
        end
        flush(stdout)
    end
end

# ------------------------------------------------------------------------------------------------ G: C7b
# any production row (any quality) of the campaign CSV: (domain, bead, basis, f) => (gain, phase, quality, wall_echo_pred)
function read_ref_any()
    lines = readlines(REF_CSV)
    hdr = split(lines[1], ',')
    ref = Dict{Tuple{String,String,String,Float64},Tuple{Float64,Float64,String,Float64}}()
    for ln in lines[2:end]
        isempty(strip(ln)) && continue
        v = split(ln, ',')
        row = Dict(hdr[i] => v[i] for i in eachindex(hdr))
        row["role"] == "production" || continue
        ref[(row["domain"], row["bead"], row["radius_basis"], parse(Float64, row["freq_Hz"]))] =
            (parse(Float64, row["gain_lab"]), parse(Float64, row["phase_lag_lab_deg"]), row["quality"], parse(Float64, row["wall_echo_pred"]))
    end
    return ref
end

function part_C7b(outdir, t_start, budget_min)
    println("\n=== C7b (informational): al32 rc, graded (order 2, ppw 32, q 1.04) at several R, next to the uniform 25x200 and 25x700 rows ===")
    println("expectation (HYPOTHESIS): graded R=100 a matches the 25x200 row, R=350 a the 25x700 row (0.3 % / 0.3 deg), R >= 1000 a converges")
    ref = read_ref_any()
    cfg = (2, 32.0, 1.04)
    for f in (52.0, 60.0, 70.0)
        println("\n--- al32 rc, ", f, " Hz ---")
        for dom in ("25x200", "25x700")
            k = ("$dom", BEAD_CSV_NAME["al32"], "rc", f)
            if haskey(ref, k)
                g, ph, ql, ec = ref[k]
                @printf("uniform %-7s: gain %.5f phase %.4f quality %-10s wall_echo_pred %.4f\n", dom, g, ph, ql, ec)
            else
                println("uniform ", dom, ": no row in the CSV")
            end
        end
        for Rdisk in (100.0, 200.0, 350.0, 1000.0, 3000.0)
            (time() - t_start) / 60 > budget_min && (println("BUDGET REACHED: R = ", Rdisk, " not started"); continue)
            bd = 2 * Rdisk
            try
                r = run_graded_case("25x$(Int(bd))", "al32", "rc", f, cfg, joinpath(outdir, "runs_c7b"); bd = bd)
                line = @sprintf("graded R=%-5g a: nr=%-5d gain %.5f phase %.4f quality %-10s wall_echo_pred %.4f", Rdisk, r.nr, r.gain, r.phase, r.quality, r.echo)
                for (dom, Rd) in (("25x200", 100.0), ("25x700", 350.0))
                    if Rdisk == Rd && haskey(ref, (dom, BEAD_CSV_NAME["al32"], "rc", f))
                        g, ph, _, _ = ref[(dom, BEAD_CSV_NAME["al32"], "rc", f)]
                        line *= @sprintf("  | vs %s: dgain %+.3f %%, dphase %+.3f deg (limits 0.3 / 0.3)", dom, 100 * (r.gain / g - 1), wrap180(r.phase - ph))
                    end
                end
                println(line)
            catch err
                println("graded R=", Rdisk, " ERROR ", sprint(showerror, err))
            end
            flush(stdout)
        end
    end
end

function pilot_main()
    length(ARGS) >= 1 || error("usage: graded_pilot.jl <outdir> [parts=ABCD] [budget_min=95]")
    outdir = ARGS[1]
    kv = Dict(split(a, '=')[1] => split(a, '=')[2] for a in ARGS[2:end])
    parts = get(kv, "parts", "ABCD")
    budget_min = parse(Float64, get(kv, "budget_min", "95"))
    mkpath(outdir)
    global_logger(setup_logging(name = "graded_pilot", jsonlines = true))
    t_start = time()
    println("graded_pilot: parts=$parts budget_min=$budget_min threads=$(Threads.nthreads()) host=$(gethostname()) started $(Dates.now())")
    occursin('A', parts) && part_C1()
    occursin('B', parts) && part_C2()
    occursin('C', parts) && part_C3(outdir, t_start, budget_min)
    occursin('D', parts) && part_C6(outdir)
    occursin('E', parts) && part_C3b(outdir, t_start, budget_min)
    occursin('F', parts) && part_C7(outdir, t_start, budget_min)
    occursin('G', parts) && part_C7b(outdir, t_start, budget_min)
    println("graded_pilot done in $(round((time() - t_start) / 60; digits = 1)) min")
end

pilot_main()
