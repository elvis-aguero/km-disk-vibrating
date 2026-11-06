#!/usr/bin/env python3
"""Build the deliverable CSV (BDF2 only) from a finished bead campaign.

Usage: python3 -I julia/cluster/beads/build_deliverable.py <campaign_dir>
Writes <campaign_dir>/deliverables/bead_gain_phase_bdf2.csv and README.txt.

Row order: production (04_sweep_R, 05_sweep_rc), gamma_scan (03_accel, 03b_accel_gt1),
numerics_bdf2 (06_numerics/bdf2_tres30). Euler rows are never written. A case that repeats an
earlier row (same bead, basis, steps/period, gamma, freq) is dropped.
Phase lag is wrapped to (-180, 180]; the post-processor (gain_phase.jl) gives mod 360.
"""
import csv, math, os, re, sys

NAMES = {"ptfe": "PTFE", "al32": "Al_1/32in", "al16": "Al_1/16in", "nylon": "Nylon_1/16in"}
COLS = ["role", "bead", "radius_basis", "diskRadius_cm", "scheme", "steps_per_period", "gamma", "freq_Hz",
        "gain_lab", "phase_lag_lab_deg", "gain_shaker", "phase_lag_shaker_deg",
        "wall_echo_pred", "wall_echo_negligible", "quality"]


def wrap(x):
    return (float(x) + 180.0) % 360.0 - 180.0


def fmt(x, n=5):
    return "" if x != x else ("%." + str(n) + "g") % x


def main(camp):
    sources = [("production", "04_sweep_R", 60), ("production", "05_sweep_rc", 60),
               ("gamma_scan", "03_accel", 60), ("gamma_scan", "03b_accel_gt1", 60),
               ("numerics_bdf2", "06_numerics/bdf2_tres30", 30)]
    rows, seen = [], set()
    for role, d, tres in sources:
        path = os.path.join(camp, d, "gain_phase_all.csv")
        if not os.path.isfile(path):
            continue
        for r in csv.DictReader(open(path)):
            if r["scheme"] != "bdf2":
                continue
            key = (r["bead"], r["basis"], tres, float(r["gamma"]), float(r["freq_Hz"]))
            if key in seen:
                continue
            seen.add(key)
            g = lambda c: float(r[c])
            has = g("gain_lab") == g("gain_lab")
            e = g("echo_pred")
            rows.append({"role": role, "bead": NAMES[r["bead"]], "radius_basis": r["basis"],
                         "diskRadius_cm": r["diskRadius_cm"], "scheme": "bdf2", "steps_per_period": tres,
                         "gamma": "%g" % g("gamma"), "freq_Hz": "%g" % g("freq_Hz"),
                         "gain_lab": fmt(g("gain_lab")),
                         "phase_lag_lab_deg": fmt(wrap(r["phase_lag_lab_deg"]) if has else float("nan")),
                         "gain_shaker": fmt(g("gain_shaker")),
                         "phase_lag_shaker_deg": fmt(wrap(r["phase_lag_shaker_deg"]) if has else float("nan")),
                         "wall_echo_pred": "%.3g" % e, "wall_echo_negligible": "true" if e < 0.01 else "false",
                         "quality": r["quality"]})
    order = {"production": 0, "gamma_scan": 1, "numerics_bdf2": 2}
    rows.sort(key=lambda r: order[r["role"]])
    out = os.path.join(camp, "deliverables")
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "bead_gain_phase_bdf2.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS)
        w.writeheader()
        w.writerows(rows)
    commit = "unknown"
    cfg = os.path.join(camp, "04_sweep_R", "run_config.toml")
    if os.path.isfile(cfg):
        m = re.search(r'git_commit\s*=\s*"([^"]+)"', open(cfg).read())
        commit = m.group(1) if m else commit
    readme = f"""bead_gain_phase_bdf2.csv  (BDF2 only; {len(rows)} rows). Campaign: {os.path.basename(camp.rstrip('/'))}, code commit {commit}
(includes the disk force-balance fix of PR 7: contact-line force counted once). Values from earlier campaigns are NOT comparable.
Row order: production (Gamma 0.5, 60 steps/period, 4 beads x R,rc x 34 freqs), gamma_scan (PTFE R, Gamma 0.1-2.0), numerics_bdf2 (PTFE R, 30 steps/period).
A case that repeats an earlier row is dropped.
Model: FaradayDisk.jl, bath-driven, BDF2, domain 25 x 200 (DTN), nr 2500. Fluid: nu 0.079 St, rho 1.148 g/cm3, gamma 68 dyn/cm, g 981. Contact angle 180 deg (sphere mass).
radius_basis: R = actual bead radius, rc = contact radius.
gain = fitted amplitude / (Gamma g / w^2) over the last 3 periods. lab = shaker trace + bath displacement.
phase_lag_*_deg = lag of the bead behind the bath displacement, signed, range (-180, 180].
wall_echo_pred: a-priori estimate (not measured) of the amplitude fraction of a wave that travels to the outer wall and back,
  exp(-2 R_domain / L_e). wall_echo_negligible = (wall_echo_pred < 0.01); where false, the finite domain may contaminate the result.
quality: good | subharmonic (a response at f/2 grew; gain and phase not reliable) | error (no result).
Time-step convergence on this code version has not been re-measured (an earlier study on the pre-fix code gave BDF2/60 errors <= 0.35% lab gain).
Definitions: julia/cluster/beads/gain_phase.jl and fit_oscillation in julia/src/convergence.jl.
"""
    open(os.path.join(out, "README.txt"), "w").write(readme)
    print(f"{len(rows)} rows -> {out}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
