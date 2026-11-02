#!/usr/bin/env python3
"""Gate between the acceleration test and the main sweeps (standard library only).

Question: does the PTFE gain / phase-lag curve depend on the acceleration Gamma?

Inputs (a campaign directory made by submit_campaign.sh):
  03_accel/gain_phase_all.csv                      Gamma <= 1 (the decision uses these)
  03b_accel_gt1/gain_phase_all.csv                 Gamma > 1  (optional, information only)
  06_numerics/euler_tres60/gain_phase_all.csv      Euler, 60 steps/period   (optional)
  06_numerics/bdf2_tres30/gain_phase_all.csv       BDF2, 30 steps/period    (optional)

Rule (all thresholds are choices of the analyst, not of the colleague):
  - Use only rows with quality == "good" (converged, no wall echo). A frequency enters the
    comparison only if EVERY Gamma <= 1 has a good row there.
  - Spread at a frequency, lab frame: (max - min) / mean of the gain over Gamma;
    largest pairwise circular difference of the phase lag.
  - Noise floor: how much the same case changes with the numerics (Euler vs BDF2, and 30
    vs 60 steps), at the reference Gamma. Threshold = max(noise floor, convergence
    tolerance), where the convergence tolerance is the code's own (ConvergenceOptions):
    2% in amplitude, 2 degrees in phase. A fitted value cannot be trusted better than that.
  - A frequency is flagged when its spread exceeds the threshold. The curve is DEPENDENT on
    Gamma when more than 10% of the compared frequencies are flagged.
  - Shaker-frame numbers and Gamma > 1 are reported, but do not decide.

Exit code: 0 = independent (the main sweeps may run), 3 = dependent (stop, check in with the
colleague), 4 = not enough good data to decide (stop).
"""
import argparse
import csv
import math
import os
import sys

GAIN_FLOOR = 0.02        # relative, = ConvergenceOptions.amplitudeRelTol
PHASE_FLOOR = 2.0        # degrees, = ConvergenceOptions.phaseAbsTolDeg
FLAG_FRACTION = 0.10
MIN_FREQS_ABS = 5
MIN_FREQS_FRAC = 0.5


def read_rows(path):
    if not os.path.isfile(path):
        return None
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def cdiff(a, b):
    return (a - b + 180.0) % 360.0 - 180.0


def good(rows, bead):
    out = {}
    for r in rows:
        if r["bead"] != bead or r["quality"] != "good":
            continue
        if not all(math.isfinite(num(r[c])) for c in ("gain_lab", "phase_lag_lab_deg", "gain_shaker", "phase_lag_shaker_deg")):
            continue
        out[(round(num(r["gamma"]), 4), round(num(r["freq_Hz"]), 3))] = r
    return out


def spread(values_gain, values_phase):
    mean = sum(values_gain) / len(values_gain)
    rel = (max(values_gain) - min(values_gain)) / mean if mean else float("inf")
    ph = max(abs(cdiff(a, b)) for a in values_phase for b in values_phase)
    return rel, ph


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("campaign")
    ap.add_argument("--bead", default="ptfe")
    ap.add_argument("--ref-gamma", type=float, default=0.5)
    ap.add_argument("--main-gamma", type=float, default=0.5)
    a = ap.parse_args()

    out_dir = os.path.join(a.campaign, "gate")
    os.makedirs(out_dir, exist_ok=True)
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    def finish(code, verdict):
        with open(os.path.join(out_dir, "gate_report.md"), "w") as f:
            f.write("\n".join(lines) + "\n")
        with open(os.path.join(out_dir, "gate_decision.txt"), "w") as f:
            f.write(f"{verdict}\nexit_code={code}\nmain_gamma={a.main_gamma}\n")
        sys.exit(code)

    rows = read_rows(os.path.join(a.campaign, "03_accel", "gain_phase_all.csv"))
    say("# Acceleration gate")
    say(f"Campaign: {a.campaign}   bead: {a.bead}   reference Gamma: {a.ref_gamma}")
    if rows is None:
        say("No 03_accel/gain_phase_all.csv: cannot decide.")
        finish(4, "INSUFFICIENT_DATA")

    table = good(rows, a.bead)
    all_rows_bead = [r for r in rows if r["bead"] == a.bead]
    gammas = sorted({round(num(r["gamma"]), 4) for r in all_rows_bead if num(r["gamma"]) <= 1.0})
    freqs = sorted({round(num(r["freq_Hz"]), 3) for r in all_rows_bead})
    counts = {}
    for r in all_rows_bead:
        counts[r["quality"]] = counts.get(r["quality"], 0) + 1
    say(f"Gamma values (<= 1): {gammas}; frequencies in the grid: {len(freqs)}; quality counts: {counts}")

    f_ok = [f for f in freqs if all((g, f) in table for g in gammas)]
    need = max(MIN_FREQS_ABS, math.ceil(MIN_FREQS_FRAC * len(freqs)))
    say(f"Frequencies where every Gamma has a good row: {len(f_ok)} (need at least {need})")
    if len(gammas) < 2 or len(f_ok) < need:
        say("VERDICT: INSUFFICIENT DATA. Too few frequencies with good rows at every Gamma.")
        finish(4, "INSUFFICIENT_DATA")
    if a.main_gamma not in gammas:
        say(f"VERDICT: INSUFFICIENT DATA. The main-sweep Gamma {a.main_gamma} is not in the tested list.")
        finish(4, "INSUFFICIENT_DATA")

    # --- noise floor from the numerics runs, at the reference Gamma ---
    deltas = {}
    for name in ("euler_tres60", "bdf2_tres30"):
        nrows = read_rows(os.path.join(a.campaign, "06_numerics", name, "gain_phase_all.csv"))
        if nrows is None:
            say(f"Numerics run {name}: missing (floor falls back to the convergence tolerance).")
            continue
        var = good(nrows, a.bead)
        dg, dp, n = 0.0, 0.0, 0
        for (g, f), r in var.items():
            ref = table.get((round(a.ref_gamma, 4), f))
            if ref is None or abs(g - a.ref_gamma) > 1e-6:
                continue
            dg = max(dg, abs(num(r["gain_lab"]) - num(ref["gain_lab"])) / num(ref["gain_lab"]))
            dp = max(dp, abs(cdiff(num(r["phase_lag_lab_deg"]), num(ref["phase_lag_lab_deg"]))))
            n += 1
        if n == 0:
            say(f"Numerics run {name}: no good rows in common with the reference (ignored).")
            continue
        deltas[name] = (dg, dp, n)
        say(f"Numerics run {name}: {n} frequencies compared; max gain change {dg:.4f} (relative), max phase change {dp:.2f} deg")
    thr_gain = max([GAIN_FLOOR] + [d[0] for d in deltas.values()])
    thr_phase = max([PHASE_FLOOR] + [d[1] for d in deltas.values()])
    say(f"Thresholds: gain spread > {thr_gain:.4f} (relative), phase spread > {thr_phase:.2f} deg")

    # --- spread across Gamma, per frequency ---
    flagged, table_lines = [], []
    sh_flag = 0
    for f in f_ok:
        rr = [table[(g, f)] for g in gammas]
        sg, sp = spread([num(r["gain_lab"]) for r in rr], [num(r["phase_lag_lab_deg"]) for r in rr])
        sgs, sps = spread([num(r["gain_shaker"]) for r in rr], [num(r["phase_lag_shaker_deg"]) for r in rr])
        bad = sg > thr_gain or sp > thr_phase
        sh_flag += (sgs > thr_gain or sps > thr_phase)
        if bad:
            flagged.append(f)
        table_lines.append(f"| {f:g} | {sg:.4f} | {sp:.2f} | {sgs:.4f} | {sps:.2f} | {'FLAG' if bad else ''} |")
    say()
    say("| f (Hz) | lab gain spread | lab phase spread (deg) | shaker gain spread | shaker phase spread (deg) | decision |")
    say("|---|---|---|---|---|---|")
    lines.extend(table_lines)
    frac = len(flagged) / len(f_ok)
    say()
    say(f"Flagged (lab frame): {len(flagged)} of {len(f_ok)} frequencies = {frac:.0%}  {flagged}")
    say(f"Informational: shaker frame would flag {sh_flag} of {len(f_ok)}.")

    # --- information only: Gamma > 1 ---
    hi = read_rows(os.path.join(a.campaign, "03b_accel_gt1", "gain_phase_all.csv"))
    if hi is None:
        say("Gamma > 1: no data (that job did not produce a file).")
    else:
        t_hi = good(hi, a.bead)
        g_hi = sorted({g for (g, _) in t_hi})
        n_flag_hi, n_cmp = 0, 0
        for f in f_ok:
            if all((g, f) in t_hi for g in g_hi) and g_hi:
                rr = [table[(g, f)] for g in gammas] + [t_hi[(g, f)] for g in g_hi]
                sg, sp = spread([num(r["gain_lab"]) for r in rr], [num(r["phase_lag_lab_deg"]) for r in rr])
                n_cmp += 1
                n_flag_hi += (sg > thr_gain or sp > thr_phase)
        say(f"Informational, adding Gamma {g_hi}: {n_flag_hi} of {n_cmp} comparable frequencies flagged.")

    if frac > FLAG_FRACTION:
        say("VERDICT: DEPENDENT on Gamma. Stop and check in with the colleague.")
        finish(3, "DEPENDENT")
    say(f"VERDICT: INDEPENDENT of Gamma within the thresholds. The main sweeps run at Gamma = {a.main_gamma}.")
    finish(0, "INDEPENDENT")


if __name__ == "__main__":
    main()
