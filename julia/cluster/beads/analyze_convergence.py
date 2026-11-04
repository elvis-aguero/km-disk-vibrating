#!/usr/bin/env python3
"""Time-step convergence analysis (Richardson extrapolation, observed order, GCI).

Usage:  python3 -I julia/cluster/beads/analyze_convergence.py <run_dir> [<run_dir> ...]
Each run_dir holds <scheme>_tres<N>/ptfe/gain_phase.csv (from 07_convergence.sh, 07b_...).
Writes convergence_table.csv and convergence_report.md into the first run_dir.

Method, per (scheme, gamma, frequency, metric), levels with step size h = 1/tres, ratio r = 2:
  observed order  p = log2( (q_h - q_h/2) / (q_h/2 - q_h/4) )   for each consecutive triple
  extrapolated    q* = q_fine + (q_fine - q_prev) / (r^p_theory - 1)   (p_theory: BDF2 2, Euler 1)
  error of level  e_L = q_L - q*          (gain in %, phase lag in deg)
  GCI_fine = 1.25 * |q_fine - q_prev| / (r^p_theory - 1)   (Roache; three-grid safety factor 1.25)
Only levels whose quality is `good` are used. The extrapolation is valid only if the observed
order is close to the theoretical order (the solutions are then in the asymptotic range); the
report says when it is not.
Target (set before looking at the results): error of the production level <= 1 % in gain, <= 1 deg in phase.
"""
import csv, glob, math, os, re, sys

THEORY = {"bdf2": 2, "euler": 1}
METRICS = [("gain_lab", "%"), ("phase_lag_lab_deg", "deg"), ("gain_shaker", "%"), ("phase_lag_shaker_deg", "deg")]
TOL = {"%": 1.0, "deg": 1.0}
PRODUCTION = ("bdf2", 60)


def wrap180(x):
    return (x + 180.0) % 360.0 - 180.0


def load(run_dirs):
    data = {}   # (scheme, tres, gamma, f) -> row
    for rd in run_dirs:
        for p in glob.glob(os.path.join(rd, "*_tres*", "ptfe", "gain_phase.csv")):
            m = re.search(r"([a-z0-9]+)_tres(\d+)", p)
            scheme, tres = m.group(1), int(m.group(2))
            for r in csv.DictReader(open(p)):
                data[(scheme, tres, float(r["gamma"]), float(r["freq_Hz"]))] = r
    return data


def val(row, metric):
    return float(row[metric])


def diff(a, b, unit, metric):
    d = a - b
    return wrap180(d) if unit == "deg" else d


def main(run_dirs):
    data = load(run_dirs)
    if not data:
        sys.exit("no gain_phase.csv files found")
    keys = sorted({(s, g, f) for (s, _, g, f) in data})
    rows, md = [], []
    for scheme, g, f in keys:
        levels = sorted(t for (s, t, gg, ff) in data if s == scheme and gg == g and ff == f)
        good = [t for t in levels if data[(scheme, t, g, f)]["quality"] == "good"]
        for metric, unit in METRICS:
            # relative values for gain (in %), absolute for phase (deg)
            q = {t: val(data[(scheme, t, g, f)], metric) for t in good}
            ref0 = q[good[-1]] if good else float("nan")
            scale = (lambda d: 100.0 * d / ref0) if unit == "%" else (lambda d: d)
            p_obs = []
            for a, b, c in zip(good, good[1:], good[2:]):
                if not (b == 2 * a and c == 2 * b):
                    continue
                d1, d2 = diff(q[a], q[b], unit, metric), diff(q[b], q[c], unit, metric)
                # q is a function of tres: a=coarse, c=fine
                p = math.log2(d1 / d2) if d1 * d2 > 0 and d2 != 0 else float("nan")
                p_obs.append((a, p))
            pt = THEORY[scheme]
            r = {"scheme": scheme, "gamma": g, "freq_Hz": f, "metric": metric, "unit": unit,
                 "levels_good": " ".join(map(str, good)), "levels_all": " ".join(map(str, levels)),
                 "p_theory": pt}
            for a, p in p_obs:
                r[f"p_obs_{a}_{2*a}_{4*a}"] = p
            if len(good) >= 2 and good[-1] == 2 * good[-2]:
                fine, prev = good[-1], good[-2]
                dstar = diff(q[fine], q[prev], unit, metric) / (2 ** pt - 1)
                r["gci_fine"] = 1.25 * abs(scale(dstar)) if unit == "%" else 1.25 * abs(dstar)
                for t in good:
                    r[f"err_{t}"] = scale(diff(q[t], q[fine], unit, metric) - dstar)   # q_t - q*, q* = q_fine + dstar
            rows.append(r)

    cols = sorted({k for r in rows for k in r}, key=lambda k: (k not in ("scheme", "gamma", "freq_Hz", "metric", "unit"), k))
    out_dir = run_dirs[0]
    with open(os.path.join(out_dir, "convergence_table.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: ("%.4g" % v if isinstance(v, float) else v) for k, v in r.items()})

    md.append("# Time-step convergence (PTFE, R, BDF2 and Euler)\n")
    md.append("Observed order p (theory: BDF2 2, Euler 1), estimated error of each level against the "
              "Richardson-extrapolated value, and GCI of the finest level. Gain errors in %, phase errors in deg.\n")
    for scheme in ("bdf2", "euler"):
        md.append(f"\n## {scheme}\n")
        md.append("| gamma | f | metric | p_obs | err by tres | GCI fine |")
        md.append("|---|---|---|---|---|---|")
        for r in rows:
            if r["scheme"] != scheme:
                continue
            ps = ", ".join("%s: %.2f" % (k.replace("p_obs_", ""), v) for k, v in r.items() if k.startswith("p_obs_"))
            es = ", ".join("%s: %+.3g" % (k[4:], v) for k, v in r.items() if k.startswith("err_"))
            gci = r.get("gci_fine", float("nan"))
            md.append(f"| {r['gamma']} | {r['freq_Hz']:g} | {r['metric']} | {ps or 'n/a'} | {es or 'n/a'} | {gci:.3g} {r['unit']} |")
    s, t = PRODUCTION
    md.append(f"\n## Verdict on the production level ({s}, {t} steps/period)\n")
    for r in rows:
        if r["scheme"] != s or f"err_{t}" not in r:
            continue
        ok = abs(r[f"err_{t}"]) <= TOL[r["unit"]]
        md.append(f"- gamma {r['gamma']}, {r['freq_Hz']:g} Hz, {r['metric']}: error {r[f'err_{t}']:+.3g} {r['unit']} "
                  f"({'within' if ok else 'ABOVE'} target {TOL[r['unit']]} {r['unit']})")
    open(os.path.join(out_dir, "convergence_report.md"), "w").write("\n".join(md) + "\n")
    print("\n".join(md))


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1:])
