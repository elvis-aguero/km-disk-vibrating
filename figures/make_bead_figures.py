#!/usr/bin/env python3
"""Figures for the bead campaign (see figstyle.py for the shared encodings).

Usage: python3 figures/make_bead_figures.py <campaign dir>
Reads  <campaign>/{03_accel,03b_accel_gt1,04_sweep_R,05_sweep_rc,06_numerics/*}/gain_phase_all.csv
Writes <campaign>/figures/fig{1..4}_*.{png,pdf} and captions.md
Needs matplotlib and numpy (on Oscar: `module load py-matplotlib/3.10.3-2hcg`).
"""
import csv
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import figstyle as fs  # noqa: E402

fs.apply()
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

NUMERIC = ["gamma", "freq_Hz", "A_ref_cm", "diskRadius_cm", "periods_total", "periods_used", "echo_pred",
           "amp_rel_change", "phase_change_deg", "gain_lab", "phase_lag_lab_deg", "gain_shaker",
           "phase_lag_shaker_deg", "max_rel_motion_over_a", "subharmonic_ratio"]
BEADS = ["ptfe", "al32", "al16", "nylon"]


def load(path):
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    for r in rows:
        for k in NUMERIC:
            r[k] = float(r[k])
        r["echo_ok"] = r["echo_ok"] == "true"
    return rows


def pick(rows, bead=None, gamma=None, basis=None):
    out = [r for r in rows
           if (bead is None or r["bead"] == bead)
           and (gamma is None or abs(r["gamma"] - gamma) < 1e-9)
           and (basis is None or r["basis"] == basis)]
    return sorted(out, key=lambda r: r["freq_Hz"])


def cols(rows, key):
    return [r[key] for r in rows]


def wrap180(d):
    return (d + 180.0) % 360.0 - 180.0


def ycol(rows, key):
    """Column for plotting. Phase lags are wrapped to (-180, 180] so that 359 deg shows as -1 deg."""
    if key.startswith("phase"):
        return [wrap180(r[key]) for r in rows]
    return [r[key] for r in rows]


def counts(rows):
    c = {}
    for r in rows:
        c[r["quality"]] = c.get(r["quality"], 0) + 1
    return c


def fmt_counts(c):
    return ", ".join(f"{k} {c[k]}" for k in fs.QUALITY_ORDER if k in c)


def frame_axes(figsize):
    fig, axes = plt.subplots(2, 2, figsize=figsize, sharex=True)
    for j, t in enumerate(["Lab frame", "Shaker frame"]):
        axes[0][j].set_title(t, fontweight="bold", fontsize=9)
    axes[0][0].set_ylabel("gain (bead / bath amplitude)")
    axes[1][0].set_ylabel("phase lag behind bath (deg)")
    for j in range(2):
        axes[1][j].set_xlabel("bath frequency (Hz)")
        axes[0][j].set_xlim(0, 95)
    axes[0][0].axhline(1.0, color="0.65", lw=0.6, ls=":", zorder=0)   # gain 1 = bead moves with the bath
    return fig, axes


def metric_keys(frame_col, row):
    return [["gain_lab", "gain_shaker"], ["phase_lag_lab_deg", "phase_lag_shaker_deg"]][row][frame_col]


def fig1_gamma(acc, hi, out):
    rows = acc + hi
    gammas = sorted({r["gamma"] for r in rows})
    fig, axes = frame_axes((fs.DOUBLE_COLUMN_IN, 4.6))
    for g in gammas:
        s = pick(rows, bead="ptfe", gamma=g)
        ls, lw = ("-", 1.3) if g <= 1.0 else (":", 1.0)
        for i in range(2):
            for j in range(2):
                fs.draw_series(axes[i][j], cols(s, "freq_Hz"), ycol(s, metric_keys(j, i)), cols(s, "quality"),
                               cols(s, "echo_ok"), fs.gamma_color(g), "o", linestyle=ls, linewidth=lw, markersize=2.8,
                               label=(f"Γ = {g:g}" if (i, j) == (0, 0) else None))
    h, l = axes[0][0].get_legend_handles_labels()
    axes[0][1].legend(h, l, frameon=False, loc="lower right", ncol=2, handlelength=1.6)   # lower right of the shaker-gain panel: no data there
    fig.legend(handles=fs.key_quality_echo(), loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.03))
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fs.save(fig, os.path.join(out, "fig1_gamma_independence"))


def fig2_beads(R, rc, out):
    fig, axes = frame_axes((fs.DOUBLE_COLUMN_IN, 4.9))
    for basis, rows in (("R", R), ("rc", rc)):
        for b in BEADS:
            s = pick(rows, bead=b, gamma=0.5)
            for i in range(2):
                for j in range(2):
                    fs.draw_series(axes[i][j], cols(s, "freq_Hz"), ycol(s, metric_keys(j, i)), cols(s, "quality"),
                                   cols(s, "echo_ok"), fs.BEAD_COLOR[b], fs.BEAD_MARKER[b],
                                   linestyle=fs.BASIS_LINESTYLE[basis], linewidth=fs.BASIS_LINEWIDTH[basis], markersize=2.3)
    handles = [Line2D([], [], color=fs.BEAD_COLOR[b], marker=fs.BEAD_MARKER[b], ms=4, lw=1.3, label=fs.BEAD_LABEL[b])
               for b in BEADS]
    handles += [Line2D([], [], color="k", lw=fs.BASIS_LINEWIDTH[k], ls=fs.BASIS_LINESTYLE[k], label=fs.BASIS_LABEL[k])
                for k in ("R", "rc")]
    handles += fs.key_quality_echo()
    fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False, bbox_to_anchor=(0.5, -0.06))
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fs.save(fig, os.path.join(out, "fig2_beads_R_vs_rc"))


def fig3_numerics(acc, eul, b30, out):
    ref = {r["freq_Hz"]: r for r in pick(acc, bead="ptfe", gamma=0.5)}
    fig, axes = plt.subplots(1, 2, figsize=(fs.DOUBLE_COLUMN_IN, 2.5))
    spec = [(eul, "Euler, 60 steps/period", "#000000", "s"), (b30, "BDF2, 30 steps/period", "#888888", "o")]
    excluded = {}
    for rows, label, color, marker in spec:
        f, dg, dp = [], [], []
        for r in pick(rows, bead="ptfe", gamma=0.5):
            ref_r = ref.get(r["freq_Hz"])
            # Only compare cases that are good in BOTH runs: a flagged case (for example a grown
            # f/2 mode) would show the artefact, not the time-stepping error.
            if ref_r is None or r["quality"] != "good" or ref_r["quality"] != "good":
                excluded.setdefault(label, []).append(f'{r["freq_Hz"]:g} Hz')
                continue
            f.append(r["freq_Hz"])
            dg.append(100.0 * (r["gain_lab"] - ref_r["gain_lab"]) / ref_r["gain_lab"])
            dp.append(wrap180(r["phase_lag_lab_deg"] - ref_r["phase_lag_lab_deg"]))
        axes[0].plot(f, dg, color=color, marker=marker, ms=3.6, lw=1.0, label=label)
        axes[1].plot(f, dp, color=color, marker=marker, ms=3.6, lw=1.0)
    for ax, tol, yl in ((axes[0], 2.0, "Δ gain vs reference (%)"), (axes[1], 2.0, "Δ phase lag vs reference (deg)")):
        ax.axhline(0, color="0.65", lw=0.6)
        for s in (-tol, tol):
            ax.axhline(s, color="0.65", lw=0.6, ls=":")
        ax.set_xlabel("bath frequency (Hz)")
        ax.set_ylabel(yl)
        ax.set_xlim(0, 95)
    axes[0].legend(frameon=False, loc="best")
    fig.tight_layout()
    fs.save(fig, os.path.join(out, "fig3_numerics"))
    return excluded


def fig3_quality(acc, hi, R, rc, out):
    freqs = sorted({r["freq_Hz"] for r in R})
    xi = {f: i for i, f in enumerate(freqs)}
    lines = []
    for g in sorted({r["gamma"] for r in acc + hi}):
        lines.append((f"PTFE, R, Γ = {g:g} (acceleration test)", pick(acc + hi, bead="ptfe", gamma=g)))
    for basis, rows in (("R", R), ("rc", rc)):
        for b in BEADS:
            lines.append((f"{fs.BEAD_LABEL[b]}, {basis}, Γ = 0.5", pick(rows, bead=b, gamma=0.5)))
    fig, ax = plt.subplots(figsize=(fs.DOUBLE_COLUMN_IN, 4.9))
    for yi, (label, rows) in enumerate(lines):
        for r in rows:
            q = r["quality"]
            if r["freq_Hz"] not in xi:
                continue
            ax.plot([xi[r["freq_Hz"]]], [yi], ls="none", marker=fs.QUALITY_MARKER[q], ms=4.2, mew=0.9,
                    color=fs.QUALITY_COLOR[q], mfc=(fs.QUALITY_COLOR[q] if q == "good" else "none") if q != "error" else None,
                    alpha=(fs.ALPHA_ECHO_OK if r["echo_ok"] else fs.ALPHA_ECHO_NOT_OK))
    ax.set_yticks(range(len(lines)))
    ax.set_yticklabels([l for l, _ in lines])
    ax.invert_yaxis()
    ticks = [f for f in freqs if f in (5, 10, 20, 30, 40, 50, 60, 70, 80, 90)]
    ax.set_xticks([xi[f] for f in ticks])
    ax.set_xticklabels([f"{f:g}" for f in ticks])
    ax.set_xlabel("bath frequency (Hz); 1 Hz steps between 50 and 70 Hz")
    handles = [Line2D([], [], ls="none", marker=fs.QUALITY_MARKER[q], color=fs.QUALITY_COLOR[q],
                      mfc=(fs.QUALITY_COLOR[q] if q == "good" else "none"), ms=4.5, label=q)
               for q in fs.QUALITY_ORDER if q in {r["quality"] for _, rows in lines for r in rows}]   # only labels that occur
    handles.append(Line2D([], [], color="k", lw=1.2, alpha=fs.ALPHA_ECHO_NOT_OK, label="faded: waves can reach the wall"))
    fig.legend(handles=handles, loc="lower center", ncol=6, frameon=False, bbox_to_anchor=(0.5, -0.04))
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fs.save(fig, os.path.join(out, "fig3_quality_map"))


def main():
    camp = sys.argv[1]
    out = os.path.join(camp, "figures")
    os.makedirs(out, exist_ok=True)
    p = lambda *a: os.path.join(camp, *a, "gain_phase_all.csv")  # noqa: E731
    acc, hi = load(p("03_accel")), load(p("03b_accel_gt1"))
    R, rc = load(p("04_sweep_R")), load(p("05_sweep_rc"))

    fig1_gamma(acc, hi, out)
    fig2_beads(R, rc, out)
    fig3_quality(acc, hi, R, rc, out)

    cap = f"""# Captions (campaign {os.path.basename(camp.rstrip('/'))})

**Figure 1**
*Gain and phase lag of the PTFE bead (radius R) do not change with the bath acceleration.*
Note. Bath-driven simulation, 5-90 Hz, fluid 7.9 cSt, BDF2, 60 steps per period. Colour: acceleration
Γ = 0.1 to 2.0 (viridis); dotted lines: Γ > 1, where the effective gravity changes sign and the model
may not be valid. Filled markers: converged without artefact; hollow: flagged (unconverged, wall echo,
subharmonic growth); faded: waves can reach the outer wall. Cases that diverged have no point.
Phase lag is wrapped to (-180, 180] degrees.
Quality counts, Γ ≤ 1: {fmt_counts(counts(acc))}; Γ > 1: {fmt_counts(counts(hi))}. The gain axis does not start at zero.

**Figure 2**
*Gain and phase lag of four beads at Γ = 0.5, with the actual radius R (solid) and the contact radius rc (dashed).*
Note. Same encodings as Figure 1; colour and marker: bead. Quality counts, R: {fmt_counts(counts(pick(R, gamma=0.5)))};
rc: {fmt_counts(counts(pick(rc, gamma=0.5)))}. The gain axis does not start at zero.

**Figure 3**
*Quality of every simulated case.*
Note. One row per sweep, one marker per frequency. Marker shape and colour give the quality label
(good, unconverged, wave_cut, subharmonic = the f/2 component exceeds 10 % of the f component, error =
numerical divergence); faded markers: waves can reach the outer wall.
"""
    with open(os.path.join(out, "captions.md"), "w") as f:
        f.write(cap)
    print("wrote", out)
    print(os.popen(f"ls {out}").read())


if __name__ == "__main__":
    main()
