"""Shared figure style and visual encodings for every figure in this repository.

Rules (from the project owner's standing preferences, and the visualization guide of the
academic-research-skills repository, academic-paper/agents/visualization_agent.md):
  - Minimal: no top/right spines, no grid, no titles inside the figure (captions carry them),
    no chart junk.
  - Information is carried by the encodings below, and the SAME encoding means the SAME thing
    in every figure of the repository. Change it here, never in a single script.
  - Sans-serif; axis labels 9 pt; ticks, legend and annotations 8 pt (nothing below 8 pt);
    300 dpi; double column = 6.9 in wide.
  - Colour-blind safe palette (Paul Tol "bright"); colour is always paired with a marker shape.

ENCODINGS
  bead           colour + marker      PTFE blue circle, Al 1/32 orange square, Al 1/16 teal
                                      triangle, nylon red diamond
  radius used    line style + width   actual radius R: solid, 1.5 pt;  contact radius rc:
                                      dashed, 1.0 pt
  acceleration   viridis colour       0.1 ... 2.0 g; Gamma > 1 (model may be invalid) dotted
  data quality   marker fill          `good`: filled.  Any other label: hollow.
  domain echo    transparency         echo_ok == False (waves can reach the outer wall):
                                      alpha 0.45;  echo_ok == True: alpha 1
  frame          panel column         lab frame left, shaker frame right
  numerics       grey / black         scheme or step-size variants of one reference
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.lines import Line2D

DOUBLE_COLUMN_IN = 6.9
DPI = 300

BEAD_COLOR = {"ptfe": "#0077BB", "al32": "#EE7733", "al16": "#009988", "nylon": "#CC3311"}
BEAD_MARKER = {"ptfe": "o", "al32": "s", "al16": "^", "nylon": "D"}
BEAD_LABEL = {"ptfe": "PTFE", "al32": "Al 1/32\"", "al16": "Al 1/16\"", "nylon": "Nylon 1/16\""}
BASIS_LINESTYLE = {"R": "-", "rc": "--"}
BASIS_LINEWIDTH = {"R": 1.5, "rc": 1.0}
BASIS_LABEL = {"R": "radius R", "rc": "contact radius rc"}

_GAMMA_POS = {0.1: 0.00, 0.25: 0.20, 0.5: 0.40, 1.0: 0.62, 1.5: 0.80, 2.0: 0.95}
QUALITY_ORDER = ["good", "unconverged", "wave_cut", "subharmonic", "error"]
QUALITY_COLOR = {"good": "#009988", "unconverged": "#EE7733", "wave_cut": "#33BBEE",
                 "subharmonic": "#EE3377", "error": "#CC3311"}
QUALITY_MARKER = {"good": "o", "unconverged": "s", "wave_cut": "^", "subharmonic": "D", "error": "x"}

ALPHA_ECHO_OK = 1.0
ALPHA_ECHO_NOT_OK = 0.45


def apply():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "font.size": 8, "axes.labelsize": 9, "axes.titlesize": 9,
        "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
        "axes.spines.top": False, "axes.spines.right": False, "axes.grid": False,
        "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "figure.dpi": 100, "savefig.dpi": DPI, "savefig.bbox": "tight", "pdf.fonttype": 42,
    })


def gamma_color(g):
    key = min(_GAMMA_POS, key=lambda k: abs(k - g))
    return cm.viridis(_GAMMA_POS[key])


def draw_series(ax, f, y, quality, echo_ok, color, marker, linestyle="-", linewidth=1.2, markersize=3.6, label=None):
    """Line through all finite points; a segment is faded if either end is not echo_ok;
    each point is a filled marker when its quality is `good`, hollow otherwise."""
    pts = [(fi, yi, q, e) for fi, yi, q, e in zip(f, y, quality, echo_ok) if yi == yi]   # drop NaN
    for (f0, y0, _, e0), (f1, y1, _, e1) in zip(pts[:-1], pts[1:]):
        a = ALPHA_ECHO_OK if (e0 and e1) else ALPHA_ECHO_NOT_OK
        ax.plot([f0, f1], [y0, y1], color=color, ls=linestyle, lw=linewidth, alpha=a, solid_capstyle="round")
    for fi, yi, q, e in pts:
        a = ALPHA_ECHO_OK if e else ALPHA_ECHO_NOT_OK
        ax.plot([fi], [yi], ls="none", marker=marker, ms=markersize, mew=0.8, color=color, alpha=a,
                mfc=(color if q == "good" else "none"), mec=color)
    if label is not None:
        ax.plot([], [], color=color, ls=linestyle, lw=linewidth, marker=marker, ms=markersize,
                mfc=color, mec=color, label=label)


def key_quality_echo():
    """Legend handles for the two universal encodings: fill = quality, transparency = domain echo."""
    return [
        Line2D([], [], ls="none", marker="o", mfc="k", mec="k", ms=4, label="converged, no artefact"),
        Line2D([], [], ls="none", marker="o", mfc="none", mec="k", ms=4, label="flagged (see Note)"),
        Line2D([], [], color="k", lw=1.2, alpha=ALPHA_ECHO_NOT_OK, label="faded: waves can reach the wall"),
    ]


def save(fig, path_no_ext):
    fig.savefig(path_no_ext + ".png")
    fig.savefig(path_no_ext + ".pdf")
    plt.close(fig)
