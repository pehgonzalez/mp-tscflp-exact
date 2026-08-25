"""Regenerates fig_ranking, the synthesis figure of the paper.

Mean re-costed optimality gap by size set for the four engines that the paper
puts side by side: the direct solve, the direct solve strengthened by the
facility-linking family F_l, and the two decomposition engines. The point of
the figure is that the ranking of the four is not stable across the four size
sets, so no engine dominates and the strengthened direct solve is the only one
whose gap stays flat as the instances grow.

Each bar carries the twenty-five per-instance gaps of its cell as grey dots, so
that the reader sees the dispersion the mean is drawn from and not only the
mean. The y axis is truncated at a ceiling; a bar or a dot that runs past it is
marked with a triangle and continues off the panel.

Usage: python make_fig_ranking.py [results_merged.csv] [bks_reference.csv] [outdir]
Requires matplotlib and numpy.
"""
import csv
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RES = sys.argv[1] if len(sys.argv) > 1 else "results_merged.csv"
BKS = sys.argv[2] if len(sys.argv) > 2 else "bks_reference.csv"
OUT = (sys.argv[3] if len(sys.argv) > 3 else "figures").rstrip("/")

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Nimbus Roman", "Times New Roman", "DejaVu Serif"],
    "mathtext.fontset": "stix", "font.size": 12, "axes.titlesize": 11.5,
    "axes.labelsize": 12, "xtick.labelsize": 11, "ytick.labelsize": 11,
    "legend.fontsize": 11, "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.9, "axes.edgecolor": "#444", "figure.dpi": 170})
COL = {"MIP": "#E69F00", "MIP+F": "#E69F00", "BBC": "#0072B2",
       "L-BBC": "#009E73"}
DOT = "#4D4D4D"
NAME = {"MIP": "MIP (direct solve)", "MIP+F": r"MIP + $F_\ell$",
        "BBC": "BBC", "L-BBC": "L-BBC"}
# The strengthened direct solve is the same engine as the direct solve with one
# extra family of constraints, so it is drawn in the same colour and set apart
# by the hatch, the convention fig_inversion already uses for a kinship.
HATCH = {"MIP+F": "///"}
SIZES = ["50-5", "50-10", "100-5", "100-10"]
TL = 1800
SLAB = [r"$50\times5$" + "\n1800 s", r"$50\times10$" + "\n1800 s",
        r"$100\times5$" + "\n1800 s", r"$100\times10$" + "\n1800 s"]
# Engine tag, admitted mode tags, method column. Both mode tags of the exact
# campaign are admitted for the decomposition engines, as in make_figures.py,
# so that a cell rerun under exact-pap supersedes its earlier recording at the
# same budget; on this data the two engines select the same runs either way.
ENGINES = [("MIP", ("exact",), "0"),
           ("MIP+F", ("exact-fl",), "0"),
           ("BBC", ("exact", "exact-pap"), "1"),
           ("L-BBC", ("exact", "exact-pap"), "2")]
CEIL = 5.0


def size_of(inst):
    p = inst.replace(".txt", "").split("-")
    return p[-2] + "-" + p[-1]


rows = list(csv.DictReader(open(RES)))


def canonical(modes, method):
    """Canonical run of every instance for one engine: seed 0, the uniform
    1800 s budget, and the latest recording if the same cell was run more than
    once. On the direct solve this keeps the July batch on 50-5, 100-5 and
    100-10, where nothing was rerun, and the August batch on 50-10, whose
    earlier recording was taken at the retired 900 s budget."""
    c = {}
    for r in rows:
        if r["seed"] != "0" or r["mode"] not in modes or r["method"] != method:
            continue
        if int(float(r["time_limit_s"])) != TL:
            continue
        inst = r["instance"].replace(".txt", "")
        if size_of(inst) not in SIZES:
            continue
        if inst not in c or r["datetime"] > c[inst]["datetime"]:
            c[inst] = r
    return c


def UB(r):
    """Re-costed upper bound of a run, the same definition make_tables.py and
    make_figures.py use: the verified re-costing when it succeeded, otherwise
    the better of the solver incumbent and the heuristic."""
    if r["verified_ok"] == "1" and float(r["verified_cost"]) > 0:
        return float(r["verified_cost"])
    o = 0.0 if r["status"] == "NOTFOUND" else float(r["obj"])
    h = float(r["heuristic_cost"])
    return o if o > 0 else (h if h > 0 else np.nan)


def gap(r):
    u = UB(r)
    b = float(r["bound"])
    return max(0.0, (u - b) / u * 100) if u > 0 and b > 0 else np.nan


G = np.full((4, len(ENGINES)), np.nan)     # mean gap
NC = np.zeros((4, len(ENGINES)), dtype=int)  # certificates
PTS = {}                                    # per-instance gaps
for j, (tag, modes, method) in enumerate(ENGINES):
    c = canonical(modes, method)
    for a, sz in enumerate(SIZES):
        sel = [i for i in sorted(c) if size_of(i) == sz]
        g = [x for x in (gap(c[i]) for i in sel) if x == x]
        PTS[(a, j)] = g
        if g:
            G[a, j] = np.mean(g)
        NC[a, j] = sum(c[i]["status"] == "OPTIMAL" for i in sel)

fig, ax = plt.subplots(figsize=(6.9, 3.75))
x = np.arange(4)
w = 0.2
rng = np.random.default_rng(0)
for j, (tag, modes, method) in enumerate(ENGINES):
    xc = x + (j - 1.5) * w
    bars = ax.bar(xc, np.minimum(G[:, j], CEIL), w, color=COL[tag],
                  edgecolor="white", hatch=HATCH.get(tag), linewidth=0.6,
                  label=NAME[tag], zorder=2)
    for a in range(4):
        v = G[a, j]
        if v != v:
            continue
        # The per-instance gaps of the cell, jittered inside the bar width.
        g = np.array(PTS[(a, j)])
        jx = xc[a] + (rng.random(len(g)) - 0.5) * (w * 0.62)
        inside = g <= CEIL
        ax.plot(jx[inside], g[inside], "o", ms=2.5, mfc=DOT, mec="white",
                mew=0.25, alpha=0.62, zorder=4, linestyle="none")
        ax.plot(jx[~inside], np.full((~inside).sum(), CEIL * 0.985), "^",
                ms=3.0, mfc=DOT, mec="none", alpha=0.62, zorder=4,
                linestyle="none")
        # A bar that runs past the ceiling is marked and labelled outside.
        over = v > CEIL
        if over:
            ax.plot([xc[a]], [CEIL + 0.07], "^", ms=6.0, mfc=COL[tag],
                    mec="white", mew=0.5, clip_on=False, zorder=5)
            ty, va = CEIL + 0.26, "bottom"
        else:
            ty, va = v + 0.07, "bottom"
        # The dots of the cell reach up to the label, so it is set on a patch
        # of background rather than left to overlap them.
        ax.text(xc[a], ty, f"{v:.2f}", ha="center", va=va, fontsize=8.4,
                color="#333", clip_on=False, zorder=6,
                bbox=dict(boxstyle="square,pad=0.10", fc="white", ec="none",
                          alpha=0.82))
        if NC[a, j]:
            ax.text(xc[a], ty + 0.36, f"{NC[a, j]}\nclosed", ha="center",
                    va="bottom", fontsize=5.8, color="#333", linespacing=0.95,
                    clip_on=False, zorder=6,
                    bbox=dict(boxstyle="square,pad=0.10", fc="white",
                              ec="none", alpha=0.82))
ax.set_xticks(x)
ax.set_xticklabels(SLAB, fontsize=10)
ax.set_xlabel(r"Plants $\times$ products")
ax.set_ylabel("optimality gap (%)")
ax.set_ylim(0, CEIL)
ax.set_yticks([0, 1, 2, 3, 4, 5])
ax.set_xlim(-0.5, 3.5)
ax.yaxis.grid(True, color="#e9e9e9", linewidth=0.8)
ax.set_axisbelow(True)
ax.tick_params(length=0)
h, l = ax.get_legend_handles_labels()
fig.legend(h, l, ncol=4, loc="lower center", bbox_to_anchor=(0.55, -0.035),
           frameon=False, handlelength=1.1, columnspacing=1.3, fontsize=10.5)
fig.tight_layout(rect=[0, 0.02, 1, 1])
fig.savefig(f"{OUT}/fig_ranking.pdf", bbox_inches="tight")
fig.savefig(f"{OUT}/fig_ranking.png", bbox_inches="tight", dpi=150)
plt.close(fig)

for a, sz in enumerate(SIZES):
    print(sz, " ".join(f"{ENGINES[j][0]}={G[a, j]:.4f}({NC[a, j]})"
                       for j in range(len(ENGINES))))
print(f"figure written to {OUT}/fig_ranking.pdf")
