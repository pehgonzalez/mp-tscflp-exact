"""Regenerates the three figures of the paper from results/results.csv.

Figure 1 (fig_difficulty): mean optimality gap by generator class and size,
    one panel per engine, with the number of instances each engine closed at
    the uniform budget printed in the cell.
Figure 2 (fig_inversion): which engine attains the best incumbent and the best
    dual bound, by size, counting a tie between the two decomposition engines
    as a decomposition win and only a tie involving the direct solve as a tie.
Figure 3 (fig_gap_scale): mean optimality gap by engine and size.

Every gap plotted here is the re-costed gap, the same quantity the tables and
the text report. The wall-clock budget is printed on the axis of each figure so
that the reader sees the budget the gaps were measured under; it is the same
1800 s on all four size sets.

Usage: python make_figures.py [results.csv] [bks_reference.csv] [outdir]
Requires matplotlib and numpy.
"""
import csv
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm

RES = sys.argv[1] if len(sys.argv) > 1 else "results/results.csv"
BKS = sys.argv[2] if len(sys.argv) > 2 else "data/bks_reference.csv"
OUT = (sys.argv[3] if len(sys.argv) > 3 else "figures").rstrip("/")

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Nimbus Roman", "Times New Roman", "DejaVu Serif"],
    "mathtext.fontset": "stix", "font.size": 12, "axes.titlesize": 11.5,
    "axes.labelsize": 12, "xtick.labelsize": 11, "ytick.labelsize": 11,
    "legend.fontsize": 11, "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.9, "axes.edgecolor": "#444", "figure.dpi": 170})
COL = {"MIP": "#E69F00", "BBC": "#0072B2", "L-BBC": "#009E73"}
TIE = "#BFBFBF"
NAME = {"MIP": "MIP (direct solve)", "BBC": "BBC", "L-BBC": "L-BBC"}
M2N = {"0": "MIP", "1": "BBC", "2": "L-BBC"}
TL = {"50-5": 1800, "50-10": 1800, "100-5": 1800, "100-10": 1800}
SIZES = ["50-5", "50-10", "100-5", "100-10"]
SLAB = [r"$50\times5$" + "\n1800 s", r"$50\times10$" + "\n1800 s",
        r"$100\times5$" + "\n1800 s", r"$100\times10$" + "\n1800 s"]
CLS = ["C1", "C2", "C3", "C4", "C5"]


def size_of(inst):
    p = inst.replace(".txt", "").split("-")
    return p[-2] + "-" + p[-1]


# Canonical run of an instance and engine: seed 0, the budget of its size set,
# and the latest recording if the same cell was run more than once. Both mode
# tags of the exact campaign are admitted, so that a cell rerun with core-point
# cuts under the tag exact-pap supersedes its earlier recording at the same
# budget rather than being ignored.
rows = list(csv.DictReader(open(RES)))
canon = {}
for r in rows:
    if r["seed"] != "0" or r["mode"] not in ("exact", "exact-pap"):
        continue
    inst = r["instance"].replace(".txt", "")
    sz = size_of(inst)
    if sz not in TL or int(float(r["time_limit_s"])) != TL[sz]:
        continue
    k = (inst, r["method"])
    if k not in canon or r["datetime"] > canon[k]["datetime"]:
        canon[k] = r


def UB(r):
    """Re-costed upper bound of a run: the facilities it opened are fixed and
    every product is re-routed through them by the solver-free min-cost-flow
    evaluator. For the two decomposition engines this reproduces the solver
    objective exactly; for the direct solve it recovers the cost that the
    incumbent's own configuration attains. Same definition as make_tables.py."""
    if r["verified_ok"] == "1" and float(r["verified_cost"]) > 0:
        return float(r["verified_cost"])
    o = 0.0 if r["status"] == "NOTFOUND" else float(r["obj"])
    h = float(r["heuristic_cost"])
    return o if o > 0 else (h if h > 0 else np.nan)


def gap(r):
    u = UB(r)
    b = float(r["bound"])
    return max(0.0, (u - b) / u * 100) if u > 0 and b > 0 else np.nan


insts = sorted(set(i for (i, m) in canon))

# Figure 1: difficulty map, one panel per engine, with certificates per cell.
# The direct solve is shown beside the two decomposition engines rather than
# folded into a best of three, because a virtual best would present the
# baseline's gaps as the reach of the method on the sets where the baseline
# wins.
H = np.full((3, 4, 5), np.nan)
NC = np.zeros((3, 4, 5), dtype=int)
for j, m in enumerate("012"):
    for a, sz in enumerate(SIZES):
        for c, cl in enumerate(CLS):
            sel = [i for i in insts if size_of(i) == sz
                   and i.split("-")[1] == cl and (i, m) in canon]
            g = [x for x in (gap(canon[(i, m)]) for i in sel) if x == x]
            if g:
                H[j, a, c] = np.mean(g)
            NC[j, a, c] = sum(canon[(i, m)]["status"] == "OPTIMAL" for i in sel)

norm = PowerNorm(gamma=0.4, vmin=0.0, vmax=float(np.nanmax(H)))
fig, axes = plt.subplots(1, 3, figsize=(7.3, 2.95), sharey=True)
for j, ax in enumerate(axes):
    im = ax.imshow(H[j], cmap="YlOrRd", norm=norm, aspect="auto")
    ax.set_xticks(range(5))
    ax.set_xticklabels(CLS, fontsize=10)
    ax.set_title(NAME[M2N[str(j)]], pad=5, fontsize=11)
    ax.set_xlabel("Generator class", fontsize=11)
    for a in range(4):
        for c in range(5):
            v = H[j, a, c]
            if v != v:
                continue
            col = "white" if norm(v) > 0.62 else "#1a1a1a"
            dy = -0.15 if NC[j, a, c] else 0.0
            ax.text(c, a + dy, f"{v:.1f}", ha="center", va="center",
                    color=col, fontsize=8.6)
            if NC[j, a, c]:
                # Rodada 9. O rotulo "N closed" era mais largo que a celula e
                # transbordava para fora do painel na primeira coluna, num
                # corpo ilegivel. Vira uma contagem entre parenteses, que cabe
                # na celula e e definida na legenda da figura.
                ax.text(c, a + 0.26, f"({NC[j, a, c]})", ha="center",
                        va="center", color=col, fontsize=7.8)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.tick_params(length=0)
axes[0].set_yticks(range(4))
axes[0].set_yticklabels(SLAB, fontsize=9.5)
axes[0].set_ylabel(r"Plants $\times$ products", fontsize=11)
cb = fig.colorbar(im, ax=axes, fraction=0.028, pad=0.02,
                  ticks=[0, 0.5, 1, 2, 5, 10, 20])
cb.set_label("mean optimality gap (%)", fontsize=10)
cb.ax.set_yticklabels(["0", "0.5", "1", "2", "5", "10", "20"], fontsize=9)
cb.ax.tick_params(length=2)
cb.outline.set_linewidth(0.6)
fig.savefig(f"{OUT}/fig_difficulty.pdf", bbox_inches="tight")
fig.savefig(f"{OUT}/fig_difficulty.png", bbox_inches="tight", dpi=150)
plt.close(fig)


# Figure 2: who attains the best incumbent and the best dual bound, by size.
# An instance on which the two decomposition engines tie above the direct solve
# counts as a decomposition win, since the comparison the text draws from this
# figure is decomposition against direct solve. Only a tie that involves the
# direct solve is a tie. The tolerance is half a unit, which on integral data
# is exact equality.
KEYS = ["MIP", "BBC", "L-BBC", "both", "tie"]
FILL = {"MIP": COL["MIP"], "BBC": COL["BBC"], "L-BBC": COL["L-BBC"],
        "both": COL["BBC"], "tie": TIE}
HATCH = {"both": "////"}
LEG = {"MIP": NAME["MIP"], "BBC": "BBC", "L-BBC": "L-BBC",
       "both": "BBC and L-BBC tied", "tie": "tie with MIP"}


def shares(metric):
    Mx = np.zeros((4, len(KEYS)))
    for a, sz in enumerate(SIZES):
        for i in [x for x in insts if size_of(x) == sz]:
            v = {m: (UB(canon[(i, m)]) if metric == "ub"
                     else float(canon[(i, m)]["bound"]))
                 for m in "012" if (i, m) in canon}
            if not v:
                continue
            best = min(v.values()) if metric == "ub" else max(v.values())
            win = [m for m in "012" if m in v and abs(v[m] - best) < 0.5]
            if "0" in win and len(win) > 1:
                k = "tie"
            elif len(win) == 1:
                k = M2N[win[0]]
            else:
                k = "both"
            Mx[a, KEYS.index(k)] += 1
    return Mx


fig, axes = plt.subplots(1, 2, figsize=(6.9, 3.6), sharey=True)
for ax, metric, ttl in zip(axes, ["ub", "bound"],
                           ["Best incumbent", "Best dual bound"]):
    Mx = shares(metric)
    P = Mx / Mx.sum(1, keepdims=True) * 100
    bottom = np.zeros(4)
    y = np.arange(4)
    for idx, k in enumerate(KEYS):
        ax.barh(y, P[:, idx], left=bottom, color=FILL[k], edgecolor="white",
                hatch=HATCH.get(k), linewidth=0.8, label=LEG[k], height=0.72)
        for a in range(4):
            cnt = int(round(Mx[a, idx]))
            if cnt > 0:
                ax.text(bottom[a] + P[a, idx] / 2, a, str(cnt), ha="center",
                        va="center", fontsize=7.6, fontweight="bold",
                        color=("#333" if k == "tie" else "white"))
        bottom += P[:, idx]
    ax.set_xlim(0, 100)
    ax.set_xlabel("share of instances (%)", fontsize=11)
    ax.set_title(ttl, pad=6, fontsize=11.5)
    ax.tick_params(length=0)
    ax.spines["left"].set_visible(False)
axes[0].set_yticks(range(4))
axes[0].set_yticklabels(SLAB, fontsize=9.5)
axes[0].set_ylabel(r"Plants $\times$ products", fontsize=11)
axes[0].invert_yaxis()
h, l = axes[0].get_legend_handles_labels()
fig.legend(h, l, ncol=3, loc="lower center", bbox_to_anchor=(0.5, -0.10),
           frameon=False, handlelength=1.1, columnspacing=1.4, fontsize=10)
fig.tight_layout(rect=[0, 0.055, 1, 1])
fig.savefig(f"{OUT}/fig_inversion.pdf", bbox_inches="tight")
fig.savefig(f"{OUT}/fig_inversion.png", bbox_inches="tight", dpi=150)
plt.close(fig)

# Figure 3: mean gap by engine and size.
G = np.full((4, 3), np.nan)
for a, sz in enumerate(SIZES):
    for j, m in enumerate("012"):
        g = [gap(canon[(i, m)]) for i in insts
             if size_of(i) == sz and (i, m) in canon]
        g = [x for x in g if x == x]
        if g:
            G[a, j] = np.mean(g)
fig, ax = plt.subplots(figsize=(5.9, 3.3))
x = np.arange(4)
w = 0.26
for j, m in enumerate("012"):
    nm = M2N[m]
    bars = ax.bar(x + (j - 1) * w, G[:, j], w, color=COL[nm], edgecolor="white",
                  linewidth=0.6, label=NAME[nm])
    for b, v in zip(bars, G[:, j]):
        if v == v:
            ax.text(b.get_x() + b.get_width() / 2, v + 0.13, f"{v:.1f}",
                    ha="center", va="bottom", fontsize=9, color="#333")
ax.set_xticks(x)
ax.set_xticklabels(SLAB, fontsize=10)
ax.set_xlabel(r"Plants $\times$ products")
ax.set_ylabel("mean optimality gap (%)")
ax.set_ylim(0, float(np.nanmax(G)) * 1.13)
ax.yaxis.grid(True, color="#e9e9e9", linewidth=0.8)
ax.set_axisbelow(True)
ax.tick_params(length=0)
h, l = ax.get_legend_handles_labels()
fig.legend(h, l, ncol=3, loc="lower center", bbox_to_anchor=(0.55, -0.035),
           frameon=False, handlelength=1.1, columnspacing=1.5, fontsize=10.5)
fig.tight_layout(rect=[0, 0.015, 1, 1])
fig.savefig(f"{OUT}/fig_gap_scale.pdf", bbox_inches="tight")
fig.savefig(f"{OUT}/fig_gap_scale.png", bbox_inches="tight", dpi=150)
plt.close(fig)
print(f"figures written to {OUT}/")
