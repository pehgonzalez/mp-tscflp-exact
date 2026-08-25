"""Regenerates fig_factorial, the component study of the Lagrangian engine.

Four arms of L-BBC on the 50-5 set, each measured against the BBC reference of
the same batch: the final dual bound of the arm minus the final dual bound of
BBC on the same instance, as a percentage of the best known solution. Positive
is a stronger bound than BBC.

The arms switch off one component each. Withholding the warm start, or
withholding the row that carries the Lagrangian bound, returns the engine to
BBC; withholding the incumbent handed to the solver costs nothing. So the gain
rests on the warm start and the bound row together, and the figure is arranged
to make that visible at a glance: the two arms that keep both components sit
away from zero, the two that lose one of them sit on it.

The five arms are taken from one batch, so that the comparison is not a
comparison across campaigns: every run plotted here, the reference included,
was recorded between 2026-08-05 and 2026-08-07 at the same 1800 s budget.

Usage: python make_fig_factorial.py [results_merged.csv] [bks_reference.csv] [outdir]
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
COL = {"MIP": "#E69F00", "BBC": "#0072B2", "L-BBC": "#009E73"}
DOT = "#4D4D4D"
SIZE = "50-5"
TL = 1800
# Arm label and its mode tag, ordered so that the two arms which keep the gain
# are adjacent at the top and the two which lose it are adjacent at the bottom.
ARMS = [("full", "exact"),
        ("start withheld", "exact-nostart"),
        ("bound row withheld", "exact-noldb"),
        ("cold start", "exact-nowarm")]
REF = ("exact", "1")   # BBC reference, same batch


def size_of(inst):
    p = inst.replace(".txt", "").split("-")
    return p[-2] + "-" + p[-1]


rows = list(csv.DictReader(open(RES)))
bks = {r["instance"]: float(r["bks"])
       for r in csv.DictReader(open(BKS))}


def canonical(mode, method):
    """Latest recording of every 50-5 instance for one arm at the 1800 s
    budget. Each arm was recorded once in the August batch, so the latest
    recording is that batch: on exact-noldb it discards an earlier August 1
    run, and on the BBC reference and the full arm it discards the July runs."""
    c = {}
    for r in rows:
        if r["seed"] != "0" or r["mode"] != mode or r["method"] != method:
            continue
        if int(float(r["time_limit_s"])) != TL:
            continue
        inst = r["instance"].replace(".txt", "")
        if size_of(inst) != SIZE:
            continue
        if inst not in c or r["datetime"] > c[inst]["datetime"]:
            c[inst] = r
    return c


ref = canonical(*REF)
D = np.full(len(ARMS), np.nan)
PTS = {}
BATCH = set()
for j, (label, mode) in enumerate(ARMS):
    arm = canonical(mode, "2")
    v = []
    for i in sorted(arm):
        if i not in ref or i not in bks:
            continue
        v.append((float(arm[i]["bound"]) - float(ref[i]["bound"])) / bks[i] * 100)
        BATCH.add(arm[i]["datetime"][:8])
        BATCH.add(ref[i]["datetime"][:8])
    PTS[j] = np.array(v)
    if v:
        D[j] = np.mean(v)

fig, ax = plt.subplots(figsize=(6.4, 3.0))
y = np.arange(len(ARMS))
h = 0.58
rng = np.random.default_rng(0)
# The level the full engine reaches, so that the arm which matches it is seen
# to match it and not merely to be long.
ax.axvline(D[0], color=COL["L-BBC"], linewidth=0.8, linestyle=(0, (1, 2.5)),
           alpha=0.75, zorder=1)
ax.barh(y, D, h, color=COL["L-BBC"], edgecolor="white", linewidth=0.6,
        zorder=2)
for j in range(len(ARMS)):
    g = PTS[j]
    jy = y[j] + (rng.random(len(g)) - 0.5) * (h * 0.70)
    ax.plot(g, jy, "o", ms=2.6, mfc=DOT, mec="white", mew=0.25, alpha=0.62,
            linestyle="none", zorder=4)
    lx = max(D[j], 0.0) + 0.07
    ax.text(lx, y[j], f"{D[j]:+.4f}", ha="left", va="center", fontsize=8.6,
            color="#333", zorder=6,
            bbox=dict(boxstyle="square,pad=0.12", fc="white", ec="none",
                      alpha=0.82))
ax.axvline(0, color="#444", linewidth=0.9, zorder=3)
ax.set_yticks(y)
ax.set_yticklabels([a for a, _ in ARMS], fontsize=10.5)
ax.set_ylim(len(ARMS) - 0.5, -0.5)
ax.set_xlim(-0.5, 3.35)
ax.set_xlabel("final dual bound minus BBC (% of BKS)")
ax.set_title(r"$50\times5$" + ", 1800 s, single batch", pad=6)
ax.xaxis.grid(True, color="#e9e9e9", linewidth=0.8)
ax.set_axisbelow(True)
ax.tick_params(length=0)
ax.spines["left"].set_visible(False)
fig.tight_layout()
fig.savefig(f"{OUT}/fig_factorial.pdf", bbox_inches="tight")
fig.savefig(f"{OUT}/fig_factorial.png", bbox_inches="tight", dpi=150)
plt.close(fig)

for j, (label, mode) in enumerate(ARMS):
    print(f"{label:20s} {mode:14s} n={len(PTS[j]):2d} mean={D[j]:+.4f}")
print("batch dates:", sorted(BATCH))
print(f"figure written to {OUT}/fig_factorial.pdf")
