"""Emits the four per-instance sidewaystable bodies of the paper.

Canonical run selection: mode 'exact', seed 0, a 1,800 s budget on all four
size sets and the most recent row per (instance, method), so a rerun of a set
replaces the earlier campaign.
Extended-budget rows (14,400 s) are excluded from these tables, which report a
uniform budget; they feed Tables 2 and 3 only.

Each engine occupies four columns, the best upper bound, the absolute dual
bound, the gap in percent and the wall-clock time, so a table is fourteen
columns wide with the instance name and the best-known value.

Usage: python make_full_tables.py [results.csv] [bks_reference.csv]
"""
import csv
import sys

RES = sys.argv[1] if len(sys.argv) > 1 else "results/results.csv"
BKSF = sys.argv[2] if len(sys.argv) > 2 else "data/bks_reference.csv"
TL = {"50-5": "1800", "50-10": "1800", "100-5": "1800", "100-10": "1800"}
SETS = [("50-5", "tab:full55", "five-product set", "1{,}800"),
        ("50-10", "tab:full510", "ten-product set", "1{,}800"),
        ("100-5", "tab:full1005", "100-plant five-product set", "1{,}800"),
        ("100-10", "tab:full10010", "100-plant ten-product set", "1{,}800")]

# Four columns per engine (upper bound, dual bound, gap, time) plus the
# instance name and the best-known value: fourteen columns in all.
HEADER = [
    "\\begin{tabular}{l rrrr rrrr rrrr r}",
    "\\toprule",
    "& \\multicolumn{4}{c}{\\textsc{MIP}} & \\multicolumn{4}{c}{\\textsc{BBC}}"
    " & \\multicolumn{4}{c}{\\textsc{L-BBC}} & \\\\",
    "\\cmidrule(lr){2-5}\\cmidrule(lr){6-9}\\cmidrule(lr){10-13}",
    "Instance & UB & LB & Gap\\% & Time & UB & LB & Gap\\% & Time"
    " & UB & LB & Gap\\% & Time & BKS \\\\",
    "\\midrule",
]


def size_of(inst):
    p = inst.split("-")
    return p[-2] + "-" + p[-1]


def num(x):
    s = f"{round(x):,}"
    return s.replace(",", "{,}")


def load():
    rows = list(csv.DictReader(open(RES)))
    bks = {r["instance"].replace(".txt", ""): float(r["bks"])
           for r in csv.DictReader(open(BKSF))}
    canon = {}
    for r in rows:
        inst = r["instance"].replace(".txt", "")
        sz = size_of(inst)
        if r["mode"] != "exact" or r["seed"] != "0":
            continue
        if sz not in TL or r["time_limit_s"] != TL[sz]:
            continue
        k = (inst, r["method"])
        if k not in canon or r["datetime"] > canon[k]["datetime"]:
            canon[k] = r
    return rows, bks, canon


def cell(r):
    """Upper bound, dual bound, gap in percent and wall-clock seconds of a run.

    The upper bound is the re-costed value of the reported solution, obtained by
    re-routing every product optimally through the facilities the run opened.
    For the two decomposition engines this reproduces the solver objective
    exactly, because their subproblem already routes optimally; for the direct
    solve it recovers the cost the incumbent's own facility configuration
    attains, which the solver may leave above the optimum of that configuration.
    The gap is recomputed against the same value, so column and gap agree.
    The dual bound is the run's own global lower bound, reported absolutely and
    unrounded here; the caller rounds it to the unit like the upper bound.
    """
    if r is None:
        return None
    if r["verified_ok"] == "1" and float(r["verified_cost"]) > 0:
        u = float(r["verified_cost"])
    else:
        o, h = float(r["obj"]), float(r["heuristic_cost"])
        u = o if (r["status"] != "NOTFOUND" and o > 0) else (h if h > 0 else None)
    b = float(r["bound"])
    if u is None:
        return None, (b if b > 0 else None), None, min(float(r["total_wall_s"]), float(r["time_limit_s"]))
    g = 100.0 * (u - b) / u if b > 0 else 100.0 * float(r["gap"])
    return u, (b if b > 0 else None), max(0.0, g), min(float(r["total_wall_s"]), float(r["time_limit_s"]))


def main():
    rows, bks, canon = load()
    for sz, label, name, budget in SETS:
        insts = sorted({i for i, _ in canon if size_of(i) == sz},
                       key=lambda i: (i.split("-")[1], i.split("-")[0]))
        print(f"%%% {label}  ({sz}, budget {budget} s, {name})")
        for h in HEADER:
            print(h)
        star = []
        prev = None
        for i in insts:
            cls = i.split("-")[1]
            if prev is not None and cls != prev:
                print("\\addlinespace")
            prev = cls
            opt = any(canon.get((i, m)) is not None
                      and canon[(i, m)]["status"] == "OPTIMAL" for m in "012")
            if opt:
                star.append(i)
            line = [i + ("$^{\\star}$" if opt else "")]
            r_of = canon
            for m in "012":
                c = cell(canon.get((i, m)))
                if c is None:
                    line += ["---", "---", "---", "---"]
                    continue
                u, b, g, t = c
                if u is None:
                    line += ["---", "---" if b is None else num(b), "---",
                             f"{round(t):d}" if t is not None else "---"]
                    continue
                v = num(u)
                if u < bks[i] - 0.5:
                    v = "\\textbf{" + v + "}"
                # Rodada 8, m1 do parecer. A convencao declarada nas legendas
                # e que um gap positivo de um run que nao fechou nunca imprime
                # 0.00, para nao ser lido como fechamento. A regra vive aqui,
                # no emissor, e nao na digitacao.
                closed = r_of[(i, m)]["status"] == "OPTIMAL"
                gtxt = f"{g:.2f}"
                if not closed and gtxt == "0.00" and g > 0:
                    gtxt = "$<$0.01"
                line += [v, "---" if b is None else num(b),
                         gtxt, f"{round(t):d}"]
            line.append(num(bks[i]))
            print(" & ".join(line) + " \\\\")
        print("\\bottomrule")
        print("\\end{tabular}")
        print(f"%%% stars: {len(star)} {star}")
        nb = sum(1 for i in insts
                 for m in "012"
                 if (c := cell(canon.get((i, m)))) and c[0] is not None
                 and c[0] < bks[i] - 0.5)
        ni = sum(1 for i in insts
                 if any((c := cell(canon.get((i, m)))) and c[0] is not None
                        and c[0] < bks[i] - 0.5 for m in "012"))
        print(f"%%% bold values: {nb}, instances with at least one: {ni}")
        miss = [(i, m) for i in insts for m in "012" if (i, m) not in canon]
        print(f"%%% missing cells: {miss}")
        print()


if __name__ == "__main__":
    main()
