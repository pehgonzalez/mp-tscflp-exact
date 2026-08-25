"""Regenerates the tables and headline numbers of the paper from results/results.csv.

Canonical run selection: mode 'exact', seed 0, one time limit per size set
(1800 s except 900 s on the fifty-plant ten-product set), and the most recent
row per (instance, method). The two extended-budget closures are read from the
14400 s rows. Engine codes in the CSV: method 0 = MIP, 1 = BBC, 2 = L-BBC.

Usage: python make_tables.py [results.csv] [bks_reference.csv]
Prints the closed-instances table, the improved best-known-solutions table,
and the headline counts quoted in the paper, so every number is traceable.
"""
import csv
import os
import sys
from collections import defaultdict


def _clamp_walls(rows):
    # Section 5.1 reports the time of a run that reached its budget at the
    # budget itself. The recorded wall can pass the limit because model
    # construction and the final solver-free re-costing sit outside the
    # solver's clock, so the clamp lives at the loader and every table,
    # figure and test inherits the declared convention.
    for r in rows:
        try:
            w, tl = float(r["total_wall_s"]), float(r["time_limit_s"])
            if w > tl:
                r["total_wall_s"] = "%.2f" % tl
        except (KeyError, TypeError, ValueError):
            pass
    return rows


RES = sys.argv[1] if len(sys.argv) > 1 else "results/results.csv"
BKS = sys.argv[2] if len(sys.argv) > 2 else "data/bks_reference.csv"
TL = {"50-5": "1800", "50-10": "1800", "100-5": "1800", "100-10": "1800"}
ENG = {"0": "MIP", "1": "BBC", "2": "L-BBC"}
# Desempate unico das duas tabelas de atribuicao, menor relogio de parede e,
# em empate exato, a ordem fixa BBC, L-BBC, MIP.
TIE = {"1": 0, "2": 1, "0": 2}


def size_of(inst):
    parts = inst.split("-")
    return parts[-2] + "-" + parts[-1]


def load():
    rows = _clamp_walls(list(csv.DictReader(open(RES))))
    bks = {r["instance"].replace(".txt", ""): float(r["bks"])
           for r in csv.DictReader(open(BKS))}
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
    # The extended pass ran at 7200 s and at 14400 s. Every one of its runs
    # belongs in the pool that attributes best upper bounds, whether or not it
    # closed, because a run that ends open can still lower the incumbent. Only
    # the subset that ends OPTIMAL may enter the certificates table, so the two
    # sets are kept apart here rather than conflated by a status filter. Rows
    # of the certification modes are excluded, since those repeat instances
    # already closed under the uniform budget and would double-count them.
    ext_rows = [r for r in rows
                if r["mode"] == "exact" and r["time_limit_s"] in ("7200", "14400")]
    # Dentro do passe estendido um certificado ja registado nao expira quando
    # uma repeticao posterior estoura o orcamento, a mesma convencao que a
    # tabela de fechamentos declara, entao a linha OPTIMAL mais recente ganha
    # do resto e so na sua ausencia vale a linha mais recente. Sem isso as duas
    # tabelas atribuiriam engines diferentes ao mesmo run.
    ext_best = {}
    for r in ext_rows:
        k = (r["instance"].replace(".txt", ""), r["method"])
        cur = ext_best.get(k)
        if cur is None:
            ext_best[k] = r
            continue
        rank = lambda x: (x["status"] == "OPTIMAL", x["datetime"])
        if rank(r) > rank(cur):
            ext_best[k] = r
    extended = list(ext_best.values())
    extended_closed = [r for r in rows
                       if r["mode"].startswith("exact") and r["seed"] == "0"
                       and r["mode"] not in ("exact-mauri",)
                       and r["status"] == "OPTIMAL"
                       and r["method"] in ("1", "2")
                       and float(r["time_limit_s"]) > 1800]
    # Pool for best-UB attribution: canonical rows plus the whole extended pass.
    pool = defaultdict(list)
    for r in list(canon.values()) + extended:
        pool[r["instance"].replace(".txt", "")].append(r)
    return rows, bks, canon, extended, extended_closed, pool


def ub(r):
    # The paper reports the solver-free re-costing of the incumbent, not the
    # solver objective, and the two differ on 55 of the 100 direct solves. The
    # verified value is therefore the primary source and the solver objective
    # only the fallback for rows written before the evaluator existed.
    v = float(r.get("verified_cost", -1) or -1)
    if v > 0:
        return v
    o = float(r["obj"])
    h = float(r["heuristic_cost"])
    return o if o > 0 else (h if h > 0 else float("inf"))


# Emissor da tabela de fechamentos (M1 da rodada 3). Cada linha declara a
# procedencia do tempo que imprime. A ordem de precedencia e estrita, primeiro
# o lote canonico homogeneo, depois o passe estendido do mesmo modo, e so
# quando nenhum dos dois fecha a instancia e que se aceita um fechamento de
# outra configuracao, sempre com nota explicita. Nenhum tempo e minimo sobre
# repeticoes de uma mesma configuracao.
PROV = {
    "exact": "an earlier first-campaign run at the same budget, which carried "
             "the core-point cuts the revision campaign retired",
    "exact-pap": "the retired core-point setting",
    "exact-var1": "a repetition of the variance study",
    "exact-var2": "a repetition of the variance study",
    "exact-var3": "a repetition of the variance study",
    "exact-cert": "the certification pass",
    "exact-cert2": "the certification pass",
    "exact-agg": "the aggregated-master ablation",
    "exact-noldb": "the withheld-bound arm of the factorial",
    "exact-nowarm": "the cold arm of the factorial",
    "exact-nostart": "the withheld-start arm of the factorial",
}


def closure_table(rows, bks, canon, insts):
    """Returns (rows_for_latex, notes, disclosure) for the certificates table."""
    # Os tres engines entram, nao so a decomposicao, para que a regra de
    # atribuicao seja literalmente a mesma da tab:bks (menor relogio de parede,
    # empate exato na ordem BBC, L-BBC, MIP) sobre o mesmo conjunto de runs.
    def opt_rows(pred):
        return [r for r in rows
                if r["mode"].startswith("exact") and r["seed"] == "0"
                and r["status"] == "OPTIMAL" and r["mode"] != "exact-mauri"
                and pred(r)]

    # Lote canonico, uma linha por (instancia, metodo), a mais recente.
    canon_opt = {}
    for (i, m), r in canon.items():
        if r["status"] == "OPTIMAL":
            canon_opt.setdefault(i, []).append(r)
    # Passe estendido, modo exact apenas, a linha mais recente por par.
    ext = {}
    for r in opt_rows(lambda r: r["mode"] == "exact"
                      and float(r["time_limit_s"]) > 1800):
        k = (r["instance"].replace(".txt", ""), r["method"])
        if k not in ext or r["datetime"] > ext[k]["datetime"]:
            ext[k] = r
    ext_opt = {}
    for (i, _), r in ext.items():
        ext_opt.setdefault(i, []).append(r)
    # Qualquer outro fechamento registado, para as notas e a divulgacao.
    other = {}
    for r in opt_rows(lambda r: True):
        other.setdefault(r["instance"].replace(".txt", ""), []).append(r)

    out, notes, disclosure = [], [], []
    for i in insts:
        pick, kind = None, None
        key = lambda r: (float(r["total_wall_s"]), TIE[r["method"]])
        if i in canon_opt:
            pick = min(canon_opt[i], key=key)
            kind = "canonical"
        elif i in ext_opt:
            pick = min(ext_opt[i], key=key)
            kind = "extended"
        elif i in other:
            pick = min(other[i], key=key)
            kind = "outside"
        if pick is None:
            continue
        out.append((i, ub(pick), float(pick["total_wall_s"]),
                    ENG[pick["method"]], kind, pick["mode"],
                    pick["time_limit_s"], pick["datetime"]))
        if kind == "outside":
            notes.append((i, PROV.get(pick["mode"], pick["mode"]),
                          ENG[pick["method"]], float(pick["total_wall_s"]),
                          pick["mode"], pick["datetime"]))
        # Divulgacao, fechamentos mais rapidos registados sob outra
        # configuracao do que a linha imprime.
        if kind != "outside":
            faster = [r for r in other.get(i, [])
                      if r["mode"] != "exact"
                      and r["time_limit_s"] == pick["time_limit_s"]
                      and float(r["total_wall_s"]) < float(pick["total_wall_s"]) - 1]
            if faster:
                f = min(faster, key=lambda r: float(r["total_wall_s"]))
                disclosure.append((i, PROV.get(f["mode"], f["mode"]),
                                   ENG[f["method"]], float(f["total_wall_s"]),
                                   float(pick["total_wall_s"])))
    order = {"50-5": 0, "50-10": 1, "100-5": 2, "100-10": 3}
    out.sort(key=lambda t: (order[size_of(t[0])], t[0]))
    disclosure.sort(key=lambda t: (order[size_of(t[0])], t[0]))
    return out, notes, disclosure


def main():
    rows, bks, canon, extended, extended_closed, pool = load()
    insts = sorted(set(i for i, _ in canon))


    def numtex(x):
        return f"{round(x):,}".replace(",", "{,}")

    print("=== Closed instances (certificates) ===")
    closed, notes, disclosure = closure_table(rows, bks, canon, insts)
    for i, opt, w, eng, kind, mode, tl, dt in closed:
        print(f"  {i:16} optimum {opt:>12,.0f}  wall {w:>7.0f}s  {eng:6}"
              f"  {kind:9} mode {mode:15} TL {tl:>5}  {dt}  BKS {bks[i]:>12,.0f}")
    print(f"  total closed: {len(closed)}")
    print("  provenance notes:")
    for i, label, eng, w, mode, dt in notes:
        print(f"    {i:16} only closure recorded under {label} ({mode}), "
              f"{eng} in {w:.0f} s on {dt}")
    print("  faster closures recorded under other configurations (disclosure):")
    for i, label, eng, w, printed in disclosure:
        print(f"    {i:16} {label:38} {eng:6} {w:6.0f} s   against {printed:6.0f} s printed")

    # --- emissao LaTeX do corpo de tab:closed ---
    # Rodada 4. A tabela ganha duas colunas, o desfecho do refechamento sob o
    # build que grava o master final e o resultado da verificacao racional
    # sobre o registro que esse refechamento deixou. As duas saem de reclose(),
    # que le o results file e rational_report.txt, nunca de digitacao.
    rc = reclose_outcome(rows)
    outdir = "paper" if os.path.isdir("paper") else "."
    with open(os.path.join(outdir, "tab_closed_body.tex"), "w") as f:
        f.write("\\begin{tabular}{lrrlrll}\n\\toprule\n"
                "Instance & Optimum & Time (s) & Engine & BKS & Batch & "
                "Re-closure \\\\\n\\midrule\n")
        for i, opt, w, eng, kind, mode, tl, dt in closed:
            mark = {"canonical": "", "extended": "$^{\\ddagger}$",
                    "outside": "$^{\\S}$"}[kind]
            batch = {"canonical": "canonical", "extended": "extended",
                     "outside": "see note"}[kind]
            f.write(f"{i} & {numtex(opt)} & {w:.0f}{mark} & \\textsc{{{eng}}} & "
                    f"{numtex(bks[i])} & {batch} & {rc.get(i, 'not rerun')} \\\\\n")
        f.write("\\bottomrule\n\\end{tabular}\n")
    print("tab_closed_body.tex:", len(closed), "linhas,", len(notes), "notas")
    print("  re-closure column:", rc)

    # --- emissao LaTeX do corpo de tab:repro (rodada 6) ---
    reps, rtot, rtopt = repetition_note(rows, [c[0] for c in closed])
    with open(os.path.join(outdir, "tab_repro_body.tex"), "w") as f:
        f.write("\\begin{tabular}{lrrll}\n\\toprule\n"
                "Instance & Runs & Closures & Engines that closed it & "
                "Gap range of the open runs (\\%) \\\\\n\\midrule\n")
        for i, n, k, who, rng in reps:
            f.write(f"{i} & {n} & {k} & {who} & {rng} \\\\\n")
        f.write("\\midrule\n"
                f"All twelve & {rtot} & {rtopt} & & \\\\\n"
                "\\bottomrule\n\\end{tabular}\n")
    print("tab_repro_body.tex:", len(reps), "linhas, runs", rtot,
          "fechamentos", rtopt)

    print("\n=== Improved best-known solutions ===")
    # Atribuicao unificada com a tabela de fechamentos (rodada 3, item 2 do
    # gerente). Onde mais de um engine atinge o valor reportado, credita-se o
    # run de menor relogio de parede, e um empate exato de relogio cai na ordem
    # fixa BBC, L-BBC, MIP. A mesma regra decide o engine da tab:closed, onde
    # nao ha empate exato de relogio, entao as duas tabelas nunca discordam.
    improved = []
    for i in insts:
        cand = [(ub(r), float(r["total_wall_s"]), r["method"]) for r in pool[i]
                if ub(r) < float("inf")]
        if not cand:
            continue
        u = min(c[0] for c in cand)
        tied = [c for c in cand if c[0] < u + 0.5]
        pick = min(tied, key=lambda c: (c[1], TIE[c[2]]))
        if u < bks[i] - 0.5:
            improved.append((i, pick[0], ENG[pick[2]], bks[i] - pick[0]))
    for i, u, eng, d in sorted(improved, key=lambda t: (size_of(t[0]), t[0])):
        print(f"  {i:16} new UB {u:>12,.0f}  by {eng:5}  improvement {d:>9,.0f}")
    by_eng = defaultdict(int)
    for _, _, eng, _ in improved:
        by_eng[eng] += 1
    print(f"  total improved: {len(improved)}  by engine: {dict(by_eng)}")

    # --- emissao LaTeX da tab:bks (gerada pela reescrita de 2026-08-17) ---
    def numtex(x):
        return f"{round(x):,}".replace(",", "{,}")
    order={"50-5":0,"50-10":1,"100-5":2,"100-10":3}
    lat=[]
    for i in insts:
        cand=[(ub(r), float(r["total_wall_s"]), r["method"], float(r["time_limit_s"])>1800) for r in pool[i] if ub(r)<float("inf")]
        if not cand: continue
        u=min(c[0] for c in cand)
        tied=[c for c in cand if c[0]<u+0.5]
        pick=min(tied, key=lambda c:(c[1], TIE[c[2]]))
        if u<bks[i]-0.5:
            # O marcador diz de qual passe vem o run creditado, nao se algum
            # outro run tambem bateria o BKS.
            dag="$^{\\ddagger}$" if pick[3] else ""
            lat.append((order[size_of(i)],i,pick[0],ENG[pick[2]],dag))
    # Duas metades lado a lado, a tabela tem 34 linhas e ocupava pagina inteira.
    rows_tex = [i+" & "+numtex(u)+" & \\textsc{"+eng+"}"+dag+" & "
                +numtex(bks[i])+" & "+numtex(bks[i]-u)+" \\\\\n"
                for _,i,u,eng,dag in sorted(lat)]
    half = (len(rows_tex)+1)//2
    head = ("\\begin{tabular}[t]{lrlrr}\n\\toprule\n"
            "Instance & New UB & Engine & BKS & Impr. \\\\\n\\midrule\n")
    with open(os.path.join("paper" if os.path.isdir("paper") else ".",
                           "tab_bks_body.tex"), "w") as f:
        f.write(head)
        f.writelines(rows_tex[:half])
        f.write("\\bottomrule\n\\end{tabular}\\hfill\n")
        f.write(head)
        f.writelines(rows_tex[half:])
        f.write("\\bottomrule\n\\end{tabular}\n")
    print("tab_bks_body.tex:", len(lat), "linhas,",
          sum(1 for t in lat if t[4]), "marcadas com ddagger")


    print("\n=== Headline counts (hundred-plant sets) ===")
    hundred = [i for i in insts if size_of(i) in ("100-5", "100-10")]
    dec_inc = mip_inc = 0
    mip_bound = 0
    lbbc_sole = 0
    dec_cert = set()
    mip_cert = set()
    for i in hundred:
        ubs = {m: ub(canon[(i, m)]) for m in "012" if (i, m) in canon}
        bounds = {m: float(canon[(i, m)]["bound"]) for m in "012"
                  if (i, m) in canon and float(canon[(i, m)]["bound"]) > 0}
        if "0" in ubs and any(m in ubs for m in "12"):
            dec = min(ubs[m] for m in "12" if m in ubs)
            if dec < ubs["0"] - 0.5:
                dec_inc += 1
            elif ubs["0"] < dec - 0.5:
                mip_inc += 1
        if "0" in bounds and any(m in bounds for m in "12"):
            if bounds["0"] > max(bounds[m] for m in "12" if m in bounds) + 0.5:
                mip_bound += 1
        best = min(ubs.values())
        win = [m for m in ubs if abs(ubs[m] - best) < 0.5]
        if win == ["2"]:
            lbbc_sole += 1
        for m in "012":
            r = canon.get((i, m))
            if r and r["status"] == "OPTIMAL":
                (mip_cert if m == "0" else dec_cert).add(i)
    print(f"  decomposition better incumbent: {dec_inc} / {len(hundred)}   MIP better: {mip_inc}")
    print(f"  MIP strictly tighter dual bound: {mip_bound} / {len(hundred)}")
    print(f"  certificates: decomposition {len(dec_cert)}, MIP {len(mip_cert)}")
    print(f"  L-BBC sole best incumbent: {lbbc_sole}")



def repetition_note(rows, closed):
    """Rodada 6, required 7 do parecer. Reconstroi do master mesclado a
    reproducibilidade de cada instancia fechada, contando TODAS as repeticoes
    registradas, com sucesso e falha, e a faixa de gap final das que nao
    fecharam. Emite o corpo de uma tabela, nunca texto digitado."""
    eng = {"0": "\\textsc{MIP}", "1": "\\textsc{BBC}", "2": "\\textsc{L-BBC}"}
    out, tot, topt = [], 0, 0
    for i in closed:
        rs = [r for r in rows if r["instance"] == i and _pos(r.get("obj"))]
        opt = [r for r in rs if r["status"] == "OPTIMAL"]
        gaps = sorted(float(r["gap"]) for r in rs
                      if r["status"] != "OPTIMAL" and r.get("gap") not in (None, ""))
        who = ", ".join(eng[m] for m in sorted({r["method"] for r in opt}))
        rng = ("%.3f--%.2f" % (gaps[0] * 100, gaps[-1] * 100)) if gaps else "---"
        out.append((i, len(rs), len(opt), who, rng))
        tot += len(rs)
        topt += len(opt)
    return out, tot, topt


def _pos(x):
    try:
        return float(x) > 0
    except (TypeError, ValueError):
        return False

def reclose_outcome(rows):
    """Outcome of the re-closure pass per instance, for the last column of
    tab:closed. Status, wall clock and final gap come from the merged results
    file, the two rational checks from rational_report.txt, which is the raw
    output of the verifier. Nothing here is typed by hand."""
    rep = {}
    path = os.path.join(os.path.dirname(os.path.abspath(RES)), "rational_report.txt")
    if os.path.exists(path):
        cur = None
        for line in open(path, encoding="utf-8"):
            line = line.strip()
            if line.startswith("== ") and line.endswith(" =="):
                cur = line[3:-3].strip()
                rep[cur] = {"c1": False, "c2": False, "slack": None}
            elif cur and line.startswith("CHECK 1 PASS"):
                rep[cur]["c1"] = True
            elif cur and line.startswith("CHECK 2 PASS"):
                rep[cur]["c2"] = True
            elif cur and "custo verificado" in line:
                rep[cur]["slack"] = line.replace("=", " ").split()[-1]
    out = {}
    print("\n=== Re-closure pass under the recording build ===")
    for r in sorted([x for x in rows if x["mode"] == "exact-cert3"],
                    key=lambda x: x["instance"]):
        i = r["instance"].replace(".txt", "")
        if r["status"].upper().startswith("OPT"):
            d = rep.get(i)
            # Rodada 8, M2 do parecer. O engine do run auditado nao coincide
            # com o engine creditado pelo relogio na coluna Engine, e o leitor
            # precisa dele para ler a divisao 3-e-3 do audit: os mestres BBC
            # nao rodam a fase e por isso nao carregam a linha Lagrangiana.
            eng = "\\textsc{%s}" % ENG[r["method"]]
            if d and d["c1"] and d["c2"]:
                out[i] = "audited, %s, slack %s" % (eng, d["slack"])
            else:
                out[i] = "closed under %s, no record" % eng
        else:
            out[i] = "expired at %.3f\\%%" % (float(r["gap"]) * 100)
        print("  %-16s %-6s %-9s %s"
              % (i, ENG[r["method"]], r["status"], out[i]))
    return out



if __name__ == "__main__":
    main()
