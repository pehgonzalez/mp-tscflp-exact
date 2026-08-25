#!/usr/bin/env python3
"""Diagnostico por classe (M9/R10) a partir do lote canonico em results_merged.csv.

Lote canonico: mode 'exact', seed 0, TL 1800, linha mais recente por (instancia, metodo).
UB = verified_cost se verified_ok=1, senao min(obj>0, heuristic_cost>0).
Todas as medias sao aritmeticas simples sobre as instancias do grupo.
"""
import csv, re
from collections import defaultdict

ROWS = list(csv.DictReader(open('results_merged.csv')))

def num(x):
    try:
        v = float(x)
        return v
    except (TypeError, ValueError):
        return None

# lote canonico
canon = {}
for r in ROWS:
    if r['mode'] != 'exact': continue
    if r['seed'] != '0': continue
    if num(r['time_limit_s']) != 1800: continue
    key = (r['instance'], r['method'])
    if key not in canon or r['datetime'] > canon[key]['datetime']:
        canon[key] = r

def ub_of(r):
    if r['verified_ok'] == '1':
        v = num(r['verified_cost'])
        if v and v > 0: return v
    cands = [num(r['obj']), num(r['heuristic_cost'])]
    cands = [c for c in cands if c and c > 0]
    return min(cands) if cands else None

# melhor UB entre metodos por instancia (referencia de gap)
best_ub = {}
for (inst, m), r in canon.items():
    u = ub_of(r)
    if u and (inst not in best_ub or u < best_ub[inst]):
        best_ub[inst] = u

pat = re.compile(r'PSC(\d)-C(\d)-(\d+)-(\d+)')
def parse(inst):
    m = pat.match(inst)
    return int(m.group(1)), int(m.group(2)), f"{m.group(3)}-{m.group(4)}"

def mean(xs):
    xs = [x for x in xs if x is not None]
    return sum(xs)/len(xs) if xs else None

def fmt(x, d=3):
    return f"{x:.{d}f}" if x is not None else "  -  "

# agregadores: por (conjunto, classe C) e por (conjunto, gerador PSC)
def collect(group_key):
    agg = defaultdict(lambda: defaultdict(list))
    insts = sorted(best_ub)
    for inst in insts:
        psc, cclass, size = parse(inst)
        gk = group_key(psc, cclass, size)
        ub = best_ub[inst]
        r0 = canon.get((inst, '0')); r1 = canon.get((inst, '1')); r2 = canon.get((inst, '2'))
        if not (r0 and r1 and r2): continue
        b0, b1, b2 = num(r0['bound']), num(r1['bound']), num(r2['bound'])
        lag = num(r2['lag_lb'])
        # gaps proprios (UB do metodo vs bound do metodo)
        for tag, r in (('MIP', r0), ('BBC', r1), ('LBBC', r2)):
            u = ub_of(r); b = num(r['bound'])
            if u and b: agg[gk][f'gap_{tag}'].append(100*(u-b)/u)
        # qualidade dual relativa ao melhor UB
        if b0: agg[gk]['dgap_MIP'].append(100*(ub-b0)/ub)
        if b1: agg[gk]['dgap_BBC'].append(100*(ub-b1)/ub)
        if b2: agg[gk]['dgap_LBBC'].append(100*(ub-b2)/ub)
        if lag and lag > 0: agg[gk]['dgap_vLD'].append(100*(ub-lag)/ub)
        # vantagem Lagrangiana no bound do ramo-e-corte
        if b1 and b2: agg[gk]['adv_LBBC_BBC'].append(100*(b2-b1)/ub)
        # v_LD vs raiz/final do compacto
        if lag and lag > 0 and b0: agg[gk]['vLD_minus_MIPbound'].append(100*(lag-b0)/ub)
        agg[gk]['n'].append(1)
    return agg

def report(title, agg, keys):
    print(f"\n== {title} ==")
    header = f"{'grupo':>12} {'n':>3} " + " ".join(f"{k:>10}" for k in keys)
    print(header)
    for gk in sorted(agg):
        vals = [mean(agg[gk][k]) if k != 'n' else len(agg[gk]['n']) for k in keys]
        cells = []
        for k, v in zip(keys, vals):
            cells.append(f"{v:>10}" if k == 'n' else f"{fmt(v):>10}")
        print(f"{str(gk):>12} " + " ".join(cells))

K = ['n', 'gap_MIP', 'gap_BBC', 'gap_LBBC', 'dgap_vLD', 'adv_LBBC_BBC', 'vLD_minus_MIPbound']

report("Por conjunto x classe C", collect(lambda p, c, s: (s, f"C{c}")), K)
report("Por conjunto x gerador PSC", collect(lambda p, c, s: (s, f"PSC{p}")), K)
report("Por conjunto (decaimento com o tamanho)", collect(lambda p, c, s: s), K)
report("Por classe C (agregado)", collect(lambda p, c, s: f"C{c}"), K)
