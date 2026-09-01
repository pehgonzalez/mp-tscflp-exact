"""Root-bound measurements for one MP-TSCFLP instance.

The study reports the hierarchy v_BR = v_LP <= v_LD <= v* as a proposition and
never measures it. This script measures the three quantities that are missing,
one instance per invocation, appending a single row to logs/bounds.csv:

  v_lp      LP relaxation of the compact model (5)-(10)
  v_br      LP value of the Benders reformulation, root cut loop run to
            convergence with continuous y and z
  root_mip  dual bound at the end of the root node of the compact MIP
  root_bbc  dual bound at the end of the root node of the Benders master

v_LD is deliberately absent. Every canonical L-BBC row of logs/results.csv
already carries it in lag_lb, and recomputing it here, under a second
implementation, would put two numbers for the same quantity in one paper.

Two choices deserve a note.

The cut used for v_br carries every coefficient, including the ones the C++
build folds into the constant when they fall below its sparsification
threshold. That threshold is numerical hygiene inside a branch-and-bound run,
not part of the definition of the reformulation, and folding a coefficient away
can only weaken a cut, so keeping all of them is what makes the measured v_br
comparable with v_lp.

The core point is off by default on every class, including the two
hundred-plant sets where the published campaign runs with it on. The table
compares bound strength across classes, so a setting that varies with the class
would confound exactly the comparison the table exists to make. The column
core_point records what was used.

Nothing here writes to logs/results.csv and nothing here touches the solver
binary, so the campaign the paper reports stays on the build that produced it.

Usage:
  python measure_bounds.py <instance.txt> [--cap 300] [--threads 16] [--seed 0]
                           [--out ../logs/bounds.csv] [--papadakos 0]
                           [--measures v_lp,v_br,root_mip,root_bbc]
  python measure_bounds.py --selftest
"""
import argparse
import datetime
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import gurobipy as gp
from gurobipy import GRB

from model_gurobipy import load_instance, build_model
from benders_gurobipy import Subproblem

EPS = 1e-6

FIELDS = ['datetime', 'instance', 'I', 'J', 'K', 'L', 'seed', 'threads', 'cap_s',
          'core_point',
          'v_lp', 'v_lp_s', 'v_lp_ok',
          'v_br', 'v_br_s', 'v_br_ok', 'v_br_rounds', 'v_br_cuts',
          'root_mip', 'root_mip_s', 'root_mip_ok',
          'root_bbc', 'root_bbc_s', 'root_bbc_ok', 'root_bbc_cuts',
          'gurobi_version']


def _bound_of(m):
    # ObjBound is undefined before the solver has produced one, and reading it
    # then raises rather than returning a sentinel.
    try:
        return float(m.ObjBound)
    except (gp.GurobiError, AttributeError):
        return None


def _theta_floor(inst, sp, l):
    # v_l with every plant and every depot open, the initial bound the C++
    # master gives theta_l.
    r = sp[l].solve([1] * inst['I'], [1] * inst['J'])
    if r is None:
        raise RuntimeError('the all-open point is infeasible for product %d, '
                           'which the benchmark is supposed to rule out' % l)
    return r[0]


def measure_v_lp(inst, env, args):
    t0 = time.time()
    m, _y, _z, _x, _w = build_model(inst, env)
    # relax() copies whatever the solver has already been told about, and
    # gurobipy batches additions until an update. Without this call the copy
    # comes back empty and its optimum is zero.
    m.update()
    r = m.relax()
    r.Params.Threads = args.threads
    r.Params.Seed = args.seed
    r.Params.TimeLimit = args.cap
    r.optimize()
    ok = r.Status == GRB.OPTIMAL
    val = float(r.ObjVal) if ok else _bound_of(r)
    dt = time.time() - t0
    m.dispose()
    r.dispose()
    return val, dt, ok


def measure_root_mip(inst, env, args):
    t0 = time.time()
    m, _y, _z, _x, _w = build_model(inst, env)
    m.Params.Threads = args.threads
    m.Params.Seed = args.seed
    m.Params.MIPGap = 0.0
    m.Params.MIPGapAbs = 0.9999
    m.Params.NodefileStart = 4.0
    # NodeLimit = 1, not 0. A limit of zero stops the solver before the root
    # cut loop has run, so what comes back is the plain LP relaxation again
    # rather than a root bound. Measured on four random instances, a limit of
    # zero reproduced v_lp to eight digits while a limit of one landed strictly
    # between v_lp and the optimum.
    m.Params.NodeLimit = 1
    m.Params.TimeLimit = args.cap
    m.optimize()
    val = _bound_of(m)
    ok = m.Status in (GRB.NODE_LIMIT, GRB.OPTIMAL) and val is not None
    dt = time.time() - t0
    m.dispose()
    return val, dt, ok


def measure_v_br(inst, env, sp, args, max_rounds=1000):
    # Cutting-plane loop on the master with y and z continuous. The LP value
    # rises monotonically as cuts arrive, so the number this returns is a valid
    # lower bound on v_BR even when the cap stops the loop early.
    I, J, L = inst['I'], inst['J'], inst['L']
    t0 = time.time()
    m = gp.Model('vbr', env=env)
    y = m.addVars(I, lb=0.0, ub=1.0)
    z = m.addVars(J, lb=0.0, ub=1.0)
    th = m.addVars(L, lb=0.0)
    m.setObjective(gp.quicksum(inst['f'][i] * y[i] for i in range(I))
                   + gp.quicksum(inst['g'][j] * z[j] for j in range(J))
                   + th.sum(), GRB.MINIMIZE)
    for l in range(L):
        D = sum(inst['q'][k][l] for k in range(inst['K']))
        m.addConstr(gp.quicksum(inst['b'][i][l] * y[i] for i in range(I)) >= D)
        m.addConstr(gp.quicksum(inst['p'][j][l] * z[j] for j in range(J)) >= D)
        th[l].LB = _theta_floor(inst, sp, l)
    m.Params.Threads = args.threads
    m.Params.Seed = args.seed

    rounds, cuts, val = 0, 0, None
    while True:
        left = args.cap - (time.time() - t0)
        if left <= 0.0:
            break
        m.Params.TimeLimit = left
        m.optimize()
        if m.Status != GRB.OPTIMAL:
            break
        rounds += 1
        val = float(m.ObjVal)
        yv = [y[i].X for i in range(I)]
        zv = [z[j].X for j in range(J)]
        tv = [th[l].X for l in range(L)]
        added = 0
        for l in range(L):
            r = sp[l].solve(yv, zv)
            if r is None:
                continue
            v, const, fac, ware = r
            tol = min(EPS * max(1.0, abs(v)), 0.4 / L)
            if tv[l] < v - tol:
                m.addConstr(th[l] >= const
                            + gp.quicksum(fac[i] * y[i] for i in range(I))
                            + gp.quicksum(ware[j] * z[j] for j in range(J)))
                added += 1
        cuts += added
        if added == 0:
            dt = time.time() - t0
            m.dispose()
            return val, dt, True, rounds, cuts
        if rounds >= max_rounds:
            break
    dt = time.time() - t0
    m.dispose()
    return val, dt, False, rounds, cuts


def measure_root_bbc(inst, env, args):
    from benders_gurobipy import BendersSolver
    t0 = time.time()
    solver = BendersSolver(inst, env=env, papadakos=bool(args.papadakos),
                           root_cuts=True)
    solver.m.Params.Threads = args.threads
    solver.m.Params.Seed = args.seed
    solver.m.Params.MIPGap = 0.0
    solver.m.Params.MIPGapAbs = 0.5
    # Same reason as in measure_root_mip, with a sharper consequence here. The
    # master carries no optimality cut until the callback has fired, so a limit
    # of zero can return a bound below the compact LP value, which is what
    # happened on all eight instances of the randomized cross-check.
    solver.m.Params.NodeLimit = 1
    solver.solve(time_limit=args.cap, output=False)
    val = _bound_of(solver.m)
    ok = solver.m.Status in (GRB.NODE_LIMIT, GRB.OPTIMAL) and val is not None
    dt = time.time() - t0
    cuts = solver.ncuts
    solver.m.dispose()
    return val, dt, ok, cuts


def fmt(v, nd=4):
    return '' if v is None else ('%.*f' % (nd, v))


def run_one(path, args):
    inst = load_instance(path)
    stem = os.path.splitext(os.path.basename(path))[0]
    env = gp.Env(empty=True)
    env.setParam('OutputFlag', 0)
    env.start()

    want = set(s.strip() for s in args.measures.split(',') if s.strip())
    row = {f: '' for f in FIELDS}
    row.update({
        'datetime': datetime.datetime.now().strftime('%Y%m%d-%H%M%S'),
        'instance': stem,
        'I': inst['I'], 'J': inst['J'], 'K': inst['K'], 'L': inst['L'],
        'seed': args.seed, 'threads': args.threads, 'cap_s': args.cap,
        'core_point': int(bool(args.papadakos)),
        'gurobi_version': '%d.%d.%d' % gp.gurobi.version(),
    })

    if 'v_lp' in want:
        v, dt, ok = measure_v_lp(inst, env, args)
        row.update(v_lp=fmt(v), v_lp_s='%.2f' % dt, v_lp_ok=int(ok))
        print('  v_lp      %s  %.1fs  ok=%d' % (fmt(v, 2), dt, ok), flush=True)

    if 'v_br' in want or 'root_bbc' in want:
        sp = [Subproblem(inst, l, env) for l in range(inst['L'])]
    if 'v_br' in want:
        v, dt, ok, rounds, cuts = measure_v_br(inst, env, sp, args)
        row.update(v_br=fmt(v), v_br_s='%.2f' % dt, v_br_ok=int(ok),
                   v_br_rounds=rounds, v_br_cuts=cuts)
        print('  v_br      %s  %.1fs  ok=%d  rounds=%d  cuts=%d'
              % (fmt(v, 2), dt, ok, rounds, cuts), flush=True)

    if 'root_mip' in want:
        v, dt, ok = measure_root_mip(inst, env, args)
        row.update(root_mip=fmt(v), root_mip_s='%.2f' % dt, root_mip_ok=int(ok))
        print('  root_mip  %s  %.1fs  ok=%d' % (fmt(v, 2), dt, ok), flush=True)

    if 'root_bbc' in want:
        v, dt, ok, cuts = measure_root_bbc(inst, env, args)
        row.update(root_bbc=fmt(v), root_bbc_s='%.2f' % dt, root_bbc_ok=int(ok),
                   root_bbc_cuts=cuts)
        print('  root_bbc  %s  %.1fs  ok=%d  cuts=%d'
              % (fmt(v, 2), dt, ok, cuts), flush=True)

    env.dispose()
    return row


def append_row(out, row):
    # The row is built in memory and written once, at the end, so an
    # interruption anywhere above leaves the file with whole rows only.
    out = os.path.abspath(out)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fresh = not os.path.exists(out) or os.path.getsize(out) == 0
    line = ','.join(str(row[f]) for f in FIELDS)
    with open(out, 'a', newline='') as fh:
        if fresh:
            fh.write(','.join(FIELDS) + '\n')
        fh.write(line + '\n')
        fh.flush()
        os.fsync(fh.fileno())


def write_selftest_instance(path):
    # Small enough for the size-limited licence that ships with the pip wheel,
    # which is what makes this runnable off the campaign machine.
    I, J, K, L = 3, 2, 4, 2
    out = ['%d %d %d %d' % (I, J, K, L)]
    for _k in range(K):
        out.append(' '.join(['10'] * L))
    for i in range(I):
        out.append(' '.join(['60'] * L))
        out.append('%d' % (100 + 10 * i))
    for l in range(L):
        for i in range(I):
            out.append(' '.join(str(1 + (i + j + l) % 4) for j in range(J)))
    for j in range(J):
        out.append(' '.join(['80'] * L))
        out.append('%d' % (50 + 5 * j))
    for l in range(L):
        for j in range(J):
            out.append(' '.join(str(2 + (j + k + l) % 3) for k in range(K)))
    with open(path, 'w') as fh:
        fh.write('\n'.join(out) + '\n')
    return path


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('instance', nargs='?')
    ap.add_argument('--cap', type=float, default=300.0,
                    help='wall-clock cap per measurement, in seconds')
    ap.add_argument('--threads', type=int, default=16)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--papadakos', type=int, default=0)
    ap.add_argument('--measures', default='v_lp,v_br,root_mip,root_bbc')
    ap.add_argument('--out', default=None,
                    help='defaults to ../logs/bounds.csv next to this script')
    ap.add_argument('--selftest', action='store_true')
    args = ap.parse_args()

    if args.out is None:
        here = os.path.dirname(os.path.abspath(__file__))
        args.out = os.path.join(here, '..', 'logs', 'bounds.csv')

    if args.selftest:
        import tempfile
        tmp = os.path.join(tempfile.mkdtemp(), 'SELFTEST-C0-3-2.txt')
        write_selftest_instance(tmp)
        args.cap = 60.0
        args.threads = 1
        args.out = os.path.join(os.path.dirname(tmp), 'bounds.csv')
        print('selftest on %s' % tmp, flush=True)
        row = run_one(tmp, args)
        append_row(args.out, row)
        lp = float(row['v_lp']) if row['v_lp'] else None
        br = float(row['v_br']) if row['v_br'] else None
        if lp is None or br is None:
            print('SELFTEST FAILED: a measurement returned nothing')
            return 1
        if abs(lp - br) > 1e-4 * max(1.0, abs(lp)):
            print('SELFTEST FAILED: v_lp %.6f and v_br %.6f disagree, and '
                  'Proposition 5 says they must not' % (lp, br))
            return 1
        tol = 1e-4 * max(1.0, abs(lp))
        rm = float(row['root_mip']) if row['root_mip'] else None
        if rm is not None and rm < lp - tol:
            print('SELFTEST FAILED: root_mip %.6f is below v_lp %.6f, which '
                  'means the root cut loop did not run' % (rm, lp))
            return 1
        rb = float(row['root_bbc']) if row['root_bbc'] else None
        if rb is not None and rb < lp - tol:
            print('SELFTEST WARNING: root_bbc %.6f is below v_lp %.6f. The '
                  'master root cut loop stopped early, so the row is usable '
                  'only if root_bbc_ok is read together with it' % (rb, lp))
        print('selftest ok, v_lp = v_br = %.6f' % lp)
        return 0

    if not args.instance:
        ap.error('an instance file is required unless --selftest is given')
    print('%s' % os.path.basename(args.instance), flush=True)
    row = run_one(args.instance, args)
    append_row(args.out, row)
    print('row appended to %s' % os.path.abspath(args.out), flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
