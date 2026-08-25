"""Independent correctness check of the MIP against exhaustive enumeration.

Generates small random MP-TSCFLP instances (seeded), solves each with:
  (a) gurobipy MIP (model_gurobipy.build_model)  — fits the size-limited pip license;
  (b) brute force: enumerate all (y,z), solve the routing LP per product with Gurobi
      (an independent formulation: pure min-cost flow, no location variables).
Asserts both optima agree. This validates the model logic that the C++ port mirrors.
"""
import itertools
import random
import sys

import gurobipy as gp
from gurobipy import GRB

from model_gurobipy import build_model


def gen_instance(seed, I=2, J=3, K=4, L=2):
    rng = random.Random(seed)
    q = [[rng.randint(5, 20) for _ in range(L)] for _ in range(K)]
    tot = [sum(q[k][l] for k in range(K)) for l in range(L)]
    b = [[rng.randint(tot[l] // I, tot[l]) for l in range(L)] for _ in range(I)]
    p = [[rng.randint(tot[l] // J, tot[l]) for l in range(L)] for _ in range(J)]
    f = [rng.randint(200, 800) for _ in range(I)]
    g = [rng.randint(100, 500) for _ in range(J)]
    c = [[[rng.randint(1, 30) for _ in range(J)] for _ in range(I)] for _ in range(L)]
    d = [[[rng.randint(1, 30) for _ in range(K)] for _ in range(J)] for _ in range(L)]
    return dict(I=I, J=J, K=K, L=L, q=q, b=b, f=f, c=c, p=p, g=g, d=d)


def routing_cost(inst, yv, zv, env):
    """Min transport cost with facilities fixed; None if infeasible."""
    I, J, K, L = inst["I"], inst["J"], inst["K"], inst["L"]
    m = gp.Model("route", env=env)
    m.Params.OutputFlag = 0
    x = m.addVars(I, J, L, lb=0.0)
    w = m.addVars(J, K, L, lb=0.0)
    m.setObjective(
        gp.quicksum(inst["c"][l][i][j] * x[i, j, l]
                    for i in range(I) for j in range(J) for l in range(L))
        + gp.quicksum(inst["d"][l][j][k] * w[j, k, l]
                      for j in range(J) for k in range(K) for l in range(L)))
    m.addConstrs(w.sum("*", k, l) >= inst["q"][k][l] for k in range(K) for l in range(L))
    m.addConstrs(x.sum("*", j, l) >= w.sum(j, "*", l) for j in range(J) for l in range(L))
    m.addConstrs(x.sum(i, "*", l) <= inst["b"][i][l] * yv[i] for i in range(I) for l in range(L))
    m.addConstrs(w.sum(j, "*", l) <= inst["p"][j][l] * zv[j] for j in range(J) for l in range(L))
    m.optimize()
    return m.ObjVal if m.Status == GRB.OPTIMAL else None


def brute_force(inst, env):
    best = None
    for yv in itertools.product((0, 1), repeat=inst["I"]):
        for zv in itertools.product((0, 1), repeat=inst["J"]):
            rc = routing_cost(inst, yv, zv, env)
            if rc is None:
                continue
            tot = rc + sum(f * v for f, v in zip(inst["f"], yv)) \
                     + sum(g * v for g, v in zip(inst["g"], zv))
            if best is None or tot < best:
                best = tot
    return best


def main():
    nseeds = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    env = gp.Env(params={"OutputFlag": 0})
    for seed in range(nseeds):
        inst = gen_instance(seed)
        m, *_ = build_model(inst, env=env)
        m.Params.OutputFlag = 0
        m.optimize()
        assert m.Status == GRB.OPTIMAL, f"seed {seed}: MIP not optimal"
        bf = brute_force(inst, env)
        assert bf is not None, f"seed {seed}: brute force infeasible?"
        assert abs(m.ObjVal - bf) < 1e-6, f"seed {seed}: MIP={m.ObjVal} != BF={bf}"
        print(f"seed {seed}: MIP == brute force == {m.ObjVal:.2f}  OK")
    print("VALIDATION PASSED")


if __name__ == "__main__":
    main()
