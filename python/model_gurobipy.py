"""MP-TSCFLP model (1)-(10) of Mauri et al. (2021) in gurobipy.

Mirror of src/exact_model.cpp, used for cross-validation of the C++ port and for
quick experiments. Run:  python model_gurobipy.py <instance> [time_limit]
"""
import sys

import gurobipy as gp
from gurobipy import GRB


def load_instance(path):
    with open(path) as f:
        tok = f.read().split()
    it = iter(tok)
    n = lambda: int(next(it))
    v = lambda: float(next(it))
    I, J, K, L = n(), n(), n(), n()
    q = [[v() for _ in range(L)] for _ in range(K)]
    b, fy = [], []
    for _ in range(I):
        b.append([v() for _ in range(L)])
        fy.append(v())
    c = [[[v() for _ in range(J)] for _ in range(I)] for _ in range(L)]
    p, gz = [], []
    for _ in range(J):
        p.append([v() for _ in range(L)])
        gz.append(v())
    d = [[[v() for _ in range(K)] for _ in range(J)] for _ in range(L)]
    return dict(I=I, J=J, K=K, L=L, q=q, b=b, f=fy, c=c, p=p, g=gz, d=d)


def build_model(inst, env=None):
    I, J, K, L = inst["I"], inst["J"], inst["K"], inst["L"]
    m = gp.Model("mptscfl", env=env) if env else gp.Model("mptscfl")
    y = m.addVars(I, vtype=GRB.BINARY, name="y")
    z = m.addVars(J, vtype=GRB.BINARY, name="z")
    x = m.addVars(I, J, L, lb=0.0, name="x")
    w = m.addVars(J, K, L, lb=0.0, name="w")
    for i in range(I):
        y[i].BranchPriority = 10
    for j in range(J):
        z[j].BranchPriority = 10
    m.setObjective(
        gp.quicksum(inst["f"][i] * y[i] for i in range(I))
        + gp.quicksum(inst["g"][j] * z[j] for j in range(J))
        + gp.quicksum(inst["c"][l][i][j] * x[i, j, l]
                      for i in range(I) for j in range(J) for l in range(L))
        + gp.quicksum(inst["d"][l][j][k] * w[j, k, l]
                      for j in range(J) for k in range(K) for l in range(L)),
        GRB.MINIMIZE)
    m.addConstrs((w.sum("*", k, l) >= inst["q"][k][l]
                  for k in range(K) for l in range(L)), name="ClientDemand")
    m.addConstrs((x.sum("*", j, l) >= w.sum(j, "*", l)
                  for j in range(J) for l in range(L)), name="FlowConservation")
    m.addConstrs((x.sum(i, "*", l) <= inst["b"][i][l] * y[i]
                  for i in range(I) for l in range(L)), name="FCapacity")
    m.addConstrs((w.sum(j, "*", l) <= inst["p"][j][l] * z[j]
                  for j in range(J) for l in range(L)), name="WCapacity")
    return m, y, z, x, w


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    inst = load_instance(sys.argv[1])
    m, *_ = build_model(inst)
    m.Params.TimeLimit = float(sys.argv[2]) if len(sys.argv) > 2 else 3600.0
    m.Params.NodefileStart = 4.0
    m.optimize()
    if m.SolCount:
        print(f"obj={m.ObjVal:.2f} bound={m.ObjBound:.2f} gap={m.MIPGap:.4%} "
              f"time={m.Runtime:.1f}s")


if __name__ == "__main__":
    main()
