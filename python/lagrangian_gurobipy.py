"""Lagrangian relaxation of the MP-TSCFLP (gurobipy reference implementation).

Dualizes flow conservation with multipliers lambda_jl >= 0.
  A(lmb): closed form over plants.
  B(lmb): single-stage multiproduct CFLP over depots, solved as MIP (gurobipy).
Subgradient with Polyak step + Lagrangian primal heuristic (repair via routing LP).

Every L(lmb) with B solved to optimality is a valid lower bound on v*.
"""
import math

import gurobipy as gp
from gurobipy import GRB


class SubproblemB:
    """min sum_j g_j z_j + sum_jkl (d_jkl + lmb_jl) w_jkl  s.t. demand, depot capacity."""

    def __init__(self, inst, env):
        J, K, L = inst["J"], inst["K"], inst["L"]
        self.inst, self.J, self.K, self.L = inst, J, K, L
        m = gp.Model("B", env=env)
        m.Params.OutputFlag = 0
        self.z = m.addVars(J, vtype=GRB.BINARY)
        self.w = m.addVars(J, K, L, lb=0.0)
        m.addConstrs(self.w.sum("*", k, l) >= inst["q"][k][l]
                     for k in range(K) for l in range(L))
        m.addConstrs(self.w.sum(j, "*", l) <= inst["p"][j][l] * self.z[j]
                     for j in range(J) for l in range(L))
        self.m = m

    def solve(self, lmb, time_limit=None):
        """MIPGapAbs = 0.25 plus a per-solve TimeLimit. ObjBound is dual
        feasible, so the returned bound is valid at any interruption. Returns
        None when no incumbent exists within the limit."""
        inst = self.inst
        self.m.setObjective(
            gp.quicksum(inst["g"][j] * self.z[j] for j in range(self.J))
            + gp.quicksum((inst["d"][l][j][k] + lmb[j][l]) * self.w[j, k, l]
                          for j in range(self.J) for k in range(self.K)
                          for l in range(self.L)), GRB.MINIMIZE)
        self.m.Params.MIPGap = 0.0
        self.m.Params.MIPGapAbs = 0.25
        if time_limit is not None:
            self.m.Params.TimeLimit = max(1.0, time_limit)
        self.m.optimize()
        if self.m.SolCount == 0:
            return None
        zbar = [int(self.z[j].X > 0.5) for j in range(self.J)]
        wout = [[sum(self.w[j, k, l].X for k in range(self.K)) for l in range(self.L)]
                for j in range(self.J)]
        return self.m.ObjBound, zbar, wout


def solve_A(inst, lmb):
    """Closed form of the plant block. Returns (value, ybar, xin[j][l] = stage-1 inflow)."""
    I, J, L = inst["I"], inst["J"], inst["L"]
    val = 0.0
    ybar = [0] * I
    xin = [[0.0] * L for _ in range(J)]
    plan = []
    for i in range(I):
        tot = inst["f"][i]
        moves = []
        for l in range(L):
            best_j = min(range(J), key=lambda j: inst["c"][l][i][j] - lmb[j][l])
            rc = inst["c"][l][i][best_j] - lmb[best_j][l]
            if rc < 0:
                tot += inst["b"][i][l] * rc
                moves.append((l, best_j, inst["b"][i][l]))
        if tot < 0:
            val += tot
            ybar[i] = 1
            plan.append((i, moves))
            for (l, j, amt) in moves:
                xin[j][l] += amt
    return val, ybar, xin


def routing_value(inst, yv, zv, env):
    """Exact transport cost given (y,z); None if infeasible. Independent check LP."""
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


def repair(inst, ybar, zbar, env):
    """Lagrangian primal heuristic: greedy (F1)/(F2) cover + exact routing."""
    I, J, L = inst["I"], inst["J"], inst["L"]
    y, z = list(ybar), list(zbar)
    for l in range(L):
        D = sum(inst["q"][k][l] for k in range(inst["K"]))
        capf = sum(inst["b"][i][l] * y[i] for i in range(I))
        for i in sorted(range(I), key=lambda i: inst["f"][i] / max(1, inst["b"][i][l])):
            if capf >= D: break
            if not y[i]:
                y[i] = 1
                capf += inst["b"][i][l]
        capw = sum(inst["p"][j][l] * z[j] for j in range(J))
        for j in sorted(range(J), key=lambda j: inst["g"][j] / max(1, inst["p"][j][l])):
            if capw >= D: break
            if not z[j]:
                z[j] = 1
                capw += inst["p"][j][l]
    rc = routing_value(inst, y, z, env)
    if rc is None:
        return None, math.inf
    total = rc + sum(f * v for f, v in zip(inst["f"], y)) \
              + sum(g * v for g, v in zip(inst["g"], z))
    return (y, z), total


def lp_dual_warmstart(inst, env):
    """Warm start: lambda0 = LP duals of flow conservation, so the first evaluation already reaches v_LP."""
    from model_gurobipy import build_model
    m, *_ = build_model(inst, env=env)
    m.Params.OutputFlag = 0
    m.update()
    r = m.relax()
    r.Params.OutputFlag = 0
    r.optimize()
    if r.Status != GRB.OPTIMAL:
        return None, None
    lmb = [[0.0] * inst["L"] for _ in range(inst["J"])]
    for c in r.getConstrs():
        if c.ConstrName.startswith("FlowConservation"):
            j, l = map(int, c.ConstrName.split("[")[1].rstrip("]").split(","))
            lmb[j][l] = max(0.0, c.Pi)
    return lmb, r.ObjVal


def lagrangian_dual(inst, env, iters=300, mu0=2.0, verbose=False, lmb0=None,
                    time_limit=None):
    """Subgradient ascent. Returns dict with best_lb, best_ub, best (y,z), lambdas."""
    import time as _time
    _t0 = _time.time()
    J, L = inst["J"], inst["L"]
    lmb = [row[:] for row in lmb0] if lmb0 else [[0.0] * L for _ in range(J)]
    spB = SubproblemB(inst, env)
    best_lb, best_ub, best_sol, best_lmb = -math.inf, math.inf, None, None
    mu = mu0
    stall = 0
    for it in range(iters):
        if time_limit is not None and _time.time() - _t0 > time_limit:
            break
        vA, ybar, xin = solve_A(inst, lmb)
        rem = None if time_limit is None else max(2.0, time_limit - (_time.time() - _t0))
        rB = spB.solve(lmb, min(rem, 60.0) if rem else 60.0)
        if rB is None:
            break  # no incumbent within the per-solve cap
        vB, zbar, wout = rB
        # L(lmb) = A + B  (the -lambda term is embedded in the reduced costs)
        lb = vA + vB
        if lb > best_lb + 1e-9:
            best_lb, best_lmb, stall = lb, [row[:] for row in lmb], 0
        else:
            stall += 1
            if stall >= 20:
                mu = max(mu / 2, 1e-3)
                stall = 0
        sol, ub = repair(inst, ybar, zbar, env)
        if ub < best_ub:
            best_ub, best_sol = ub, sol
        # subgradient of L at lmb: s_jl = sum_k w_jkl - sum_i x_ijl
        s = [[wout[j][l] - xin[j][l] for l in range(L)] for j in range(J)]
        norm2 = sum(s[j][l] ** 2 for j in range(J) for l in range(L))
        if norm2 < 1e-12:
            break  # dualized constraints satisfied tightly: L(lmb) is optimal dual value
        # Guard against best_ub = inf (infeasible repair) blowing up the step.
        target = best_ub - lb if math.isfinite(best_ub) else max(abs(lb), 1.0)
        t = mu * max(target, 1e-6) / norm2
        for j in range(J):
            for l in range(L):
                lmb[j][l] = max(0.0, lmb[j][l] + t * s[j][l])
        if verbose and it % 20 == 0:
            print(f"  it={it:4d} L={lb:12.2f} best_lb={best_lb:12.2f} "
                  f"ub={best_ub:12.2f} mu={mu:.3f}")
    return dict(best_lb=best_lb, best_ub=best_ub, best_sol=best_sol, lmb=best_lmb)
