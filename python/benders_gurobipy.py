"""Branch-and-Benders-cut for the MP-TSCFLP (gurobipy reference implementation)."""
import sys

import gurobipy as gp
from gurobipy import GRB

from model_gurobipy import load_instance

EPS = 1e-6


class Subproblem:
    def __init__(self, inst, l, env):
        I, J, K = inst["I"], inst["J"], inst["K"]
        self.I, self.J, self.K, self.l = I, J, K, l
        m = gp.Model(f"sp{l}", env=env)
        m.Params.OutputFlag = 0
        x = m.addVars(I, J, lb=0.0)
        w = m.addVars(J, K, lb=0.0)
        m.setObjective(
            gp.quicksum(inst["c"][l][i][j] * x[i, j] for i in range(I) for j in range(J))
            + gp.quicksum(inst["d"][l][j][k] * w[j, k] for j in range(J) for k in range(K)))
        self.dem = [m.addConstr(w.sum("*", k) >= inst["q"][k][l]) for k in range(K)]
        self.cons = [m.addConstr(x.sum("*", j) - w.sum(j, "*") >= 0.0) for j in range(J)]
        self.fcap = [m.addConstr(x.sum(i, "*") <= 0.0) for i in range(I)]
        self.wcap = [m.addConstr(w.sum(j, "*") <= 0.0) for j in range(J)]
        self.b = [inst["b"][i][l] for i in range(I)]
        self.p = [inst["p"][j][l] for j in range(J)]
        self.q = [inst["q"][k][l] for k in range(K)]
        self.m = m

    def solve(self, yv, zv):
        for i in range(self.I):
            self.fcap[i].RHS = self.b[i] * yv[i]
        for j in range(self.J):
            self.wcap[j].RHS = self.p[j] * zv[j]
        self.m.optimize()
        if self.m.Status != GRB.OPTIMAL:
            return None
        const = sum(self.q[k] * self.dem[k].Pi for k in range(self.K))
        fac = [self.b[i] * self.fcap[i].Pi for i in range(self.I)]
        ware = [self.p[j] * self.wcap[j].Pi for j in range(self.J)]
        return self.m.ObjVal, const, fac, ware


class BendersSolver:
    def __init__(self, inst, env=None, papadakos=False, root_cuts=True):
        self.inst = inst
        self.env = env or gp.Env()
        self.papadakos = papadakos
        self.root_cuts = root_cuts
        I, J, L = inst["I"], inst["J"], inst["L"]
        self.sp = [Subproblem(inst, l, self.env) for l in range(L)]
        m = gp.Model("benders_master", env=self.env)
        self.y = m.addVars(I, vtype=GRB.BINARY, name="y")
        self.z = m.addVars(J, vtype=GRB.BINARY, name="z")
        self.theta = m.addVars(L, lb=0.0, name="theta")
        m.setObjective(
            gp.quicksum(inst["f"][i] * self.y[i] for i in range(I))
            + gp.quicksum(inst["g"][j] * self.z[j] for j in range(J))
            + self.theta.sum(), GRB.MINIMIZE)
        for l in range(L):
            D = sum(inst["q"][k][l] for k in range(inst["K"]))
            m.addConstr(gp.quicksum(inst["b"][i][l] * self.y[i] for i in range(I)) >= D)
            m.addConstr(gp.quicksum(inst["p"][j][l] * self.z[j] for j in range(J)) >= D)
        for l in range(L):
            r = self.sp[l].solve([1] * I, [1] * J)
            assert r is not None, "benchmark guarantees all-open feasibility"
            self.theta[l].LB = r[0]
        m.Params.LazyConstraints = 1
        if root_cuts:
            m.Params.PreCrush = 1
        self.m = m
        self.core_y = [1.0] * I
        self.core_z = [1.0] * J
        self.ncuts = 0

    def _add_cut(self, model, l, const, fac, ware, lazy, where=None):
        # Dropped coefficients (<= 0) move into the constant (weakens only).
        dropped = sum(fac[i] for i in range(self.inst["I"]) if abs(fac[i]) <= EPS) \
                + sum(ware[j] for j in range(self.inst["J"]) if abs(ware[j]) <= EPS)
        expr = (const + dropped) \
             + gp.quicksum(fac[i] * self.y[i] for i in range(self.inst["I"]) if abs(fac[i]) > EPS) \
             + gp.quicksum(ware[j] * self.z[j] for j in range(self.inst["J"]) if abs(ware[j]) > EPS)
        if lazy:
            model.cbLazy(self.theta[l] >= expr)
        elif where is not None:
            model.cbCut(self.theta[l] >= expr)
        else:
            model.addConstr(self.theta[l] >= expr)
        self.ncuts += 1

    def _callback(self, model, where):
        if where == GRB.Callback.MIPSOL:
            yv = model.cbGetSolution(self.y)
            zv = model.cbGetSolution(self.z)
            tv = model.cbGetSolution(self.theta)
            yb = [round(yv[i]) for i in range(self.inst["I"])]
            zb = [round(zv[j]) for j in range(self.inst["J"])]
            cut_added = False
            for l in range(self.inst["L"]):
                r = self.sp[l].solve(yb, zb)
                if r is None:
                    model.cbLazy(
                        gp.quicksum(self.y[i] for i in range(self.inst["I"]) if yb[i] == 0)
                        + gp.quicksum(self.z[j] for j in range(self.inst["J"]) if zb[j] == 0)
                        >= 1)
                    cut_added = True
                    continue
                v, const, fac, ware = r
                tol = min(EPS * max(1.0, abs(v)), 0.4 / self.inst["L"])
                if tv[l] < v - tol:
                    cut_added = True
                    self._add_cut(model, l, const, fac, ware, lazy=True)
                    if self.papadakos:
                        rc = self.sp[l].solve(self.core_y, self.core_z)
                        if rc is not None:
                            self._add_cut(model, l, rc[1], rc[2], rc[3], lazy=True)
            # Core point updated once per round, not per product.
            if self.papadakos and cut_added:
                self.core_y = [0.5 * (c + v) for c, v in zip(self.core_y, yb)]
                self.core_z = [0.5 * (c + v) for c, v in zip(self.core_z, zb)]
        elif where == GRB.Callback.MIPNODE and self.root_cuts:
            if model.cbGet(GRB.Callback.MIPNODE_STATUS) != GRB.OPTIMAL:
                return
            if model.cbGet(GRB.Callback.MIPNODE_NODCNT) > 0:
                return
            yv = model.cbGetNodeRel(self.y)
            zv = model.cbGetNodeRel(self.z)
            tv = model.cbGetNodeRel(self.theta)
            for l in range(self.inst["L"]):
                r = self.sp[l].solve([yv[i] for i in range(self.inst["I"])],
                                     [zv[j] for j in range(self.inst["J"])])
                if r is None:
                    continue
                v, const, fac, ware = r
                tol = min(EPS * max(1.0, abs(v)), 0.4 / self.inst["L"])
                if tv[l] < v - tol:
                    self._add_cut(model, l, const, fac, ware, lazy=False, where=where)

    def solve(self, time_limit=3600.0, output=True):
        self.m.Params.TimeLimit = time_limit
        self.m.Params.OutputFlag = 1 if output else 0
        self.m.optimize(lambda model, where: self._callback(model, where))
        res = dict(status=self.m.Status, ncuts=self.ncuts, runtime=self.m.Runtime)
        if self.m.SolCount:
            res.update(obj=self.m.ObjVal, bound=self.m.ObjBound, gap=self.m.MIPGap,
                       y=[int(self.y[i].X > 0.5) for i in range(self.inst["I"])],
                       z=[int(self.z[j].X > 0.5) for j in range(self.inst["J"])])
        return res


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    inst = load_instance(sys.argv[1])
    tl = float(sys.argv[2]) if len(sys.argv) > 2 else 3600.0
    solver = BendersSolver(inst)
    r = solver.solve(tl)
    if "obj" in r:
        print(f"obj={r['obj']:.2f} bound={r['bound']:.2f} gap={r['gap']:.4%} "
              f"cuts={r['ncuts']} time={r['runtime']:.1f}s")


if __name__ == "__main__":
    main()
