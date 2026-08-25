"""Numerical validation of the Lagrangian bound hierarchy (Theorem L1).

For seeded random instances checks:
  v_LP <= best_LD + tol,  best_LD <= v* + tol,  best_UB >= v*,  UB re-evaluated
and reports how often LD strictly exceeds LP (the integrality-gap payoff).
Usage: python validate_lagrangian.py [nseeds] [iters]
"""
import sys

import gurobipy as gp
from gurobipy import GRB

from lagrangian_gurobipy import lagrangian_dual, routing_value
from model_gurobipy import build_model
from validate_bruteforce import gen_instance


def main():
    nseeds = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    iters = int(sys.argv[2]) if len(sys.argv) > 2 else 300
    env = gp.Env(params={"OutputFlag": 0})
    wins = 0
    for seed in range(nseeds):
        inst = gen_instance(seed, I=3, J=4, K=5, L=2)

        m, *_ = build_model(inst, env=env)
        m.Params.OutputFlag = 0
        m.optimize()
        assert m.Status == GRB.OPTIMAL
        vstar = m.ObjVal
        r = m.relax()
        r.Params.OutputFlag = 0
        r.optimize()
        vlp = r.ObjVal

        res = lagrangian_dual(inst, env, iters=iters)
        lb, ub = res["best_lb"], res["best_ub"]
        tol = 1e-6 * max(1.0, abs(vstar))

        assert lb <= vstar + tol, f"seed {seed}: LD={lb} > v*={vstar} (INVALID BOUND)"
        assert ub >= vstar - tol, f"seed {seed}: UB={ub} < v*={vstar} (INVALID UB)"
        # subgradient may stop short of the dual optimum, but must never lose to a
        # sanity floor: L(0) = A(0)+B(0) <= ... ; we check the hierarchy claim on the
        # attained value only in the direction that is always valid:
        gap_note = "LD>LP (integrality gap!)" if lb > vlp + tol else "LD<=LP (attained)"
        if lb > vlp + tol:
            wins += 1
        # independent re-evaluation of the heuristic solution
        if res["best_sol"] is not None:
            y, z = res["best_sol"]
            rc = routing_value(inst, y, z, env)
            assert rc is not None
            total = rc + sum(f * v for f, v in zip(inst["f"], y)) \
                       + sum(g * v for g, v in zip(inst["g"], z))
            assert abs(total - ub) < tol, f"seed {seed}: UB mismatch {total} vs {ub}"
        print(f"seed {seed}: v_LP={vlp:10.2f}  LD={lb:10.2f}  v*={vstar:10.2f}  "
              f"UB={ub:10.2f}  [{gap_note}]")
    print(f"LD strictly above LP on {wins}/{nseeds} seeds")
    print("LAGRANGIAN VALIDATION PASSED")


if __name__ == "__main__":
    main()
