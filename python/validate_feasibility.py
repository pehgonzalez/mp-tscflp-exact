"""Validation of the a priori feasibility characterization.

For seeded random instances and every opening vector (y, z), checks that the
routing subproblem of every product is feasible if and only if, for each
product l, the total open plant capacity and the total open depot capacity
both cover the total demand D_l. Feasibility is decided by an independent
routing LP, so the equivalence is tested against a formulation that knows
nothing about the characterization.
"""
import itertools
import sys

import gurobipy as gp

from validate_bruteforce import gen_instance, routing_cost


def condition_F(inst, yv, zv):
    for l in range(inst["L"]):
        D = sum(inst["q"][k][l] for k in range(inst["K"]))
        if sum(inst["b"][i][l] * yv[i] for i in range(inst["I"])) < D:
            return False
        if sum(inst["p"][j][l] * zv[j] for j in range(inst["J"])) < D:
            return False
    return True


def main():
    nseeds = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    env = gp.Env(params={"OutputFlag": 0})
    for seed in range(nseeds):
        inst = gen_instance(seed)
        agree = 0
        for yv in itertools.product((0, 1), repeat=inst["I"]):
            for zv in itertools.product((0, 1), repeat=inst["J"]):
                lp_feasible = routing_cost(inst, yv, zv, env) is not None
                cond = condition_F(inst, yv, zv)
                assert lp_feasible == cond, \
                    f"seed {seed}: F-condition {cond} but LP feasible {lp_feasible} at y={yv} z={zv}"
                agree += 1
        print(f"seed {seed}: {agree} opening vectors, characterization exact  OK")
    print("FEASIBILITY VALIDATION PASSED")


if __name__ == "__main__":
    main()
