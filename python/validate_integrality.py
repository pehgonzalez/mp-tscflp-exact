"""Validation of the integrality of the optimal value (integer data).

For seeded random instances with integer data, enumerates all (y, z),
solves the routing LP per opening vector, and asserts that every finite
routing cost and the overall optimal value are integers up to LP tolerance.
This is the property that turns an absolute optimality gap below one into
an exact certificate.
"""
import sys

import gurobipy as gp

from validate_bruteforce import gen_instance, routing_cost


def main():
    nseeds = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    env = gp.Env(params={"OutputFlag": 0})
    import itertools
    for seed in range(nseeds):
        inst = gen_instance(seed)
        best = None
        checked = 0
        for yv in itertools.product((0, 1), repeat=inst["I"]):
            for zv in itertools.product((0, 1), repeat=inst["J"]):
                rc = routing_cost(inst, yv, zv, env)
                if rc is None:
                    continue
                assert abs(rc - round(rc)) < 1e-6, \
                    f"seed {seed}: non-integral routing cost {rc} at y={yv} z={zv}"
                checked += 1
                tot = rc + sum(f * v for f, v in zip(inst["f"], yv)) \
                         + sum(g * v for g, v in zip(inst["g"], zv))
                if best is None or tot < best:
                    best = tot
        assert best is not None and abs(best - round(best)) < 1e-6, \
            f"seed {seed}: non-integral optimum {best}"
        print(f"seed {seed}: {checked} routing values integral, optimum {best:.0f} integral  OK")
    print("INTEGRALITY VALIDATION PASSED")


if __name__ == "__main__":
    main()
