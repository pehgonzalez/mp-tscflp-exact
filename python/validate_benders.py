import sys
import gurobipy as gp
from gurobipy import GRB
from benders_gurobipy import BendersSolver
from model_gurobipy import build_model
from validate_bruteforce import brute_force, gen_instance, routing_cost


def main():
    nseeds = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    env = gp.Env(params={"OutputFlag": 0})
    for seed in range(nseeds):
        inst = gen_instance(seed, I=3, J=4, K=5, L=2)
        m, *_ = build_model(inst, env=env)
        m.Params.OutputFlag = 0
        m.optimize()
        assert m.Status == GRB.OPTIMAL
        mip = m.ObjVal
        bf = brute_force(inst, env)
        ncuts = {}
        for pap in (False, True):
            b = BendersSolver(inst, env=env, papadakos=pap)
            r = b.solve(time_limit=120, output=False)
            ncuts[pap] = r["ncuts"]
            assert r["status"] == GRB.OPTIMAL, f"seed {seed} pap={pap}: not optimal"
            assert abs(r["obj"] - mip) < 1e-5 * max(1, abs(mip)), \
                f"seed {seed} pap={pap}: Benders={r['obj']} != MIP={mip}"
            rc = routing_cost(inst, r["y"], r["z"], env)
            assert rc is not None
            total = rc + sum(f * v for f, v in zip(inst["f"], r["y"])) \
                       + sum(g * v for g, v in zip(inst["g"], r["z"]))
            assert abs(total - r["obj"]) < 1e-5 * max(1, abs(total)), \
                f"seed {seed} pap={pap}: reported obj != re-evaluated cost"
        assert abs(mip - bf) < 1e-6 * max(1, abs(mip))
        print(f"seed {seed}: Benders == MIP == BF == {mip:.2f}  "
              f"(cuts std={ncuts[False]}, pap={ncuts[True]})  OK")
    print("BENDERS VALIDATION PASSED")


if __name__ == "__main__":
    main()
