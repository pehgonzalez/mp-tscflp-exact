// Lagrangian relaxation of the MP-TSCFLP (reference implementation validated
// in python/lagrangian_gurobipy.py and python/validate_lagrangian.py).
//
// Dualizes flow conservation (4) with lambda_jl >= 0.
//   A(lmb): closed form over plants, O(IJL).
//   B(lmb): depot-side CFLP (z binary), solved by Gurobi; ObjBound keeps every
//           reported L(lmb) a valid lower bound even with MIPGap > 0.
// Subgradient with Polyak step + Lagrangian primal heuristic (greedy (F1)/(F2)
// cover + exact routing via the solver-free FlowEvaluator).
#ifndef MPTSCFL_LAGRANGIAN_HPP
#define MPTSCFL_LAGRANGIAN_HPP

#include <vector>

#include "gurobi_c++.h"
#include "flow_evaluator.hpp"
#include "instance.hpp"

namespace mptscfl {

struct LagrangianResult {
    double best_lb = -1e100;       // valid dual bound v_LD attained
    double best_ub = 1e100;        // certified primal cost (FlowEvaluator-routed)
    bool has_solution = false;
    std::vector<int> y, z;         // best primal solution
    std::vector<std::vector<double>> lmb; // multipliers of the best bound [j][l]
    int iterations = 0;
    double runtime = 0.0;
};

class LagrangianSolver {
public:
    explicit LagrangianSolver(const Instance& inst);

    // Subgradient ascent; stops at iters, time_limit (s) or vanishing subgradient.
    // warm=false starts the multipliers at zero instead of the LP duals, the
    // ablation arm of the rebuilt warm-start study.
    LagrangianResult solve(int iters, double time_limit, bool verbose = true,
                           bool warm = true);

private:
    struct AResult {
        double value = 0.0;
        std::vector<int> ybar;
        std::vector<std::vector<double>> xin; // [j][l] stage-1 inflow at A's optimum
    };
    AResult solve_A(const std::vector<std::vector<double>>& lmb) const;

    // Warm start: lambda0 = LP-relaxation duals of the flow-conservation rows
    // makes the FIRST evaluation A+B_MIP already >= v_LP. Returns lmb (empty on
    // failure/timeout; caller falls back to zeros).
    std::vector<std::vector<double>> lp_dual_warmstart(double time_limit) const;

    // Greedy (F1)/(F2) cover + exact routing; returns certified cost (or +inf).
    double repair(std::vector<int> y, std::vector<int> z,
                  std::vector<int>& out_y, std::vector<int>& out_z) const;

    void build_B();

    const Instance& p_;
    FlowEvaluator eval_;
    GRBEnv env_;
    GRBModel B_;
    std::vector<GRBVar> zB_;
    std::vector<std::vector<std::vector<GRBVar>>> wB_; // [j][k][l]
};

} // namespace mptscfl
#endif
