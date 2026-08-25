// Two-steps decomposition heuristic + exact warm start.
//   per product l: solve stage 2 (warehouses->customers) as a 1-stage FLP; the flow
//   leaving each warehouse becomes the demand of stage 1 (factories->warehouses);
//   union of opened facilities across products -> route all products (min-cost flow)
//   -> close unused facilities -> heuristic solution = MIP start + cutoff for the exact model.
#ifndef MPTSCFL_TWO_STEPS_SOLVER_HPP
#define MPTSCFL_TWO_STEPS_SOLVER_HPP

#include <vector>

#include "exact_model.hpp"
#include "flow_evaluator.hpp"
#include "instance.hpp"
#include "single_stage_model.hpp"

namespace mptscfl {

struct TwoStepsResult {
    bool heuristic_feasible = false;
    double heuristic_cost = 0.0;
    std::vector<int> y, z;   // heuristic (merged, CNUF-filtered) solution
    ExactResult exact;       // result of the exact phase warm-started with it
};

class TwoStepsSolver {
public:
    explicit TwoStepsSolver(const Instance& inst) : p_(inst), eval_(inst) {}

    // time_limit applies to the exact phase; method as in ExactModel::run.
    // logfile, if non-empty, receives the Gurobi log of the exact phase.
    TwoStepsResult solve(double time_limit, int method = 0, bool run_exact = true,
                         const std::string& logfile = "");

private:
    const Instance& p_;
    FlowEvaluator eval_;
};

} // namespace mptscfl
#endif
