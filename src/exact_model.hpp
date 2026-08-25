// Exact MIP for the MP-TSCFLP, the compact model of Mauri et al. (2021),
// solved directly with the Gurobi C++ API (branching priorities, cutoff,
// node files, IIS support). The manual branch-and-Benders-cut lives in
// benders_model.*.
#ifndef MPTSCFL_EXACT_MODEL_HPP
#define MPTSCFL_EXACT_MODEL_HPP

#include <string>
#include <vector>

#include "gurobi_c++.h"
#include "instance.hpp"

namespace mptscfl {

struct ExactResult {
    Status status = Status::NotFound;
    double obj = 0.0;
    double bound = 0.0;
    double gap = 0.0;
    double runtime = 0.0;
    std::vector<int> y; // open factories
    std::vector<int> z; // open warehouses
};

class ExactModel {
public:
    explicit ExactModel(const Instance& inst);

    // MIP start (heuristic solution) for y/z; flows are completed by the solver.
    void set_start(const std::vector<int>& ybar, const std::vector<int>& zbar);

    // Full Gurobi log to file (survives console truncation).
    void set_log(const std::string& path);
    void set_seed(int seed);

    // Proof mode (default true): MIPGap 0, MIPGapAbs 0.9999, the integral-value
    // certificate regime. false restores the solver defaults, MIPGap 1e-4 and
    // MIPGapAbs 1e-10, the protocol that mode "exact-mauri" reproduces.
    void set_proof_mode(bool on) { proof_mode_ = on; }

    // Adds the two aggregate feasibility rows (F_l) of Proposition 2 per
    // product, valid for the compact model; the "exact-fl" baseline of the
    // revision campaign.
    void add_feasibility_rows();

    // Adds the constraint objective >= lb, the Lagrangian bound handed to the
    // compact model in the "exact-mipguided" baseline.
    void add_objective_lower_bound(double lb);

    // method: 0 = plain branch-and-bound (Gurobi default cuts/heuristics on).
    //         1 = reserved for manual Benders (Phase 1); currently falls back to 0
    //             with a warning, since Gurobi has no automatic Benders.
    ExactResult run(double time_limit, double cutoff = -1.0, int method = 0, int threads = 0);

    // Local branching around (ybar, zbar): adds the two Hamming-ball constraints of
    // Fischetti & Lodi (2003) exactly as in Exato::localbranching, solves, removes them.
    ExactResult local_branching(const std::vector<int>& ybar, const std::vector<int>& zbar,
                                int delta1, int delta2, double time_limit,
                                double cutoff = -1.0, int method = 0);

    GRBModel& model() { return model_; }

    // The location variables, so an external GRBCallback (repair_callback.hpp)
    // can read them in MIPSOL and write them back in MIPNODE. Header-only
    // accessors: exact_model.cpp is untouched. Valid from construction on,
    // since build() runs in the constructor.
    const std::vector<GRBVar>& factory_vars() const { return y_; }
    const std::vector<GRBVar>& warehouse_vars() const { return z_; }

private:
    void build();

    const Instance& p_;
    GRBEnv env_;
    GRBModel model_;
    bool proof_mode_ = true;
    std::vector<GRBVar> y_, z_;
    std::vector<std::vector<std::vector<GRBVar>>> x_; // [i][j][l]
    std::vector<std::vector<std::vector<GRBVar>>> w_; // [j][k][l]
};

} // namespace mptscfl
#endif
