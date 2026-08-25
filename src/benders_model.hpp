// Branch-and-Benders-cut for the MP-TSCFLP (reference implementation validated
// by brute force: python/benders_gurobipy.py and python/validate_benders.py).
//
// Master in (y, z, theta_l), one transportation-LP subproblem per product,
// disaggregated optimality cuts as lazy constraints, a-priori aggregate capacity
// inequalities (no feasibility cuts needed), theta_l >= v_l(all-open)
// initial bounds, optional Papadakos core-point cuts and root user cuts.
#ifndef MPTSCFL_BENDERS_MODEL_HPP
#define MPTSCFL_BENDERS_MODEL_HPP

#include <memory>
#include <mutex>
#include <string>
#include <vector>

#include "gurobi_c++.h"
#include "exact_model.hpp" // ExactResult, Status
#include "instance.hpp"

// Defined once BendersOptions carries theta_integer; main.cpp tests it so the
// same main.cpp compiles before and after this patch.
#define MPTSCFL_HAS_THETA_INTEGER 1

namespace mptscfl {

struct BendersCut {
    double value = 0.0;             // v_l(y,z)
    double constant = 0.0;          // sum_k q_kl * alpha_k
    std::vector<double> fac_coef;   // b_il * Pi_fcap_i   (<= 0)
    std::vector<double> ware_coef;  // p_jl * Pi_wcap_j   (<= 0)
    bool feasible = false;
};

// Um corte de otimalidade guardado para o despejo do certificado. É a
// BendersCut mais o produto a que ela pertence (-1 no master agregado).
struct CertCut {
    int l = -1;
    double constant = 0.0;
    std::vector<double> fac_coef;
    std::vector<double> ware_coef;
};

// Transportation LP of one product; capacities enter as mutable RHS, so the same
// model is re-solved (dual simplex warm start) at every separation call.
class ProductSubproblem {
public:
    ProductSubproblem(const Instance& inst, int l, GRBEnv& env);
    BendersCut solve(const std::vector<double>& yv, const std::vector<double>& zv);

private:
    const Instance& p_;
    int l_;
    GRBModel lp_;
    std::vector<GRBConstr> dem_, fcap_, wcap_;
};

struct BendersOptions {
    bool papadakos = false;  // extra core-point cuts (validity unconditional, sec. 5)
    bool root_cuts = true;   // user cuts at the root LP
    bool aggregated = false; // ablation: a single theta for the whole routing stage
                             // and one summed cut per round (classic single-cut
                             // master), so the disaggregation claim of Section
                             // 4.1.2 is measured rather than asserted
    bool theta_integer = false; // integral-theta master (mode "exact-thetaint"):
                             // theta_l is declared GRB_INTEGER, the master runs
                             // at MIPGapAbs 0.9999 and the Lagrangian bound is
                             // transferred as ceil(v_LD) - 1e-6. Sound because
                             // every v_l(y,z) is a min-cost flow over integral
                             // data, hence integral at every feasible (y,z)
    int threads = 0;
    double eps = 1e-6;       // relative cut-violation tolerance
};

class BendersModel {
public:
    BendersModel(const Instance& inst, BendersOptions opt = {});

    void set_start(const std::vector<int>& ybar, const std::vector<int>& zbar);
    void set_log(const std::string& path);
    void set_seed(int seed);

    // Valid by weak duality (Theorem L1): fixed costs + sum_l theta_l >= v_LD.
    // With opt_.theta_integer the right-hand side is strengthened to
    // ceil(v_LD) - 1e-6, valid because the objective is integral.
    void add_global_lower_bound(double v_ld);

    ExactResult run(double time_limit, double cutoff = -1.0);

    long long cuts_added() const { return ncuts_; }

private:
    friend class BendersCallback;

    // Certificado racional (item R2). Só liga com MPTSCFL_CERT_DUMP setada.
    bool cert_dump_ = false;          // = getenv("MPTSCFL_CERT_DUMP") != nullptr
    std::vector<CertCut> cert_cuts_;  // cortes na ordem em que foram criados
    std::mutex cert_mtx_;             // callbacks do Gurobi podem vir de threads
    double v_ld_ = -1.0;              // guardado por add_global_lower_bound
    void cert_record(int l, const BendersCut& c);
    void cert_write(const ExactResult& res) const;

    const Instance& p_;
    BendersOptions opt_;
    GRBEnv env_;
    GRBModel master_;
    std::vector<GRBVar> y_, z_, theta_;
    std::vector<std::unique_ptr<ProductSubproblem>> sp_;
    std::vector<double> core_y_, core_z_; // Papadakos core point
    long long ncuts_ = 0;

    // theta-repair: incumbent accepted with theta_l > v_l gets re-injected at the
    // next MIPNODE with theta_l := v_l exactly (reported obj = true cost).
    bool has_pending_ = false;
    std::vector<double> pend_y_, pend_z_, pend_theta_;

    // set by the callback on any exception: poisons the result
    bool callback_failed_ = false;
};

} // namespace mptscfl
#endif
