#include "benders_model.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

#include "flow_evaluator.hpp"
#include "param_inject.hpp"

namespace mptscfl {

// ---------------------------------------------------------------- subproblem
ProductSubproblem::ProductSubproblem(const Instance& inst, int l, GRBEnv& env)
    : p_(inst), l_(l), lp_(env) {
    const int I = p_.nfactories, J = p_.nwarehouses, K = p_.ncustomers;
    lp_.set(GRB_IntParam_OutputFlag, 0);
    apply_env_params(lp_); // certificate block: the cut coefficients come from
                           // these duals, so the tightened regime reaches them too

    std::vector<std::vector<GRBVar>> x(I, std::vector<GRBVar>(J));
    std::vector<std::vector<GRBVar>> w(J, std::vector<GRBVar>(K));
    for (int i = 0; i < I; ++i)
        for (int j = 0; j < J; ++j)
            x[i][j] = lp_.addVar(0, GRB_INFINITY, p_.flowcost_fw[l_][i][j], GRB_CONTINUOUS,
                                 "x");
    for (int j = 0; j < J; ++j)
        for (int k = 0; k < K; ++k)
            w[j][k] = lp_.addVar(0, GRB_INFINITY, p_.flowcost_wc[l_][j][k], GRB_CONTINUOUS,
                                 "w");
    lp_.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

    dem_.reserve(K);
    for (int k = 0; k < K; ++k) {
        GRBLinExpr e;
        for (int j = 0; j < J; ++j) e += w[j][k];
        dem_.push_back(lp_.addConstr(e >= p_.customer_demand[k][l_], "dem"));
    }
    for (int j = 0; j < J; ++j) {
        GRBLinExpr in, out;
        for (int i = 0; i < I; ++i) in += x[i][j];
        for (int k = 0; k < K; ++k) out += w[j][k];
        lp_.addConstr(in - out >= 0.0, "cons");
    }
    fcap_.reserve(I);
    for (int i = 0; i < I; ++i) {
        GRBLinExpr e;
        for (int j = 0; j < J; ++j) e += x[i][j];
        fcap_.push_back(lp_.addConstr(e <= 0.0, "fcap")); // RHS set per call
    }
    wcap_.reserve(J);
    for (int j = 0; j < J; ++j) {
        GRBLinExpr e;
        for (int k = 0; k < K; ++k) e += w[j][k];
        wcap_.push_back(lp_.addConstr(e <= 0.0, "wcap"));
    }
}

BendersCut ProductSubproblem::solve(const std::vector<double>& yv,
                                    const std::vector<double>& zv) {
    const int I = p_.nfactories, J = p_.nwarehouses, K = p_.ncustomers;
    for (int i = 0; i < I; ++i)
        fcap_[i].set(GRB_DoubleAttr_RHS, p_.factory_capacity[i][l_] * yv[i]);
    for (int j = 0; j < J; ++j)
        wcap_[j].set(GRB_DoubleAttr_RHS, p_.warehouse_capacity[j][l_] * zv[j]);
    lp_.optimize();

    BendersCut cut;
    if (lp_.get(GRB_IntAttr_Status) != GRB_OPTIMAL) return cut; // infeasible safeguard
    cut.feasible = true;
    cut.value = lp_.get(GRB_DoubleAttr_ObjVal);
    for (int k = 0; k < K; ++k)
        cut.constant += p_.customer_demand[k][l_] * dem_[k].get(GRB_DoubleAttr_Pi);
    cut.fac_coef.resize(I);
    for (int i = 0; i < I; ++i)
        cut.fac_coef[i] = p_.factory_capacity[i][l_] * fcap_[i].get(GRB_DoubleAttr_Pi);
    cut.ware_coef.resize(J);
    for (int j = 0; j < J; ++j)
        cut.ware_coef[j] = p_.warehouse_capacity[j][l_] * wcap_[j].get(GRB_DoubleAttr_Pi);
    return cut;
}

// ----------------------------------------------------------------- callback
class BendersCallback : public GRBCallback {
public:
    explicit BendersCallback(BendersModel& o) : o_(o) {}

protected:
    void callback() override {
        // A swallowed exception would silently ACCEPT the incumbent (no lazy
        // cut added). Any failure must abort the solve and poison the result.
        try {
            if (where == GRB_CB_MIPSOL) {
                separate(/*integer=*/true);
            } else if (where == GRB_CB_MIPNODE &&
                       getIntInfo(GRB_CB_MIPNODE_STATUS) == GRB_OPTIMAL) {
                inject_pending(); // theta-repair (any node)
                if (o_.opt_.root_cuts && getDoubleInfo(GRB_CB_MIPNODE_NODCNT) < 0.5)
                    separate(/*integer=*/false);
            }
        } catch (GRBException& e) {
            std::cerr << "[Benders callback] FATAL: " << e.getMessage() << "\n";
            o_.callback_failed_ = true;
            abort();
        } catch (...) {
            std::cerr << "[Benders callback] FATAL: unknown exception\n";
            o_.callback_failed_ = true;
            abort();
        }
    }

private:
    // Sparsified cut with validity-preserving compensation:
    // coefficients are <= 0, so moving a dropped term into the constant (its value
    // at y=1) only weakens the cut: coef*y >= coef for y in [0,1].
    GRBLinExpr cut_expr(const BendersCut& c) {
        double cst = c.constant;
        GRBLinExpr e;
        for (int i = 0; i < (int)c.fac_coef.size(); ++i) {
            if (std::abs(c.fac_coef[i]) > o_.opt_.eps) e += c.fac_coef[i] * o_.y_[i];
            else cst += c.fac_coef[i];
        }
        for (int j = 0; j < (int)c.ware_coef.size(); ++j) {
            if (std::abs(c.ware_coef[j]) > o_.opt_.eps) e += c.ware_coef[j] * o_.z_[j];
            else cst += c.ware_coef[j];
        }
        return e + cst;
    }

    // Certificate bookkeeping: the dump must carry the cut exactly AS ADDED.
    // cut_expr() folds every coefficient below eps into the constant, and the
    // aggregated master posts the sum of the per-product expressions as a
    // single cut, so the same folding and the same summation happen here.
    void accumulate_as_added(const BendersCut& c, BendersCut& acc) {
        if (acc.fac_coef.size() < c.fac_coef.size())
            acc.fac_coef.resize(c.fac_coef.size(), 0.0);
        if (acc.ware_coef.size() < c.ware_coef.size())
            acc.ware_coef.resize(c.ware_coef.size(), 0.0);
        acc.constant += c.constant;
        for (int i = 0; i < (int)c.fac_coef.size(); ++i) {
            if (std::abs(c.fac_coef[i]) > o_.opt_.eps) acc.fac_coef[i] += c.fac_coef[i];
            else acc.constant += c.fac_coef[i];
        }
        for (int j = 0; j < (int)c.ware_coef.size(); ++j) {
            if (std::abs(c.ware_coef[j]) > o_.opt_.eps) acc.ware_coef[j] += c.ware_coef[j];
            else acc.constant += c.ware_coef[j];
        }
    }

    // theta-repair injection: post the stored solution with theta_l = v_l exact.
    void inject_pending() {
        if (!o_.has_pending_) return;
        const int I = o_.p_.nfactories, J = o_.p_.nwarehouses;
        for (int i = 0; i < I; ++i) setSolution(o_.y_[i], o_.pend_y_[i]);
        for (int j = 0; j < J; ++j) setSolution(o_.z_[j], o_.pend_z_[j]);
        for (int l = 0; l < (int)o_.theta_.size(); ++l)
            setSolution(o_.theta_[l], o_.pend_theta_[l]);
        useSolution();
        o_.has_pending_ = false;
    }

    // The aggregated master carries one theta for the whole routing stage, so
    // a round solves every product, sums the values and the cut expressions,
    // and adds at most one optimality cut. The tolerance aggregates the same
    // way, the per-product cap 0.4/L summing to at most 0.4, so the
    // certificate arithmetic of Section 4.1.3 is unchanged.
    void separate_aggregated(bool integer) {
        const int I = o_.p_.nfactories, J = o_.p_.nwarehouses, L = o_.p_.ncommodities;
        std::vector<double> yv(I), zv(J);
        for (int i = 0; i < I; ++i)
            yv[i] = integer ? std::round(getSolution(o_.y_[i])) : getNodeRel(o_.y_[i]);
        for (int j = 0; j < J; ++j)
            zv[j] = integer ? std::round(getSolution(o_.z_[j])) : getNodeRel(o_.z_[j]);
        const double tva = integer ? getSolution(o_.theta_[0]) : getNodeRel(o_.theta_[0]);

        std::vector<BendersCut> cs(L);
        double total = 0.0;
        for (int l = 0; l < L; ++l) {
            cs[l] = o_.sp_[l]->solve(yv, zv);
            if (!cs[l].feasible) { // impossible under (F1),(F2); combinatorial safeguard
                std::cerr << "[benders] combinatorial safeguard fired for product "
                          << l << ": run not covered by the exact certificate\n";
                if (integer) {
                    GRBLinExpr e;
                    for (int i = 0; i < I; ++i) if (yv[i] < 0.5) e += o_.y_[i];
                    for (int j = 0; j < J; ++j) if (zv[j] < 0.5) e += o_.z_[j];
                    addLazy(e >= 1.0);
                    ++o_.ncuts_;
                }
                return; // no aggregate value exists this round
            }
            total += cs[l].value;
        }
        const double tol = std::min(o_.opt_.eps * std::max(1.0, total), 0.4);
        if (tva < total - tol) {
            GRBLinExpr agg;
            for (int l = 0; l < L; ++l) agg += cut_expr(cs[l]);
            if (integer) addLazy(o_.theta_[0] >= agg);
            else addCut(o_.theta_[0] >= agg);
            ++o_.ncuts_;
            if (o_.cert_dump_) { // one CUT record with l = -1: the summed cut
                BendersCut rec;
                for (int l = 0; l < L; ++l) accumulate_as_added(cs[l], rec);
                o_.cert_record(-1, rec);
            }
            if (integer && o_.opt_.papadakos) {
                GRBLinExpr core;
                BendersCut core_rec;
                bool ok = true;
                for (int l = 0; l < L; ++l) {
                    BendersCut pc = o_.sp_[l]->solve(o_.core_y_, o_.core_z_);
                    if (!pc.feasible) { ok = false; break; }
                    core += cut_expr(pc);
                    if (o_.cert_dump_) accumulate_as_added(pc, core_rec);
                }
                if (ok) { addLazy(o_.theta_[0] >= core); ++o_.ncuts_; }
                if (ok && o_.cert_dump_) o_.cert_record(-1, core_rec);
                for (int i = 0; i < I; ++i) o_.core_y_[i] = 0.5 * (o_.core_y_[i] + yv[i]);
                for (int j = 0; j < J; ++j) o_.core_z_[j] = 0.5 * (o_.core_z_[j] + zv[j]);
            }
            return;
        }
        // theta-repair, aggregated: the incumbent stands but theta overstates
        // the routing cost, so the stored solution is re-posted with the sum.
        if (integer && tva - total > o_.opt_.eps) {
            o_.pend_y_.assign(yv.begin(), yv.end());
            o_.pend_z_.assign(zv.begin(), zv.end());
            o_.pend_theta_.assign(1, total);
            o_.has_pending_ = true;
        }
    }

    void separate(bool integer) {
        if (o_.opt_.aggregated) { separate_aggregated(integer); return; }
        const int I = o_.p_.nfactories, J = o_.p_.nwarehouses, L = o_.p_.ncommodities;
        std::vector<double> yv(I), zv(J), tv(L);
        for (int i = 0; i < I; ++i)
            yv[i] = integer ? std::round(getSolution(o_.y_[i])) : getNodeRel(o_.y_[i]);
        for (int j = 0; j < J; ++j)
            zv[j] = integer ? std::round(getSolution(o_.z_[j])) : getNodeRel(o_.z_[j]);
        for (int l = 0; l < L; ++l)
            tv[l] = integer ? getSolution(o_.theta_[l]) : getNodeRel(o_.theta_[l]);

        std::vector<double> vl(L, -1.0); // exact subproblem values (for theta-repair)
        bool cut_added = false;
        for (int l = 0; l < L; ++l) {
            BendersCut c = o_.sp_[l]->solve(yv, zv);
            if (!c.feasible) { // impossible under (F1),(F2); combinatorial safeguard
                std::cerr << "[benders] combinatorial safeguard fired for product "
                          << l << ": run not covered by the exact certificate\n";
                if (integer) {
                    GRBLinExpr e;
                    for (int i = 0; i < I; ++i) if (yv[i] < 0.5) e += o_.y_[i];
                    for (int j = 0; j < J; ++j) if (zv[j] < 0.5) e += o_.z_[j];
                    addLazy(e >= 1.0);
                    ++o_.ncuts_;
                    cut_added = true; // never schedule theta-repair here
                }
                continue;
            }
            vl[l] = c.value;
            // A relative tolerance alone would allow theta to sit up to
            // eps*|v_l| BELOW v_l; summed over products this can exceed the
            // absolute-gap certificate. Cap the tolerance so L*tol < 0.5.
            const double tol =
                std::min(o_.opt_.eps * std::max(1.0, std::abs(c.value)), 0.4 / L);
            if (tv[l] < c.value - tol) {
                cut_added = true;
                if (integer) addLazy(o_.theta_[l] >= cut_expr(c));
                else addCut(o_.theta_[l] >= cut_expr(c));
                ++o_.ncuts_;
                if (o_.cert_dump_) {
                    BendersCut rec;
                    accumulate_as_added(c, rec);
                    o_.cert_record(l, rec);
                }
                if (integer && o_.opt_.papadakos) {
                    BendersCut pc = o_.sp_[l]->solve(o_.core_y_, o_.core_z_);
                    if (pc.feasible) {
                        addLazy(o_.theta_[l] >= cut_expr(pc));
                        ++o_.ncuts_;
                        if (o_.cert_dump_) {
                            BendersCut rec;
                            accumulate_as_added(pc, rec);
                            o_.cert_record(l, rec);
                        }
                    }
                }
            }
        }
        // The Papadakos core point is updated once per round, not per product.
        if (integer && o_.opt_.papadakos && cut_added) {
            for (int i = 0; i < I; ++i) o_.core_y_[i] = 0.5 * (o_.core_y_[i] + yv[i]);
            for (int j = 0; j < J; ++j) o_.core_z_[j] = 0.5 * (o_.core_z_[j] + zv[j]);
        }
        // theta-repair: incumbent will be accepted (no cut) but some theta_l
        // overestimates v_l -> schedule re-injection with exact values.
        if (integer && !cut_added) {
            double slack = 0.0;
            for (int l = 0; l < L; ++l)
                if (vl[l] >= 0.0) slack += std::max(0.0, tv[l] - vl[l]);
            if (slack > o_.opt_.eps) {
                o_.pend_y_.assign(yv.begin(), yv.end());
                o_.pend_z_.assign(zv.begin(), zv.end());
                o_.pend_theta_.assign(vl.begin(), vl.end());
                o_.has_pending_ = true;
            }
        }
    }

    BendersModel& o_;
};

// -------------------------------------------------------------------- master
BendersModel::BendersModel(const Instance& inst, BendersOptions opt)
    : p_(inst), opt_(opt), env_(true), master_((env_.start(), env_)) {
    const int I = p_.nfactories, J = p_.nwarehouses, L = p_.ncommodities;

    sp_.reserve(L);
    for (int l = 0; l < L; ++l) sp_.push_back(std::make_unique<ProductSubproblem>(p_, l, env_));

    for (int i = 0; i < I; ++i) {
        std::ostringstream n; n << "y[" << i << "]";
        y_.push_back(master_.addVar(0, 1, p_.fixedcost_factory[i], GRB_BINARY, n.str()));
        y_.back().set(GRB_IntAttr_BranchPriority, 10);
    }
    for (int j = 0; j < J; ++j) {
        std::ostringstream n; n << "z[" << j << "]";
        z_.push_back(master_.addVar(0, 1, p_.fixedcost_warehouse[j], GRB_BINARY, n.str()));
        z_.back().set(GRB_IntAttr_BranchPriority, 10);
    }
    // Aggregated ablation: one theta stands for the whole routing stage. The
    // rest of the class indexes theta_ by its own size, so every downstream
    // loop reads correctly in either shape.
    const int nth = opt_.aggregated ? 1 : L;
    for (int l = 0; l < nth; ++l) {
        std::ostringstream n; n << "theta[" << l << "]";
        // Integral-theta master (opt_.theta_integer): v_l(y,z) is a min-cost flow
        // over integral data, so theta_l = v_l at every optimum and restricting
        // theta_l to Z cuts no optimal solution while rounding up each node LP.
        theta_.push_back(master_.addVar(0, GRB_INFINITY, 1.0,
                                        opt_.theta_integer ? GRB_INTEGER : GRB_CONTINUOUS,
                                        n.str()));
    }
    master_.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

    // (F1_l),(F2_l): a-priori feasibility.
    for (int l = 0; l < L; ++l) {
        const double D = p_.total_demand(l);
        GRBLinExpr ef, ew;
        for (int i = 0; i < I; ++i) ef += p_.factory_capacity[i][l] * y_[i];
        for (int j = 0; j < J; ++j) ew += p_.warehouse_capacity[j][l] * z_[j];
        std::ostringstream nf, nw;
        nf << "F1[" << l << "]"; nw << "F2[" << l << "]";
        master_.addConstr(ef >= D, nf.str());
        master_.addConstr(ew >= D, nw.str());
    }
    // theta_l >= v_l(all-open): valid initial bounds by capacity monotonicity.
    // Aggregated: the sum of the same values bounds the single theta, and one
    // infeasible product leaves the bound off exactly as it leaves that
    // product's bound off in the disaggregated master.
    std::vector<double> ones_y(I, 1.0), ones_z(J, 1.0);
    // Integral costs make the true theta integral, so rounding the bound up keeps
    // it valid; the -1e-6 protects a bound that is already integral but came out
    // of the LP represented from above.
    auto theta_lb = [&](double lb) { return opt_.theta_integer ? std::ceil(lb - 1e-6) : lb; };
    if (opt_.aggregated) {
        double s = 0.0;
        bool all_ok = true;
        for (int l = 0; l < L; ++l) {
            BendersCut c = sp_[l]->solve(ones_y, ones_z);
            if (c.feasible) s += c.value; else all_ok = false;
        }
        if (all_ok) theta_[0].set(GRB_DoubleAttr_LB, theta_lb(s));
    } else {
        for (int l = 0; l < L; ++l) {
            BendersCut c = sp_[l]->solve(ones_y, ones_z);
            if (c.feasible) theta_[l].set(GRB_DoubleAttr_LB, theta_lb(c.value));
        }
    }

    core_y_.assign(I, 1.0);
    core_z_.assign(J, 1.0);
    master_.set(GRB_IntParam_LazyConstraints, 1);
    if (opt_.root_cuts) master_.set(GRB_IntParam_PreCrush, 1);
    if (opt_.threads > 0) master_.set(GRB_IntParam_Threads, opt_.threads);
    apply_env_params(master_);

    const char* cd = std::getenv("MPTSCFL_CERT_DUMP");
    cert_dump_ = (cd != nullptr && *cd != '\0');
    if (cert_dump_)
        std::cout << "[cert] certificate dump enabled, target directory " << cd << "\n";
}

void BendersModel::set_log(const std::string& path) {
    master_.set(GRB_StringParam_LogFile, path);
    // The full solver log lives in the file; the console keeps only the
    // campaign's own progress lines.
    master_.set(GRB_IntParam_LogToConsole, 0);
}

void BendersModel::set_seed(int seed) {
    master_.set(GRB_IntParam_Seed, seed);
    // subproblem LPs are deterministic re-solves; master seed governs the search
}

void BendersModel::add_global_lower_bound(double v_ld) {
    v_ld_ = v_ld;
    // Integral objective: the Lagrangian bound may be lifted to the next integer.
    // The 1e-6 margin absorbs the floating error of v_LD itself, which comes out
    // of a subgradient loop and can sit a hair above the true value.
    const double rhs = opt_.theta_integer ? std::ceil(v_ld) - 1e-6 : v_ld;
    GRBLinExpr e;
    for (int i = 0; i < p_.nfactories; ++i) e += p_.fixedcost_factory[i] * y_[i];
    for (int j = 0; j < p_.nwarehouses; ++j) e += p_.fixedcost_warehouse[j] * z_[j];
    // Over theta_ by its own size, not by product count: the aggregated master
    // carries a single theta standing for the whole sum, and the inequality
    // reads identically in either shape.
    for (size_t l = 0; l < theta_.size(); ++l) e += theta_[l];
    master_.addConstr(e >= rhs, "LD_bound");
}

void BendersModel::set_start(const std::vector<int>& ybar, const std::vector<int>& zbar) {
    for (int i = 0; i < p_.nfactories; ++i) y_[i].set(GRB_DoubleAttr_Start, ybar[i]);
    for (int j = 0; j < p_.nwarehouses; ++j) z_[j].set(GRB_DoubleAttr_Start, zbar[j]);
}

void BendersModel::cert_record(int l, const BendersCut& c) {
    if (!cert_dump_) return;
    std::lock_guard<std::mutex> lk(cert_mtx_);
    CertCut cc;
    cc.l = l;
    cc.constant = c.constant;
    cc.fac_coef = c.fac_coef;
    cc.ware_coef = c.ware_coef;
    cert_cuts_.push_back(std::move(cc));
}

void BendersModel::cert_write(const ExactResult& res) const {
    const char* dir = std::getenv("MPTSCFL_CERT_DUMP");
    if (dir == nullptr || *dir == '\0') return;
    try {
        std::filesystem::create_directories(dir);
        // id = stem da instancia + carimbo de tempo: unico por run, e casa com
        // o nome do log do Gurobi na inspecao manual.
        const std::string stem = std::filesystem::path(p_.file_name).stem().string();
        std::time_t tt = std::time(nullptr);
        char ts[32];
        std::strftime(ts, sizeof(ts), "%Y%m%d-%H%M%S", std::localtime(&tt));
        const std::string path =
            (std::filesystem::path(dir) / (stem + "_" + ts + ".certdump")).string();
        std::ofstream o(path);
        if (!o) { std::cout << "[cert] cannot open " << path << "\n"; return; }

        auto num = [](double v) {           // %.17g, round-trip exato do double
            char b[40];
            std::snprintf(b, sizeof(b), "%.17g", v);
            return std::string(b);
        };

        o << "# mptscfl certdump v1\n"
          << "# CUT l const F <fac coefs> W <ware coefs> means\n"
          << "#   theta_l >= const + sum_i fac_i y_i + sum_j ware_j z_j\n"
          << "# FL l FAC rhs means sum_i b_il y_i >= rhs (idem WARE with p_jl z_j)\n"
          << "# every real is printed with %.17g; l = -1 marks the aggregated theta\n";
        o << "FORMAT 1\n";
        o << "INSTANCE " << p_.file_name << "\n";
        o << "DIMS " << p_.nfactories << " " << p_.nwarehouses << " "
          << p_.ncustomers << " " << p_.ncommodities << "\n";
        o << "RUN " << 1 << " " << opt_.threads << " " << (opt_.papadakos ? 1 : 0)
          << " " << (opt_.aggregated ? 1 : 0) << " "
          << (opt_.theta_integer ? 1 : 0) << "\n";
        o << "STATUS OPTIMAL\n";
        o << "OBJ " << num(res.obj) << "\n";
        o << "BOUND " << num(res.bound) << "\n";
        o << "VLD " << (v_ld_ > 0 ? num(v_ld_) : std::string("NONE")) << "\n";

        // Custo re-roteado, o mesmo numero que o main.cpp grava em verified_cost.
        std::string verified = "NONE";
        if (!res.y.empty() && !res.z.empty()) {
            FlowEvaluator ev(p_);
            FlowResult fr = ev.evaluate(res.y, res.z);
            if (fr.feasible) verified = num(ev.solution_cost(res.y, res.z, fr));
        }
        o << "VERIFIED " << verified << "\n";

        o << "Y"; for (int v : res.y) o << " " << v; o << "\n";
        o << "Z"; for (int v : res.z) o << " " << v; o << "\n";

        o << "NCUTS " << cert_cuts_.size() << "\n";
        for (const CertCut& c : cert_cuts_) {
            o << "CUT " << c.l << " " << num(c.constant) << " F";
            for (double v : c.fac_coef) o << " " << num(v);
            o << " W";
            for (double v : c.ware_coef) o << " " << num(v);
            o << "\n";
        }

        // As desigualdades (F_l) agregadas: sum_i b_il y_i >= Q_l e
        // sum_j p_jl z_j >= Q_l, com Q_l a demanda total do produto l. Sao as
        // mesmas linhas que o construtor do master ja adiciona a priori.
        for (int l = 0; l < p_.ncommodities; ++l) {
            const std::string q = num(p_.total_demand(l));
            o << "FL " << l << " FAC " << q << "\n";
            o << "FL " << l << " WARE " << q << "\n";
        }
        o << "END\n";
        std::cout << "[cert] certificate written to " << path << " ("
                  << cert_cuts_.size() << " cuts)\n";
    } catch (const std::exception& e) {
        // Um certificado que falha nunca pode derrubar o run.
        std::cout << "[cert] dump failed: " << e.what() << "\n";
    }
}

ExactResult BendersModel::run(double time_limit, double cutoff) {
    master_.set(GRB_DoubleParam_TimeLimit, time_limit);
    // Proof mode: the optimal value is integral. In the Benders master
    // obj = fixed + sum theta, and acceptance allows theta_l >= v_l - tol_l with
    // sum_l tol_l <= 0.4. MIPGapAbs = 0.5 then certifies: bound > obj - 0.5 >=
    // true - 0.9, and integrality closes the argument.
    master_.set(GRB_DoubleParam_MIPGap, 0.0);
    master_.set(GRB_DoubleParam_MIPGapAbs, 0.5);
    // Integral theta makes the whole objective integral: a gap below 1 is a proof.
    if (opt_.theta_integer) master_.set(GRB_DoubleParam_MIPGapAbs, 0.9999);
    master_.set(GRB_DoubleParam_Cutoff, cutoff > 0 ? cutoff : GRB_INFINITY);

    BendersCallback cb(*this);
    master_.setCallback(&cb);
    master_.optimize();

    ExactResult r;
    if (callback_failed_) { // never trust a run whose separation failed
        std::cerr << "[Benders] result POISONED by callback failure; discarding.\n";
        r.status = Status::NotFound;
        return r;
    }
    const int st = master_.get(GRB_IntAttr_Status);
    if (st == GRB_INFEASIBLE) { r.status = Status::Infeasible; return r; }
    if (master_.get(GRB_IntAttr_SolCount) == 0) {
        r.status = Status::NotFound;
        // On CUTOFF (or timeout) without incumbent the dual bound is
        // still valid and lets the caller certify a heuristic UB.
        try { r.bound = master_.get(GRB_DoubleAttr_ObjBound); } catch (GRBException&) {}
        return r;
    }
    r.status = (st == GRB_OPTIMAL) ? Status::OptimalFound : Status::SolutionFound;
    r.obj = master_.get(GRB_DoubleAttr_ObjVal);
    r.bound = master_.get(GRB_DoubleAttr_ObjBound);
    r.gap = master_.get(GRB_DoubleAttr_MIPGap);
    r.runtime = master_.get(GRB_DoubleAttr_Runtime);
    r.y.resize(p_.nfactories);
    r.z.resize(p_.nwarehouses);
    for (int i = 0; i < p_.nfactories; ++i) r.y[i] = y_[i].get(GRB_DoubleAttr_X) > 0.5;
    for (int j = 0; j < p_.nwarehouses; ++j) r.z[j] = z_[j].get(GRB_DoubleAttr_X) > 0.5;
    std::cout << "[Benders] cuts added: " << ncuts_ << "\n";
    if (cert_dump_ && r.status == Status::OptimalFound) cert_write(r);
    return r;
}

} // namespace mptscfl
