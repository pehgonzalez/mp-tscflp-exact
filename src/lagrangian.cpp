#include "lagrangian.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <numeric>

#include "param_inject.hpp"

namespace mptscfl {

LagrangianSolver::LagrangianSolver(const Instance& inst)
    : p_(inst), eval_(inst), env_(true), B_((env_.start(), env_)) {
    build_B();
}

void LagrangianSolver::build_B() {
    const int J = p_.nwarehouses, K = p_.ncustomers, L = p_.ncommodities;
    B_.set(GRB_IntParam_OutputFlag, 0);
    apply_env_params(B_);
    zB_.reserve(J);
    for (int j = 0; j < J; ++j)
        zB_.push_back(B_.addVar(0, 1, p_.fixedcost_warehouse[j], GRB_BINARY, "z"));
    wB_.assign(J, std::vector<std::vector<GRBVar>>(K, std::vector<GRBVar>(L)));
    for (int j = 0; j < J; ++j)
        for (int k = 0; k < K; ++k)
            for (int l = 0; l < L; ++l) // objective coefficient set per iteration
                wB_[j][k][l] = B_.addVar(0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "w");
    B_.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
    for (int k = 0; k < K; ++k)
        for (int l = 0; l < L; ++l) {
            GRBLinExpr e;
            for (int j = 0; j < J; ++j) e += wB_[j][k][l];
            B_.addConstr(e >= p_.customer_demand[k][l], "dem");
        }
    for (int j = 0; j < J; ++j)
        for (int l = 0; l < L; ++l) {
            GRBLinExpr e;
            for (int k = 0; k < K; ++k) e += wB_[j][k][l];
            B_.addConstr(e <= p_.warehouse_capacity[j][l] * zB_[j], "cap");
        }
}

LagrangianSolver::AResult
LagrangianSolver::solve_A(const std::vector<std::vector<double>>& lmb) const {
    const int I = p_.nfactories, J = p_.nwarehouses, L = p_.ncommodities;
    AResult r;
    r.ybar.assign(I, 0);
    r.xin.assign(J, std::vector<double>(L, 0.0));
    for (int i = 0; i < I; ++i) {
        double tot = p_.fixedcost_factory[i];
        std::vector<std::pair<int, int>> moves; // (l, best_j)
        for (int l = 0; l < L; ++l) {
            int bj = 0;
            double brc = p_.flowcost_fw[l][i][0] - lmb[0][l];
            for (int j = 1; j < J; ++j) {
                double rc = p_.flowcost_fw[l][i][j] - lmb[j][l];
                if (rc < brc) { brc = rc; bj = j; }
            }
            if (brc < 0.0) {
                tot += p_.factory_capacity[i][l] * brc;
                moves.push_back({l, bj});
            }
        }
        if (tot < 0.0) {
            r.value += tot;
            r.ybar[i] = 1;
            for (auto [l, j] : moves) r.xin[j][l] += p_.factory_capacity[i][l];
        }
    }
    return r;
}

double LagrangianSolver::repair(std::vector<int> y, std::vector<int> z,
                                std::vector<int>& out_y, std::vector<int>& out_z) const {
    const int I = p_.nfactories, J = p_.nwarehouses, L = p_.ncommodities;
    // Greedy cover of (F1_l), (F2_l) by fixed-cost/capacity ratio.
    for (int l = 0; l < L; ++l) {
        const double D = p_.total_demand(l);
        double capf = 0.0, capw = 0.0;
        for (int i = 0; i < I; ++i) capf += p_.factory_capacity[i][l] * y[i];
        for (int j = 0; j < J; ++j) capw += p_.warehouse_capacity[j][l] * z[j];
        std::vector<int> idx(I);
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), [&](int a, int b) {
            return p_.fixedcost_factory[a] * std::max(1.0, p_.factory_capacity[b][l]) <
                   p_.fixedcost_factory[b] * std::max(1.0, p_.factory_capacity[a][l]);
        });
        for (int i : idx) {
            if (capf >= D) break;
            if (!y[i]) { y[i] = 1; capf += p_.factory_capacity[i][l]; }
        }
        std::vector<int> jdx(J);
        std::iota(jdx.begin(), jdx.end(), 0);
        std::sort(jdx.begin(), jdx.end(), [&](int a, int b) {
            return p_.fixedcost_warehouse[a] * std::max(1.0, p_.warehouse_capacity[b][l]) <
                   p_.fixedcost_warehouse[b] * std::max(1.0, p_.warehouse_capacity[a][l]);
        });
        for (int j : jdx) {
            if (capw >= D) break;
            if (!z[j]) { z[j] = 1; capw += p_.warehouse_capacity[j][l]; }
        }
    }
    FlowResult fr = eval_.evaluate(y, z);
    if (!fr.feasible) return 1e100;
    // CNUF: closing unused facilities only removes fixed charges (routing unchanged).
    FlowResult fr2 = eval_.evaluate(fr.used_factories, fr.used_warehouses);
    if (fr2.feasible) {
        out_y = fr.used_factories;
        out_z = fr.used_warehouses;
        return eval_.solution_cost(out_y, out_z, fr2);
    }
    out_y = y;
    out_z = z;
    return eval_.solution_cost(y, z, fr);
}

std::vector<std::vector<double>>
LagrangianSolver::lp_dual_warmstart(double time_limit) const {
    const int I = p_.nfactories, J = p_.nwarehouses, K = p_.ncustomers,
              L = p_.ncommodities;
    std::vector<std::vector<double>> lmb;
    try {
        GRBModel lp(const_cast<GRBEnv&>(env_));
        lp.set(GRB_IntParam_OutputFlag, 0);
        apply_env_params(lp);
        lp.set(GRB_DoubleParam_TimeLimit, time_limit);
        std::vector<GRBVar> y(I), z(J);
        for (int i = 0; i < I; ++i)
            y[i] = lp.addVar(0, 1, p_.fixedcost_factory[i], GRB_CONTINUOUS, "y");
        for (int j = 0; j < J; ++j)
            z[j] = lp.addVar(0, 1, p_.fixedcost_warehouse[j], GRB_CONTINUOUS, "z");
        std::vector<std::vector<std::vector<GRBVar>>> x(
            I, std::vector<std::vector<GRBVar>>(J, std::vector<GRBVar>(L)));
        std::vector<std::vector<std::vector<GRBVar>>> w(
            J, std::vector<std::vector<GRBVar>>(K, std::vector<GRBVar>(L)));
        for (int i = 0; i < I; ++i)
            for (int j = 0; j < J; ++j)
                for (int l = 0; l < L; ++l)
                    x[i][j][l] = lp.addVar(0, GRB_INFINITY, p_.flowcost_fw[l][i][j],
                                           GRB_CONTINUOUS, "x");
        for (int j = 0; j < J; ++j)
            for (int k = 0; k < K; ++k)
                for (int l = 0; l < L; ++l)
                    w[j][k][l] = lp.addVar(0, GRB_INFINITY, p_.flowcost_wc[l][j][k],
                                           GRB_CONTINUOUS, "w");
        lp.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
        for (int k = 0; k < K; ++k)
            for (int l = 0; l < L; ++l) {
                GRBLinExpr e;
                for (int j = 0; j < J; ++j) e += w[j][k][l];
                lp.addConstr(e >= p_.customer_demand[k][l], "dem");
            }
        std::vector<std::vector<GRBConstr>> cons(J, std::vector<GRBConstr>(L));
        for (int j = 0; j < J; ++j)
            for (int l = 0; l < L; ++l) {
                GRBLinExpr in, out;
                for (int i = 0; i < I; ++i) in += x[i][j][l];
                for (int k = 0; k < K; ++k) out += w[j][k][l];
                cons[j][l] = lp.addConstr(in - out >= 0.0, "cons");
            }
        for (int i = 0; i < I; ++i)
            for (int l = 0; l < L; ++l) {
                GRBLinExpr e;
                for (int j = 0; j < J; ++j) e += x[i][j][l];
                lp.addConstr(e <= p_.factory_capacity[i][l] * y[i], "fcap");
            }
        for (int j = 0; j < J; ++j)
            for (int l = 0; l < L; ++l) {
                GRBLinExpr e;
                for (int k = 0; k < K; ++k) e += w[j][k][l];
                lp.addConstr(e <= p_.warehouse_capacity[j][l] * z[j], "wcap");
            }
        lp.optimize();
        if (lp.get(GRB_IntAttr_Status) != GRB_OPTIMAL) return lmb; // timeout: fallback
        lmb.assign(J, std::vector<double>(L, 0.0));
        for (int j = 0; j < J; ++j)
            for (int l = 0; l < L; ++l)
                lmb[j][l] = std::max(0.0, cons[j][l].get(GRB_DoubleAttr_Pi));
        std::cout << "[Lagrangian] LP-dual warm start: v_LP = "
                  << lp.get(GRB_DoubleAttr_ObjVal) << "\n";
    } catch (GRBException& e) {
        std::cerr << "[Lagrangian] LP warm start failed: " << e.getMessage() << "\n";
        lmb.clear();
    }
    return lmb;
}

LagrangianResult LagrangianSolver::solve(int iters, double time_limit, bool verbose,
                                         bool warm) {
    const int J = p_.nwarehouses, K = p_.ncustomers, L = p_.ncommodities;
    LagrangianResult res;
    const auto t0 = std::chrono::steady_clock::now();
    auto elapsed = [&] {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    };
    // Start from the LP-relaxation duals: the first evaluation already reaches v_LP.
    // With warm=false the multipliers start at zero, the cold arm of the ablation.
    std::vector<std::vector<double>> lmb;
    if (warm) {
        lmb = lp_dual_warmstart(std::min(0.4 * time_limit, 120.0));
    } else if (verbose) {
        std::cout << "[Lagrangian] warm start disabled: multipliers start at zero\n";
    }
    if (lmb.empty()) lmb.assign(J, std::vector<double>(L, 0.0)); // fallback / cold
    double mu = 2.0;
    int stall = 0;

    // Seed best_ub with the certified all-open cost so the Polyak step can
    // never blow up on an early infeasible repair. Infeasible all-open = bad instance.
    {
        std::vector<int> ones_y(p_.nfactories, 1), ones_z(J, 1);
        FlowResult fr = eval_.evaluate(ones_y, ones_z);
        if (!fr.feasible) {
            std::cerr << "[Lagrangian] all-open routing infeasible: aborting.\n";
            return res;
        }
        FlowResult fr2 = eval_.evaluate(fr.used_factories, fr.used_warehouses);
        res.best_ub = fr2.feasible
                          ? eval_.solution_cost(fr.used_factories, fr.used_warehouses, fr2)
                          : eval_.solution_cost(ones_y, ones_z, fr);
        res.y = fr.used_factories;
        res.z = fr.used_warehouses;
        res.has_solution = true;
    }
    // Tight subproblem bounds and a per-iteration time cap keep every
    // reported bound valid at any interruption.
    B_.set(GRB_DoubleParam_MIPGap, 0.0);
    B_.set(GRB_DoubleParam_MIPGapAbs, 0.25);

    // --- instrumentação R3 (somente contagem e impressão) ---
    long long instr_evals = 0;   // avaliações de A+B completadas
    long long instr_trunc = 0;   // avaliações em que o MIP de B parou no tempo
    long long instr_stall = 0;   // iterações sem melhora de best_lb
    int       instr_iters = 0;   // iterações efetivamente executadas

    for (int it = 0; it < iters && elapsed() < time_limit; ++it) {
        B_.set(GRB_DoubleParam_TimeLimit, std::max(1.0, time_limit - elapsed()));
        AResult A = solve_A(lmb);

        for (int j = 0; j < J; ++j)
            for (int k = 0; k < K; ++k)
                for (int l = 0; l < L; ++l)
                    wB_[j][k][l].set(GRB_DoubleAttr_Obj, p_.flowcost_wc[l][j][k] + lmb[j][l]);
        const auto instr_tB0 = std::chrono::steady_clock::now();
        B_.optimize();
        const double instr_tB =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - instr_tB0).count();
        const bool instr_truncB = (B_.get(GRB_IntAttr_Status) != GRB_OPTIMAL);
        ++instr_evals;
        if (instr_truncB) ++instr_trunc;
        if (B_.get(GRB_IntAttr_Status) != GRB_OPTIMAL &&
            B_.get(GRB_IntAttr_SolCount) == 0)
            break; // should not happen; bail out conservatively
        const double vB_bound = B_.get(GRB_DoubleAttr_ObjBound); // valid even w/ gap
        const double lb = A.value + vB_bound;

        if (lb > res.best_lb + 1e-9 * std::max(1.0, std::abs(res.best_lb))) {
            res.best_lb = lb;
            res.lmb = lmb;
            stall = 0;
        } else if (++stall >= 20) {
            mu = std::max(mu / 2.0, 1e-3);
            stall = 0;
        }

        instr_iters = it + 1;
        if (lb <= res.best_lb - 1e-9) ++instr_stall;   // não melhorou nesta iteração
        if (verbose) {
            std::cout << "[lag] it=" << (it + 1)
                      << " A+B=" << lb
                      << " best=" << res.best_lb
                      << " tB=" << instr_tB << "s"
                      << " truncB=" << (instr_truncB ? 1 : 0) << "\n";
        }

        std::vector<int> zbar(J);
        std::vector<std::vector<double>> wout(J, std::vector<double>(L, 0.0));
        for (int j = 0; j < J; ++j) {
            zbar[j] = zB_[j].get(GRB_DoubleAttr_X) > 0.5;
            for (int l = 0; l < L; ++l)
                for (int k = 0; k < K; ++k)
                    wout[j][l] += wB_[j][k][l].get(GRB_DoubleAttr_X);
        }
        std::vector<int> ry, rz;
        double ub = repair(A.ybar, zbar, ry, rz);
        if (ub < res.best_ub) {
            res.best_ub = ub;
            res.y = ry;
            res.z = rz;
            res.has_solution = true;
        }

        // subgradient s_jl = sum_k w_jkl - sum_i x_ijl ; Polyak step
        double norm2 = 0.0;
        for (int j = 0; j < J; ++j)
            for (int l = 0; l < L; ++l) {
                const double s = wout[j][l] - A.xin[j][l];
                norm2 += s * s;
            }
        if (norm2 < 1e-12) break; // dualized constraints tight: dual optimum reached
        const double t = mu * std::max(res.best_ub - lb, 1e-6) / norm2;
        for (int j = 0; j < J; ++j)
            for (int l = 0; l < L; ++l) {
                const double s = wout[j][l] - A.xin[j][l];
                lmb[j][l] = std::max(0.0, lmb[j][l] + t * s);
            }
        res.iterations = it + 1;
        if (verbose && it % 10 == 0)
            std::cout << "[Lagrangian] it=" << it << " L=" << lb
                      << " best_lb=" << res.best_lb << " best_ub=" << res.best_ub
                      << " t=" << elapsed() << "s\n";
    }
    res.runtime = elapsed();
    if (verbose)
        std::cout << "[Lagrangian] done: v_LD=" << res.best_lb << " UB=" << res.best_ub
                  << " iters=" << res.iterations << " time=" << res.runtime << "s\n";
    if (verbose) {
        std::cout << "[lag] iters=" << instr_iters
                  << " evals=" << instr_evals
                  << " trunc=" << instr_trunc
                  << " stalls=" << instr_stall << "\n";
    }
    return res;
}

} // namespace mptscfl
