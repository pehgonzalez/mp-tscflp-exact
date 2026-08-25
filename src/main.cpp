// CLI: mptscfl <instance> <lb_method> <time_limit> [mode] [seed] [threads]
//   lb_method: 0 = B&B, 1 = branch-and-Benders-cut, 2 = Lagrangian-guided Benders
//   mode:      "two-steps" (default) | "exact"
//   seed:      Gurobi Seed (default 0); threads: Gurobi Threads (default 0 = auto)
//
// Logging (append-only schema, nothing lost if console truncates):
//   logs/<inst>_m<method>_<mode>_s<seed>_<timestamp>.gurobi.log — full solver log
//   logs/results.csv — one row per run; legacy-schema files are renamed aside.
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>

#include "benders_model.hpp"
#include "exact_model.hpp"
#include "flow_evaluator.hpp"
#include "instance.hpp"
#include "lagrangian.hpp"
#include "repair_callback.hpp"
#include "two_steps_solver.hpp"

using namespace mptscfl;

static const char* CSV_HEADER =
    "datetime,instance,mode,method,seed,threads,time_limit_s,status,obj,bound,gap,"
    "solver_time_s,total_wall_s,lag_lb,lag_ub,lag_time_s,benders_cuts,"
    "heuristic_cost,verified_cost,verified_ok,gurobi_version,gurobi_log";

static void usage(const char* prog) {
    std::cout << "Usage: " << prog
              << " <instance> <lb_method> <time_limit> [mode] [seed] [threads]\n"
              << "  lb_method: 0 = B&B, 1 = branch-and-Benders-cut,\n"
              << "             2 = Lagrangian-guided Benders\n"
              << "  mode: two-steps (default) | exact;  seed: default 0;  threads: 0=auto\n";
}

static std::string status_str(Status s) {
    switch (s) {
        case Status::OptimalFound: return "OPTIMAL";
        case Status::SolutionFound: return "FEASIBLE";
        case Status::Infeasible: return "INFEASIBLE";
        default: return "NOTFOUND";
    }
}

int main(int argc, char* argv[]) {
    if (argc < 4) { usage(argv[0]); return 1; }
    const std::string datafile = argv[1];
    const int method = std::atoi(argv[2]);
    const double time_limit = std::atof(argv[3]);
    const std::string mode = (argc > 4) ? argv[4] : "two-steps";
    const int seed = (argc > 5) ? std::atoi(argv[5]) : 0;
    const int threads = (argc > 6) ? std::atoi(argv[6]) : 0;
    const bool papadakos = (argc > 7) && std::atoi(argv[7]) != 0; // ablation flag
    // Any mode starting with "exact" is exact (labels such as "exact-pap"
    // distinguish ablation configurations in the CSV without a new column).
    const bool exact_mode = mode.rfind("exact", 0) == 0;
    // Two further ablation switches ride on the mode label the same way.
    // "agg" collapses the per-product thetas into one (single-cut master);
    // "noldb" withholds the Lagrangian bound inequality from the L-BBC master
    // while keeping its start and cutoff, which is the channel-separation arm
    // the revision runs. Neither label collides with exact, exact-pap or the
    // variance tags exact-var1..3.
    const bool aggregated = mode.find("agg") != std::string::npos;
    const bool no_ld_bound = mode.find("noldb") != std::string::npos;
    // Revision-campaign modes, matched exactly so no legacy label collides.
    // "exact-nowarm"    method 2, subgradient starts at zero multipliers;
    // "exact-nostart"   method 2, cutoff installed, master start withheld,
    //                   the reproduction of the defective earlier build;
    // "exact-fl"        method 0, the rows (F_l) added to the compact model;
    // "exact-mipguided" method 0, v_LD as objective bound, repaired solution
    //                   as MIP start, cutoff, the second baseline;
    // "exact-mauri"     method 0, solver default tolerances, the protocol of
    //                   Mauri et al. (2021);
    // "exact-flguided"  method 0, both baselines combined, the rows (F_l)
    //                   together with v_LD, the repaired start and the cutoff.
    // Third-round modes (R3, R5, R12 of the referee list):
    // "exact-oneeval"   method 2, the Lagrangian phase reduced to ONE evaluation
    //                   of L(lambda) at the LP duals: warm start kept, zero
    //                   subgradient steps (iters=1; the loop of
    //                   LagrangianSolver::solve is `for it < iters`, so iters=1
    //                   is exactly one evaluation of A+B and no multiplier
    //                   update is ever applied to a second one). Measures how
    //                   much of the transfer is the warm start alone;
    // "exact-thetaint"  method 2 and method 1, the master declares theta_l
    //                   integer, runs at MIPGapAbs 0.9999 and receives
    //                   ceil(v_LD) - 1e-6 instead of v_LD. Needs the
    //                   benders_model patch (see patch_thetaint.md): the flag
    //                   travels in BendersOptions::theta_integer, which only
    //                   exists once the patch is applied. Until then the mode
    //                   refuses to run instead of silently behaving as "exact";
    // "exact-repaircb"  method 0, the compact MIP with the incumbent-improvement
    //                   callback of repair_callback.hpp (MIPSOL re-routes the
    //                   candidate with the FlowEvaluator, MIPNODE injects the
    //                   cheaper repaired solution);
    // "exact-flgcb"     method 0, exact-flguided plus that same callback.
    const bool no_warm    = mode == "exact-nowarm";
    const bool no_start   = mode == "exact-nostart";
    const bool mip_fl     = mode == "exact-fl" || mode == "exact-flguided" ||
                            mode == "exact-flgcb";
    const bool mip_guided = mode == "exact-mipguided" || mode == "exact-flguided" ||
                            mode == "exact-flgcb";
    const bool mauri_prot = mode == "exact-mauri";
    const bool one_eval   = mode == "exact-oneeval";
    const bool theta_int  = mode == "exact-thetaint";
    const bool repair_cb  = mode == "exact-repaircb" || mode == "exact-flgcb";
    if (const char* pp = std::getenv("MPTSCFL_GRB_PARAMS"))
        std::cout << "Gurobi params injected into every model: " << pp << "\n";

    Instance inst;
    try {
        inst.load_file(datafile);
    } catch (const std::exception& e) {
        std::cerr << e.what() << "\n";
        return 1;
    }
    std::cout << "Instance loaded - " << datafile << " (" << inst.nfactories << "x"
              << inst.nwarehouses << "x" << inst.ncustomers << ", L=" << inst.ncommodities
              << ") seed=" << seed << " threads=" << threads << "\n";

    const std::string stem = std::filesystem::path(datafile).stem().string();
    std::time_t tt = std::time(nullptr);
    char ts[32];
    std::strftime(ts, sizeof(ts), "%Y%m%d-%H%M%S", std::localtime(&tt));
    std::filesystem::create_directories("logs");
    std::ostringstream basen;
    basen << "logs/" << stem << "_m" << method << "_" << mode << "_s" << seed << "_" << ts;
    const std::string gurobilog = basen.str() + ".gurobi.log";
    std::cout << "Gurobi log: " << gurobilog << "\n";

    std::ostringstream gv;
    gv << GRB_VERSION_MAJOR << "." << GRB_VERSION_MINOR << "." << GRB_VERSION_TECHNICAL;

    const auto wall0 = std::chrono::steady_clock::now();
    // Remaining TOTAL budget (fix for observed overrun: model-construction time — the
    // Lagrangian B model, the Benders subproblem LPs and the all-open theta bounds —
    // was not discounted from any phase; now every solve gets wall-clock remainder).
    auto remaining = [&](double floor_s = 1.0) {
        const double spent =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - wall0).count();
        return std::max(floor_s, time_limit - spent);
    };
    try {
        ExactResult res;
        double heuristic_cost = -1.0;
        long long cuts = 0;
        double lag_lb = -1.0, lag_ub = -1.0, lag_time = -1.0;
        std::vector<int> heur_y, heur_z; // for CUTOFF promotion

        if (method == 2) {
            const double lag_budget = std::min(0.15 * time_limit, 600.0);
            LagrangianSolver lag(inst);
            // iters=1 is the single-evaluation arm: the ascent loop runs once,
            // evaluates A+B at the warm-started multipliers and stops before any
            // subgradient step could produce a second evaluation.
            const int lag_iters = one_eval ? 1 : 1000;
            if (one_eval)
                std::cout << "[main] single-evaluation ablation: one L(lambda) at the "
                          << "LP duals, no subgradient iteration\n";
            LagrangianResult ld = lag.solve(/*iters=*/lag_iters, lag_budget,
                                            /*verbose=*/true, /*warm=*/!no_warm);
            lag_lb = ld.best_lb > -1e99 ? ld.best_lb : -1.0;
            lag_ub = ld.has_solution ? ld.best_ub : -1.0;
            lag_time = ld.runtime;
            BendersOptions bo;
            bo.threads = threads;
            bo.papadakos = papadakos;
            bo.aggregated = aggregated;
#ifdef MPTSCFL_HAS_THETA_INTEGER
            if (theta_int) {
                bo.theta_integer = true;
                std::cout << "[main] integral-theta master: theta_l in Z, MIPGapAbs "
                          << "0.9999, bound transfer at ceil(v_LD) - 1e-6\n";
            }
#else
            if (theta_int) {
                std::cerr << "[main] mode exact-thetaint requires the benders_model "
                             "patch of patch_thetaint.md (BendersOptions::theta_integer "
                             "is absent in this build); refusing to run so the CSV never "
                             "records a plain exact run under the thetaint label\n";
                return 3;
            }
#endif
            BendersModel benders(inst, bo);
            benders.set_log(gurobilog);
            benders.set_seed(seed);
            if (ld.best_lb > -1e99 && !no_ld_bound) benders.add_global_lower_bound(ld.best_lb);
            if (no_ld_bound)
                std::cout << "[main] channel ablation: Lagrangian bound inequality withheld, "
                          << "start and cutoff kept\n";
            double cutoff = -1.0;
            if (ld.has_solution) {
                heuristic_cost = ld.best_ub;
                heur_y = ld.y;
                heur_z = ld.z;
                if (no_start) {
                    std::cout << "[main] transfer ablation: master start withheld, "
                              << "cutoff kept, the defective-build reproduction\n";
                } else {
                    benders.set_start(ld.y, ld.z);
                }
                cutoff = ld.best_ub + 0.499; // integral optimum: noise-robust
            }
            res = benders.run(remaining(), cutoff);
            cuts = benders.cuts_added();
        } else if (method == 1) {
            BendersOptions bo;
            bo.threads = threads;
            bo.papadakos = papadakos;
            bo.aggregated = aggregated;
#ifdef MPTSCFL_HAS_THETA_INTEGER
            if (theta_int) {
                bo.theta_integer = true;
                std::cout << "[main] integral-theta master: theta_l in Z, MIPGapAbs "
                          << "0.9999 (no v_LD in method 1)\n";
            }
#else
            if (theta_int) {
                std::cerr << "[main] mode exact-thetaint requires the benders_model "
                             "patch of patch_thetaint.md; refusing to run\n";
                return 3;
            }
#endif
            BendersModel benders(inst, bo);
            benders.set_log(gurobilog);
            benders.set_seed(seed);
            if (!exact_mode) {
                TwoStepsSolver tss(inst);
                auto h = tss.solve(time_limit, /*method=*/0, /*run_exact=*/false);
                if (h.heuristic_feasible) {
                    heuristic_cost = h.heuristic_cost;
                    heur_y = h.y;
                    heur_z = h.z;
                    benders.set_start(h.y, h.z);
                    res = benders.run(remaining(), h.heuristic_cost + 0.499);
                } else {
                    res = benders.run(remaining());
                }
            } else {
                res = benders.run(remaining());
            }
            cuts = benders.cuts_added();
        } else if (exact_mode) {
            ExactModel exact(inst);
            exact.set_log(gurobilog);
            exact.set_seed(seed);
            if (mauri_prot) {
                exact.set_proof_mode(false);
                std::cout << "[main] Mauri et al. protocol: solver default tolerances, "
                          << "no proof mode\n";
            }
            if (mip_fl) exact.add_feasibility_rows();
            double mip_cutoff = -1.0;
            if (mip_guided) {
                const double lag_budget = std::min(0.15 * time_limit, 600.0);
                LagrangianSolver lag(inst);
                LagrangianResult ld = lag.solve(/*iters=*/1000, lag_budget);
                lag_lb = ld.best_lb > -1e99 ? ld.best_lb : -1.0;
                lag_ub = ld.has_solution ? ld.best_ub : -1.0;
                lag_time = ld.runtime;
                if (ld.best_lb > -1e99) exact.add_objective_lower_bound(ld.best_lb);
                if (ld.has_solution) {
                    heuristic_cost = ld.best_ub;
                    heur_y = ld.y;
                    heur_z = ld.z;
                    exact.set_start(ld.y, ld.z);
                    mip_cutoff = ld.best_ub + 0.499;
                }
                std::cout << "[main] baseline ablation: v_LD and the repaired "
                          << "solution handed to the compact model\n";
            }
            // The improvement callback has to be installed BEFORE run(), which
            // calls optimize() internally. It stays alive for the whole solve.
            std::unique_ptr<RepairCallback> rcb;
            if (repair_cb) {
                rcb.reset(new RepairCallback(inst, exact.factory_vars(),
                                             exact.warehouse_vars()));
                exact.model().setCallback(rcb.get());
                std::cout << "[main] repair callback installed: MIPSOL re-routes every "
                          << "incumbent candidate, MIPNODE injects the cheaper one\n";
            }
            res = exact.run(remaining(), mip_cutoff, method, threads);
            if (rcb)
                std::cout << "[main] repair callback: " << rcb->repairs_found()
                          << " candidates improved, " << rcb->injections()
                          << " injected, " << rcb->failures() << " callback errors\n";
        } else {
            TwoStepsSolver tss(inst);
            auto h = tss.solve(time_limit, method, /*run_exact=*/true, gurobilog);
            heuristic_cost = h.heuristic_feasible ? h.heuristic_cost : -1.0;
            if (h.heuristic_feasible) { heur_y = h.y; heur_z = h.z; }
            res = h.exact;
        }

        // CUTOFF/timeout without incumbent but with a certified heuristic UB
        // and bound >= UB - 0.5 proves the heuristic solution optimal (integral value).
        if (res.status == Status::NotFound && heuristic_cost > 0 && !heur_y.empty() &&
            res.bound >= heuristic_cost - 0.5 && res.bound > 0) {
            res.status = Status::OptimalFound;
            res.obj = heuristic_cost;
            res.y = heur_y;
            res.z = heur_z;
            std::cout << "[main] heuristic UB certified optimal by dual bound "
                      << res.bound << "\n";
        } else if (res.status == Status::NotFound && heuristic_cost > 0 &&
                   !heur_y.empty()) {
            // With an active cutoff the master may end without
            // any incumbent even though the heuristic solution is the best known.
            // Report it as the incumbent instead of swallowing it as NOTFOUND.
            res.status = Status::SolutionFound;
            res.obj = heuristic_cost;
            res.y = heur_y;
            res.z = heur_z;
            if (res.bound > 0) res.gap = (res.obj - res.bound) / std::abs(res.obj);
            std::cout << "[main] no improving incumbent; reporting heuristic UB "
                      << heuristic_cost << " (bound " << res.bound << ")\n";
        }

        // Independent certification (solver-free re-routing of (y,z)).
        double verified = -1.0;
        bool verified_ok = false;
        if (!res.y.empty()) {
            FlowEvaluator ev(inst);
            FlowResult fr = ev.evaluate(res.y, res.z);
            if (fr.feasible) {
                verified = ev.solution_cost(res.y, res.z, fr);
                verified_ok = true;
            }
        }
        const double total_wall =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - wall0).count();

        std::cout << std::fixed << std::setprecision(2)
                  << "status=" << status_str(res.status) << " obj=" << res.obj
                  << " bound=" << res.bound << " gap=" << std::setprecision(6) << res.gap
                  << std::setprecision(2) << " solver_time=" << res.runtime
                  << " total_wall=" << total_wall << "\n";
        if (!res.y.empty())
            std::cout << "verified_cost="
                      << (verified_ok ? std::to_string(verified) : "INFEASIBLE (bug!)")
                      << "\n";

        // Summary row. A legacy-schema CSV is moved aside once.
        const std::string csvpath = "logs/results.csv";
        if (std::filesystem::exists(csvpath)) {
            std::ifstream in(csvpath);
            std::string first;
            std::getline(in, first);
            in.close();
            if (first != CSV_HEADER)
                std::filesystem::rename(csvpath, "logs/results_legacy.csv");
        }
        const bool fresh = !std::filesystem::exists(csvpath);
        std::ofstream csv(csvpath, std::ios::app);
        if (fresh) csv << CSV_HEADER << "\n";
        csv << ts << "," << stem << "," << mode << "," << method << "," << seed << ","
            << threads << "," << time_limit << "," << status_str(res.status) << ","
            << std::fixed << std::setprecision(2) << res.obj << "," << res.bound << ","
            << std::setprecision(8) << res.gap << std::setprecision(2) << ","
            << res.runtime << "," << total_wall << "," << lag_lb << "," << lag_ub << ","
            << lag_time << "," << cuts << "," << heuristic_cost << ","
            << (verified_ok ? verified : -1.0) << "," << (verified_ok ? 1 : 0) << ","
            << gv.str() << "," << gurobilog << "\n";
        std::cout << "Summary row appended to " << csvpath << "\n";
    } catch (GRBException& e) {
        std::cerr << "Gurobi error " << e.getErrorCode() << ": " << e.getMessage() << "\n";
        return 2;
    }
    return 0;
}
