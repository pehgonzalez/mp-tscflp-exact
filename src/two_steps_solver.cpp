#include "two_steps_solver.hpp"

#include <chrono>
#include <iostream>

namespace mptscfl {

TwoStepsResult TwoStepsSolver::solve(double time_limit, int method, bool run_exact,
                                     const std::string& logfile) {
    // The heuristic time must come out of the total budget, as done for the
    // Lagrangian phase in lb_method 2.
    const auto t0 = std::chrono::steady_clock::now();
    TwoStepsResult r;
    std::vector<int> used_fac(p_.nfactories, 0), used_ware(p_.nwarehouses, 0);

    // Steps 1-2 per product, as in TwoStepsSolver::solve + merge() of the original code.
    for (int l = 0; l < p_.ncommodities; ++l) {
        // Stage 2: warehouses (facilities) -> customers.
        SingleStageModel n1(p_.nwarehouses, p_.ncustomers, &p_.t_customer_demand[l],
                            p_.t_warehouse_capacity[l], p_.fixedcost_warehouse,
                            p_.flowcost_wc[l]);
        if (n1.solve() == Status::Infeasible) {
            std::cerr << "[TwoSteps] stage-2 infeasible for product " << l << "\n";
            return r;
        }
        // Stage 1: factories (facilities) -> warehouses, demand = warehouse usage.
        SingleStageModel n2(p_.nfactories, p_.nwarehouses, nullptr,
                            p_.t_factory_capacity[l], p_.fixedcost_factory,
                            p_.flowcost_fw[l]);
        n2.set_demand(n1.facility_usage());
        if (n2.solve() == Status::Infeasible) {
            std::cerr << "[TwoSteps] stage-1 infeasible for product " << l << "\n";
            return r;
        }
        for (int j = 0; j < p_.nwarehouses; ++j) used_ware[j] |= n1.opened()[j];
        for (int i = 0; i < p_.nfactories; ++i) used_fac[i] |= n2.opened()[i];
    }

    // Merged evaluation: route every product through the union, then CNUF.
    FlowResult fr = eval_.evaluate(used_fac, used_ware);
    if (!fr.feasible) {
        std::cerr << "[TwoSteps] merged routing infeasible (should not happen with union)\n";
        return r;
    }
    // Close-not-used-facilities, then re-evaluate on the filtered set (cost can only drop
    // in fixed charges; the routing on the filtered set is identical by construction).
    FlowResult fr2 = eval_.evaluate(fr.used_factories, fr.used_warehouses);
    if (fr2.feasible && eval_.solution_cost(fr.used_factories, fr.used_warehouses, fr2) <=
                            eval_.solution_cost(used_fac, used_ware, fr)) {
        r.y = fr.used_factories;
        r.z = fr.used_warehouses;
        r.heuristic_cost = eval_.solution_cost(r.y, r.z, fr2);
    } else {
        r.y = used_fac;
        r.z = used_ware;
        r.heuristic_cost = eval_.solution_cost(r.y, r.z, fr);
    }
    r.heuristic_feasible = true;
    std::cout << "[TwoSteps] heuristic cost: " << r.heuristic_cost << "\n";

    if (run_exact) {
        const double spent =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        ExactModel exact(p_);
        if (!logfile.empty()) exact.set_log(logfile);
        exact.set_start(r.y, r.z);
        // Integral optimum: UB + 0.499 is a noise-robust cutoff that still
        // admits solutions of value UB.
        r.exact = exact.run(std::max(1.0, time_limit - spent),
                            r.heuristic_cost + 0.499, method);
    }
    return r;
}

} // namespace mptscfl
