// Two-stage min-cost flow evaluator: given open factories/warehouses, computes the
// cheapest routing of every product and the facilities actually used.
// Backend: built-in successive-shortest-paths solver (mcmf.hpp) behind a
// solver-independent interface.
#ifndef MPTSCFL_FLOW_EVALUATOR_HPP
#define MPTSCFL_FLOW_EVALUATOR_HPP

#include <cmath>
#include <limits>
#include <vector>

#include "instance.hpp"
#include "mcmf.hpp"

namespace mptscfl {

struct FlowResult {
    bool feasible = false;
    double transport_cost = 0.0;           // sum over products of routing cost
    std::vector<int> used_factories;        // 0/1, flow > 0
    std::vector<int> used_warehouses;       // 0/1, flow > 0
};

class FlowEvaluator {
public:
    explicit FlowEvaluator(const Instance& inst) : p_(inst) {}

    // Routes all products through the open facilities. Node layout per product:
    // S -> factory i (cap b_il) -> warehouse j (in/out split, cap p_jl) -> customer k -> T (cap q_kl).
    FlowResult evaluate(const std::vector<int>& open_fac, const std::vector<int>& open_ware) const {
        FlowResult r;
        r.used_factories.assign(p_.nfactories, 0);
        r.used_warehouses.assign(p_.nwarehouses, 0);
        r.transport_cost = 0.0;
        for (int l = 0; l < p_.ncommodities; ++l) {
            if (!route_product(l, open_fac, open_ware, r)) return r; // infeasible
        }
        r.feasible = true;
        return r;
    }

    // Full objective for a location vector: fixed costs of open facilities + transport.
    // Facilities left unused by the flow still pay their fixed cost here (caller may CNUF first).
    double solution_cost(const std::vector<int>& open_fac, const std::vector<int>& open_ware,
                         const FlowResult& r) const {
        if (!r.feasible) return std::numeric_limits<double>::infinity();
        double c = r.transport_cost;
        for (int i = 0; i < p_.nfactories; ++i) c += open_fac[i] * p_.fixedcost_factory[i];
        for (int j = 0; j < p_.nwarehouses; ++j) c += open_ware[j] * p_.fixedcost_warehouse[j];
        return c;
    }

private:
    using ll = long long;
    static ll as_int(double v) { return static_cast<ll>(std::llround(v)); }

    bool route_product(int l, const std::vector<int>& open_fac, const std::vector<int>& open_ware,
                       FlowResult& r) const {
        const int I = p_.nfactories, J = p_.nwarehouses, K = p_.ncustomers;
        const int S = 0, T = 1 + I + 2 * J + K;
        auto fac = [&](int i) { return 1 + i; };
        auto win = [&](int j) { return 1 + I + j; };
        auto wout = [&](int j) { return 1 + I + J + j; };
        auto cus = [&](int k) { return 1 + I + 2 * J + k; };

        ll demand = as_int(p_.total_demand(l));
        MinCostMaxFlow g(T + 1);
        std::vector<int> arc_sf(I, -1), arc_ww(J, -1);

        for (int i = 0; i < I; ++i)
            if (open_fac[i])
                arc_sf[i] = g.add_arc(S, fac(i), as_int(p_.factory_capacity[i][l]), 0);
        for (int j = 0; j < J; ++j)
            if (open_ware[j])
                arc_ww[j] = g.add_arc(win(j), wout(j), as_int(p_.warehouse_capacity[j][l]), 0);
        for (int i = 0; i < I; ++i)
            if (open_fac[i])
                for (int j = 0; j < J; ++j)
                    if (open_ware[j])
                        g.add_arc(fac(i), win(j), demand, as_int(p_.flowcost_fw[l][i][j]));
        for (int j = 0; j < J; ++j)
            if (open_ware[j])
                for (int k = 0; k < K; ++k)
                    g.add_arc(wout(j), cus(k), demand, as_int(p_.flowcost_wc[l][j][k]));
        for (int k = 0; k < K; ++k)
            g.add_arc(cus(k), T, as_int(p_.customer_demand[k][l]), 0);

        auto [flow, cost] = g.run(S, T, demand);
        if (flow < demand) return false;
        r.transport_cost += static_cast<double>(cost);
        for (int i = 0; i < I; ++i)
            if (arc_sf[i] >= 0 && g.flow(arc_sf[i]) > 0) r.used_factories[i] = 1;
        for (int j = 0; j < J; ++j)
            if (arc_ww[j] >= 0 && g.flow(arc_ww[j]) > 0) r.used_warehouses[j] = 1;
        return true;
    }

    const Instance& p_;
};

} // namespace mptscfl
#endif
