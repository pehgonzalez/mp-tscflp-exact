// Incumbent-improvement callback for the compact MIP (modes "exact-repaircb"
// and "exact-flgcb" of the third-round campaign, referee item R5).
//
// The compact model pays the fixed cost of every open facility, and the solver
// regularly accepts an incumbent that opens a facility the flow never uses (the
// LP relaxation has no incentive to close it once the flows are fixed). The
// FlowEvaluator re-routes (y, z) exactly, by successive shortest paths over
// integer data, and reports which facilities actually carry flow. Closing the
// idle ones keeps the routing feasible and strictly lowers the objective by
// their fixed cost, so the repaired point is a valid, cheaper solution of the
// same MIP.
//
// Division of labour between the two callback points, deliberately conservative:
//   MIPSOL  read y/z of the candidate, re-route, CNUF-filter, cost it. If the
//           exact cost beats the candidate objective by more than 0.5 (the data
//           is integral, so anything smaller is solver noise) the repaired point
//           is REMEMBERED. Nothing is written at MIPSOL: setSolution is only
//           legal at MIPNODE.
//   MIPNODE inject the remembered point once, with setSolution + useSolution,
//           and forget it. Only the location variables are written; the flow
//           variables x/w are left for the solver to complete, exactly as the
//           MIP start of ExactModel::set_start does.
//
// The callback never aborts the solve and never throws: any Gurobi or evaluator
// exception is counted and swallowed, since a callback that escapes with an
// exception poisons the whole optimization.
#ifndef MPTSCFL_REPAIR_CALLBACK_HPP
#define MPTSCFL_REPAIR_CALLBACK_HPP

#include <cmath>
#include <vector>

#include "gurobi_c++.h"
#include "flow_evaluator.hpp"
#include "instance.hpp"

namespace mptscfl {

class RepairCallback : public GRBCallback {
public:
    // yv/zv must be the location variables of the model this callback is
    // attached to; the callback keeps references to the model's vectors.
    RepairCallback(const Instance& inst, const std::vector<GRBVar>& yv,
                   const std::vector<GRBVar>& zv)
        : p_(inst), eval_(inst), y_(yv), z_(zv) {}

    long long repairs_found() const { return nfound_; }
    long long injections() const { return ninject_; }
    long long failures() const { return nfail_; }

protected:
    void callback() override {
        try {
            if (where == GRB_CB_MIPSOL) {
                capture();
            } else if (where == GRB_CB_MIPNODE) {
                inject_pending();
            }
        } catch (const GRBException&) {
            ++nfail_;   // never let it escape: it would poison the optimization
        } catch (const std::exception&) {
            ++nfail_;
        }
    }

private:
    // MIPSOL: re-route the candidate and remember a strictly cheaper repair.
    void capture() {
        const int I = p_.nfactories, J = p_.nwarehouses;
        std::vector<int> yy(I, 0), zz(J, 0);
        for (int i = 0; i < I; ++i) yy[i] = getSolution(y_[i]) > 0.5 ? 1 : 0;
        for (int j = 0; j < J; ++j) zz[j] = getSolution(z_[j]) > 0.5 ? 1 : 0;

        FlowResult fr = eval_.evaluate(yy, zz);
        if (!fr.feasible) return;  // solver handed us something we cannot route

        // CNUF: keep only the facilities the optimal routing actually uses.
        std::vector<int> ry = fr.used_factories, rz = fr.used_warehouses;
        FlowResult rr = eval_.evaluate(ry, rz);
        if (!rr.feasible) return;  // cannot happen with a min-cost routing; be safe
        const double repaired = eval_.solution_cost(ry, rz, rr);

        const double cand_obj = getDoubleInfo(GRB_CB_MIPSOL_OBJ);
        if (!(repaired < cand_obj - 0.5)) return;  // no integral improvement

        ++nfound_;
        // Um pendente ainda nao injetado so e substituido por algo melhor.
        if (has_pending_ && repaired >= pend_cost_) return;
        pend_y_ = ry;
        pend_z_ = rz;
        pend_cost_ = repaired;
        has_pending_ = true;
    }

    // MIPNODE: hand the remembered point to the solver, once.
    void inject_pending() {
        if (!has_pending_) return;
        if (getIntInfo(GRB_CB_MIPNODE_STATUS) != GRB_OPTIMAL) return;  // node LP unusable
        const int I = p_.nfactories, J = p_.nwarehouses;
        if (static_cast<int>(pend_y_.size()) != I ||
            static_cast<int>(pend_z_.size()) != J) {
            has_pending_ = false;
            return;
        }
        for (int i = 0; i < I; ++i) setSolution(y_[i], static_cast<double>(pend_y_[i]));
        for (int j = 0; j < J; ++j) setSolution(z_[j], static_cast<double>(pend_z_[j]));
        // Flows stay unset on purpose: the solver completes them by solving the
        // remaining transportation LP, which is exactly the routing we costed.
        useSolution();
        ++ninject_;
        has_pending_ = false;
    }

    const Instance& p_;
    FlowEvaluator eval_;
    const std::vector<GRBVar>& y_;
    const std::vector<GRBVar>& z_;

    bool has_pending_ = false;
    std::vector<int> pend_y_, pend_z_;
    double pend_cost_ = 0.0;

    long long nfound_ = 0, ninject_ = 0, nfail_ = 0;
};

} // namespace mptscfl
#endif
