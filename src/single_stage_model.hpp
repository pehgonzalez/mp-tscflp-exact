// Single-stage capacitated FLP (one product).
// Used by TwoStepsSolver: stage 2 (warehouses->customers) then stage 1 (factories->warehouses).
//   min  sum_i f_i y_i + sum_ij c_ij x_ij
//   s.t. sum_i x_ij >= demand_j          (each "customer" j)
//        sum_j x_ij <= cap_i * y_i       (each "facility" i)
#ifndef MPTSCFL_SINGLE_STAGE_MODEL_HPP
#define MPTSCFL_SINGLE_STAGE_MODEL_HPP

#include <vector>

#include "gurobi_c++.h"
#include "instance.hpp"

namespace mptscfl {

class SingleStageModel {
public:
    SingleStageModel(int nfac, int ncus,
                     const std::vector<double>* demand,        // may be null; set later
                     const std::vector<double>& capacity,
                     const std::vector<double>& fixedcost,
                     const std::vector<std::vector<double>>& flowcost)
        : nfac_(nfac), ncus_(ncus), capacity_(capacity), fixedcost_(fixedcost),
          flowcost_(flowcost), env_(true), model_((env_.start(), env_)) {
        if (demand) demand_ = *demand;
        build();
    }

    // Mirrors D1solver::set_new_customer_demand (demands become RHS at solve time).
    void set_demand(const std::vector<double>& d) { demand_ = d; }

    Status solve(double time_limit = 3600.0);

    const std::vector<int>& opened() const { return opened_; }
    // Amount shipped out of each facility (D1solver::get_new_demand).
    const std::vector<double>& facility_usage() const { return usage_; }
    double objective() const { return obj_; }

private:
    void build();

    int nfac_, ncus_;
    std::vector<double> demand_, capacity_, fixedcost_;
    std::vector<std::vector<double>> flowcost_; // [i][j]
    GRBEnv env_;
    GRBModel model_;
    std::vector<GRBVar> y_;
    std::vector<std::vector<GRBVar>> x_;
    std::vector<GRBConstr> demand_constr_;
    std::vector<int> opened_;
    std::vector<double> usage_;
    double obj_ = 0.0;
};

} // namespace mptscfl
#endif
