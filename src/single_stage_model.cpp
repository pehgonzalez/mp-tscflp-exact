#include "single_stage_model.hpp"

#include <sstream>

namespace mptscfl {

void SingleStageModel::build() {
    y_.reserve(nfac_);
    for (int i = 0; i < nfac_; ++i) {
        std::ostringstream n; n << "y[" << i << "]";
        y_.push_back(model_.addVar(0, 1, fixedcost_[i], GRB_BINARY, n.str()));
    }
    x_.assign(nfac_, std::vector<GRBVar>(ncus_));
    for (int i = 0; i < nfac_; ++i)
        for (int j = 0; j < ncus_; ++j) {
            std::ostringstream n; n << "x[" << i << "," << j << "]";
            x_[i][j] = model_.addVar(0, GRB_INFINITY, flowcost_[i][j], GRB_CONTINUOUS, n.str());
        }
    model_.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);

    demand_constr_.reserve(ncus_);
    for (int j = 0; j < ncus_; ++j) {
        GRBLinExpr e;
        for (int i = 0; i < nfac_; ++i) e += x_[i][j];
        std::ostringstream n; n << "Demand[" << j << "]";
        demand_constr_.push_back(model_.addConstr(e >= 0.0, n.str())); // RHS set in solve()
    }
    for (int i = 0; i < nfac_; ++i) {
        GRBLinExpr e;
        for (int j = 0; j < ncus_; ++j) e += x_[i][j];
        std::ostringstream n; n << "Capacity[" << i << "]";
        model_.addConstr(e <= capacity_[i] * y_[i], n.str());
    }
}

Status SingleStageModel::solve(double time_limit) {
    for (int j = 0; j < ncus_; ++j)
        demand_constr_[j].set(GRB_DoubleAttr_RHS, demand_[j]);
    model_.set(GRB_DoubleParam_TimeLimit, time_limit);
    model_.set(GRB_IntParam_OutputFlag, 0);
    model_.optimize();

    const int st = model_.get(GRB_IntAttr_Status);
    if (st == GRB_INFEASIBLE) return Status::Infeasible;
    if (model_.get(GRB_IntAttr_SolCount) == 0) return Status::NotFound;

    obj_ = model_.get(GRB_DoubleAttr_ObjVal);
    opened_.assign(nfac_, 0);
    usage_.assign(nfac_, 0.0);
    for (int i = 0; i < nfac_; ++i) {
        opened_[i] = y_[i].get(GRB_DoubleAttr_X) > 0.5;
        double s = 0.0;
        for (int j = 0; j < ncus_; ++j) s += x_[i][j].get(GRB_DoubleAttr_X);
        usage_[i] = s;
    }
    return (st == GRB_OPTIMAL) ? Status::OptimalFound : Status::SolutionFound;
}

} // namespace mptscfl
