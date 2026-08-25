// Injection of Gurobi parameters from the environment, for the certificate
// verification block of the revision. MPTSCFL_GRB_PARAMS holds pairs written
// "Name=Value;Name=Value", applied to every model this binary builds, after
// the code's own settings of record and before any solve, so a rerun under
// NumericFocus=3 and tightened tolerances differs from the campaign run in
// nothing else. Unset, the variable changes nothing at all, which is what
// keeps the campaign binary and the verification binary the same binary.
//
// An unknown parameter name throws GRBException out of the run, by design.
// A verification whose requested regime was silently not applied would be
// worse than no verification.
#ifndef MPTSCFL_PARAM_INJECT_HPP
#define MPTSCFL_PARAM_INJECT_HPP

#include <cstdlib>
#include <sstream>
#include <string>

#include "gurobi_c++.h"

namespace mptscfl {

inline void apply_env_params(GRBModel& m) {
    const char* raw = std::getenv("MPTSCFL_GRB_PARAMS");
    if (!raw || !*raw) return;
    std::stringstream ss(raw);
    std::string pair;
    while (std::getline(ss, pair, ';')) {
        const auto eq = pair.find('=');
        if (eq == std::string::npos || eq == 0) continue;
        m.getEnv().set(pair.substr(0, eq), pair.substr(eq + 1));
    }
}

} // namespace mptscfl
#endif
