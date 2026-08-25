// Solver-free unit tests: instance parser, MCMF, flow evaluator.
// Run: ./test_core <path-to-PSC1-C1-50-5.txt>
#include <cassert>
#include <cmath>
#include <iostream>

#include "../src/flow_evaluator.hpp"
#include "../src/instance.hpp"
#include "../src/mcmf.hpp"

using namespace mptscfl;

static void test_mcmf_known() {
    // 2 sources -> 2 sinks classic transportation: supplies {10,10}, demands {10,10},
    // costs [[1,4],[3,2]] -> optimal = 10*1 + 10*2 = 30.
    MinCostMaxFlow g(6); // 0=S,1,2=sources,3,4=sinks,5=T
    g.add_arc(0, 1, 10, 0);
    g.add_arc(0, 2, 10, 0);
    int a11 = g.add_arc(1, 3, 20, 1);
    g.add_arc(1, 4, 20, 4);
    g.add_arc(2, 3, 20, 3);
    int a22 = g.add_arc(2, 4, 20, 2);
    g.add_arc(3, 5, 10, 0);
    g.add_arc(4, 5, 10, 0);
    auto [f, c] = g.run(0, 5);
    assert(f == 20 && c == 30);
    assert(g.flow(a11) == 10 && g.flow(a22) == 10);

    // Infeasibility detection: demand exceeds capacity.
    MinCostMaxFlow h(3);
    h.add_arc(0, 1, 5, 1);
    h.add_arc(1, 2, 3, 0);
    auto [f2, c2] = h.run(0, 2, 5);
    assert(f2 == 3 && c2 == 3);
    std::cout << "test_mcmf_known: OK\n";
}

static void test_tiny_evaluator() {
    // Hand-built instance: 1 factory, 1 warehouse, 1 customer, 1 product.
    Instance p;
    p.nfactories = 1; p.nwarehouses = 1; p.ncustomers = 1; p.ncommodities = 1;
    p.customer_demand = {{7}};
    p.factory_capacity = {{10}};
    p.warehouse_capacity = {{10}};
    p.fixedcost_factory = {100};
    p.fixedcost_warehouse = {50};
    p.flowcost_fw = {{{2}}};
    p.flowcost_wc = {{{3}}};
    FlowEvaluator ev(p);
    auto r = ev.evaluate({1}, {1});
    assert(r.feasible);
    assert(std::abs(r.transport_cost - 7 * (2 + 3)) < 1e-9);
    assert(std::abs(ev.solution_cost({1}, {1}, r) - (35 + 150)) < 1e-9);
    auto r2 = ev.evaluate({0}, {1});
    assert(!r2.feasible);
    std::cout << "test_tiny_evaluator: OK\n";
}

static void test_parse_and_route(const char* path) {
    Instance p;
    p.load_file(path);
    assert(p.nfactories == 50 && p.nwarehouses == 100 && p.ncustomers == 200 &&
           p.ncommodities == 5);
    // Sanity: capacities read in the same order as the legacy parser.
    assert(p.TotalCapFac > 0 && p.TotalCapWare > 0);
    for (int l = 0; l < p.ncommodities; ++l) assert(p.total_demand(l) > 0);

    // All facilities open must be feasible (benchmark guarantees capacity >= demand).
    FlowEvaluator ev(p);
    std::vector<int> all_f(p.nfactories, 1), all_w(p.nwarehouses, 1);
    auto r = ev.evaluate(all_f, all_w);
    assert(r.feasible);
    assert(r.transport_cost > 0);
    double full = ev.solution_cost(all_f, all_w, r);
    // CNUF can only reduce total cost.
    auto r2 = ev.evaluate(r.used_factories, r.used_warehouses);
    assert(r2.feasible);
    double filtered = ev.solution_cost(r.used_factories, r.used_warehouses, r2);
    assert(filtered <= full + 1e-6);
    std::cout << "test_parse_and_route: OK (transport=" << r.transport_cost
              << ", all-open=" << full << ", CNUF=" << filtered << ")\n";
}

int main(int argc, char* argv[]) {
    test_mcmf_known();
    test_tiny_evaluator();
    if (argc > 1) test_parse_and_route(argv[1]);
    else std::cout << "(skipping instance test: no path given)\n";
    std::cout << "ALL TESTS PASSED\n";
    return 0;
}
