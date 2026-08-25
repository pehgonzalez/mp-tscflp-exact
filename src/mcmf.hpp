// Integer min-cost max-flow via successive shortest paths (Dijkstra + Johnson potentials).
// All instance data is integral (verified at parse time), so long long
// arithmetic is exact.
#ifndef MPTSCFL_MCMF_HPP
#define MPTSCFL_MCMF_HPP

#include <cstdint>
#include <limits>
#include <queue>
#include <vector>

namespace mptscfl {

class MinCostMaxFlow {
public:
    using ll = long long;
    static constexpr ll INF = std::numeric_limits<ll>::max() / 4;

    explicit MinCostMaxFlow(int n) : n_(n), head_(n, -1), pot_(n, 0), dist_(n), preve_(n) {}

    // Adds directed arc u->v; returns arc id (use it to query flow()).
    int add_arc(int u, int v, ll cap, ll cost) {
        arcs_.push_back({v, head_[u], cap, cost});
        head_[u] = static_cast<int>(arcs_.size()) - 1;
        arcs_.push_back({u, head_[v], 0, -cost}); // residual
        head_[v] = static_cast<int>(arcs_.size()) - 1;
        return static_cast<int>(arcs_.size()) - 2;
    }

    // Sends min-cost flow s->t up to `want`. Returns {flow_sent, total_cost}.
    // Requires all original arc costs >= 0 (true for the PSC data).
    std::pair<ll, ll> run(int s, int t, ll want = INF) {
        ll flow = 0, cost = 0;
        std::fill(pot_.begin(), pot_.end(), 0);
        while (flow < want && dijkstra(s, t)) {
            for (int v = 0; v < n_; ++v)
                if (dist_[v] < INF) pot_[v] += dist_[v];
            ll push = want - flow;
            for (int e = preve_[t]; e != -1; e = preve_[arcs_[e ^ 1].to])
                push = std::min(push, arcs_[e].cap);
            for (int e = preve_[t]; e != -1; e = preve_[arcs_[e ^ 1].to]) {
                arcs_[e].cap -= push;
                arcs_[e ^ 1].cap += push;
                cost += push * arcs_[e].cost;
            }
            flow += push;
        }
        return {flow, cost};
    }

    // Flow on original arc `id` (as returned by add_arc).
    ll flow(int id) const { return arcs_[id ^ 1].cap; }

private:
    struct Arc { int to, next; ll cap, cost; };

    bool dijkstra(int s, int t) {
        std::fill(dist_.begin(), dist_.end(), INF);
        std::fill(preve_.begin(), preve_.end(), -1);
        using QE = std::pair<ll, int>;
        std::priority_queue<QE, std::vector<QE>, std::greater<>> pq;
        dist_[s] = 0;
        pq.push({0, s});
        while (!pq.empty()) {
            auto [d, u] = pq.top(); pq.pop();
            if (d > dist_[u]) continue;
            for (int e = head_[u]; e != -1; e = arcs_[e].next) {
                const Arc& a = arcs_[e];
                if (a.cap <= 0) continue;
                ll nd = d + a.cost + pot_[u] - pot_[a.to];
                if (nd < dist_[a.to]) {
                    dist_[a.to] = nd;
                    preve_[a.to] = e;
                    pq.push({nd, a.to});
                }
            }
        }
        return dist_[t] < INF;
    }

    int n_;
    std::vector<int> head_;
    std::vector<Arc> arcs_;
    std::vector<ll> pot_, dist_;
    std::vector<int> preve_;
};

} // namespace mptscfl
#endif
