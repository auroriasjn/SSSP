//
// Created by Jeremy Ng on 2/23/26.
//

#include "../graph.h"
#include "parallel_dijkstra.h"

// CPSC 424. Adaptation of https://github.com/cmuparlay/parlaylib/blob/master/examples/bucketed_dijkstra.h
// Assuming a *real positive integral weight*
//
// Created by Jeremy Ng on 2/23/26.
//

#include "../graph.h"
#include "../types.h"
#include "../parallel_types.h"

#include <queue>
#include <vector>

void ParallelDijkstraSolver::solve(const Graph& g, Vertex source) {
    const auto n = g.num_vertices();

    // Atomic distance array so parallel relaxations are safe.
    dist = DistSeq(n);
    parlay::parallel_for(0, n, [&](std::size_t i) {
        dist[i].store(INF, std::memory_order_relaxed);
    });
    dist[source].store(0, std::memory_order_relaxed);

    std::priority_queue<PQNode, std::vector<PQNode>, std::greater<>> pq;
    pq.emplace(0, source);

    while (!pq.empty()) {
        const auto [du, u] = pq.top();
        pq.pop();

        // Standard stale-entry check.
        if (du != dist[u].load(std::memory_order_relaxed)) continue;

        const auto& nbrs = g.neighbors(u);
        const std::size_t m = nbrs.size();

        // Mark whether each relaxation succeeded, and what to push.
        parlay::sequence<uint8_t> improved(m, 0);
        parlay::sequence<PQNode> updates(m);

        parlay::parallel_for(0, m, [&](std::size_t i) {
            const auto& e = nbrs[i];
            const Vertex v = e.to;
            const Distance nd = du + static_cast<Distance>(e.weight);

            if (parlay::write_min(&dist[v], nd, std::less<Distance>())) {
                improved[i] = 1;
                updates[i] = {nd, v};
            }
        });

        // Sequentially push successful decreases.
        for (std::size_t i = 0; i < m; ++i) {
            if (improved[i]) {
                pq.push(updates[i]);
            }
        }
    }
}