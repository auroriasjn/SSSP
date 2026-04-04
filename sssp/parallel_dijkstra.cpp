//
// Created by Jeremy Ng on 2/23/26.
//

#include "../graph.h"
#include "../utils.h"
#include "parallel_dijkstra.h"

#include <queue>
#include <vector>

void ParallelDijkstraSolver::solve(const Graph& g, Vertex source) {
    const auto n = g.num_vertices();

    dist = DistSeq(n);
    reinsertions = SizeSeq(n);

    parlay::parallel_for(0, n, [&](std::size_t i) {
        dist[i].store(INF, std::memory_order_relaxed);
        reinsertions[i].store(0, std::memory_order_relaxed);
    });
    dist[source].store(static_cast<Distance>(0), std::memory_order_relaxed);

    std::priority_queue<PQNode, std::vector<PQNode>, std::greater<>> pq;
    pq.emplace(0, source);

    // Threshold for spawning threads.
    const std::size_t PARALLEL_THRESHOLD = 512;
    while (!pq.empty()) {
        const auto [du, u] = pq.top();
        pq.pop();

        if (du != dist[u].load(std::memory_order_relaxed)) continue;

        const auto& nbrs = g.neighbors(u);
        const std::size_t m = nbrs.size();

        // Sequential fallback for low-degree vertices.
        if (m < PARALLEL_THRESHOLD) {
            for (std::size_t i = 0; i < m; ++i) {
                const auto& e = nbrs[i];
                const Vertex v = e.to;
                const Distance nd = du + static_cast<Distance>(e.weight);

                Distance old;
                if (write_min_old(dist[v], nd, old)) {
                    if (old != INF) {
                        reinsertions[v].fetch_add(1, std::memory_order_relaxed);
                    }
                    pq.emplace(nd, v);
                }
            }
        }
        else {
            parlay::sequence<uint8_t> improved(m, 0);
            parlay::sequence<PQNode> updates(m);

            parlay::parallel_for(0, m, [&](std::size_t i) {
                const auto& e = nbrs[i];
                const Vertex v = e.to;
                const Distance nd = du + static_cast<Distance>(e.weight);

                Distance old;
                if (write_min_old(dist[v], nd, old)) {
                    if (old != INF) {
                        reinsertions[v].fetch_add(1, std::memory_order_relaxed);
                    }
                    improved[i] = 1;
                    updates[i] = {nd, v};
                }
            });

            for (std::size_t i = 0; i < m; ++i) {
                if (improved[i]) {
                    pq.push(updates[i]);
                }
            }
        }
    }
}