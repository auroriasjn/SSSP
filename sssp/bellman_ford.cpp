//
// Created by Jeremy Ng on 3/30/26.
//

#include "../graph.h"
#include "bellman_ford.h"

//
// Created by Jeremy Ng on 3/30/26.
//

#include "../graph.h"
#include "bellman_ford.h"

void BellmanFordSolver::solve(const Graph& g, Vertex source) {
    const size_t n = g.num_vertices();
    dist.assign(n, INF);
    dist[source] = 0;

    reinsertions.clear();
    reinsertions.reserve(n > 0 ? n - 1 : 0);

    VertexList frontier{source};
    VertexList next_frontier;
    std::vector<uint8_t> in_next(n, 0);
    std::vector<uint8_t> seen_before(n, 0);

    // source is already in a frontier initially
    seen_before[source] = 1;

    for (size_t pass = 0; pass + 1 < n && !frontier.empty(); ++pass) {
        next_frontier.clear();
        std::fill(in_next.begin(), in_next.end(), 0);

        uint64_t pass_reinsertions = 0;

        for (Vertex u : frontier) {
            if (dist[u] == INF) continue;

            for (const auto& e : g.neighbors(u)) {
                Vertex v = e.to;
                Distance nd = dist[u] + e.weight;

                if (nd < dist[v]) {
                    dist[v] = nd;

                    if (!in_next[v]) {
                        in_next[v] = 1;
                        next_frontier.push_back(v);

                        if (seen_before[v]) {
                            ++pass_reinsertions;
                        } else {
                            seen_before[v] = 1;
                        }
                    }
                }
            }
        }

        reinsertions.push_back(pass_reinsertions);
        frontier.swap(next_frontier);
    }
}