#include "dijkstra.h"
#include "sssp_solver.h"
#include "../graph.h"

void DijkstraSolver::solve(const Graph& g, Vertex source) {
    dist.assign(g.num_vertices(), INF);
    dist[source] = 0.0;

    std::priority_queue<PQNode, std::vector<PQNode>, std::greater<>> pq;
    pq.emplace(0.0, source);

    while (!pq.empty()) {
        auto [du, u] = pq.top();
        pq.pop();

        if (du > dist[u]) continue;

        for (const auto& edge : g.neighbors(u)) {
            auto v = edge.to;
            auto w = edge.weight;

            auto nd = du + w;
            if (nd < dist[v]) {
                dist[v] = nd;
                pq.emplace(nd, v);
            }
        }
    }
}