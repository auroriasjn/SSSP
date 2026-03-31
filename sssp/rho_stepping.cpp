//
// Created by Jeremy Ng on 3/2/26.
//

#include "../graph.h"
#include "../utils.h"
#include "rho_stepping.h"

#ifndef NDEBUG
#define DBG(x) do { x; } while (0)
#else
#define DBG(x) do {} while (0)
#endif

void RhoSteppingSolver::solve(const Graph& g, Vertex source) {
    const size_t n = g.num_vertices();

    // Authoritative tentative distances
    dist = DistSeq(n);
    parlay::parallel_for(0, n, [&](std::size_t i) {
        dist[i].store(INF, std::memory_order_relaxed);
    });
    dist[source].store(static_cast<Distance>(0), std::memory_order_relaxed);

    // Array-based LaB-PQ
    ArrayLaBPQ pq(dist);
    pq.update(source);

    DBG(std::cerr << "INF = " << INF << "\n";);
    DBG(std::cerr << "dist[source] = " << dist[source].load() << "\n";);

    while (!pq.empty()) {
        // Snapshot current active records
        auto active = pq.active_vertices();
        if (active.empty()) {
            DBG(std::cerr << "[solve] Active vertices empty. \n";);
            break;
        }

        // Compute rho-threshold from current active frontier
        Distance theta = pq.get_threshold(active, dist, rho_);
        DBG(std::cerr << "[solve] Current threshold: " << theta << "\n";);

        // Extract everything with key <= theta
        auto frontier = pq.extract(theta);
        if (frontier.empty()) {
            // Defensive: should not usually happen, but avoids infinite loops
            DBG(std::cerr << "[solve] frontier empty. Breaking..." << "\n");
            break;
        }

        // Relax outgoing edges from all extracted vertices
        DBG(std::cerr << "[solve] frontier size is currently " << frontier.size() << "\n";);
        DBG({
            std::cerr << "[solve] frontier currently contains:\n  ";
            for (auto v : frontier) {
                std::cerr << v;
            }
            std::cerr << std::endl;
        });

        parlay::parallel_for(0, frontier.size(), [&](std::size_t i) {
            Vertex u = frontier[i];
            Distance du = dist[u].load(std::memory_order_relaxed);

            if (du == INF) {
                DBG(std::cerr << "[solve] [parallel-for] infinite distance detected at vertex " << u << "\n";);
                return;
            }

            const auto& neigh = g.neighbors(u);
            parlay::parallel_for(0, neigh.size(), [&](std::size_t j) {
                const auto& e = neigh[j];
                Vertex v = e.to;
                Distance nd = du + e.weight;

                if (write_min(dist[v], nd)) {
                    pq.update(v);
                }
            });
        });
    }
}
