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
    reinsertions = SizeSeq(n);

    parlay::parallel_for(0, n, [&](std::size_t i) {
        dist[i].store(INF, std::memory_order_relaxed);
        reinsertions[i].store(0, std::memory_order_relaxed);
    });
    dist[source].store(static_cast<Distance>(0), std::memory_order_relaxed);

    // Array-based LaB-PQ
    ArrayLaBPQ pq(dist);
    pq.update(source);

    DBG(std::cerr << "INF = " << INF << "\n";);
    DBG(std::cerr << "dist[source] = " << dist[source].load() << "\n";);

    // Threshold for spawning inner threads on high-degree vertices.
    // Tune this based on graph density: 512 or 1024 are good defaults.
    const std::size_t PARALLEL_THRESHOLD = 1024;

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
            DBG(std::cerr << "[solve] frontier empty. Breaking..." << "\n");
            break;
        }

        DBG(std::cerr << "[solve] frontier size is currently " << frontier.size() << "\n";);
        DBG({
                std::cerr << "[solve] frontier currently contains:\n  ";
                for (auto v : frontier) {
                    std::cerr << v;
                }
                std::cerr << std::endl;
            });

        // The outer parallel_for provides massive parallelism across the frontier
        parlay::parallel_for(0, frontier.size(), [&](std::size_t i) {
            Vertex u = frontier[i];
            Distance du = dist[u].load(std::memory_order_relaxed);

            if (du == INF) {
                DBG(std::cerr << "[solve] [parallel-for] infinite distance detected at vertex " << u << "\n";);
                return;
            }

            const auto& neigh = g.neighbors(u);
            const std::size_t m = neigh.size();

            // OPTIMIZATION: Fast sequential fallback for low-degree vertices.
            // This prevents nested parallel overhead on sparse graphs.
            if (m < PARALLEL_THRESHOLD) {
                for (std::size_t j = 0; j < m; ++j) {
                    const auto& e = neigh[j];
                    Vertex v = e.to;
                    Distance nd = du + e.weight;

                    Distance old;
                    if (write_min_old(dist[v], nd, old)) {
                        if (old != INF) {
                            reinsertions[v].fetch_add(1, std::memory_order_relaxed);
                        }
                        pq.update(v);
                    }
                }
            }
                // OPTIMIZATION: Nested parallelism kicks in ONLY for massive hub vertices.
            else {
                parlay::parallel_for(0, m, [&](std::size_t j) {
                    const auto& e = neigh[j];
                    Vertex v = e.to;
                    Distance nd = du + e.weight;

                    Distance old;
                    if (write_min_old(dist[v], nd, old)) {
                        if (old != INF) {
                            reinsertions[v].fetch_add(1, std::memory_order_relaxed);
                        }
                        pq.update(v);
                    }
                });
            }
        });
    }
}