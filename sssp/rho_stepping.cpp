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

Distance RhoSteppingSolver::get_threshold(const VertexSeq& frontier, const DistSeq& dist_a) {
    constexpr std::size_t SSSP_SAMPLES = 1000;
    const std::size_t frontier_size = frontier.size();

    if (frontier_size == 0) {
        return INF;
    }

    if (frontier_size <= rho_) {
        DBG(std::cerr << "[get_threshold] frontier size less than rho.\n";);
        auto frontier_dist = parlay::delayed_tabulate(frontier_size, [&](std::size_t i) {
            return dist_a[frontier[i]].load(std::memory_order_relaxed);
        });
        return *parlay::max_element(frontier_dist);
    }

    std::array<Distance, SSSP_SAMPLES + 1> sample_dist{};
    for (std::size_t i = 0; i <= SSSP_SAMPLES; ++i) {
        Vertex v = frontier[hash_value(seed_ + static_cast<uint32_t>(i)) % frontier_size];
        sample_dist[i] = dist_a[v].load(std::memory_order_relaxed);
    }

    seed_ += static_cast<uint32_t>(SSSP_SAMPLES + 1);

    std::size_t id = static_cast<std::size_t>(
            (static_cast<double>(rho_) / static_cast<double>(frontier_size)) * SSSP_SAMPLES
    );
    if (id > SSSP_SAMPLES) id = SSSP_SAMPLES;

    std::sort(sample_dist.begin(), sample_dist.end());
    return sample_dist[id];
}

void RhoSteppingSolver::solve(const Graph& g, Vertex source) {
    const size_t n = g.num_vertices();

    // Authoritative tentative distances
    DistSeq dist_a(n);
    parlay::parallel_for(0, n, [&](std::size_t i) {
        dist_a[i].store(INF, std::memory_order_relaxed);
    });
    dist_a[source].store(static_cast<Distance>(0), std::memory_order_relaxed);

    // Array-based LaB-PQ
    ArrayLaBPQ pq(dist_a);
    pq.update(source);

    DBG(std::cerr << "INF = " << INF << "\n";);
    DBG(std::cerr << "dist[source] = " << dist_a[source].load() << "\n";);

    while (!pq.empty()) {
        // Snapshot current active records
        auto active = pq.active_vertices();
        if (active.empty()) {
            DBG(std::cerr << "[solve] Active vertices empty. \n";);
            break;
        }

        // Compute rho-threshold from current active frontier
        Distance theta = get_threshold(active, dist_a);
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
            Distance du = dist_a[u].load(std::memory_order_relaxed);

            if (du == INF) {
                DBG(std::cerr << "[solve] [parallel-for] infinite distance detected at vertex " << u << "\n";);
                return;
            }

            const auto& neigh = g.neighbors(u);
            parlay::parallel_for(0, neigh.size(), [&](std::size_t j) {
                const auto& e = neigh[j];
                Vertex v = e.to;
                Distance nd = du + e.weight;

                if (write_min(dist_a[v], nd)) {
                    pq.update(v);
                }
            });
        });
    }

    // Copy out final distances for the solver API
    dist.assign(n, INF);
    parlay::parallel_for(0, n, [&](std::size_t i) {
        dist[i] = dist_a[i].load(std::memory_order_relaxed);
    });
}
