//
// Created by Jeremy Ng on 2/23/26.
//

#include "../graph.h"
#include "parallel_dijkstra.h"

// CPSC 424. Adaptation of https://github.com/cmuparlay/parlaylib/blob/master/examples/bucketed_dijkstra.h
// Assuming a *real positive integral weight*
void ParallelDijkstraSolver::solve(const Graph& g, Vertex source) {
    const auto n = static_cast<std::size_t>(g.num_vertices());

    // Allocate space.
    dist.assign(n, INF);

    parlay::sequence<std::atomic<Distance>> dist_a(n);
    parlay::parallel_for(0, n, [&](std::size_t i) {
        dist_a[i].store(INF, std::memory_order_relaxed);
    });
    dist_a[source].store(static_cast<Distance>(0), std::memory_order_relaxed);

    parlay::sequence<NestV> buckets(1);
    buckets[0] = NestV(1, SeqV(1, source));

    Bucket max_bucket = 0;
    for (Bucket b = 0; b <= max_bucket; ++b) {
        auto frontier = parlay::filter(parlay::flatten(buckets[b]), [&](Vertex v) {
            return dist_a[v].load(std::memory_order_relaxed) == static_cast<Distance>(b);
        });
        if (frontier.empty()) continue;

        // Candidates are (bucket_index, vertex)
        auto candidates = delayed::flatten(parlay::map(frontier, [&](Vertex u) {
            const Distance du = dist_a[u].load(std::memory_order_relaxed);
            return delayed::map(g.neighbors(u), [=](const Graph::Edge& e) {
                const Distance nd = du + static_cast<Distance>(e.weight);
                const auto k = static_cast<Bucket>(nd);
                return std::pair<Bucket, Vertex>(k, e.to);
            });
        }));

        // Keep only successful decreases (atomic)
        auto kept = delayed::to_sequence(delayed::filter(candidates, [&](const auto& kv) {
            const Bucket k = kv.first;
            const Vertex v = kv.second;
            const auto nd = static_cast<Distance>(k);

            return parlay::write_min(&dist_a[v], nd, std::less<Distance>());
        }));
        if (kept.empty()) continue;

        const Bucket new_max = parlay::reduce(
                delayed::map(kept, [](const auto& kv) { return kv.first; }),
                parlay::maximum<Bucket>()
        );

        if (new_max >= buckets.size()) buckets.resize(new_max + 1);

        auto grouped = parlay::group_by_index(kept, new_max + 1);
        parlay::parallel_for((Bucket)0, new_max + 1, [&](Bucket i) {
            if (grouped[i].empty()) return;
            buckets[i].push_back(std::move(grouped[i]));
        });

        if (new_max > max_bucket) max_bucket = new_max;
    }

    // Copy atomic results back into dist for the API
    parlay::parallel_for(0, n, [&](std::size_t i) {
        dist[i] = dist_a[i].load(std::memory_order_relaxed);
    });
}
