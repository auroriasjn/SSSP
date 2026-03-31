//
// Created by Jeremy Ng on 3/23/26.
//

#ifndef SENIORTHESIS_PARALLEL_BUNDLE_DIJKSTRA_H
#define SENIORTHESIS_PARALLEL_BUNDLE_DIJKSTRA_H

#include <algorithm>
#include <queue>
#include <unordered_set>
#include <vector>
#include <random>
#include <cstdint>
#include <optional>

#include <parlay/primitives.h>
#include <parlay/sequence.h>

#include "../parallel_types.h"
#include "../labpq/array_labpq.h"
#include "sssp_solver.h"

class ParallelBundleDijkstraSolver : public SSSPSolver {
private:
    using BoolSeq = parlay::sequence<bool>;

    // Final SSSP result
    DistSeq dist;

    // Bundle construction
    VertexSeq R_vertices;                 // compact list of representatives
    parlay::sequence<uint8_t> in_R;       // membership bitmap

    NestDist local_dist;                  // local_dist[v] = [(u, d(v,u)), ...]
    NestV bundle;
    NestV ball;
    VertexSeq b;

public:
    // Bundle Construction and Bundle Relaxation
    ~ParallelBundleDijkstraSolver() override = default;

    void construct(const Graph& g, Vertex source);

    void relax(Vertex v, Distance cand, ArrayLaBPQ& pq);

    void solve(const Graph& g, Vertex source) override;

    Distance distance(Vertex v) const override {
        return dist[v].load(std::memory_order_relaxed);
    }

    size_t num_vertices() const override {
        return dist.size();
    }

    const char* name() const override {
        return "Parallel Bundle Dijkstra";
    }

    std::optional<Distance> find_dist(Vertex from, Vertex to) const {
        for (const auto& [u, d] : local_dist[from]) {
            if (u == to) {
                return d;
            }
        }
        return std::nullopt;
    }
};

#endif //SENIORTHESIS_PARALLEL_BUNDLE_DIJKSTRA_H
