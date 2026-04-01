//
// Created by Jeremy Ng on 2/23/26.
//

#ifndef SENIORTHESIS_PARALLEL_DIJKSTRA_H
#define SENIORTHESIS_PARALLEL_DIJKSTRA_H


#include <vector>
#include <queue>
#include <limits>
#include <utility>
#include <atomic>

#include <parlay/primitives.h>
#include <parlay/sequence.h>
#include <parlay/delayed.h>

#include "sssp_solver.h"
#include "../parallel_types.h"

namespace delayed = parlay::delayed;

class Graph;

class ParallelDijkstraSolver : public SSSPSolver {
private:
    DistSeq dist;
    SizeSeq reinsertions;

public:
    void solve(const Graph& g, Vertex source) override;

    Distance distance(Vertex v) const override {
        return dist[v].load(std::memory_order_relaxed);
    }

    size_t reinserts(Vertex v) const override {
        return reinsertions[v].load(std::memory_order_relaxed);;
    }

    size_t num_vertices() const override {
        return dist.size();
    }

    const char* name() const override {
        return "Parallel Dijkstra";
    }
};

#endif //SENIORTHESIS_PARALLEL_DIJKSTRA_H
