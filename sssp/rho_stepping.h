//
// Created by Jeremy Ng on 3/2/26.
//

#ifndef SENIORTHESIS_RHO_STEPPING_H
#define SENIORTHESIS_RHO_STEPPING_H

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
#include "../labpq/array_labpq.h"

class Graph;

class RhoSteppingSolver : public SSSPSolver {
private:
    DistSeq dist;
    std::size_t rho_ = 1 << 20;
    uint32_t seed_ = 0;

public:
    explicit RhoSteppingSolver(std::size_t rho = (1u << 20)) : rho_(rho) {}

    void solve(const Graph& g, Vertex source) override;

    Distance distance(Vertex v) const override {
        return dist[v].load(std::memory_order_relaxed);
    }

    size_t num_vertices() const override {
        return dist.size();
    }

    const char* name() const override {
        return "Rho Stepping";
    }
};

#endif //SENIORTHESIS_RHO_STEPPING_H
