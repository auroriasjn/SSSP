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
#include "../labpq/array_labpq.h"

class Graph;

class RhoSteppingSolver : public SSSPSolver {
private:
    using VertexSeq = parlay::sequence<Vertex>;
    using DistSeq = parlay::sequence<std::atomic<Distance>>;

    std::vector<Distance> dist;
    std::size_t rho_ = 1 << 20;
    uint32_t seed_ = 0;

    Distance get_threshold(const VertexSeq& frontier, const DistSeq& dist_a);

public:
    explicit RhoSteppingSolver(std::size_t rho = (1u << 20)) : rho_(rho) {}

    void solve(const Graph& g, Vertex source) override;

    const char* name() const override {
        return "Rho Stepping";
    }

    const std::vector<Distance>& distances() const override {
        return dist;
    }
};

#endif //SENIORTHESIS_RHO_STEPPING_H
