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

namespace delayed = parlay::delayed;

class Graph;

class ParallelDijkstraSolver : public SSSPSolver {
private:
    std::vector<Distance> dist;

public:
    void solve(const Graph& g, Vertex source) override;
    const std::vector<Distance>& distances() const override {
        return dist;
    }

    const char* name() const override {
        return "Parallel Dijkstra";
    }
};

#endif //SENIORTHESIS_PARALLEL_DIJKSTRA_H
