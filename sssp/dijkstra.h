//
// Created by Jeremy Ng on 2/23/26.
//

#ifndef SENIORTHESIS_DIJKSTRA_H
#define SENIORTHESIS_DIJKSTRA_H

#include <vector>
#include <queue>
#include <limits>
#include <utility>

#include "sssp_solver.h"

class Graph;

class DijkstraSolver : public SSSPSolver {

private:
    std::vector<Distance> dist;
    std::vector<std::size_t> reinsertions;

public:
    void solve(const Graph& g, Vertex source) override;

    Distance distance(Vertex v) const override {
        return dist[v];
    }

    std::size_t num_vertices() const override {
        return dist.size();
    }

    const char* name() const override {
        return "Dijkstra";
    }

    std::size_t reinserts(Vertex v) const override {
        return reinsertions[v];
    }
};

#endif //SENIORTHESIS_DIJKSTRA_H
