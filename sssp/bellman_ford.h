//
// Created by Jeremy Ng on 3/30/26.
//

#ifndef SENIORTHESIS_BELLMAN_FORD_H
#define SENIORTHESIS_BELLMAN_FORD_H


#include "sssp_solver.h"

class Graph;

class BellmanFordSolver : public SSSPSolver {

private:
    std::vector<Distance> dist;
    std::vector<size_t> reinsertions;

public:
    void solve(const Graph& g, Vertex source) override;

    Distance distance(Vertex v) const override {
        return dist[v];
    }

    std::size_t num_vertices() const override {
        return dist.size();
    }

    const char* name() const override {
        return "Bellman-Ford";
    }

    size_t reinserts(Vertex v) const override {
        return reinsertions[v];
    }
};


#endif //SENIORTHESIS_BELLMAN_FORD_H
