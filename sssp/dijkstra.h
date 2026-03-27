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
public:
    void solve(const Graph& g, Vertex source) override;

    const std::vector<Distance>& distances() const override {
        return dist;
    }

    const char* name() const override {
        return "Dijkstra";
    }
};

#endif //SENIORTHESIS_DIJKSTRA_H
