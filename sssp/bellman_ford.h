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
public:
    void solve(const Graph& g, Vertex source) override;

    const std::vector<Distance>& distances() const override {
        return dist;
    }

    const char* name() const override {
        return "Bellman-Ford";
    }
};


#endif //SENIORTHESIS_BELLMAN_FORD_H
