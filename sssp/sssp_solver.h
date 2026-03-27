#ifndef SENIORTHESIS_SSSPALGORITHM_H
#define SENIORTHESIS_SSSPALGORITHM_H

#include "../types.h"

class Graph;

class SSSPSolver {
public:
    virtual ~SSSPSolver() = default;

    // Solve SSSP from source
    virtual void solve(const Graph& g, Vertex source) = 0;

    // Access results
    virtual const std::vector<Distance>& distances() const = 0;

    // Optional: solver name for benchmarking output
    virtual const char* name() const = 0;
};

#endif