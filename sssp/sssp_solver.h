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
    virtual size_t num_vertices() const = 0;
    virtual Distance distance(Vertex v) const = 0;

    // Optional: solver name for benchmarking output
    virtual const char* name() const = 0;

    // Relaxations
    virtual size_t reinserts(Vertex v) const = 0;
};

#endif