//
// Created by Jeremy Ng on 3/2/26.
//

#ifndef SENIORTHESIS_BUNDLE_DIJKSTRA_H
#define SENIORTHESIS_BUNDLE_DIJKSTRA_H

#include <algorithm>
#include <queue>
#include <unordered_set>
#include <vector>
#include <random>
#include <cstdint>

#include "sssp_solver.h"

class BundleDijkstraSolver : public SSSPSolver {
private:
    // Final SSSP result
    std::vector<Distance> dist_s;
    std::vector<size_t> reinsertions;

    // Bundle construction
    VertexSet R;
    std::vector<std::unordered_map<Vertex, Distance>> local_dist;

    std::vector<VertexList> bundle;
    std::vector<VertexList> ball;
    std::vector<Vertex> b;

public:
    // Bundle Construction and Bundle Relaxation
    ~BundleDijkstraSolver() = default;

    void construct(const Graph& g, Vertex source);
    void relax(Vertex v, Distance d, std::priority_queue<PQNode, std::vector<PQNode>, std::greater<>> &pq);

    void solve(const Graph& g, Vertex source) override;

    Distance distance(Vertex v) const override {
        return dist_s[v];
    }

    std::size_t num_vertices() const override {
        return dist_s.size();
    }

    const char* name() const override {
        return "Bundle Dijkstra";
    }

    size_t reinserts(Vertex v) const override {
        return reinsertions[v];
    }
};

#endif //SENIORTHESIS_BUNDLE_DIJKSTRA_H
