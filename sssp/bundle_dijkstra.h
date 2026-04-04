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
    std::vector<Distance> dist_s;
    std::vector<size_t> reinsertions;

    std::vector<uint8_t> in_R;

    // NEW: Caches the distance to the bundled vertex for O(1) relaxation
    std::vector<Distance> dist_to_b;

    std::vector<std::vector<std::pair<Vertex, Distance>>> local_dist;
    std::vector<VertexList> bundle;
    std::vector<VertexList> ball;
    std::vector<Vertex> b;

public:
    ~BundleDijkstraSolver() = default;

    void construct(const Graph& g, Vertex source);

    // FIXED: Signature now correctly matches the .cpp file's vector container
    void relax(Vertex v, Distance d, std::vector<PQNode>& pq_container);

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