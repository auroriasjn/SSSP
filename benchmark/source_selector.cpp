//
// Created by Jeremy Ng on 2/23/26.
//

#include <algorithm>
#include "source_selector.h"
#include "../graph.h"

// Helper method to choose first k vertices to get shortest paths for.
std::vector<Vertex> SourceSelector::first_k(const Graph& g, size_t k) {
    std::vector<Vertex> result;
    result.reserve(k);

    for (auto v : g.vertices()) {
        if (result.size() == k) break;
        result.push_back(v);
    }
    return result;
}

std::vector<Vertex> SourceSelector::random_subset(const Graph& g, size_t k, uint64_t seed) {
    const size_t n = g.num_vertices();
    if (k >= n) {
        std::vector<Vertex> all;
        all.reserve(n);
        for (Vertex i = 0; i < n; ++i) all.push_back(i);
        return all;
    }

    std::vector<Vertex> verts;
    verts.reserve(n);
    for (Vertex i = 0; i < n; ++i) verts.push_back(i);

    std::mt19937_64 rng(seed);
    std::shuffle(verts.begin(), verts.end(), rng);
    verts.resize(k);
    return verts;
}