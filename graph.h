//
// Created by Jeremy Ng on 2/23/26.
//
#ifndef SENIORTHESIS_GRAPH_H
#define SENIORTHESIS_GRAPH_H

#include <unordered_map>
#include <vector>
#include <stdexcept>

#include "types.h"

class Graph {
public:
    struct Edge {
        Vertex to;
        Weight weight;
    };

private:
    bool directed_ = false;

    // external <-> internal mapping
    std::unordered_map<ExtVertex, Vertex> ext_to_int_;
    std::vector<ExtVertex> int_to_ext_;  // index -> external id

    std::vector<std::vector<Edge>> adj_;
    std::uint64_t n_edges_ = 0; // logical edges (undirected counted once)

public:
    explicit Graph(bool is_directed = false) : directed_(is_directed) {}

    bool directed() const { return directed_; }

    std::uint64_t num_vertices() const { return adj_.size(); }
    std::uint64_t num_edges() const { return n_edges_; }

    // Ensure a vertex exists internally; returns internal index.
    Vertex intern(ExtVertex id) {
        auto it = ext_to_int_.find(id);
        if (it != ext_to_int_.end()) return it->second;

        auto idx = static_cast<Vertex>(adj_.size());
        ext_to_int_.emplace(id, idx);
        int_to_ext_.push_back(id);
        adj_.emplace_back(); // new adjacency list
        return idx;
    }

    // Add by external IDs (what your parser reads)
    void add_edge(ExtVertex u_id, ExtVertex v_id, Weight w = 1.0) {
        Vertex u = intern(u_id);
        Vertex v = intern(v_id);

        adj_[u].push_back({v, w});
        if (!directed_) {
            adj_[v].push_back({u, w});
        }
        ++n_edges_;
    }

    // Core adjacency access for solvers (internal vertices)
    const std::vector<Edge>& neighbors(Vertex u) const {
        return adj_[u];
    }

    // Optional helpers if you need to translate in/out
    ExtVertex external_id(Vertex u) const { return int_to_ext_[u]; }
    Vertex internal_id(ExtVertex id) const {
        auto it = ext_to_int_.find(id);
        if (it == ext_to_int_.end()) throw std::out_of_range("Unknown external vertex id");
        return it->second;
    }

    std::vector<Vertex> vertices() const {
        std::vector<Vertex> v;
        v.reserve(adj_.size());
        for (Vertex i = 0; i < adj_.size(); ++i) v.push_back(i);
        return v;
    }
};

#endif // SENIORTHESIS_GRAPH_H