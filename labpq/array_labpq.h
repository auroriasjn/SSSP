#pragma once

#include <atomic>
#include <algorithm>
#include <limits>
#include <vector>

#include <parlay/primitives.h>
#include <parlay/sequence.h>

#include "../graph.h"
#include "../sssp/sssp_solver.h"

#ifndef NDEBUG
#define DBG(x) do { x; } while (0)
#else
#define DBG(x) do {} while (0)
#endif

class ArrayLaBPQ {
public:
    using DistAtomSeq = parlay::sequence<std::atomic<Distance>>;
    using VertexSeq   = parlay::sequence<Vertex>;

private:
    DistAtomSeq* delta_;                          // authoritative distances
    parlay::sequence<std::atomic<bool>> in_q_;   // membership flags
    std::atomic<std::size_t> active_count_{0};

public:
    explicit ArrayLaBPQ(DistAtomSeq& delta)
            : delta_(&delta), in_q_(delta.size()) {
        parlay::parallel_for(0, in_q_.size(), [&](std::size_t i) {
            in_q_[i].store(false, std::memory_order_relaxed);
        });
    }

    void clear() {
        parlay::parallel_for(0, in_q_.size(), [&](std::size_t i) {
            in_q_[i].store(false, std::memory_order_relaxed);
        });
        active_count_.store(0, std::memory_order_relaxed);
    }

    bool empty() const {
        return active_count_.load(std::memory_order_acquire) == 0;
    }

    std::size_t size() const {
        return active_count_.load(std::memory_order_acquire);
    }

    void update(Vertex v) {
        bool expected = false;
        if (in_q_[v].compare_exchange_strong(
                expected, true,
                std::memory_order_acq_rel,
                std::memory_order_relaxed)) {
            active_count_.fetch_add(1, std::memory_order_acq_rel);
            DBG(std::cerr << "[LaB-PQ] updating vertex " << v << std::endl;);
        }
    }

    // Returns and removes all active vertices with delta[v] <= theta.
    VertexSeq extract(Distance theta) {
        auto ids = parlay::tabulate(in_q_.size(), [](std::size_t i) {
            return static_cast<Vertex>(i);
        });

        auto out = parlay::filter(ids, [&](Vertex v) {
            if (!in_q_[v].load(std::memory_order_acquire)) return false;
            Distance dv = (*delta_)[v].load(std::memory_order_relaxed);
            return dv <= theta;
        });

        parlay::parallel_for(0, out.size(), [&](std::size_t i) {
            Vertex v = out[i];
            bool was_present = in_q_[v].exchange(false, std::memory_order_acq_rel);
            if (was_present) {
                active_count_.fetch_sub(1, std::memory_order_acq_rel);
            }
        });

        return out;
    }

    // Snapshot all currently active vertices without removing them.
    VertexSeq active_vertices() const {
        auto ids = parlay::tabulate(in_q_.size(), [](std::size_t i) {
            return static_cast<Vertex>(i);
        });

        return parlay::filter(ids, [&](Vertex v) {
            return in_q_[v].load(std::memory_order_acquire);
        });
    }
};