#pragma once

#include <atomic>
#include <algorithm>
#include <limits>
#include <vector>
#include <array>

#include <parlay/primitives.h>
#include <parlay/sequence.h>

#include "../parallel_types.h"
#include "../graph.h"
#include "../sssp/sssp_solver.h"
#include "../utils.h"

#ifndef NDEBUG
#define DBG(x) do { x; } while (0)
#else
#define DBG(x) do {} while (0)
#endif

// Let's see if this is slightly quicker.
class ArrayLaBPQ {
private:
    DistSeq* delta_;
    parlay::sequence<std::atomic<bool>> in_q_;

    parlay::sequence<Vertex> active_;
    std::atomic<std::size_t> active_size_{0};

    uint32_t seed_ = 0;

public:
    explicit ArrayLaBPQ(DistSeq& delta)
            : delta_(&delta),
              in_q_(delta.size()),
              active_(delta.size()) {
        parlay::parallel_for(0, in_q_.size(), [&](std::size_t i) {
            in_q_[i].store(false, std::memory_order_relaxed);
        });
    }

    void clear() {
        std::size_t sz = active_size_.load(std::memory_order_acquire);
        parlay::parallel_for(0, sz, [&](std::size_t i) {
            in_q_[active_[i]].store(false, std::memory_order_relaxed);
        });
        active_size_.store(0, std::memory_order_release);
    }

    bool empty() const {
        return active_size_.load(std::memory_order_acquire) == 0;
    }

    std::size_t size() const {
        return active_size_.load(std::memory_order_acquire);
    }

    void update(Vertex v) {
        bool expected = false;
        if (in_q_[v].compare_exchange_strong(
                expected, true,
                std::memory_order_acq_rel,
                std::memory_order_relaxed)) {
            std::size_t pos = active_size_.fetch_add(1, std::memory_order_acq_rel);
            active_[pos] = v;
        }
    }

    VertexSeq extract(Distance theta) {
        std::size_t sz = active_size_.load(std::memory_order_acquire);

        auto candidates = parlay::tabulate(sz, [&](std::size_t i) {
            return active_[i];
        });

        auto ready = parlay::filter(candidates, [&](Vertex v) {
            return (*delta_)[v].load(std::memory_order_relaxed) <= theta;
        });

        auto deferred = parlay::filter(candidates, [&](Vertex v) {
            return (*delta_)[v].load(std::memory_order_relaxed) > theta;
        });

        parlay::parallel_for(0, ready.size(), [&](std::size_t i) {
            in_q_[ready[i]].store(false, std::memory_order_release);
        });

        parlay::parallel_for(0, deferred.size(), [&](std::size_t i) {
            active_[i] = deferred[i];
        });

        active_size_.store(deferred.size(), std::memory_order_release);
        return ready;
    }

    VertexSeq active_vertices() const {
        std::size_t sz = active_size_.load(std::memory_order_acquire);
        return parlay::tabulate(sz, [&](std::size_t i) {
            return active_[i];
        });
    }

    // Get the threshold
    Distance get_threshold(const VertexSeq& frontier, const DistSeq& dist_a, const size_t k) {
        constexpr std::size_t SSSP_SAMPLES = 1000;
        const std::size_t frontier_size = frontier.size();

        if (frontier_size == 0) {
            return INF;
        }

        if (frontier_size <= k) {
            DBG(std::cerr << "[get_threshold] frontier size less than rho.\n";);
            auto frontier_dist = parlay::delayed_tabulate(frontier_size, [&](std::size_t i) {
                return dist_a[frontier[i]].load(std::memory_order_relaxed);
            });
            return *parlay::max_element(frontier_dist);
        }

        std::array<Distance, SSSP_SAMPLES + 1> sample_dist{};
        for (std::size_t i = 0; i <= SSSP_SAMPLES; ++i) {
            Vertex v = frontier[hash_value(seed_ + static_cast<uint32_t>(i)) % frontier_size];
            sample_dist[i] = dist_a[v].load(std::memory_order_relaxed);
        }

        seed_ += static_cast<uint32_t>(SSSP_SAMPLES + 1);

        auto id = static_cast<std::size_t>(
                (static_cast<double>(k) / static_cast<double>(frontier_size)) * SSSP_SAMPLES
        );
        if (id > SSSP_SAMPLES) id = SSSP_SAMPLES;

        std::sort(sample_dist.begin(), sample_dist.end());
        return sample_dist[id];
    }
};