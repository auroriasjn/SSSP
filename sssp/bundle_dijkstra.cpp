#include "../graph.h"
#include "bundle_dijkstra.h"

#include <cmath>
#include <queue>
#include <random>
#include <vector>
#include <cassert>
#include <iostream>

#ifndef NDEBUG
#define DBG(x) do { x; } while (0)
#else
#define DBG(x) do {} while (0)
#endif

namespace {
    inline bool is_invalid_vertex(Vertex v) {
        return v == Vertex(-1);
    }
}

void BundleDijkstraSolver::construct(const Graph& g, Vertex source) {
    const size_t n = g.num_vertices();

    size_t k = 1;
    if (n > 2) {
        const double ln_n = std::log(static_cast<double>(n));
        const double ln_ln_n = std::log(std::max(2.0, ln_n));
        const double val = std::sqrt(ln_n / ln_ln_n);
        k = std::max<size_t>(1, static_cast<size_t>(std::ceil(val)));
    }

    std::mt19937 gen(42);
    std::bernoulli_distribution coin(1.0 / static_cast<double>(k));

    const size_t L = std::max<size_t>(
            1,
            static_cast<size_t>(std::ceil(k * std::log(std::max<size_t>(2, k))))
    );

    in_R.assign(n, 0);
    ball.assign(n, {});
    bundle.assign(n, {});
    b.assign(n, Vertex(-1));
    dist_s.clear();
    local_dist.assign(n, {});
    dist_to_b.assign(n, INF); // FIXED: properly allocated to prevent Segfault
    reinsertions.assign(n, 0); // <-- ADD THIS LINE

    std::vector<Distance> d_buf(n, INF);
    std::vector<Vertex> d_dirty;

    std::vector<std::vector<std::pair<Vertex, Distance>>> extracted(n);
    for(size_t i = 0; i < n; ++i) {
        extracted[i].reserve(L);
    }

    std::vector<uint8_t> in_R1(n, 0);
    for (auto v : g.vertices()) {
        if (v == source || coin(gen)) {
            in_R1[v] = 1;
        }
    }

    std::vector<uint8_t> in_R2(n, 0);
    std::vector<PQNode> pq_container;
    pq_container.reserve(L * 4);

    for (auto v : g.vertices()) {
        if (in_R1[v]) continue;

        for (Vertex x : d_dirty) d_buf[x] = INF;
        d_dirty.clear();

        pq_container.clear();
        d_buf[v] = static_cast<Distance>(0);
        d_dirty.push_back(v);

        pq_container.push_back({0, v});

        size_t popped = 0;
        bool hit_R1 = false;

        while (!pq_container.empty() && popped < L) {
            std::pop_heap(pq_container.begin(), pq_container.end(), std::greater<>());
            auto [du, u] = pq_container.back();
            pq_container.pop_back();

            if (du != d_buf[u]) continue;

            extracted[v].emplace_back(u, du);
            ++popped;

            if (in_R1[u]) {
                hit_R1 = true;
                break;
            }

            for (const auto& edge : g.neighbors(u)) {
                const Vertex x = edge.to;
                const Distance nd = du + edge.weight;
                if (nd < d_buf[x]) {
                    if (d_buf[x] == INF) d_dirty.push_back(x);
                    d_buf[x] = nd;

                    pq_container.push_back({nd, x});
                    std::push_heap(pq_container.begin(), pq_container.end(), std::greater<>());
                }
            }
        }

        if (!hit_R1) in_R2[v] = 1;
    }

    for (Vertex x : d_dirty) d_buf[x] = INF;
    d_dirty.clear();

    for (size_t i = 0; i < n; ++i) {
        in_R[i] = in_R1[i] || in_R2[i];
    }

    for (auto v : g.vertices()) {
        if (in_R[v]) continue;

        local_dist[v].reserve(extracted[v].size());
        ball[v].reserve(extracted[v].size());

        bool found_rep = false;
        for (auto [u, du] : extracted[v]) {
            if (in_R[u]) {
                b[v] = u;
                dist_to_b[v] = du; // FIXED: Now safely caching the distance
                local_dist[v].emplace_back(u, du);
                found_rep = true;
                break;
            }
            ball[v].push_back(u);
            local_dist[v].emplace_back(u, du);
        }
        assert(found_rep && "No representative found for vertex outside R");
    }

    for (size_t u = 0; u < n; ++u) {
        bundle[u].clear();
        if (in_R[u]) {
            bundle[u].push_back(u);
            b[u] = u;
        }
    }
    for (auto v : g.vertices()) {
        if (in_R[v]) continue;
        assert(!is_invalid_vertex(b[v]));
        bundle[b[v]].push_back(v);
    }
}

void BundleDijkstraSolver::relax(
        Vertex v,
        Distance D,
        std::vector<PQNode>& pq_container
) {
    if (D >= dist_s[v]) return;
    dist_s[v] = D;

    if (in_R[v]) {
        pq_container.push_back({D, v});
        std::push_heap(pq_container.begin(), pq_container.end(), std::greater<>());
        return;
    }

    const Vertex root = b[v];
    if (!is_invalid_vertex(root)) {
        const Distance nextD = D + dist_to_b[v];
        if (nextD < dist_s[root]) {
            dist_s[root] = nextD;
            pq_container.push_back({nextD, root});
            std::push_heap(pq_container.begin(), pq_container.end(), std::greater<>());
        }
    }
}

void BundleDijkstraSolver::solve(const Graph& g, Vertex source) {
    const size_t n = g.num_vertices();

    construct(g, source);

    dist_s.assign(n, INF);
    dist_s[source] = static_cast<Distance>(0);

    std::vector<PQNode> pq_container;
    pq_container.reserve(n);

    for (size_t u = 0; u < n; ++u) {
        if (in_R[u]) {
            pq_container.push_back({dist_s[u], u});
        }
    }
    std::make_heap(pq_container.begin(), pq_container.end(), std::greater<>());

    while (!pq_container.empty()) {
        std::pop_heap(pq_container.begin(), pq_container.end(), std::greater<>());
        auto [du, u] = pq_container.back();
        pq_container.pop_back();

        if (du != dist_s[u]) continue;
        if (du == INF) break;
        if (!in_R[u]) continue;

        for (Vertex v : bundle[u]) {
            for (const auto& [y, d_vy] : local_dist[v]) {
                if (dist_s[y] < INF) {
                    relax(v, dist_s[y] + d_vy, pq_container);
                }
            }

            auto relax_from_z2 = [&](Vertex z2, Distance dist_z2_to_v) {
                for (const auto& edge : g.neighbors(z2)) {
                    const Vertex z1 = edge.to;
                    if (dist_s[z1] < INF) {
                        relax(v, dist_s[z1] + edge.weight + dist_z2_to_v, pq_container);
                    }
                }
            };

            for (const auto& [z2, d_vz2] : local_dist[v]) {
                if (z2 != u) {
                    relax_from_z2(z2, d_vz2);
                }
            }
            relax_from_z2(v, 0);
        }

        for (Vertex x : bundle[u]) {
            if (dist_s[x] >= INF) continue;

            for (const auto& edge : g.neighbors(x)) {
                const Vertex y = edge.to;
                const Distance through_x = dist_s[x] + edge.weight;

                relax(y, through_x, pq_container);

                if (!in_R[y]) {
                    for (const auto& [z1, d_yz1] : local_dist[y]) {
                        if (z1 != b[y]) {
                            relax(z1, through_x + d_yz1, pq_container);
                        }
                    }
                }
            }
        }
    }
}