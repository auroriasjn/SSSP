//
// Created by Jeremy Ng on 3/23/26.
//

#include "../types.h"
#include "../graph.h"
#include "parallel_bundle_dijkstra.h"

#include <cmath>
#include <queue>
#include <unordered_map>
#include <vector>
#include <iostream>

namespace {
    inline bool is_invalid_vertex(Vertex v) {
        return v == Vertex(-1);
    }

    inline uint64_t mix64(uint64_t x) {
        x += 0x9e3779b97f4a7c15ULL;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        return x ^ (x >> 31);
    }

    inline size_t get_k(size_t n) {
        // Base case
        if (n <= 2) return 1;

        // k logic
        const double ln_n = std::log(static_cast<double>(n));
        const double ln_ln_n = std::log(std::max(2.0, ln_n));
        const double val = std::sqrt(ln_n / ln_ln_n);

        DBG(std::cerr << "k: " << std::ceil(val) << std::endl;);

        return std::max<size_t>(1, static_cast<size_t>(std::ceil(val)));
    }
}

void ParallelBundleDijkstraSolver::construct(const Graph& g, Vertex source) {
    const size_t n = g.num_vertices();
    auto verts = g.vertices();

    // Same k as sequential version
    size_t k = get_k(n);
    const size_t L = std::max<size_t>(
            1,
            static_cast<size_t>(std::ceil(k * std::log(std::max<size_t>(2, k))))
    );

    const uint64_t seed = 42;
    const double p = 1.0 / static_cast<double>(k);

    // Initializers
    parlay::sequence<uint8_t> in_R1(n, 0);
    parlay::sequence<uint8_t> in_R2(n, 0);
    in_R = parlay::sequence<uint8_t>(n, 0);

    ball = NestV(n);
    bundle = NestV(n);
    b = VertexSeq(n, Vertex(-1));
    local_dist = NestDist(n);

    NestV extracted(n);
    NestDist temp_dist(n);

    // ---- Phase 1: sample R1 ----
    parlay::parallel_for(0, verts.size(), [&](std::size_t i) {
        Vertex v = verts[i];
        uint64_t h = mix64(static_cast<uint64_t>(v) ^ seed);
        double x = static_cast<double>(h) /
                   static_cast<double>(std::numeric_limits<uint64_t>::max());

        if (v == source || x < p) {
            in_R1[v] = 1;
        }
    });

    // ---- Phase 2: truncated Dijkstra from each v not in R1 ----
    TIME_BLOCK("Truncated Dijkstras", {
        parlay::parallel_for(0, verts.size(), [&](std::size_t i) {
            Vertex v = verts[i];
            if (in_R1[v]) return;

            std::priority_queue<PQNode, std::vector<PQNode>, std::greater<>> pq;
            std::unordered_map<Vertex, Distance> d_v;
            VertexList ext;
            std::vector<std::pair<Vertex, Distance>> td;

            pq.emplace(static_cast<Distance>(0), v);
            d_v[v] = static_cast<Distance>(0);

            size_t popped = 0;
            bool hit_R1 = false;

            while (!pq.empty() && popped < L) {
                auto [du, u] = pq.top();
                pq.pop();

                auto it_u = d_v.find(u);
                if (it_u == d_v.end() || du != it_u->second) {
                    continue;
                }

                ext.push_back(u);
                td.emplace_back(u, du);
                ++popped;

                if (in_R1[u]) {
                    hit_R1 = true;
                    break;
                }

                for (const auto &edge: g.neighbors(u)) {
                    const Vertex x = edge.to;
                    const Distance nd = du + static_cast<Distance>(edge.weight);

                    auto it_x = d_v.find(x);
                    if (it_x == d_v.end() || nd < it_x->second) {
                        d_v[x] = nd;
                        pq.emplace(nd, x);
                    }
                }
            }

            if (!hit_R1) {
                in_R2[v] = 1;
            }

            extracted[v] = VertexSeq(ext.begin(), ext.end());
            temp_dist[v] = DistMapSeq(td.begin(), td.end());
        });
    });

    // ---- Phase 3: finalize R = R1 ∪ R2 ----
    parlay::parallel_for(0, n, [&](size_t v) {
        uint8_t is_in_R = static_cast<uint8_t>(in_R1[v] || in_R2[v]);
        in_R[v] = is_in_R;

        if (is_in_R) {
            b[v] = v;
        }
    });

    // ---- Phase 4: compute b(v), ball(v), local_dist[v] for v not in R ----
    TIME_BLOCK("Bundle Computation", {
        parlay::parallel_for(0, verts.size(), [&](std::size_t i) {
            Vertex v = verts[i];
            if (in_R[v]) {
                ball[v] = VertexSeq{};
                local_dist[v] = DistMapSeq{{v, static_cast<Distance>(0)}};
                return;
            }

            const auto &ext = extracted[v];
            const auto &td = temp_dist[v];

            VertexList bv;
            std::vector<std::pair<Vertex, Distance>> ldv;

            for (size_t j = 0; j < ext.size(); ++j) {
                Vertex u = ext[j];
                Distance du = td[j].second;

                if (in_R[u]) {
                    b[v] = u;
                    ldv.emplace_back(u, du);
                    break;
                }

                bv.push_back(u);
                ldv.emplace_back(u, du);
            }
            ball[v] = VertexSeq(bv.begin(), bv.end());
            local_dist[v] = DistMapSeq(ldv.begin(), ldv.end());
        });
    });

    // ---- Phase 5: build bundle(u) = {u} ∪ {v notin R : b(v)=u} ----
    auto assign_pairs = parlay::delayed_tabulate(n, [&](size_t v) {
        return std::pair<size_t, Vertex>(static_cast<size_t>(b[v]), static_cast<Vertex>(v));
    });
    bundle = parlay::group_by_index(assign_pairs, n);
}

void ParallelBundleDijkstraSolver::relax(Vertex v, Distance cand, ArrayLaBPQ& pq) {
    Distance old = dist[v].load(std::memory_order_relaxed);

    while (cand < old) {
        if (dist[v].compare_exchange_weak(old, cand,
                                          std::memory_order_acq_rel,
                                          std::memory_order_relaxed)) {

            successful_relax_per_vertex[v]++;

            if (in_R[v]) {
                reinsertions[v]++;
                pq.update(v);
            } else {
                Vertex root = b[v];
                if (!is_invalid_vertex(root)) {
                    recursive_root_hits[root]++;

                    pq.update(root);

                    auto off = find_dist(v, root);
                    if (off.has_value()) {
                        relax(root, cand + *off, pq);
                    }
                }
            }
            return;
        }
    }
}

void ParallelBundleDijkstraSolver::solve(const Graph& g, Vertex source) {
    const std::size_t n = g.num_vertices();

    TIME_BLOCK("construct", construct(g, source));

    TIME_BLOCK("storage", {
        dist = DistSeq(n);
        reinsertions = SizeSeq(n);
        successful_relax_per_vertex = SizeSeq(n);
        recursive_root_hits = SizeSeq(n);

        parlay::parallel_for(0, n, [&](std::size_t i) {
            dist[i].store(INF, std::memory_order_relaxed);
            reinsertions[i].store(0, std::memory_order_relaxed);
            successful_relax_per_vertex[i].store(0, std::memory_order_relaxed);
            recursive_root_hits[i].store(0, std::memory_order_relaxed);
        });
        dist[source].store(0, std::memory_order_relaxed);
    });

    ArrayLaBPQ pq(dist);
    const std::size_t batch_size = get_k(n);

    // Match reference semantics: all representatives are active initially.
    TIME_BLOCK("parallel-update", {
        parlay::parallel_for(0, n, [&](std::size_t i) {
            if (in_R[i]) {
                pq.update(i);
            }
        });
    });

    while (!pq.empty()) {
        VertexSeq active = pq.active_vertices();
        if (active.empty()) break;

        // Only one thread can call this
        const Distance theta = pq.get_threshold(active, dist, batch_size);
        VertexSeq frontier = pq.extract(theta);
        if (frontier.empty()) break;

        parlay::sequence<Task> tasks =
                parlay::flatten(parlay::tabulate(frontier.size(), [&](std::size_t i) {
                    const Vertex u = frontier[i];
                    const VertexSeq& owned = bundle[u];
                    return parlay::tabulate(owned.size(), [&](std::size_t j) {
                        return Task{u, owned[j]};
                    });
                }));

        TIME_BLOCK("Step 1", {
            parlay::parallel_for(0, tasks.size(), [&](std::size_t t) {
                const auto [u, v] = tasks[t];

                const Distance du = dist[u].load(std::memory_order_relaxed);
                if (du == INF) return;

                if (v == u) {
                    relax(v, du, pq);
                } else {
                    auto it_uv = find_dist(v, u);
                    if (it_uv.has_value()) {
                        relax(v, du + *it_uv, pq);
                    }
                }

                for (Vertex y : ball[v]) {
                    const Distance dy = dist[y].load(std::memory_order_relaxed);
                    if (dy == INF) continue;

                    auto it_vy = find_dist(v, y);
                    if (it_vy.has_value()) {
                        relax(v, dy + *it_vy, pq);
                    }
                }

                auto relax_from_z2 = [&](Vertex z2) {
                    Distance dist_z2_to_v = 0;
                    if (z2 != v) {
                        auto it = find_dist(v, z2);
                        if (!it.has_value()) return;
                        dist_z2_to_v = *it;
                    }

                    for (const auto& edge : g.neighbors(z2)) {
                        const Vertex z1 = edge.to;
                        const Distance dz1 = dist[z1].load(std::memory_order_relaxed);
                        if (dz1 == INF) continue;

                        relax(v, dz1 + edge.weight + dist_z2_to_v, pq);
                    }
                };

                for (Vertex z2 : ball[v]) {
                    relax_from_z2(z2);
                }
                relax_from_z2(v);
            });
        });

        TIME_BLOCK("Step 2", {
            parlay::parallel_for(0, tasks.size(), [&](std::size_t t) {
                const auto [u, x] = tasks[t];

                const Distance du = dist[u].load(std::memory_order_relaxed);
                if (du == INF) return;

                const Distance dx = dist[x].load(std::memory_order_relaxed);
                if (dx == INF) return;

                for (const auto &edge: g.neighbors(x)) {
                    const Vertex y = edge.to;
                    const Distance through_x = dx + edge.weight;

                    relax(y, through_x, pq);

                    if (!in_R[y]) {
                        for (Vertex z1: ball[y]) {
                            auto it_yz1 = find_dist(y, z1);
                            if (it_yz1.has_value()) {
                                relax(z1, through_x + *it_yz1, pq);
                            }
                        }
                    }
                }
            });
        });
    }

    auto print_top = [&](const char* label, const SizeSeq& arr) {
        std::vector<std::pair<Vertex, size_t>> vals;

        for (Vertex v = 0; v < arr.size(); ++v) {
            if (arr[v] > 0) vals.emplace_back(v, arr[v]);
        }

        std::sort(vals.begin(), vals.end(),
                  [](auto& a, auto& b) { return a.second > b.second; });

        std::cerr << label << "\n";
        for (size_t i = 0; i < std::min<size_t>(10, vals.size()); ++i) {
            std::cerr << "  v=" << vals[i].first
                      << " count=" << vals[i].second << "\n";
        }
    };

    DBG(std::cerr << "\n=== Debug ===\n";);
    DBG(print_top("Top reinserted vertices:", reinsertions););
    DBG(print_top("Top successful relax:", successful_relax_per_vertex););
    DBG(print_top("Top recursive root hits:", recursive_root_hits););
}