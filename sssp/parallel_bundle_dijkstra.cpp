//
// Created by Jeremy Ng on 3/23/26.
//

#include "../graph.h"
#include "parallel_bundle_dijkstra.h"

#include <cmath>
#include <queue>
#include <random>
#include <unordered_map>
#include <unordered_set>
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

    inline uint64_t mix64(uint64_t x) {
        x += 0x9e3779b97f4a7c15ULL;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        return x ^ (x >> 31);
    }
}

void ParallelBundleDijkstraSolver::construct(const Graph& g, Vertex source) {
    const size_t n = g.num_vertices();

    // Same k as sequential version
    size_t k = 1;
    if (n > 2) {
        const double ln_n = std::log(static_cast<double>(n));
        const double ln_ln_n = std::log(std::max(2.0, ln_n));
        const double val = std::sqrt(ln_n / ln_ln_n);
        k = std::max<size_t>(1, static_cast<size_t>(std::ceil(val)));
    }

    const size_t L = std::max<size_t>(
            1,
            static_cast<size_t>(std::ceil(k * std::log(std::max<size_t>(2, k))))
    );

    const uint64_t seed = 42;
    const double p = 1.0 / static_cast<double>(k);

    auto verts = g.vertices();

    // ---- initialize parallel state ----
    parlay::sequence<uint8_t> in_R1(n, 0);
    parlay::sequence<uint8_t> in_R2(n, 0);
    in_R = parlay::sequence<uint8_t>(n, 0);

    ball = NestV(n);
    bundle = NestV(n);
    b = VertexSeq(n, Vertex(-1));
    local_dist = NestDist(n);

    dist_s.assign(n, INF);
    dist_s[source] = 0;

    // Temporary per-vertex results from truncated Dijkstra
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

    VertexSeq R1_vertices = parlay::filter(verts, [&](Vertex v) {
        return in_R1[v];
    });

    // ---- Phase 2: truncated Dijkstra from each v not in R1 ----
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

            for (const auto& edge : g.neighbors(u)) {
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

    // ---- Phase 3: finalize R = R1 ∪ R2 ----
    parlay::parallel_for(0, n, [&](size_t v) {
        in_R[v] = static_cast<uint8_t>(in_R1[v] || in_R2[v]);
    });

    R_vertices = parlay::filter(verts, [&](Vertex v) {
        return in_R[v];
    });

    // Representatives are their own bundle roots
    parlay::parallel_for(0, R_vertices.size(), [&](size_t i) {
        Vertex u = R_vertices[i];
        b[u] = u;
    });

    // ---- Phase 4: compute b(v), ball(v), local_dist[v] for v not in R ----
    parlay::parallel_for(0, verts.size(), [&](std::size_t i) {
        Vertex v = verts[i];
        if (in_R[v]) {
            // Optional: define representative vertices as having singleton ball/local-dist
            ball[v] = VertexSeq{};
            local_dist[v] = DistMapSeq{{v, static_cast<Distance>(0)}};
            return;
        }

        const auto& ext = extracted[v];
        const auto& td  = temp_dist[v];

        VertexList bv;
        std::vector<std::pair<Vertex, Distance>> ldv;

        bool found_rep = false;

        for (size_t j = 0; j < ext.size(); ++j) {
            Vertex u = ext[j];
            Distance du = td[j].second;  // same traversal order as ext

            if (in_R[u]) {
                b[v] = u;
                ldv.emplace_back(u, du);
                found_rep = true;
                break;
            }

            bv.push_back(u);
            ldv.emplace_back(u, du);
        }

        if (!found_rep) {
#ifndef NDEBUG
            std::cerr << "[construct] No representative found for v=" << v << "\n";
            std::abort();
#endif
        }

        ball[v] = VertexSeq(bv.begin(), bv.end());
        local_dist[v] = DistMapSeq(ldv.begin(), ldv.end());
    });

    // ---- Phase 5: build bundle(u) = {u} ∪ {v notin R : b(v)=u} ----
    auto nonR_vertices = parlay::filter(verts, [&](Vertex v) {
        return !in_R[v];
    });

    auto assign_pairs = parlay::map(nonR_vertices, [&](Vertex v) {
        return std::pair<size_t, Vertex>(static_cast<size_t>(b[v]), v);
    });

    auto grouped = parlay::group_by_index(assign_pairs, n);

    // grouped[u] contains all v with b[v] = u
    bundle = NestV(n);
    parlay::parallel_for(0, n, [&](size_t u) {
        std::vector<Vertex> bu;

        if (in_R[u]) {
            bu.push_back(static_cast<Vertex>(u));  // keep representative in own bundle
        }

        for (Vertex v : grouped[u]) {
            bu.push_back(v);
        }

        bundle[u] = VertexSeq(bu.begin(), bu.end());
    });
}

void ParallelBundleDijkstraSolver::relax(
        Vertex v,
        Distance D,
        std::priority_queue<PQNode, std::vector<PQNode>, std::greater<>>& pq
) {
    // TODO: FIX THIS
    DBG(std::cerr << "[relax] enter v=" << v
                  << " D=" << D
                  << " dist_s[v]=" << dist_s[v] << "\n";);

    if (D >= dist_s[v]) {
        DBG(std::cerr << "[relax] skip (no improvement) v=" << v << "\n";);
        return;
    }

    dist_s[v] = D;

    DBG(std::cerr << "[relax] updated dist_s[" << v << "] = " << D << "\n";);

    // Case 1: v ∈ R
    if (in_R[v]) {
        DBG(std::cerr << "[relax] v in R -> push to pq: ("
                      << D << ", " << v << ")\n";);
        pq.emplace(D, v);
        return;
    }

    // Case 2: v ∉ R, propagate to representative/root
    const Vertex root = b[v];

    DBG(std::cerr << "[relax] v not in R, root=" << root << "\n";);

    if (is_invalid_vertex(root)) {
        DBG(std::cerr << "[relax] WARNING: invalid root for v=" << v << "\n";);
        return;
    }

    auto delta = find_dist(v, root);
    if (!delta.has_value()) {
        DBG(std::cerr << "[relax] WARNING: no local_dist entry for v="
                      << v << " root=" << root << "\n";);
        return;
    }

    const Distance nextD = D + *delta;

    DBG(std::cerr << "[relax] propagate -> root=" << root
                  << " edge=" << *delta
                  << " nextD=" << nextD << "\n";);

    relax(root, nextD, pq);
}

void ParallelBundleDijkstraSolver::solve(const Graph& g, Vertex source) {
    const size_t n = g.num_vertices();

#ifndef NDEBUG
#define DBG(x) do { x; } while (0)
#else
#define DBG(x) do {} while (0)
#endif
    // Step 1 of bundle construction
    construct(g, source);

    DBG(std::cerr << "[solve] source in R? " << R.count(source) << "\n";);
    DBG(std::cerr << "[solve] bundle[source].size()=" << bundle[source].size() << "\n";);
    DBG(std::cerr << "[solve] source_in_bundle[source]? "
                  << (std::find(bundle[source].begin(), bundle[source].end(), source) != bundle[source].end())
                  << "\n";);

    dist_s.assign(n, INF);
    dist_s[source] = static_cast<Distance>(0);

    DBG(std::cerr << "[solve] source=" << source
                  << " in R? " << R.count(source)
                  << " bundle[source].size()=" << bundle[source].size()
                  << " source_in_own_bundle? "
                  << std::count(bundle[source].begin(), bundle[source].end(), source)
                  << "\n";);

    std::priority_queue<PQNode, std::vector<PQNode>, std::greater<>> pq;

#ifndef NDEBUG
    std::vector<size_t> pop_count(n, 0);
    std::vector<size_t> relax_count(n, 0);

    DBG(std::cerr << "[solve] source=" << source
                  << " n=" << n
                  << " |R|=" << R.size() << "\n";);
#endif

    // Algorithm 1 initializes the heap with all vertices of R, keyed by d(.)
    for (auto u: R) {
        pq.emplace(dist_s[u], u);  // source gets 0, others get INF
        DBG(std::cerr << "[solve] initial pq push u=" << u
                      << " key=" << dist_s[u] << "\n";);
    }

    while (!pq.empty()) {
        auto [du, u] = pq.top();
        pq.pop();

#ifndef NDEBUG
        pop_count[u]++;
        DBG(std::cerr << "\n[solve] pop u=" << u
                      << " du=" << du
                      << " dist_s[u]=" << dist_s[u]
                      << " pq_size=" << pq.size()
                      << " pop_count=" << pop_count[u] << "\n";);
#endif

        // Skip stale heap entries.
        if (du != dist_s[u]) {
            DBG(std::cerr << "[solve] stale pop, skip u=" << u << "\n";);
            continue;
        }

        if (du == INF) {
            DBG(std::cerr << "[solve] infinite pop, skip u=" << u << "\n";);
            break;
        }

        // Only R-vertices should drive the main loop.
        if (!R.count(u)) {
            DBG(std::cerr << "[solve] WARNING: popped non-R vertex u=" << u << "\n";);
            continue;
        }

        DBG(std::cerr << "[solve] process u=" << u
                      << " bundle[u].size()=" << bundle[u].size() << "\n";);

        // -------------------------
        // Step 1
        // -------------------------
        for (Vertex v: bundle[u]) {
            DBG(std::cerr << "[solve][step1] bundle root u=" << u
                          << " visiting v=" << v
                          << " ball[v].size()=" << ball[v].size()
                          << " dist_s[v]=" << dist_s[v] << "\n";);

            // Relax(v, d(u) + dist(u,v))
            {
                auto it_uv = local_dist[v].find(u); // dist(v,u) = dist(u,v)
                if (it_uv != local_dist[v].end() && dist_s[u] < INF) {
                    Distance cand = dist_s[u] + it_uv->second;
#ifndef NDEBUG
                    relax_count[v]++;
#endif
                    DBG(std::cerr << "[solve][step1] direct from u: relax("
                                  << v << ", " << cand << ")"
                                  << " via local_dist[" << v << "][" << u << "]="
                                  << it_uv->second << "\n";);
                    relax(v, cand, pq);
                } else {
                    DBG(std::cerr << "[solve][step1] skip direct from u for v=" << v
                                  << " reason="
                                  << (it_uv == local_dist[v].end() ? "missing local_dist " : "")
                                  << (dist_s[u] >= INF ? "dist_s[u]=INF" : "")
                                  << "\n";);
                }
            }

            // For y in Ball(v): Relax(v, d(y) + dist(y,v))
            for (Vertex y: ball[v]) {
                auto it_vy = local_dist[v].find(y);
                if (it_vy != local_dist[v].end() && dist_s[y] < INF) {
                    Distance cand = dist_s[y] + it_vy->second;
#ifndef NDEBUG
                    relax_count[v]++;
#endif
                    DBG(std::cerr << "[solve][step1] from ball y=" << y
                                  << " -> relax(" << v << ", " << cand << ")"
                                  << " with dist_s[y]=" << dist_s[y]
                                  << " local_dist[" << v << "][" << y << "]="
                                  << it_vy->second << "\n";);
                    relax(v, cand, pq);
                } else {
                    DBG(std::cerr << "[solve][step1] skip ball y=" << y
                                  << " for v=" << v
                                  << " reason="
                                  << (it_vy == local_dist[v].end() ? "missing local_dist " : "")
                                  << (dist_s[y] >= INF ? "dist_s[y]=INF" : "")
                                  << "\n";);
                }
            }

            // For z2 in Ball(v) U {v}, for z1 in N(z2):
            // Relax(v, d(z1) + w(z1,z2) + dist(z2,v))
            auto relax_from_z2 = [&](Vertex z2) {
                auto dist_z2_to_v = static_cast<Distance>(0);

                if (z2 != v) {
                    auto it = local_dist[v].find(z2);
                    if (it == local_dist[v].end()) {
                        DBG(std::cerr << "[solve][step1] WARNING: missing local_dist["
                                      << v << "][" << z2 << "]\n";);
                        return;
                    }
                    dist_z2_to_v = it->second;
                }

                DBG(std::cerr << "[solve][step1] relax_from_z2 z2=" << z2
                              << " -> v=" << v
                              << " dist_z2_to_v=" << dist_z2_to_v
                              << " deg(z2)=" << g.neighbors(z2).size() << "\n";);

                for (const auto &edge: g.neighbors(z2)) {
                    const Vertex z1 = edge.to;
                    if (dist_s[z1] < INF) {
                        Distance cand = dist_s[z1] + edge.weight + dist_z2_to_v;
#ifndef NDEBUG
                        relax_count[v]++;
#endif
                        DBG(std::cerr << "[solve][step1] z1=" << z1
                                      << " z2=" << z2
                                      << " -> relax(" << v << ", " << cand << ")"
                                      << " with dist_s[z1]=" << dist_s[z1]
                                      << " w=" << edge.weight
                                      << " dist(z2,v)=" << dist_z2_to_v << "\n";);
                        relax(v, cand, pq);
                    } else {
                        DBG(std::cerr << "[solve][step1] skip z1=" << z1
                                      << " because dist_s[z1]=INF\n";);
                    }
                }
            };

            for (Vertex z2: ball[v]) {
                relax_from_z2(z2);
            }
            relax_from_z2(v);
        }

        // -------------------------
        // Step 2
        // -------------------------
        for (Vertex x: bundle[u]) {
            if (dist_s[x] >= INF) {
                DBG(std::cerr << "[solve][step2] skip x=" << x
                              << " because dist_s[x]=INF\n";);
                continue;
            }

            DBG(std::cerr << "[solve][step2] expand x=" << x
                          << " dist_s[x]=" << dist_s[x]
                          << " deg(x)=" << g.neighbors(x).size() << "\n";);

            for (const auto &edge: g.neighbors(x)) {
                const Vertex y = edge.to;
                const Distance through_x = dist_s[x] + edge.weight;

#ifndef NDEBUG
                relax_count[y]++;
#endif
                DBG(std::cerr << "[solve][step2] x=" << x
                              << " -> y=" << y
                              << " w=" << edge.weight
                              << " through_x=" << through_x
                              << " dist_s[y]=" << dist_s[y] << "\n";);

                relax(y, through_x, pq);

                // Ball(y) is defined only when y not in R.
                if (!R.count(y)) {
                    DBG(std::cerr << "[solve][step2] y=" << y
                                  << " not in R, ball[y].size()=" << ball[y].size() << "\n";);

                    for (Vertex z1: ball[y]) {
                        auto it_yz1 = local_dist[y].find(z1);
                        if (it_yz1 != local_dist[y].end()) {
                            Distance cand = through_x + it_yz1->second;
#ifndef NDEBUG
                            relax_count[z1]++;
#endif
                            DBG(std::cerr << "[solve][step2] via ball(y): y=" << y
                                          << " z1=" << z1
                                          << " -> relax(" << z1 << ", " << cand << ")"
                                          << " local_dist[" << y << "][" << z1 << "]="
                                          << it_yz1->second << "\n";);
                            relax(z1, cand, pq);
                        } else {
                            DBG(std::cerr << "[solve][step2] WARNING: missing local_dist["
                                          << y << "][" << z1 << "]\n";);
                        }
                    }
                } else {
                    DBG(std::cerr << "[solve][step2] y=" << y << " in R, skip ball(y)\n";);
                }
            }
        }
    }

#ifndef NDEBUG
    DBG(std::cerr << "\n[solve] finished\n";);

    size_t expanded_vertices = 0;
    size_t repeated_pops = 0;
    for (size_t i = 0; i < n; ++i) {
        if (pop_count[i] > 0) expanded_vertices++;
        if (pop_count[i] > 1) {
            repeated_pops++;
            std::cerr << "[solve] repeated pop vertex=" << i
                      << " count=" << pop_count[i]
                      << " final_dist=" << dist_s[i] << "\n";
        }
        if (relax_count[i] > 20) {
            std::cerr << "[solve] high relax_count vertex=" << i
                      << " count=" << relax_count[i]
                      << " final_dist=" << dist_s[i] << "\n";
        }
    }

    std::cerr << "[solve] expanded_vertices=" << expanded_vertices
              << " repeated_pops=" << repeated_pops << "\n";
#endif

#undef DBG
}