//
// Created by Jeremy Ng on 3/2/26.
//

#include "../graph.h"
#include "bundle_dijkstra.h"

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

    VertexSet R1, R2;
    std::vector<std::vector<Vertex>> extracted(n);
    std::vector<std::unordered_map<Vertex, Distance>> temp_dist(n);

    R.clear();
    ball.assign(n, {});
    bundle.assign(n, {});
    b.assign(n, Vertex(-1));
    dist_s.clear();
    local_dist.assign(n, {});

    // Sample R1
    for (auto v : g.vertices()) {
        if (v == source || coin(gen)) {
            R1.insert(v);
        }
    }

    // Truncated Dijkstra for v in V \ R1
    for (auto v : g.vertices()) {
        if (R1.count(v)) {
            continue;
        }

        std::priority_queue<PQNode, std::vector<PQNode>, std::greater<>> pq;
        std::unordered_map<Vertex, Distance> d_v;

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

            extracted[v].push_back(u);
            temp_dist[v][u] = du;
            ++popped;

            if (R1.count(u)) {
                hit_R1 = true;
                break;
            }

            for (const auto& edge : g.neighbors(u)) {
                const Vertex x = edge.to;
                const Distance nd = du + edge.weight;

                auto it_x = d_v.find(x);
                if (it_x == d_v.end() || nd < it_x->second) {
                    d_v[x] = nd;
                    pq.emplace(nd, x);
                }
            }
        }

        if (!hit_R1) {
            R2.insert(v);
        }
    }

    // Finalize R = R1 ∪ R2
    R = R1;
    R.insert(R2.begin(), R2.end());

    // Compute b(v), Ball(v), local_dist[v][.]
    for (auto v : g.vertices()) {
        if (R.count(v)) {
            continue;
        }

        bool found_rep = false;
        for (Vertex u : extracted[v]) {
            if (R.count(u)) {
                b[v] = u;
                local_dist[v][u] = temp_dist[v][u];
                found_rep = true;
                break;
            }

            ball[v].push_back(u);
            local_dist[v][u] = temp_dist[v][u];
        }

        assert(found_rep && "No representative found for vertex outside R");
    }

    // Build Bundle(u) = { v in V \ R : b(v) = u }
    for (auto u : R) {
        bundle[u].clear();
        bundle[u].push_back(u);
        b[u] = u;
    }

    for (auto v : g.vertices()) {
        if (R.count(v)) {
            continue;
        }
        assert(!is_invalid_vertex(b[v]));
        bundle[b[v]].push_back(v);
    }

#ifndef NDEBUG
    {
        size_t num_R = R.size();
        size_t num_not_R = 0;
        size_t num_nonempty_bundles = 0;
        size_t total_bundle_members = 0;

        for (auto v: g.vertices()) {
            if (!R.count(v)) {
                ++num_not_R;
            }
            if (!bundle[v].empty()) {
                ++num_nonempty_bundles;
                total_bundle_members += bundle[v].size();
            }
        }

        std::cerr << "|R| = " << num_R << "\n";
        std::cerr << "|V \\ R| = " << num_not_R << "\n";
        std::cerr << "# nonempty bundles = " << num_nonempty_bundles << "\n";
        std::cerr << "# total bundle members = " << total_bundle_members << "\n";

        size_t owner_count_sum = 0;
        for (Vertex u: g.vertices()) {
            owner_count_sum += bundle[u].size();
        }

        assert(owner_count_sum == n);

        for (Vertex v: g.vertices()) {
            if (!R.count(v)) {
                assert(!is_invalid_vertex(b[v]));
                assert(std::find(bundle[b[v]].begin(), bundle[b[v]].end(), v) != bundle[b[v]].end());
                assert(local_dist[v].find(b[v]) != local_dist[v].end());
            }
        }

        for (Vertex v: g.vertices()) {
            if (!R.count(v)) {
                for (Vertex x: ball[v]) {
                    assert(local_dist[v].find(x) != local_dist[v].end());
                }
            }
        }
    }
#endif

}

void BundleDijkstraSolver::relax(
        Vertex v,
        Distance D,
        std::priority_queue<PQNode, std::vector<PQNode>, std::greater<>>& pq
) {
    DBG(std::cerr << "[relax] enter v=" << v << " D=" << D
                  << " dist_s[v]=" << dist_s[v] << "\n";);

    if (D >= dist_s[v]) {
        DBG(std::cerr << "[relax] skip (no improvement) v=" << v << "\n";);
        return;
    }

    dist_s[v] = D;

    DBG(std::cerr << "[relax] updated dist_s[" << v << "] = " << D << "\n";);

    // Case 1: v ∈ R → push to heap
    if (R.count(v)) {
        DBG(std::cerr << "[relax] v in R → push to pq: ("
                      << D << ", " << v << ")\n";);
        pq.emplace(D, v);
        return;
    }

    // Case 2: propagate to bundle root
    const Vertex root = b[v];

    DBG(std::cerr << "[relax] v not in R, root=" << root << "\n";);

    if (!is_invalid_vertex(root)) {
        auto it = local_dist[v].find(root);

        if (it != local_dist[v].end()) {
            Distance nextD = D + it->second;

            DBG(std::cerr << "[relax] propagate → root=" << root
                          << " edge=" << it->second
                          << " nextD=" << nextD << "\n";);

            relax(root, nextD, pq);
        } else {
            DBG(std::cerr << "[relax] WARNING: no local_dist entry for v="
                          << v << " root=" << root << "\n";);
        }
    } else {
        DBG(std::cerr << "[relax] WARNING: invalid root for v=" << v << "\n";);
    }

#undef DBG
}

void BundleDijkstraSolver::solve(const Graph& g, Vertex source) {
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
    for (auto u : R) {
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
        for (Vertex v : bundle[u]) {
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
            for (Vertex y : ball[v]) {
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

                for (const auto& edge : g.neighbors(z2)) {
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

            for (Vertex z2 : ball[v]) {
                relax_from_z2(z2);
            }
            relax_from_z2(v);
        }

        // -------------------------
        // Step 2
        // -------------------------
        for (Vertex x : bundle[u]) {
            if (dist_s[x] >= INF) {
                DBG(std::cerr << "[solve][step2] skip x=" << x
                              << " because dist_s[x]=INF\n";);
                continue;
            }

            DBG(std::cerr << "[solve][step2] expand x=" << x
                          << " dist_s[x]=" << dist_s[x]
                          << " deg(x)=" << g.neighbors(x).size() << "\n";);

            for (const auto& edge : g.neighbors(x)) {
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

                    for (Vertex z1 : ball[y]) {
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