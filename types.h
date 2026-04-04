//
// Created by Jeremy Ng on 3/2/26.
//

#ifndef SENIORTHESIS_TYPES_H
#define SENIORTHESIS_TYPES_H

#include <cstdint>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <limits>

#ifndef NDEBUG
#define TIME_BLOCK(name, code)                                  \
    do {                                                        \
        auto _start = std::chrono::high_resolution_clock::now();\
        code;                                                   \
        auto _end = std::chrono::high_resolution_clock::now();  \
        double _elapsed =                                       \
            std::chrono::duration<double>(_end - _start).count();\
        std::cerr << "[TIMER] " << name << ": " << _elapsed     \
                  << " seconds\n";                              \
    } while (0)
#else
#define TIME_BLOCK(name, code) code
#endif

#ifndef NDEBUG
#define DBG(x) do { x; } while (0)
#else
#define DBG(x) do {} while (0)
#endif

using ExtVertex = std::uint64_t;
using Vertex    = std::uint32_t;
using Weight    = std::uint64_t;
using Distance  = std::uint64_t;

using VertexList = std::vector<Vertex>;
using VertexSet = std::unordered_set<Vertex>;
using Task = std::pair<Vertex, Vertex>;
using PQNode = std::pair<Distance, Vertex>;

constexpr Distance INF = std::numeric_limits<Distance>::max();

#endif //SENIORTHESIS_TYPES_H
