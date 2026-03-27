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

using ExtVertex = std::uint64_t;
using Vertex    = std::uint32_t;
using Weight    = std::uint64_t;
using Distance  = std::uint64_t;

using VertexList = std::vector<Vertex>;
using VertexSet = std::unordered_set<Vertex>;
using PQNode = std::pair<Distance, Vertex>;

constexpr Distance INF = std::numeric_limits<Distance>::max();

#endif //SENIORTHESIS_TYPES_H
