//
// Created by Jeremy Ng on 3/24/26.
//

#ifndef SENIORTHESIS_PARALLEL_TYPES_H
#define SENIORTHESIS_PARALLEL_TYPES_H

#include <parlay/sequence.h>
#include <parlay/primitives.h>
#include "types.h"

using VertexSeq = parlay::sequence<Vertex>;
using DistPair = std::pair<Vertex, Distance>;
using DistMapSeq = parlay::sequence<DistPair>;
using NestV = parlay::sequence<VertexSeq>;
using NestDist = parlay::sequence<DistMapSeq>;

#endif //SENIORTHESIS_PARALLEL_TYPES_H
