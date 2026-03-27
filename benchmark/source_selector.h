//
// Created by Jeremy Ng on 2/23/26.
//

#ifndef SENIORTHESIS_SOURCE_SELECTOR_H
#define SENIORTHESIS_SOURCE_SELECTOR_H

#include <random>
#include "../types.h"

class Graph;

class SourceSelector {
public:
    static std::vector<Vertex>
    random_subset(const Graph& g, size_t k, std::uint64_t seed = 42);

    static std::vector<Vertex>
    first_k(const Graph& g, size_t k);
};

#endif //SENIORTHESIS_SOURCE_SELECTOR_H
