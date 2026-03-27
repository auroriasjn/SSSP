//
// Created by Jeremy Ng on 3/23/26.
//

#ifndef SENIORTHESIS_LABPQ_H
#define SENIORTHESIS_LABPQ_H

#include "../types.h"

#include <parlay/primitives.h>
#include <parlay/sequence.h>
#include <parlay/delayed.h>

class LaBPQ {
public:
    using VertexSeq = parlay::sequence<Vertex>;
    using DistSeq   = parlay::sequence<std::atomic<Distance>>;

    virtual ~LaBPQ() = default;

    virtual void update(Vertex v) = 0;
    virtual VertexSeq extract(Distance theta) = 0;
    virtual bool empty() const = 0;
};


#endif //SENIORTHESIS_LABPQ_H
