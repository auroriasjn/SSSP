#ifndef SENIORTHESIS_UTILS_H
#define SENIORTHESIS_UTILS_H

#pragma once

#include <atomic>
#include <cstdlib>
#include <parlay/utilities.h>

#include "types.h"

constexpr std::size_t BLOCK_SIZE = 1024;

template <typename T>
inline bool write_min(std::atomic<T>& a, T b) {
    T cur = a.load(std::memory_order_relaxed);
    while (b < cur) {
        if (a.compare_exchange_weak(
                cur, b,
                std::memory_order_acq_rel,
                std::memory_order_relaxed)) {
            return true;
        }
    }
    return false;
}

template <typename T>
inline bool write_min_old(std::atomic<T>& a, T b, T& old_val) {
    T cur = a.load(std::memory_order_relaxed);
    while (b < cur) {
        if (a.compare_exchange_weak(
                cur, b,
                std::memory_order_acq_rel,
                std::memory_order_relaxed)) {
            old_val = cur;
            return true;
        }
    }
    old_val = cur;
    return false;
}

template <typename T>
inline bool write_max(std::atomic<T>& a, T b) {
    T cur = a.load(std::memory_order_relaxed);
    while (cur < b) {
        if (a.compare_exchange_weak(
                cur, b,
                std::memory_order_acq_rel,
                std::memory_order_relaxed)) {
            return true;
        }
    }
    return false;
}

template <typename T>
inline T hash_value(T a) {
    if constexpr (sizeof(T) == 4) {
        return parlay::hash32(a);
    } else if constexpr (sizeof(T) == 8) {
        return parlay::hash64(a);
    } else {
        std::abort();
    }
}

template <typename T>
inline T hash_value_2(T a) {
    if constexpr (sizeof(T) == 4) {
        return parlay::hash32_2(a);
    } else if constexpr (sizeof(T) == 8) {
        return parlay::hash64_2(a);
    } else {
        std::abort();
    }
}

#endif