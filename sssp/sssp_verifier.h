#ifndef SENIORTHESIS_SSSP_VERIFIER_H
#define SENIORTHESIS_SSSP_VERIFIER_H

#include <cassert>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>

#include "../graph.h"
#include "sssp_solver.h"
#include "dijkstra.h"

class SSSPVerifier {
public:
    static void verify(
            const Graph& g,
            Vertex source,
            const SSSPSolver& candidate,
            const std::string& out_file
    ) {
        std::ostream* out = &std::cerr;
        std::ofstream file;

        if (!out_file.empty()) {
            file.open(out_file, std::ios::binary | std::ios::app);
            if (!file) {
                std::cerr << "Failed to open output file: "
                          << out_file << "\n";
                std::exit(EXIT_FAILURE);
            }
            out = &file;
        }

        using clock = std::chrono::high_resolution_clock;

        DijkstraSolver reference;

        auto start = clock::now();
        reference.solve(g, source);
        auto end = clock::now();

        double elapsed =
                std::chrono::duration<double>(end - start).count();

        (*out) << "  Reference Dijkstra time: "
               << elapsed << " sec\n";

        const size_t n = reference.num_vertices();
        assert(n == candidate.num_vertices());

        size_t mismatch_count = 0;

        for (size_t i = 0; i < n; ++i) {
            const Distance expected = reference.distance(static_cast<Vertex>(i));
            const Distance actual = candidate.distance(static_cast<Vertex>(i));

            if (expected != actual) {
                (*out) << "  Mismatch at vertex " << i
                       << " | expected: " << expected
                       << " | actual: " << actual
                       << "\n";
                mismatch_count++;
            }
        }

        if (mismatch_count > 0) {
            (*out) << "Total mismatches: "
                   << mismatch_count << "\n";
        }

        assert(mismatch_count == 0 && "SSSP verification failed");
    }
};

#endif // SENIORTHESIS_SSSP_VERIFIER_H