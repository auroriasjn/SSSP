//
// Created by Jeremy Ng on 2/23/26.
//

#include "benchmark_solver.h"
#include "../parser/parser.h"
#include "../sssp/sssp_solver.h"
#include "../sssp/sssp_verifier.h"

Graph BenchmarkRunner::preprocess(const std::string& file_name, bool weighted) {
    auto start = std::chrono::high_resolution_clock::now();

    Graph g = Parser::parse(file_name, weighted);

    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "Graph load time: "
              << std::chrono::duration<double>(end - start).count()
              << " sec\n";

    return g;
}

void BenchmarkRunner::run(
        const std::string& in_file,
        const std::string& out_file,
        const SolverFactory& factory,
        const VertexList& sources,
        bool weighted,
        bool verify
) {
    Graph g = preprocess(in_file, weighted);

    for (auto source : sources) {
        auto solver = factory();
        auto start = std::chrono::high_resolution_clock::now();

        solver->solve(g, source);

        auto end = std::chrono::high_resolution_clock::now();
        double elapsed = std::chrono::duration<double>(end - start).count();

        ExtVertex ext_source = g.external_id(source);
        write_result(*solver, out_file, ext_source, elapsed);

        if (verify) {
            SSSPVerifier::verify(g, source, *solver, out_file);
        }
    }
}

void BenchmarkRunner::write_result(
        const SSSPSolver& solver,
        const std::string& out_file,
        ExtVertex source_ext,
        double elapsed_seconds
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

    (*out)
            << "Solver: " << solver.name()
            << " | Source: " << source_ext
            << " | Time: " << elapsed_seconds
            << " sec\n";
}