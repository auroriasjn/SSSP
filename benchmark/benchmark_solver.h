//
// Created by Jeremy Ng on 2/23/26.
//
#ifndef SENIORTHESIS_BENCHMARK_SOLVER_H
#define SENIORTHESIS_BENCHMARK_SOLVER_H

#include <memory>
#include <functional>
#include <string>
#include <vector>
#include <cstdint>
#include <fstream>
#include <iostream>

#include "../graph.h"

class SSSPSolver;

class BenchmarkRunner {
public:
    using SolverFactory =
            std::function<std::unique_ptr<SSSPSolver>()>;

    static void run(
            const std::string& in_file,
            const std::string& out_file,
            const SolverFactory& factory,
            const std::vector<Vertex>& sources,
            bool weighted,
            bool verify
    );

private:
    static Graph preprocess(const std::string& file_name, bool weighted);
    static void write_result(
            const SSSPSolver& solver,
            const std::string& out_file,
            ExtVertex source,
            double elapsed_seconds
    );
};

#endif //SENIORTHESIS_BENCHMARK_SOLVER_H
