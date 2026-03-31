#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
#include <functional>
#include <cstdlib>
#include <charconv>

#include "benchmark/benchmark_solver.h"
#include "benchmark/source_selector.h"
#include "parser/parser.h"

#include "sssp/sssp_solver.h"
#include "sssp/dijkstra.h"
#include "sssp/parallel_dijkstra.h"
#include "sssp/parallel_bundle_dijkstra.h"
#include "sssp/bundle_dijkstra.h"
#include "sssp/rho_stepping.h"
#include "sssp/bellman_ford.h"

#include "graph.h"

int main(int argc, char* argv[]) {
    if (argc == 1) {
        std::cerr <<
                  "Usage: " << argv[0] << " -i input_file [-o output_file] [-a algorithm] [-w] [-v] [-n number]\n"
                                          "Options:\n"
                                          "\t-i  input file path (required)\n"
                                          "\t-o  output file path\n"
                                          "\t-a  algorithm: [dijkstra] [parallel-dijkstra] [bd] [parallel-bd] [rho-stepping] [bf]\n"
                                          "\t-w  weighted graph\n"
                                          "\t-v  verify result\n"
                                          "\t-n  number of sources to analyze\n";
        return EXIT_FAILURE;
    }

    // ---------------------------
    // Argument Parsing
    // ---------------------------
    std::unordered_map<std::string, std::string> options;
    std::vector<std::string> flags;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg[0] == '-') {
            // Flag without value
            if (arg == "-w" || arg == "-v") {
                flags.push_back(arg);
            }
            // Option with value
            else if (arg == "-i" || arg == "-a" || arg == "-n" || arg == "-o") {
                if (i + 1 >= argc) {
                    std::cerr << "Error: Missing value for " << arg << "\n";
                    return EXIT_FAILURE;
                }
                options[arg] = argv[++i];
            }
            else {
                std::cerr << "Error: Unknown option " << arg << "\n";
                return EXIT_FAILURE;
            }
        }
        else {
            std::cerr << "Error: Unexpected argument " << arg << "\n";
            return EXIT_FAILURE;
        }
    }

    // ---------------------------
    // Extract Options
    // ---------------------------

    if (options.find("-i") == options.end()) {
        std::cerr << "Error: -i input file is required\n";
        return EXIT_FAILURE;
    }

    std::string in_file = options["-i"];
    std::string out_file = options.count("-o")
                           ? options["-o"]
                           : "";

    std::string algorithm = options.count("-a")
                            ? options["-a"]
                            : "dijkstra";

    size_t n = 100;
    if (options.count("-n")) {
        const std::string& s = options["-n"];

        auto [ptr, ec] = std::from_chars(
                s.data(),
                s.data() + s.size(),
                n
        );
        if (ec != std::errc() || ptr != s.data() + s.size()) {
            std::cerr << "Invalid value for -n: " << s << "\n";
            std::exit(EXIT_FAILURE);
        }
    }

    bool weighted = std::find(flags.begin(), flags.end(), "-w") != flags.end();
    bool verify   = std::find(flags.begin(), flags.end(), "-v") != flags.end();

    // ---------------------------
    // Load Graph
    // ---------------------------
    Graph g = Parser::parse(in_file, weighted);

    auto sources = SourceSelector::random_subset(g, 100);
    // ---------------------------
    // Solver Factory
    // ---------------------------

    std::function<std::unique_ptr<SSSPSolver>()> factory;

    if (algorithm == "dijkstra") {
        factory = []() {
            return std::make_unique<DijkstraSolver>();
        };
    } else if (algorithm == "parallel-dijkstra") {
        factory = []() {
            return std::make_unique<ParallelDijkstraSolver>();
        };
    } else if (algorithm == "bd") {
        factory = []() {
            return std::make_unique<BundleDijkstraSolver>();
        };
    } else if (algorithm == "parallel-bd") {
        factory = []() {
            return std::make_unique<ParallelBundleDijkstraSolver>();
        };
    } else if (algorithm == "rho-stepping") {
        factory = []() {
            return std::make_unique<RhoSteppingSolver>();
        };
    } else if (algorithm == "bf") {
        factory = []() {
            return std::make_unique<BellmanFordSolver>();
        };
    } else {
        std::cerr << "Error: Unknown algorithm '" << algorithm << "'\n";
        return EXIT_FAILURE;
    }

    // ---------------------------
    // Run Benchmark
    // ---------------------------

    BenchmarkRunner::run(
            in_file,
            out_file,
            factory,
            sources,
            weighted,
            verify
    );

    return EXIT_SUCCESS;
}