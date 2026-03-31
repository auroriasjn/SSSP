cmake_minimum_required(VERSION 3.28)
project(SeniorThesis)

set(CMAKE_CXX_STANDARD 17)
find_package(Parlay REQUIRED)

add_executable(SeniorThesis sssp.cpp
        sssp/dijkstra.cpp
        sssp/sssp_solver.h
        sssp/sssp_verifier.h
        benchmark/source_selector.cpp
        benchmark/benchmark_solver.cpp
        parser/parser.cpp
        sssp/parallel_dijkstra.cpp
        sssp/rho_stepping.cpp
        sssp/rho_stepping.h
        sssp/bundle_dijkstra.cpp
        sssp/bundle_dijkstra.h
        types.h
        sssp/parallel_bundle_dijkstra.cpp
        sssp/parallel_bundle_dijkstra.h
        labpq/labpq.cpp
        labpq/labpq.h
        labpq/array_labpq.cpp
        labpq/array_labpq.h
        utils.h
        sssp/parallel_bundle_dijkstra.cpp
        parallel_types.h
        sssp/parallel_bundle_dijkstra.cpp
        sssp/bellman_ford.cpp
        sssp/bellman_ford.h
        sssp/bellman_ford.h
        sssp/parallel_bundle_dijkstra.cpp)
