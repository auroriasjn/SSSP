#include "parser.h"

#include <iostream>
#include <fstream>
#include <string>
#include <charconv>
#include <stdexcept>

// Parser makes the graph *undirected*
Graph Parser::parse(const std::string& file_name, const bool weighted) {
    std::ifstream file(file_name, std::ios::binary);
    if (!file)
        throw std::runtime_error("Failed to open file");

    Graph graph;
    for (std::string line; std::getline(file, line); ) {
        if (line.empty()) continue;

        // Expected format:
        // + U V timestamp
        if (line[0] != '+')
            continue;
        const char* ptr = line.data() + 2;

        Vertex u, v;
        unsigned int w;

        // Parse u
        auto [p1, ec1] = std::from_chars(ptr, ptr + line.size(), u);
        if (ec1 != std::errc()) continue;

        ptr = p1 + 1;

        // Parse v
        auto [p2, ec2] = std::from_chars(ptr, ptr + line.size(), v);
        if (ec2 != std::errc()) continue;

        ptr = p2 + 1;

        // Parse timestamp
        if (weighted) {
            auto [p3, ec3] = std::from_chars(ptr, ptr + line.size(), w);
            if (ec3 != std::errc()) continue;

            graph.add_edge(u, v, w);
        } else {
            // Unweighted edge add
            graph.add_edge(u, v, 1.0);
        }
    }

    return graph;
}