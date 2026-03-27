//
// Created by Jeremy Ng on 2/23/26.
//

#ifndef SENIORTHESIS_PARSER_H
#define SENIORTHESIS_PARSER_H

#include <string>
#include <cstdint>

#include "../graph.h"

class Parser {
public:
    Parser() = default;
    static Graph parse(const std::string& file_name, bool weighted);
};

#endif //SENIORTHESIS_PARSER_H
