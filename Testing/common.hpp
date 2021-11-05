#pragma once

#include <string>
#include <vector>
#include <sstream>

#include <tc/group.hpp>
#include <tc/groups.hpp>
#include <tc/solver.hpp>

std::vector<unsigned int> parse_vec(const std::string &part) {
    std::istringstream iss(part);

    std::vector<unsigned int> res;

    std::string token;
    while (std::getline(iss, token, ' ')) {
        res.push_back(std::stoul(token));
    }

    return res;
}

template<unsigned int N>
int compute(
    const tc::Group<N> &group,
    const std::vector<unsigned int> &vgens
) {
    auto table = tc::solve<N>(group, vgens);
    return table.order();
}