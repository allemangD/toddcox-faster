#pragma once

#include <string>
#include <vector>
#include <sstream>

#include <tc/group.hpp>
#include <tc/groups.hpp>
#include <tc/solver.hpp>

tc::Symbol parse_vec(const std::string &part) {
    std::istringstream iss(part);

    std::vector<unsigned int> vec;

    std::string token;
    while (std::getline(iss, token, ' ')) {
        vec.push_back(std::stoul(token));
    }

    return Eigen::Map<tc::Symbol>(vec.data(), vec.size());
}

size_t compute(
    const tc::Group &group,
    const tc::Symbol &gens
) {
    auto table = tc::solve(group, gens);
    return table.order();
}
