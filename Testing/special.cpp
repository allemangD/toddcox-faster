#include "common.hpp"

#include <iostream>

int main(int argc, char **argv) {
    auto name = std::string(argv[1]);
    auto vgens = parse_vec(argv[2]);
    auto target = std::stoul(argv[3]);

    unsigned int order;

    tc::Group group(0);

    if (name == "E6") {
        group = tc::group::E(6);
    }
    if (name == "E7") {
        group = tc::group::E(7);
    }
    if (name == "E8") {
        group = tc::group::E(8);
    }
    if (name == "B6") {
        group = tc::group::B(6);
    }
    if (name == "B7") {
        group = tc::group::B(7);
    }
    if (name == "B8") {
        group = tc::group::B(8);
    }

    order = compute(group, vgens);

    std::cout << "Order: " << order << ":" << target << std::endl;
    return order != target;
}
