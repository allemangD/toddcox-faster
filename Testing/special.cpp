#include "common.hpp"

#include <iostream>

int main(int argc, char **argv) {
    auto name = std::string(argv[1]);
    auto vgens = parse_vec(argv[2]);
    auto target = std::stoul(argv[3]);

    unsigned int order;

    if (name == "E6") {
        auto group = tc::group::E<6>();
        order = compute<6>(group, vgens);
    }
    if (name == "E7") {
        auto group = tc::group::E<7>();
        order = compute<7>(group, vgens);
    }
    if (name == "E8") {
        auto group = tc::group::E<8>();
        order = compute<8>(group, vgens);
    }
    if (name == "B6") {
        auto group = tc::group::B<6>();
        order = compute<6>(group, vgens);
    }
    if (name == "B7") {
        auto group = tc::group::B<7>();
        order = compute<7>(group, vgens);
    }
    if (name == "B8") {
        auto group = tc::group::B<8>();
        order = compute<8>(group, vgens);
    }

    std::cout << "Order: " << order << ":" << target << std::endl;
    return order != target;
}
