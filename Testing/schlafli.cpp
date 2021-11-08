#include "common.hpp"

#include <iostream>

int main(int argc, char **argv) {
    auto vsymbol = parse_vec(argv[1]);
    auto vgens = parse_vec(argv[2]);
    auto target = std::stoul(argv[3]);

    tc::Symbol symbol(vsymbol.size());
    symbol << Eigen::Map<tc::Symbol>(vsymbol.data(), vsymbol.size());
    tc::Group group = tc::schlafli(symbol);
    auto order = compute(group, vgens);

    std::cout << "Order: " << order << ":" << target << std::endl;
    return order != target;
}
