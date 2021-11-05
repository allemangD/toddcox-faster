#include "common.hpp"

#include <iostream>

int main(int argc, char **argv) {
    auto vsymbol = parse_vec(argv[1]);
    auto vgens = parse_vec(argv[2]);
    auto target = std::stoul(argv[3]);

    int order;
    switch (vsymbol.size()) {
    case 0: {
        Eigen::Vector<unsigned int, 0> symbol(vsymbol.data());
        order = compute<1>(tc::schlafli<1>(symbol), vgens);
        break;
    }
    case 1: {
        Eigen::Vector<unsigned int, 1> symbol(vsymbol.data());
        order = compute<2>(tc::schlafli<2>(symbol), vgens);
        break;
    }
    case 2: {
        Eigen::Vector<unsigned int, 2> symbol(vsymbol.data());
        order = compute<3>(tc::schlafli<3>(symbol), vgens);
        break;
    }
    case 3: {
        Eigen::Vector<unsigned int, 3> symbol(vsymbol.data());
        order = compute<4>(tc::schlafli<4>(symbol), vgens);
        break;
    }
    case 4: {
        Eigen::Vector<unsigned int, 4> symbol(vsymbol.data());
        order = compute<5>(tc::schlafli<5>(symbol), vgens);
        break;
    }
    }

    std::cout << "Order: " << order << ":" << target << std::endl;
    return order != target;
}
