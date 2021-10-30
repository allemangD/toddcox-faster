#include "tc/core.hpp"
#include "tc/groups.hpp"

#include <ctime>
#include <iostream>
#include <iomanip>

template<class G>
void test(const G &group) {
    auto s = std::clock();
    auto cosets = group.solve();
    auto e = std::clock();

    double diff = (double) (e - s) / CLOCKS_PER_SEC;
    int order = cosets.size();

    std::cout
        << std::setw(7) << group.name << ", "
        << std::setw(7) << order << ", "
        << std::fixed << std::setprecision(6) << diff << "s"
        << std::endl;
}

int main() {
    test(tc::group::H(2));
    test(tc::group::H(3));
    test(tc::group::H(4));
    test(tc::group::T(100));
    test(tc::group::T(500));
    test(tc::group::T(1000));
    test(tc::group::E(6));
    test(tc::group::E(7));
    test(tc::group::B(6));
    test(tc::group::B(7));
    test(tc::group::B(8));

    return 0;
}
