#include <ctime>
#include <iostream>
#include <iomanip>

#include <tc/complex.hpp>
#include <tc/groups.hpp>

template<unsigned int N, class G>
void test(const G &group) {
    auto s = std::clock();
    auto combos = tc::combinations(group.gens, N - 1);
    auto data = tc::merge<N>(tc::hull<N>(group, combos, {}));
    auto e = std::clock();

    double diff = (double) (e - s) / CLOCKS_PER_SEC;
    int count = data.size();

    std::cout
        << std::setw(2) << N << ", "
        << std::setw(7) << group.name << ", "
        << std::setw(7) << count << ", "
        << std::fixed << std::setprecision(6) << diff << "s"
        << std::endl;
}

int main() {
    test<4>(tc::group::H(4));
    test<4>(tc::group::B(4));
    test<4>(tc::group::B(5));
    test<4>(tc::group::B(6));
    test<4>(tc::group::E(6));

    test<3>(tc::group::H(3));
    test<3>(tc::group::B(4));
    test<3>(tc::group::B(5));
    test<3>(tc::group::B(6));
    test<3>(tc::group::E(6));

    return 0;
}
