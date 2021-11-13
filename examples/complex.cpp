#include <ctime>
#include <iostream>
#include <iomanip>

#include <tc/complex.hpp>
#include <tc/groups.hpp>

template<class G>
void test(const G &group, unsigned int N) {
    auto s = std::clock();
    auto combos = tc::combinations(group.gens, N - 1);
    auto data = tc::merge(tc::hull(group, combos));
    auto e = std::clock();

    double diff = (double) (e - s) / CLOCKS_PER_SEC;
    int count = data.cols();

    std::cout
        << std::setw(2) << N << ", "
        << std::setw(7) << group.name << ", "
        << std::setw(7) << count << ", "
        << std::fixed << std::setprecision(6) << diff << "s"
        << std::endl;
}

int main() {
    test(tc::group::H(4), 4);
    test(tc::group::B(4), 4);
    test(tc::group::B(5), 4);
    test(tc::group::B(6), 4);
    test(tc::group::E(6), 4);
    test(tc::group::H(3), 3);
    test(tc::group::B(4), 3);
    test(tc::group::B(5), 3);
    test(tc::group::B(6), 3);
    test(tc::group::E(6), 3);

    return 0;
}
