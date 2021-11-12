#include <iostream>
#include <tc/complex.hpp>
#include <tc/groups.hpp>

int main() {
    tc::Symbol symbol(3);
    symbol << 5, 3, 3;
    auto group = tc::schlafli(symbol);

    constexpr int N = 4;
    std::vector<tc::Symbol> combos = tc::combinations(group.gens, N - 1);
    auto data = tc::merge<N>(tc::hull<4>(group, combos, {}));

    std::cout << data.size() << std::endl;

    return 0;
}
