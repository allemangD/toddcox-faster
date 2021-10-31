#include "tc/core.hpp"
#include "tc/groups.hpp"

#include <iostream>

int main() {
    auto cube = tc::group::B(3);
    auto vars = cube.solve();

    std::string start;
    std::vector<std::string> names = {"a", "b", "c"};
    auto words = vars.path.walk(start, names, std::plus<>());

    for (const auto &word: words) {
        std::cout << (word.empty() ? "-" : word) << std::endl;
    }

    return 0;
}
