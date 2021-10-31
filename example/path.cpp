#include "tc/core.hpp"
#include "tc/groups.hpp"

#include <iostream>

int main() {
    auto cube = tc::group::B(3);
    auto vars = cube.solve();

    auto words = vars.path.walk<std::string, std::string>(
        "",
        {"a", "b", "c"},
        [](auto a, auto b) { return a + b; }
    );

    for (const auto &word : words) {
        std::cout << word << std::endl;
    }

    return 0;
}
