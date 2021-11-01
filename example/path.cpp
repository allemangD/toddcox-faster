#include <vector>
#include <string>
#include <iostream>

#include <tc/solver.hpp>
#include <tc/groups.hpp>

int main() {
    auto cube = tc::group::B<3>();
    auto vars = tc::solve(cube);

    std::string start;
    std::vector<std::string> names = {"a", "b", "c"};
    auto words = vars.path().walk(start, names, std::plus<>());

    for (const auto &word: words) {
        std::cout << (word.empty() ? "-" : word) << std::endl;
    }

    return 0;
}
