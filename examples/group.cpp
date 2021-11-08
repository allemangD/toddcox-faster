#include <iostream>

#include <tc/solver.hpp>
#include <tc/groups.hpp>

int main() {
    tc::Symbol symbol(4);
    symbol << 5, 3, 2, 3;
    tc::Group group = tc::schlafli(symbol);
//    tc::Group group = tc::group::E(7);

    size_t srank = 3;

    tc::SubGroups subs = tc::subgroups(group, srank);

    std::cout << "Group " << group.name << " (" << subs.size() << " subgroups)" << std::endl;
    std::cout << group << std::endl;

    for (const auto &sub: subs) {
        for (int i = 0; i < srank; ++i) {
            for (int j = 0; j < srank; ++j) {
                auto sub_mult = sub(i, j);
                auto src_mult = group(sub.gens(i), sub.gens(j));

                if (sub_mult != src_mult) {
                    std::cout << "Incorrect subgroup " << sub.name << std::endl;
                    std::cout << sub << std::endl;
                    return 1;
                }
            }
        }
    }

    return 0;
}
