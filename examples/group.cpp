#include <iostream>

#include <tc/solver.hpp>
#include <tc/groups.hpp>

int main() {
    constexpr unsigned int Rank = 5;
    constexpr unsigned int SRank = 3;

    tc::Group<Rank> group = tc::schlafli<Rank>({5, 3, 2, 3});
//    tc::Group<Rank> group = tc::group::E<Rank>();

    tc::SubGroups<Rank, SRank> subs = tc::subgroups<Rank, SRank>(group);

    std::cout << "Group " << group.name << " (" << subs.size() << " subgroups)" << std::endl;
    std::cout << group << std::endl;

    for (const auto &sub: subs) {
        for (int i = 0; i < SRank; ++i) {
            for (int j = 0; j < SRank; ++j) {
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
