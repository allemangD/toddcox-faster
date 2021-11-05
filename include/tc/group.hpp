#pragma once

#include <string>
#include <sstream>
#include <numeric>

#include <Eigen/Eigen>

namespace {
    template<class T>
    std::string stringify(const T &vec) {
        std::stringstream ss;
        ss << "[" << vec.transpose() << "]";
        return ss.str();
    }
}

namespace tc {
    template<unsigned int Rank>
    using Symbol = Eigen::Vector<unsigned int, Rank>;

    template<unsigned int Rank>
    Symbol<Rank> iota() {
        Symbol<Rank> res;
        for (int i = 0; i < Rank; ++i) {
            res(i) = i;
        }
        return res;
    }

    /// A Coxeter Matrix
    template<unsigned int Rank>
    class Group : public Eigen::Matrix<unsigned int, Rank, Rank> {
    public:
        using Base = Eigen::Matrix<unsigned int, Rank, Rank>;

        std::string name = "G";
        Symbol<Rank> gens = iota<Rank>();

        using Base::Base;
    };

    template<unsigned int Rank, unsigned int SRank>
    Group<SRank> subgroup(const Group<Rank> &group, Symbol<SRank> gens) {
        Group<SRank> res;
        res.name = group.name + ":" + stringify(gens);
        res.gens = gens;

        for (int i = 0; i < SRank; ++i) {
            for (int j = 0; j < SRank; ++j) {
                res(i, j) = group(gens[i], gens[j]);
            }
        }

        return res;
    }

    template<unsigned int Rank, unsigned int SRank>
    Symbol<Rank> inverse(Symbol<SRank> gens) {
        Symbol<Rank> res;
        res.fill(0);
        for (int i = 0; i < SRank; ++i) {
            res(gens(i)) = i;
        }
        return res;
    }

    constexpr unsigned int Factorial(unsigned int n) {
        unsigned int res = 1;
        for (int i = 1; i <= n; ++i) {
            res *= i;
        }
        return res;
    }

    constexpr unsigned int Choose(unsigned int n, unsigned int k) {
        return Factorial(n) / Factorial(k) / Factorial(n - k);
    }

    template<unsigned int Rank, unsigned int SRank>
    using SubGroups = std::array<Group<SRank>, Choose(Rank, SRank)>;

    template<unsigned int Rank, unsigned int SRank>
    SubGroups<Rank, SRank> subgroups(const Group<Rank> &group) {
        std::vector<bool> mask(Rank, false);
        std::fill_n(mask.begin(), SRank, true);

        SubGroups<Rank, SRank> res;

        size_t i = 0;
        Symbol<SRank> row;
        do {
            for (int j = 0, k = 0; j < Rank; ++j) {
                if (mask[j])
                    row(k++) = j;
            }
            res[i++] = subgroup<Rank, SRank>(group, row);
        } while (std::prev_permutation(mask.begin(), mask.end()));

        return res;
    }

    /**
     * Create a named coxeter matrix from a simplified schlafli symbol
     */
    template<unsigned int Rank>
    Group<Rank> schlafli(const Symbol<Rank - 1> &mults, const std::string &name) {
        Group<Rank> res;
        res.name = name;

        res.fill(2);
        res.diagonal().fill(1);
        res.topRightCorner(Rank - 1, Rank - 1).diagonal() << mults;
        res.bottomLeftCorner(Rank - 1, Rank - 1).diagonal() << mults;

        return res;
    }

    /**
     * Create a coxeter matrix from a simplified schlafli symbol.
     */
    template<unsigned int Rank>
    Group<Rank> schlafli(const Symbol<Rank - 1> &mults) {
        return schlafli<Rank>(mults, stringify(mults));
    }

    template<unsigned int GR, unsigned int HR>
    Group<GR + HR> product(const Group<GR> &g, const Group<HR> &h) {
        Group<GR + HR> res;
        res.name = g.name + "*" + h.name;

        res.fill(2);

        int off = 0;
        res.block(off, off, GR, GR) << g.array() + off;
        off += GR;

        res.block(off, off, HR, HR) << h.array() + off;
        off += HR;

        return res;
    }

    template<unsigned int GR, unsigned int P>
    Group<GR * P> power(const Group<GR> &g) {
        Group<GR * P> res;
        res.name = g.name + "^" + P;

        res.fill(2);

        for (int k = 0; k < P; ++k) {
            int off = k * GR;
            res.block(off, off, GR, GR) << g.array() + off;
        }

        return res;
    }
}
