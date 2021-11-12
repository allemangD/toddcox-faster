#pragma once

#include <string>
#include <sstream>
#include <numeric>

#include <Eigen/Eigen>

namespace tc {
    template<class T>
    std::string stringify(const T &vec) {
        std::stringstream ss;
        ss << "[" << vec.transpose() << "]";
        return ss.str();
    }
}

namespace tc {
    using Symbol = Eigen::Vector<unsigned int, Eigen::Dynamic>;

    using MatrixXui = Eigen::Matrix<unsigned int, Eigen::Dynamic, Eigen::Dynamic>;

    /// A Coxeter Matrix
    class Group : public MatrixXui {
    public:
        using Base = MatrixXui;

        std::string name = "G";
        Symbol gens;

        explicit Group(size_t rank) : Base(rank, rank), gens(rank) {
            for (Eigen::Index i = 0; i < rank; ++i) {
                gens(i) = i;
            }
        }

        [[nodiscard]] size_t rank() const {
            return rows();
        }
    };

    Group subgroup(const Group &group, const Symbol &gens) {
        size_t rank = group.size();
        size_t srank = gens.size();

        Group res(srank);
        res.name = group.name + ":" + stringify(gens);
        res.gens = gens;

        for (Eigen::Index i = 0; i < srank; ++i) {
            for (Eigen::Index j = 0; j < srank; ++j) {
                res(i, j) = group(gens[i], gens[j]);
            }
        }

        return res;
    }

    unsigned int factorial(unsigned int n) {
        unsigned int res = 1;
        for (int i = 1; i <= n; ++i) {
            res *= i;
        }
        return res;
    }

    unsigned int choose(unsigned int n, unsigned int k) {
        return factorial(n) / factorial(k) / factorial(n - k);
    }

    using SubGroups = std::vector<Group>;

    SubGroups subgroups(const Group &group, size_t srank) {
        size_t rank = group.rank();

        std::vector<bool> mask(rank, false);
        std::fill_n(mask.begin(), srank, true);

        SubGroups res;
        res.reserve(choose(rank, srank));

        Symbol row(srank);
        do {
            for (int j = 0, k = 0; j < rank; ++j) {
                if (mask[j])
                    row(k++) = j;
            }
            res.push_back(subgroup(group, row));
        } while (std::prev_permutation(mask.begin(), mask.end()));

        return res;
    }

    /**
     * Determine which of g_gens are the correct names for sg_gens within the current context
     */
    Symbol recontext_gens(
        size_t rank,
        Symbol g_gens,
        Symbol sg_gens
    ) {
        std::sort(g_gens.begin(),  g_gens.end());
        std::sort(sg_gens.begin(),  sg_gens.end());

        int inv_gen_map[rank];
        for (int i = 0; i < g_gens.size(); ++i) {
            inv_gen_map[g_gens[i]] = i;
        }

        Symbol s_sg_gens(sg_gens.size());
        for (int i = 0; i < sg_gens.size(); ++i) {
            s_sg_gens[i] = inv_gen_map[sg_gens[i]];
        }
        std::sort(s_sg_gens.begin(),  s_sg_gens.end());

        return s_sg_gens;
    }

    /**
     * Create a named coxeter matrix from a simplified schlafli symbol
     */
    Group schlafli(const Symbol &mults, const std::string &name) {
        size_t rank = mults.size() + 1;

        Group res(rank);
        res.name = name;

        res.fill(2);
        res.diagonal().fill(1);
        res.topRightCorner(rank - 1, rank - 1).diagonal() << mults;
        res.bottomLeftCorner(rank - 1, rank - 1).diagonal() << mults;

        return res;
    }

    /**
     * Create a coxeter matrix from a simplified schlafli symbol.
     */
    Group schlafli(const Symbol &mults) {
        return schlafli(mults, stringify(mults));
    }

    Group product(const Group &g, const Group &h) {
        Group res(g.rank() + h.rank());
        res.name = g.name + "*" + h.name;

        res.fill(2);

        Eigen::Index off = 0;
        res.block(off, off, g.rank(), g.rank()) << g.array() + off;
        off += (Eigen::Index) g.rank();

        res.block(off, off, h.rank(), h.rank()) << h.array() + off;
        off += (Eigen::Index) h.rank();

        return res;
    }

    Group power(const Group &g, size_t p) {
        Group res(g.rank() * p);

        res.name = g.name + "^" + std::to_string(p);

        res.fill(2);

        for (Eigen::Index k = 0; k < p; ++k) {
            auto off = (Eigen::Index) g.rank() * k;
            res.block(off, off, g.rank(), g.rank()) << g.array() + off;
        }

        return res;
    }
}
