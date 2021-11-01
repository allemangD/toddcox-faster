#pragma once

#include <string>
#include <sstream>

#include <Eigen/Eigen>

namespace {
    template<class T>
    std::string stringify(const T &vec) {
        std::stringstream ss;
        ss << "[" << vec << "]";
        return ss.str();
    }
}

namespace tc {
    /// A Schlafli Matrix
    template<unsigned int Rank>
    class Group : public Eigen::Matrix<unsigned int, Rank, Rank> {
    public:
        using Base = Eigen::Matrix<unsigned int, Rank, Rank>;

        std::string name;

        using Base::Base;
    };

    template<unsigned int Rank>
    using Symbol = Eigen::Vector<unsigned int, Rank>;

    template<unsigned int Rank>
    Group<Rank> schlafli(const Symbol<Rank - 1> &mults, const std::string &name) {
        Group<Rank> res;
        res.name = name;

        res.fill(2);
        res.topRightCorner(Rank - 1, Rank - 1).diagonal() << mults;
        res.bottomLeftCorner(Rank - 1, Rank - 1).diagonal() << mults;

        return res;
    }

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
