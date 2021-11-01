#pragma once

#include "core.hpp"

namespace tc::group {
    /**
     * Universal Coxeter Group
     */
    template<unsigned int Rank>
    Group <Rank> U() {
        std::stringstream ss;
        ss << "U(" << Rank << ")";

        Group<Rank> res;
        res.name = ss.str();
        res.fill(2);
        return res;
    }

    /**
     * Simplex
     */
    template<unsigned int Rank>
    Group <Rank> A() {
        std::stringstream ss;
        ss << "A(" << Rank << ")";

        if (Rank == 0) {
            Group<Rank> res;
            res.name = ss.str();
            return res;
        }

        std::array<unsigned int, Rank - 1> mults;
        mults.fill(3);

        return schlafli<Rank>(mults, ss.str());
    }

    /**
     * Cube, Orthoplex
     */
    template<unsigned int Rank>
    Group <Rank> B() {
        std::stringstream ss;
        ss << "B(" << Rank << ")";

        tc::Symbol<Rank - 1> mults;
        mults.fill(3);
        mults(0) = 4;

        return schlafli<Rank>(mults, ss.str());
    }

    /**
     * Demicube, Orthoplex
     */
    template<unsigned int Rank>
    Group <Rank> D() {
        std::stringstream ss;
        ss << "D(" << Rank << ")";

        tc::Symbol<Rank - 1> mults;
        mults.fill(3);
        mults(Rank - 2) = 2;

        Group<Rank> g = schlafli<Rank>(mults, ss.str());
        g(1, Rank - 1) = 3;
        g(Rank - 1, 1) = 3;

        return g;
    }

    /**
     * E groups
     */
    template<unsigned int Rank>
    Group <Rank> E() {
        std::stringstream ss;
        ss << "E(" << Rank << ")";

        tc::Symbol<Rank - 1> mults;
        mults.fill(3);
        mults(Rank - 2) = 2;

        Group<Rank> g = schlafli<Rank>(mults, ss.str());
        g(2, Rank - 1) = 3;
        g(Rank - 1, 2) = 3;

        return g;
    }

    /**
     * 24 Cell
     */
    Group<4> F4() {
        return schlafli<4>({3, 4, 3}, "F4");
    }

    /**
     * Hexagon
     */
    Group<2> G2() {
        return schlafli<2>(tc::Symbol<1>(6), "G2");
    }

    /**
     * Icosahedron
     */
    template<unsigned int Rank>
    Group <Rank> H() {
        std::stringstream ss;
        ss << "H(" << Rank << ")";

        tc::Symbol<Rank - 1> mults;
        mults.fill(3);
        mults(0) = 5;

        return schlafli<Rank>(mults, ss.str());
    }

    /**
     * Polygonal
     */
    Group<2> I2(unsigned int n) {
        std::stringstream ss;
        ss << "I2(" << n << ")";

        return schlafli<2>(tc::Symbol<1>(n), ss.str());
    }

    /**
     * Toroidal. I2(n) * I2(m)
     */
    Group<4> T(unsigned int n, unsigned int m) {
        std::stringstream ss;
        ss << "T(" << n << "," << m << ")";

        return schlafli<4>({n, 2, m}, ss.str());
    }

    /**
     * Toroidal. T(n, n)
     */
    Group<4> T(unsigned int n) {
        std::stringstream ss;
        ss << "T(" << n << ")";

        return schlafli<4>({n, 2, n}, ss.str());
    }
}
