#pragma once

#include "core.hpp"

namespace tc {
    /**
     * Construct a group from a (simplified) Schlafli Symbol of the form [a, b, ..., c]
     * @param mults: The sequence of multiplicites between adjacent generators.
     */
    template<unsigned int Rank>
    Group <Rank> schlafli(const std::array<unsigned int, Rank - 1> &mults, const std::string &name) {
        Group<Rank> g(name);

        for (int i = 0; i < Rank - 1; i++) {
            g(i, i + 1) = mults[i];
        }

        return g;
    }

    /**
     * Construct a group from a (simplified) Schlafli Symbol of the form [a, b, ..., c]
     * @param mults: The sequence of multiplicites between adjacent generators.
     */
    template<unsigned int Rank>
    Group <Rank> schlafli(const std::array<unsigned int, Rank - 1> &mults) {
        std::stringstream ss;
        ss << "[";
        if (Rank) {
            for (size_t i = 0; i < Rank - 2; ++i) {
                ss << mults[i] << ",";
            }
            ss << mults[Rank - 1];
        }
        ss << "]";

        return schlafli<Rank>(mults, ss.str());
    }

    namespace group {
        /**
         * Simplex
         */
        template<unsigned int Rank>
        Group <Rank> A() {
            std::stringstream ss;
            ss << "A(" << Rank << ")";

            if (Rank == 0)
                return Group<Rank>(ss.str());

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

            std::array<unsigned int, Rank - 1> mults;
            mults.fill(3);
            mults[0] = 4;

            return schlafli<Rank>(mults, ss.str());
        }

        /**
         * Demicube, Orthoplex
         */
        template<unsigned int Rank>
        Group <Rank> D() {
            std::stringstream ss;
            ss << "D(" << Rank << ")";

            std::array<unsigned int, Rank - 1> mults;
            mults.fill(3);
            mults[Rank - 2] = 2;

            Group<Rank> g = schlafli<Rank>(mults, ss.str());
            g(1, Rank - 1) = 3;

            return g;
        }

        /**
         * E groups
         */
        template<unsigned int Rank>
        Group <Rank> E() {
            std::stringstream ss;
            ss << "E(" << Rank << ")";

            std::array<unsigned int, Rank - 1> mults;
            mults.fill(3);
            mults[Rank - 2] = 2;

            Group<Rank> g = schlafli<Rank>(mults, ss.str());
            g(2, Rank - 1) = 3;

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
            return schlafli<2>({6}, "G2");
        }

        /**
         * Icosahedron
         */
        template<unsigned int Rank>
        Group <Rank> H() {
            std::stringstream ss;
            ss << "H(" << Rank << ")";

            std::array<unsigned int, Rank - 1> mults;
            mults.fill(3);
            mults[0] = 5;

            return schlafli<Rank>(mults, ss.str());
        }

        /**
         * Polygonal
         */
        Group<2> I2(unsigned int n) {
            std::stringstream ss;
            ss << "I2(" << n << ")";

            return schlafli<2>({n}, ss.str());
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
}
