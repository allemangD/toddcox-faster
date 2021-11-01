#pragma once

#include "core.hpp"

namespace tc {
    /**
     * Construct a group from a (simplified) Schlafli Symbol of the form [a, b, ..., c]
     * @param mults: The sequence of multiplicites between adjacent generators.
     */
    Group schlafli(const std::vector<int> &mults, const std::string &name) {
        int ngens = (int) mults.size() + 1;

        Group g(ngens, name);

        for (int i = 0; i < (int) mults.size(); i++) {
            g(i, i + 1) = mults[i];
        }

        return g;
    }

    /**
     * Construct a group from a (simplified) Schlafli Symbol of the form [a, b, ..., c]
     * @param mults: The sequence of multiplicites between adjacent generators.
     */
    Group schlafli(const std::vector<int> &mults) {
        std::stringstream ss;
        ss << "[";
        if (!mults.empty()) {
            for (size_t i = 0; i < mults.size() - 1; ++i) {
                ss << mults[i] << ",";
            }
            ss << mults.back();
        }
        ss << "]";

        return schlafli(mults, ss.str());
    }

    namespace group {
        /**
         * Simplex
         */
        Group A(int dim) {
            std::stringstream ss;
            ss << "A(" << dim << ")";

            if (dim == 0)
                return Group(0, ss.str());

            const std::vector<int> &mults = std::vector<int>(dim - 1, 3);

            return schlafli(mults, ss.str());
        }

        /**
         * Cube, Orthoplex
         */
        Group B(int dim) {
            std::stringstream ss;
            ss << "B(" << dim << ")";

            std::vector<int> mults(dim - 1, 3);
            mults[0] = 4;

            return schlafli(mults, ss.str());
        }

        /**
         * Demicube, Orthoplex
         */
        Group D(int dim) {
            std::stringstream ss;
            ss << "D(" << dim << ")";

            std::vector<int> mults(dim - 1, 3);
            mults[dim - 2] = 2;

            Group g = schlafli(mults, ss.str());
            g(1, dim - 1) = 3;

            return g;
        }

        /**
         * E groups
         */
        Group E(int dim) {
            std::stringstream ss;
            ss << "E(" << dim << ")";

            std::vector<int> mults(dim - 1, 3);
            mults[dim - 2] = 2;

            Group g = schlafli(mults, ss.str());
            g(2, dim - 1) = 3;

            return g;
        }

        /**
         * 24 Cell
         */
        Group F4() {
            return schlafli({3, 4, 3}, "F4");
        }

        /**
         * Hexagon
         */
        Group G2() {
            return schlafli({6}, "G2");
        }

        /**
         * Icosahedron
         */
        Group H(int dim) {
            std::stringstream ss;
            ss << "H(" << dim << ")";

            std::vector<int> mults(dim - 1, 3);
            mults[0] = 5;

            return schlafli(mults, ss.str());
        }

        /**
         * Polygonal
         */
        Group I2(int n) {
            std::stringstream ss;
            ss << "I2(" << n << ")";

            return schlafli({n}, ss.str());
        }

        /**
         * Toroidal. I2(n) * I2(m)
         */
        Group T(int n, int m) {
            std::stringstream ss;
            ss << "T(" << n << "," << m << ")";

            return schlafli({n, 2, m}, ss.str());
        }

        /**
         * Toroidal. T(n, n)
         */
        Group T(int n) {
            std::stringstream ss;
            ss << "T(" << n << ")";

            return schlafli({n, 2, n}, ss.str());
        }
    }
}
