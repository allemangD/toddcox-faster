#pragma once

#include "group.hpp"

namespace tc::group {
    /**
     * Universal Coxeter Group
     */
    Group U(size_t rank) {
        std::string name = "U(" + std::to_string(rank) + ")";

        Group res(rank);
        res.name = name;
        res.fill(2);
        return res;
    }

    /**
     * Simplex
     */
    Group A(size_t rank) {
        std::string name = "A(" + std::to_string(rank) + ")";

        if (rank == 0) {
            Group res(rank);
            res.name = name;
            return res;
        }

        tc::Symbol symbol(rank - 1);
        symbol.fill(3);

        return schlafli(symbol, name);
    }

    /**
     * Cube, Orthoplex
     */
    Group B(size_t rank) {
        std::string name = "B(" + std::to_string(rank) + ")";

        tc::Symbol symbol(rank - 1);
        symbol.fill(3);
        symbol(0) = 4;

        return schlafli(symbol, name);
    }

    /**
     * Demicube, Orthoplex
     */
    Group D(size_t rank) {
        std::string name = "D(" + std::to_string(rank) + ")";

        tc::Symbol symbol(rank - 1);
        symbol.fill(3);
        symbol((Eigen::Index) rank - 2) = 2;

        Group g = schlafli(symbol, name);
        g(1, (Eigen::Index) rank - 1) = 3;
        g((Eigen::Index) rank - 1, 1) = 3;

        return g;
    }

    /**
     * E groups
     */
    Group E(size_t rank) {
        std::string name = "E(" + std::to_string(rank) + ")";

        tc::Symbol symbol(rank - 1);
        symbol.fill(3);
        symbol((Eigen::Index) rank - 2) = 2;

        Group g = schlafli(symbol, name);
        g(2, (Eigen::Index) rank - 1) = 3;
        g((Eigen::Index) rank - 1, 2) = 3;

        return g;
    }

    /**
     * 24 Cell
     */
    Group F4() {
        tc::Symbol symbol(3);
        symbol << 3, 4, 3;
        return schlafli(symbol, "F4");
    }

    /**
     * Hexagon
     */
    Group G2() {
        tc::Symbol symbol(1);
        symbol << 6;
        return schlafli(symbol, "G2");
    }

    /**
     * Icosahedron
     */
    Group H(size_t rank) {
        std::string name = "H(" + std::to_string(rank) + ")";

        tc::Symbol symbol(rank - 1);
        symbol.fill(3);
        symbol(0) = 5;

        return schlafli(symbol, name);
    }

    /**
     * Polygonal
     */
    Group I2(unsigned int n) {
        std::string name = "I2(" + std::to_string(n) + ")";

        tc::Symbol symbol(1);
        symbol << n;
        return schlafli(symbol, name);
    }

    /**
     * Toroidal. I2(n) * I2(m)
     */
    Group T(unsigned int n, unsigned int m) {
        std::string name = "T(" + std::to_string(n) + "," + std::to_string(m) + ")";

        tc::Symbol symbol(3);
        symbol << n, 2, m;

        return schlafli(symbol, name);
    }

    /**
     * Toroidal. T(n, n)
     */
    Group T(unsigned int n) {
        std::string name = "T(" + std::to_string(n) + ")";

        tc::Symbol symbol(3);
        symbol << n, 2, n;

        return schlafli(symbol, name);
    }
}
