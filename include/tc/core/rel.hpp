#pragma once

#include <array>

struct Rel {
    std::array<int, 2> gens;
    int mult;

    Rel() = default;

    Rel(const Rel &) = default;

    Rel(int a, int b, int m)
        : gens({a, b}), mult(m) {
    }

    [[nodiscard]] Rel shift(int off) const {
        return Rel(gens[0] + off, gens[1] + off, mult);
    }
};
