#pragma once

#include <vector>

#include "path.hpp"

struct Cosets {
    int ngens;
    std::vector<int> data;
    Path path;

    Cosets(const Cosets &) = default;

    explicit Cosets(int ngens)
        : ngens(ngens) {
    }

    void add_row() {
        data.resize(data.size() + ngens, -1);
        path.add_row();
    }

    void put(int coset, int gen, int target) {
        data[coset * ngens + gen] = target;
        data[target * ngens + gen] = coset;

        if (path.get(target).from_idx == -1) {
            path.put(coset, gen, target);
        }
    }

    void put(int idx, int target) {
        // todo remove
        int coset = idx / ngens;
        int gen = idx % ngens;
        put(coset, gen, target);
    }

    [[nodiscard]] int get(int coset, int gen) const {
        return data[coset * ngens + gen];
    }

    [[nodiscard]] int get(int idx) const {
        return data[idx];
    }

    [[nodiscard]] size_t size() const {
        return path.size();
    }
};