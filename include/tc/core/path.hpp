#pragma once

#include <vector>
#include <functional>

#include "action.hpp"

struct Path {
    std::vector<Action> path;

    Path() = default;

    Path(const Path &) = default;

    void add_row() {
        path.resize(path.size() + 1);
    }

    [[nodiscard]] Action get(int to_idx) const {
        return path[to_idx];
    }

    void put(int from_idx, int gen, int to_idx) {
        path[to_idx] = Action(from_idx, gen);
    }

    [[nodiscard]] size_t size() const {
        return path.size();
    }

    template<class T, class F>
    std::vector<T> walk(const T &start, const F &op) {
        std::vector<T> res;
        res.reserve(path.size());
        res.push_back(start);

        for (size_t i = 1; i < path.size(); ++i) {
            auto &action = path[i];
            auto &from = res[action.from_idx];
            auto &gen = action.gen;

            res.push_back(op(from, gen));
        }

        return res;
    }

    template<class T, class E, class F>
    std::vector<T> walk(const T &start, const E &gens, const F &op) {
        return walk(start, [&](const T &from, const int gen) {
            return op(from, gens[gen]);
        });
    }
};
