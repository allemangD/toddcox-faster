#pragma once

#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <memory>
#include <queue>

#include "rel.hpp"
#include "cosets.hpp"

namespace tc {
    struct Group;
    struct SubGroup;

    struct Group {
        int ngens;
        std::vector<std::vector<int>> _mults;
        std::string name;

        Group(const Group &) = default;

        explicit Group(
            int ngens,
            const std::vector<Rel> &rels = {},
            std::string name = "G"
        ) : ngens(ngens), name(std::move(name)) {
            _mults.resize(ngens);

            for (auto &mult: _mults) {
                mult.resize(ngens, 2);
            }

            for (const auto &rel: rels) {
                set(rel);
            }
        }

        void set(const Rel &r) {
            _mults[r.gens[0]][r.gens[1]] = r.mult;
            _mults[r.gens[1]][r.gens[0]] = r.mult;
        }

        [[nodiscard]] int get(int a, int b) const {
            return _mults[a][b];
        }

        [[nodiscard]] std::vector<Rel> get_rels() const {
            std::vector<Rel> res;
            for (int i = 0; i < ngens - 1; ++i) {
                for (int j = i + 1; j < ngens; ++j) {
                    res.emplace_back(i, j, get(i, j));
                }
            }
            return res;
        }

        [[nodiscard]] Group product(const Group &other) const {
            std::stringstream ss;
            ss << name << "*" << other.name;

            Group g(ngens + other.ngens, get_rels(), ss.str());

            for (const auto &rel: other.get_rels()) {
                g.set(rel.shift(ngens));
            }

            return g;
        }

        [[nodiscard]] Group power(int p) const {
            std::stringstream ss;
            ss << name << "^" << p;

            Group g(ngens * p, {}, ss.str());
            for (const auto &rel: get_rels()) {
                for (int off = 0; off < g.ngens; off += ngens) {
                    g.set(rel.shift(off));
                }
            }

            return g;
        }

        [[nodiscard]] SubGroup subgroup(
            const std::vector<int> &gens
        ) const;

        [[nodiscard]] Cosets solve(
            const std::vector<int> &sub_gens = {}
        ) const;
    };

    struct SubGroup : public Group {
        std::vector<int> gen_map;
        const Group &parent;

        SubGroup(const Group &parent, std::vector<int> gen_map) : Group(gen_map.size()), parent(parent) {

            std::sort(gen_map.begin(), gen_map.end());
            this->gen_map = gen_map;

            for (size_t i = 0; i < gen_map.size(); ++i) {
                for (size_t j = 0; j < gen_map.size(); ++j) {
                    int mult = parent.get(gen_map[i], gen_map[j]);
                    set(Rel(i, j, mult));
                }
            }
        }
    };

    SubGroup Group::subgroup(const std::vector<int> &gens) const {
        return SubGroup(*this, gens);
    }

    Group operator*(const Group &g, const Group &h) {
        return g.product(h);
    }

    Group operator^(const Group &g, int p) {
        return g.power(p);
    }
}

#include "solve.hpp"
