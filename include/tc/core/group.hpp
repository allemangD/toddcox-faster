#pragma once

#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <memory>
#include <queue>

#include "rel.hpp"
#include "cosets.hpp"

#include <Eigen/Eigen>
#include <iostream>

namespace tc {
    struct Group;
    struct SubGroup;

    struct Group {
        using Matrix = Eigen::MatrixXi;

        int ngens;
        std::string name;
        Matrix _data;
        Eigen::SelfAdjointView<Matrix, Eigen::Upper> _mults;

        Group(const Group &g)
            : ngens(g.ngens),
            name(g.name),
            _data(g._data),
            _mults(_data) {
        }

        Group(Group &&g) noexcept
            : ngens(g.ngens),
            name(std::move(g.name)),
            _data(std::move(g._data)),
            _mults(_data) {
        }

        explicit Group(
            int ngens,
            std::string name = "G"
        ) : ngens(ngens),
            name(std::move(name)),
            _data(ngens, ngens),
            _mults(_data) {
            _data.fill(2);
        }

        Matrix::Scalar &operator()(int a, int b) {
            return _mults(a, b);
        }

        Matrix::Scalar operator()(int a, int b) const {
            return _mults(a, b);
        }

        [[nodiscard]] std::vector<Rel> get_rels() const {
            std::vector<Rel> res;
            for (int i = 0; i < ngens - 1; ++i) {
                for (int j = i + 1; j < ngens; ++j) {
                    res.emplace_back(i, j, operator()(i, j));
                }
            }
            return res;
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
                    int mult = parent(gen_map[i], gen_map[j]);
                    operator()(i, j) = mult;
                }
            }
        }
    };

    SubGroup Group::subgroup(const std::vector<int> &gens) const {
        return SubGroup(*this, gens);
    }

    Group product(const Group &g, const Group &h) {
        std::stringstream ss;
        ss << g.name << "*" << h.name;

        Group res(g.ngens + h.ngens, ss.str());

        int off = 0;

        for (int i = 0; i < g.ngens; ++i) {
            for (int j = i; j < g.ngens; ++j) {
                res(i + off, j + off) = g(i, j);
            }
        }
        off += g.ngens;

        for (int i = 0; i < h.ngens; ++i) {
            for (int j = i; j < h.ngens; ++j) {
                res(i + off, j + off) = h(i, j);
            }
        }

        return res;
    }

    Group power(const Group &g, int p) {
        std::stringstream ss;
        ss << g.name << "^" << p;

        Group res(g.ngens * p, ss.str());

        for (int i = 0; i < g.ngens; ++i) {
            for (int j = i; j < g.ngens; ++j) {
                for (int k = 0; k < p; ++k) {
                    int off = k * g.ngens;
                    res(i + off, j + off) = g(i, j);
                }
            }
        }

        return res;
    }

    Group operator*(const Group &g, const Group &h) {
        return product(g, h);
    }

    Group operator^(const Group &g, int p) {
        return power(g, p);
    }
}

#include "solve.hpp"
