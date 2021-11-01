#pragma once

#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <memory>
#include <queue>

#include "rel.hpp"

#include <Eigen/Eigen>
#include <iostream>

namespace tc {
    template<unsigned int Rank>
    struct Group;

    template<unsigned int Rank, unsigned int PRank>
    struct SubGroup;

    template<unsigned int Rank>
    class Group {
    public:
        using Matrix = Eigen::Matrix<unsigned int, Rank, Rank>;

    private:
    public:
        std::string _name;
        Matrix _data;
        Eigen::SelfAdjointView<Matrix, Eigen::Upper> _mults;

    public:
        Group(const Group<Rank> &g) :
            _name(g._name),
            _data(g._data),
            _mults(_data) {
        }

        Group(Group &&g) noexcept:
            _name(std::move(g._name)),
            _data(std::move(g._data)),
            _mults(_data) {
        }

        explicit Group(std::string name = "G") :
            _name(std::move(name)),
            _data(),
            _mults(_data) {
            _data.fill(2);
        }

        unsigned int rank() const {
            return _data.rows();
        }

        std::string name() const {
            return _name;
        }

        typename Matrix::Scalar &operator()(int a, int b) {
            return _mults(a, b);
        }

        typename Matrix::Scalar operator()(int a, int b) const {
            return _mults(a, b);
        }

        [[nodiscard]] std::vector<Rel> get_rels() const {
            std::vector<Rel> res;
            for (int i = 0; i < Rank - 1; ++i) {
                for (int j = i + 1; j < Rank; ++j) {
                    res.emplace_back(i, j, _mults(i, j));
                }
            }
            return res;
        }
    };

    template<unsigned int GR, unsigned int HR>
    Group<GR + HR> product(const Group<GR> &g, const Group<HR> &h) {
        std::stringstream ss;
        ss << g.name << "*" << h.name;

        Group<GR + HR> res(ss.str());

        int off = 0;
        for (int i = 0; i < GR; ++i) {
            for (int j = i; j < GR; ++j) {
                res(i + off, j + off) = g(i, j);
            }
        }
        off += GR;

        for (int i = 0; i < HR; ++i) {
            for (int j = i; j < HR; ++j) {
                res(i + off, j + off) = h(i, j);
            }
        }
        off += HR;

        return res;
    }

    template<unsigned int GR, unsigned int P>
    Group<GR * P> power(const Group<GR> &g) {
        std::stringstream ss;
        ss << g.name << "^" << P;

        Group<GR * P> res(ss.str());

        for (int k = 0; k < P; ++k) {
            int off = k * GR;

            for (int i = 0; i < GR; ++i) {
                for (int j = i; j < GR; ++j) {
                    res(i + off, j + off) = g(i, j);
                }
            }
        }

        return res;
    }

//    Group operator*(const Group &g, const Group &h) {
//        return product(g, h);
//    }
//
//    Group operator^(const Group &g, int p) {
//        return power(g, p);
//    }
}
