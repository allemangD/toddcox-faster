#pragma once

#include <vector>

namespace tc {
    template<unsigned int Rank>
    class Path;

    template<unsigned int Rank>
    class Cosets {
    private:
        std::vector<int> data;

    public:
        Cosets(const Cosets<Rank> &) = default;

        Cosets() = default;

        void add_row() {
            data.resize(data.size() + Rank, -1);
        }

        void put(int coset, int gen, int target) {
            data[coset * Rank + gen] = target;
            data[target * Rank + gen] = coset;
        }

        [[nodiscard]] int get(int coset, int gen) const {
            return data[coset * Rank + gen];
        }

        [[nodiscard]] size_t order() const {
            return data.size() / Rank;
        }

        Path<Rank> path() const;
    };

    template<unsigned int Rank>
    class Path {
    private:
        friend class Cosets<Rank>;

        std::vector<unsigned int> source;
        std::vector<unsigned int> gen;
        size_t _order;

        explicit Path(size_t order) : _order(order), source(order), gen(order) {}

    public:
        size_t order() const {
            return _order;
        }

        template<class T, class F>
        std::vector <T> walk(const T &start, const F &op) {
            std::vector <T> res;
            res.reserve(order());
            res.push_back(start);

            for (size_t i = 1; i < order(); ++i) {
                auto val = op(res[source[i]], gen[i]);
                res.push_back(val);
            }

            return res;
        }

        template<class T, class E, class F>
        std::vector <T> walk(const T &start, const E &gens, const F &op) {
            return walk(start, [&](const T &s, const int g) {
                return op(s, gens[g]);
            });
        }
    };

    template<unsigned int Rank>
    Path<Rank> Cosets<Rank>::path() const {
        Path<Rank> res(order());
        std::vector<bool> set(order());

        for (int coset = 0; coset < order(); ++coset) {
            for (int gen = 0; gen < Rank; ++gen) {
                int target = get(coset, gen);

                if (!set[target]) {
                    res.source[target] = coset;
                    res.gen[target] = gen;
                    set[target] = true;
                }
            }
        }

        return res;
    }
}
