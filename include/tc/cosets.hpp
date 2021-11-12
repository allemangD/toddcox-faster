#pragma once

#include <vector>

namespace tc {
    class Path;

    class Cosets {
    private:
        std::vector<int> data;
        size_t _rank;

    public:
        Cosets(const Cosets &) = default;

        explicit Cosets(size_t rank) : _rank(rank) {}

        void add_row() {
            data.resize(data.size() + rank(), -1);
        }

        void put(int coset, int gen, int target) {
            data[coset * rank() + gen] = target;
            data[target * rank() + gen] = coset;
        }

        [[nodiscard]] int get(int coset, int gen) const {
            return data[coset * rank() + gen];
        }

        [[nodiscard]] size_t rank() const {
            return _rank;
        }

        [[nodiscard]] size_t order() const {
            if (!_rank) return 0;
            return data.size() / _rank;
        }

        Path path() const;
    };

    class Path {
    private:
        friend class Cosets;

        std::vector<unsigned int> source;
        std::vector<unsigned int> gen;
        size_t _rank;
        size_t _order;

        explicit Path(size_t rank, size_t order) : _rank(rank), _order(order), source(order), gen(order) {}

    public:
        [[nodiscard]] size_t rank() const {
            return _rank;
        }

        [[nodiscard]] size_t order() const {
            return _order;
        }

        template<class T, class F>
        std::vector<T> walk(const T &start, const F &op) const {
            std::vector<T> res;
            res.reserve(order());
            res.push_back(start);

            for (size_t i = 1; i < order(); ++i) {
                auto val = op(res[source[i]], gen[i]);
                res.push_back(val);
            }

            return res;
        }

        template<class T, class E, class F>
        std::vector<T> walk(const T &start, const E &gens, const F &op) const {
            return walk(start, [&](const T &s, const int g) {
                return op(s, gens[g]);
            });
        }
    };

    Path Cosets::path() const {
        Path res(rank(), order());
        std::vector<bool> set(order());

        for (int coset = 0; coset < order(); ++coset) {
            for (int gen = 0; gen < rank(); ++gen) {
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
