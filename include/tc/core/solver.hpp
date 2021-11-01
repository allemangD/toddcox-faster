#pragma once

#include <memory>
#include <vector>

namespace tc {
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

        [[nodiscard]] int get(int coset, int gen) const {
            return data[coset * ngens + gen];
        }

        [[nodiscard]] size_t size() const {
            return path.size();
        }
    };
}

namespace {
    struct Row {
        int gnr;
        int *lst;
    };

    struct Table {
    private:
    public:
        Rel rel;
        std::vector<Row> rows;

    public:
        explicit Table(const Rel &rel) : rel(rel) {
        }
    };

    template<class T, size_t BlockSize = 4096>
    class BlockAllocator {
        /// 4096 seems to be the best (on my machine anway) from profiling.
    private:
        int block = 0;
        int next = 0;
        std::vector<T *> data = {build()};

        T *build() {
            T *blk = new T[BlockSize];
            std::fill_n(blk, BlockSize, 0);
            return blk;
        }

    public:
        T *operator()() {
            if (next >= BlockSize) {
                data.push_back(build());
                block++;
                next = 0;
            }

            return &data[block][next++];
        }

        ~BlockAllocator() {
            for (auto &blk: data) {
                delete[] blk;
            }
        }
    };

    template<unsigned int Rank>
    class Tables {
    public:
        static constexpr unsigned int Rels = Rank * (Rank + 1) / 2 - Rank;

    private:
        int *null_lst_ptr = new int;
        BlockAllocator<int> alloc;

        std::array<std::shared_ptr<Table>, Rels> tables;
        std::array<std::vector<std::shared_ptr<Table>>, Rank> deps;

    public:
        explicit Tables(const tc::Group<Rank> &group) {
            const auto &rels = group.get_rels();
            for (int i = 0; i < Rels; ++i) {
                const auto &rel = rels[i];
                auto table = std::make_shared<Table>(rel);
                tables[i] = table;
                deps[rel.gens[0]].push_back(table);
                deps[rel.gens[1]].push_back(table);
            }
        }

        void add_row() {
            // std::vector already does block allocation.
            for (const auto &table: tables) {
                table->rows.emplace_back();
            }
        }

        void initialize(int target, const tc::Cosets &cosets) {
            for (auto &table: tables) {
                const Rel &rel = table->rel;
                Row &row = table->rows[target];

                if (row.lst == nullptr) {
                    if (cosets.get(target, rel.gens[0]) != target and
                        cosets.get(target, rel.gens[1]) != target) {
                        row.lst = alloc();
                        row.gnr = 0;
                    } else {
                        row.lst = null_lst_ptr;
                        row.gnr = -1;
                    }
                }
            }
        }

        ~Tables() {
            delete null_lst_ptr;
        }

        void learn(int coset, int gen, int target, const tc::Cosets &cosets, std::priority_queue<int> &facts) {
            if (target == coset) {
                for (auto &table: deps[gen]) {
                    Row &target_row = table->rows[target];

                    if (target_row.lst == nullptr) {
                        target_row.gnr = -1;
                    }
                }
            }

            for (auto &table: deps[gen]) {
                Row &target_row = table->rows[target];
                Row &coset_row = table->rows[coset];
                const Rel &rel = table->rel;

                if (target_row.lst == nullptr) {
                    target_row.lst = coset_row.lst;
                    target_row.gnr = coset_row.gnr + 1;

                    if (coset_row.gnr < 0) {
                        target_row.gnr -= 2;
                    }

                    if (target_row.gnr == rel.mult) {
                        // forward learn
                        int lst = *target_row.lst;
                        int gen_ = rel.gens[rel.gens[0] == gen];
                        facts.push(lst * Rank + gen_);
                    } else if (target_row.gnr == -rel.mult) {
                        // stationary learn
                        int gen_ = rel.gens[rel.gens[0] == gen];
                        facts.push(target * Rank + gen_);
                    } else if (target_row.gnr == rel.mult - 1) {
                        // determined family
                        *target_row.lst = target;
                    }
                }
            }
        }
    };
}

namespace tc {
    template<unsigned int Rank>
    tc::Cosets solve(const Group <Rank> &g, const std::vector<int> &sub_gens = {}) {
        tc::Cosets cosets(Rank);
        cosets.add_row();

        if (Rank == 0) {
            return cosets;
        }

        for (int gen: sub_gens) {
            if (gen < Rank)
                cosets.put(0, gen, 0);
        }

        Tables<Rank> tables(g);
        tables.add_row();
        tables.initialize(0, cosets);

        std::priority_queue<int> facts;

        for (int idx = 0; idx < cosets.data.size(); idx++) {
            int coset = idx / Rank;
            int gen = idx % Rank;
            if (cosets.get(coset, gen) >= 0) continue;

            int target = cosets.size();
            cosets.add_row();
            tables.add_row();

            facts.push(idx);

            // todo nothing before the current coset will be used.
            //  delete all table rows using old cosets to free memory early.
            //  probably some unrolled linked list would be good; just drop
            //  old blocks.

            while (!facts.empty()) {
                int fact_idx = facts.top();
                facts.pop();

                int coset_ = fact_idx / Rank;
                int gen_ = fact_idx % Rank;

                if (cosets.get(coset_, gen_) != -1)
                    continue;

                cosets.put(coset_, gen_, target);

                tables.learn(coset_, gen_, target, cosets, facts);
            }

            tables.initialize(target, cosets);
        }

        return cosets;
    }
}
