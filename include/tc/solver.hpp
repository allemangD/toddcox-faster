#pragma once

#include <algorithm>
#include <array>
#include <memory>
#include <vector>
#include <queue>

#include "group.hpp"
#include "cosets.hpp"

namespace {
    struct Row {
        int gnr;
        int *lst;
    };

    struct Table {
    private:
    public:
        int i, j, mult;

        std::vector<Row> rows;

    public:
        explicit Table(int i, int j, int mult) :
            i(i), j(j), mult(mult) {
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

    class Tables {
    private:
        int *null_lst_ptr = new int;
        BlockAllocator<int> alloc;

        std::vector<std::shared_ptr<Table>> tables;
        std::vector<std::vector<std::shared_ptr<Table>>> deps;

        size_t _rank;
        size_t _rels;

    public:
        explicit Tables(const tc::Group &group) : _rank(group.rank()), _rels(rank() * (rank() + 1) / 2 - rank()) {
            deps.resize(rank());
            for (int i = 0; i < rank() - 1; ++i) {
                for (int j = i + 1; j < rank(); ++j) {
                    auto table = std::make_shared<Table>(i, j, group(i, j));
                    tables.push_back(table);
                    deps[i].push_back(table);
                    deps[j].push_back(table);
                }
            }
        }

        [[nodiscard]] size_t rank() const {
            return _rank;
        }

        [[nodiscard]] size_t rels() const {
            return _rels;
        }

        void add_row() {
            // std::vector already does block allocation.
            for (const auto &table: tables) {
                table->rows.emplace_back();
            }
        }

        void initialize(int target, const tc::Cosets &cosets) {
            for (auto &table: tables) {
                Row &row = table->rows[target];

                if (row.lst == nullptr) {
                    if (cosets.get(target, table->i) != target and
                        cosets.get(target, table->j) != target) {
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

        void learn(int coset, int gen, int target, const tc::Cosets &cosets, std::priority_queue<size_t> &facts) {
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

                if (target_row.lst == nullptr) {
                    target_row.lst = coset_row.lst;
                    target_row.gnr = coset_row.gnr + 1;

                    if (coset_row.gnr < 0) {
                        target_row.gnr -= 2;
                    }

                    if (target_row.gnr == table->mult) {
                        // forward learn
                        int lst = *target_row.lst;
                        int gen_ = (table->i == gen) ? table->j : table->i;
                        facts.push(lst * rank() + gen_);
                    } else if (target_row.gnr == -table->mult) {
                        // stationary learn
                        int gen_ = (table->i == gen) ? table->j : table->i;
                        facts.push(target * rank() + gen_);
                    } else if (target_row.gnr == table->mult - 1) {
                        // determined family
                        *target_row.lst = target;
                    }
                }
            }
        }
    };
}

namespace tc {
    /**
     * Assumes that g is a coxeter group - that is, self-adjoint and the diagonal is 2.
     */
    tc::Cosets solve(const Group &group, const Symbol &s_gens) {
        size_t rank = group.rank();

        tc::Cosets cosets(rank);
        cosets.add_row();

        if (rank == 0) {
            return cosets;
        }

        for (unsigned int gen: s_gens) {
            if (gen < rank)
                cosets.put(0, gen, 0);
        }

        Tables tables(group);
        tables.add_row();
        tables.initialize(0, cosets);

        std::priority_queue<size_t> facts;

        for (int coset = 0; coset < cosets.order(); coset++) {
            for (int gen = 0; gen < rank; ++gen) {
                if (cosets.get(coset, gen) >= 0) continue; // todo vector<bool> set

                int target = cosets.order();
                cosets.add_row();
                tables.add_row();

                facts.push(coset * rank + gen);

                // todo nothing before the current coset will be used.
                //  delete all table rows using old cosets to free memory early.
                //  probably some unrolled linked list would be good; just drop
                //  old blocks.

                while (!facts.empty()) {
                    int fact_idx = facts.top();
                    facts.pop();

                    int coset_ = fact_idx / rank;
                    int gen_ = fact_idx % rank;

                    if (cosets.get(coset_, gen_) != -1)
                        continue;

                    cosets.put(coset_, gen_, target);

                    tables.learn(coset_, gen_, target, cosets, facts);
                }

                tables.initialize(target, cosets);
            }
        }

        return cosets;
    }

/**
 * Solve the cosets generated by sg_gens within the subgroup generated by g_gens of the group context
 */
    Cosets solve(
        const Group &context,
        const Symbol &g_gens,
        const Symbol &sg_gens
    ) {
        const Symbol &proper_sg_gens = recontext_gens(context.rank(), g_gens, sg_gens);
        const Group &group = subgroup(context, g_gens);

        return solve(group, proper_sg_gens);
    }
}
