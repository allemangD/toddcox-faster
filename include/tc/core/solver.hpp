#pragma once

#include <memory>
#include <vector>

namespace {
    struct Row {
        int gnr;
        std::shared_ptr<int> lst;
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

    struct Tables {
    private:
        std::shared_ptr<int> null_lst_ptr = std::make_shared<int>();
        int ngens;

        std::vector<std::shared_ptr<Table>> tables;
        std::vector<std::vector<std::shared_ptr<Table>>> deps;

    public:
        explicit Tables(const tc::Group &group)
            : ngens(group.ngens), deps(ngens) {
            for (const auto &rel: group.get_rels()) {
                auto table = std::make_shared<Table>(rel);
                tables.push_back(table);
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

        void initialize(int target, const Cosets &cosets) {
            for (auto &table: tables) {
                const Rel &rel = table->rel;
                Row &row = table->rows[target];

                if (row.lst == nullptr) {
                    if (cosets.get(target, rel.gens[0]) != target and
                        cosets.get(target, rel.gens[1]) != target) {
                        row.lst = std::make_shared<int>();
                        row.gnr = 0;
                    } else {
                        row.lst = null_lst_ptr;
                        row.gnr = -1;
                    }
                }
            }
        }

        void learn(int coset, int gen, int target, const Cosets &cosets, std::priority_queue<int> &facts) {
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
                        facts.push(lst * ngens + gen_);
                    } else if (target_row.gnr == -rel.mult) {
                        // stationary learn
                        int gen_ = rel.gens[rel.gens[0] == gen];
                        facts.push(target * ngens + gen_);
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
    Cosets solve(const Group &g, const std::vector<int> &sub_gens = {}) {
        Cosets cosets(g.ngens);
        cosets.add_row();

        if (g.ngens == 0) {
            return cosets;
        }

        for (int gen: sub_gens) {
            if (gen < g.ngens)
                cosets.put(0, gen, 0);
        }

        Tables tables(g);
        tables.add_row();
        tables.initialize(0, cosets);

        std::priority_queue<int> facts;

        for (int idx = 0; idx < cosets.data.size(); idx++) {
            int coset = idx / g.ngens;
            int gen = idx % g.ngens;
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

                int coset_ = fact_idx / g.ngens;
                int gen_ = fact_idx % g.ngens;

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
