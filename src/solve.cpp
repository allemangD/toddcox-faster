#include "tc/core.hpp"

#include <algorithm>
#include <memory>

namespace tc {
    struct Row {
        int gnr;
        std::shared_ptr<int> lst;
    };

    struct RelTables {
    private:
        std::vector<Row> rows;
        size_t nrels;

    public:
        explicit RelTables(size_t nrels) : nrels(nrels) {
        }

        Row &operator()(size_t irel, size_t coset) {
            // todo check bounds
            size_t idx = coset * nrels + irel;
            return rows[idx];
        }

        void add_row() {
            // std::vector already does block allocation.
            rows.resize(rows.size() + nrels);
        }

        void del_rows_to(size_t coset) {
            /// strictly refers to freeing pre-allocated blocks of for *gnrs and **lst_ptrs.
            /// actual `rows` is unchanged.
        }
    };

    Cosets Group::solve(const std::vector<int> &sub_gens) const {
        Cosets cosets(ngens);
        cosets.add_row();

        if (ngens == 0) {
            return cosets;
        }

        for (int g: sub_gens) {
            if (g < ngens)
                cosets.put(0, g, 0);
        }


        const auto &rels = get_rels(); // todo move to Group member
        const auto nrels = rels.size();

        // todo encapsulate
        std::vector<std::vector<int>> gen_map(ngens);
        int rel_idx = 0;
        for (Rel m: rels) {
            gen_map[m.gens[0]].push_back(rel_idx);
            gen_map[m.gens[1]].push_back(rel_idx);
            rel_idx++;
        }

        std::shared_ptr<int> null_lst_ptr = std::make_shared<int>();

        RelTables tables(nrels);
        tables.add_row();

        for (int irel = 0; irel < nrels; irel++) {
            Row &row = tables(irel, 0);
            const Rel &rel = rels[irel];

            if (cosets.get(rel.gens[0]) == -1 && cosets.get(rel.gens[1]) == -1) {
                row.lst = std::make_shared<int>();
                row.gnr = 0;
            } else {
                row.lst = null_lst_ptr;
                row.gnr = -1;
            }
        }

        int idx = 0;
        int coset, gen, target, fact_idx, lst, gen_;
        while (true) {
            while (idx < cosets.data.size() and cosets.get(idx) >= 0)
                idx++;

            if (idx == cosets.data.size()) {
                tables.del_rows_to(idx / ngens);
                break;
            }

            target = cosets.size();
            cosets.add_row();
            tables.add_row();

            std::vector<int> facts;
            facts.push_back(idx);

            coset = idx / ngens;
            gen = idx % ngens;

            tables.del_rows_to(coset);

            while (!facts.empty()) {
                fact_idx = facts.back();
                facts.pop_back();

                if (cosets.get(fact_idx) != -1)
                    continue;

                cosets.put(fact_idx, target);

                coset = fact_idx / ngens;
                gen = fact_idx % ngens;

                if (target == coset) {
                    for (int irel: gen_map[gen]) {
                        Row &target_row = tables(irel, target);
                        if (target_row.lst == nullptr)
                            target_row.gnr = -1;
                    }
                }

                for (int irel: gen_map[gen]) {
                    Row &target_row = tables(irel, target);
                    Row &coset_row = tables(irel, coset);

                    if (target_row.lst == nullptr) {
                        const Rel &rel = rels[irel];
                        target_row.lst = coset_row.lst;
                        target_row.gnr = coset_row.gnr + 1;

                        if (coset_row.gnr < 0) {
                            target_row.gnr -= 2;
                        }

                        if (target_row.gnr == rel.mult) {
                            lst = *target_row.lst;
                            gen_ = rel.gens[rel.gens[0] == gen];
                            facts.push_back(lst * ngens + gen_);
                        } else if (target_row.gnr == -rel.mult) {
                            gen_ = rel.gens[rel.gens[0] == gen];
                            facts.push_back(target * ngens + gen_);
                        } else if (target_row.gnr == rel.mult - 1) {
                            *target_row.lst = target;
                        }
                    }
                }

                std::sort(facts.begin(), facts.end(), std::greater<>());
            }

            for (int irel = 0; irel < nrels; irel++) {
                Row &target_row = tables(irel, target);
                const Rel &rel = rels[irel];

                if (target_row.lst == nullptr) {
                    if ((cosets.get(target, rel.gens[0]) != target) and
                        (cosets.get(target, rel.gens[1]) != target)) {
                        target_row.lst = std::make_shared<int>();
                        target_row.gnr = 0;
                    } else {
                        target_row.lst = null_lst_ptr;
                        target_row.gnr = -1;
                    }
                }
            }
        }

        return cosets;
    }
}
