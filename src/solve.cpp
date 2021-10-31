#include "tc/core.hpp"

#include <algorithm>
#include <memory>
#include <queue>

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
    };

    std::vector<std::vector<int>> dependency_map(int ngens, const std::vector<Rel> &rels) {
        std::vector<std::vector<int>> deps(ngens);
        for (int irel = 0; irel < rels.size(); ++irel) {
            const Rel &rel = rels[irel];
            deps[rel.gens[0]].push_back(irel);
            deps[rel.gens[1]].push_back(irel);
        }
        return deps;
    }

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

        auto deps = dependency_map(ngens, rels);

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

        for (int idx = 0; idx < cosets.data.size(); idx++) {
            if (cosets.get(idx) >= 0) continue;

            int target = cosets.size();
            cosets.add_row();
            tables.add_row();

            std::priority_queue<int> facts;
            facts.push(idx);

            // todo nothing before the current coset will be used.
            //  delete all table rows using old cosets to free memory early.
            //  probably some unrolled linked list would be good; just drop
            //  old blocks.

            while (!facts.empty()) {
                int fact_idx = facts.top();
                facts.pop();

                if (cosets.get(fact_idx) != -1)
                    continue;

                cosets.put(fact_idx, target);

                int coset = fact_idx / ngens;
                int gen = fact_idx % ngens;

                if (target == coset) {
                    for (int irel: deps[gen]) {
                        Row &target_row = tables(irel, target);
                        if (target_row.lst == nullptr) {
                            target_row.gnr = -1;
                        }
                    }
                }

                for (int irel: deps[gen]) {
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
