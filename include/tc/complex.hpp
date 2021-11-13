#pragma once

#include <tc/group.hpp>
#include <tc/solver.hpp>
#include <cmath>
#include <optional>
#include <numeric>
#include <iostream>
#include <utility>

namespace tc {
    std::vector<Symbol> combinations(const Symbol &symbol, size_t srank) {
        size_t rank = symbol.size();

        std::vector<bool> mask(rank, false);
        std::fill_n(mask.begin(), srank, true);

        std::vector<Symbol> combos;
        combos.reserve(choose(rank, srank));

        Symbol row(srank);
        do {
            for (int j = 0, k = 0; j < rank; ++j) {
                if (mask[j]) {
                    row(k++) = symbol(j);
                }
            }
            combos.emplace_back(row);
        } while (std::prev_permutation(mask.begin(), mask.end()));

        return combos;
    }

    Symbol fan(const Symbol &prim, unsigned root) {
        Symbol res(prim.size() + 1);
        res << prim, root;
        return res;
    }

    std::vector<Symbol> fan(const std::vector<Symbol>& prims, int root) {
        std::vector<Symbol> res;
        res.reserve(prims.size());
        for (const auto &prim: prims) {
            Symbol s(prim.size() + 1);
            s << prim, root;
            res.push_back(s);
        }
        return res;
    }

    void flip(Symbol &prim) {
        if (prim.size() > 1) {
            std::swap(prim[0], prim[1]);
        }
    }

    void apply(const tc::Cosets &table, unsigned int gen, Symbol &prim) {
        for (auto &ind: prim) {
            ind = table.get(ind, gen);
        }
        flip(prim);
    }

    void apply(const tc::Cosets &table, unsigned int gen, std::vector<Symbol> &prims) {
        for (auto &prim: prims) {
            apply(table, gen, prim);
        }
    }

/**
 * Produce a list of all generators for the group context. The range [0..group.rank).
 */
    std::vector<unsigned int> generators(const tc::Group &context) {
        std::vector<unsigned int> g_gens(context.rank());
        std::iota(g_gens.begin(), g_gens.end(), 0);
        return g_gens;
    }

/**
 * Determine whether the orientation of the group sg_gens is reversed from the group g_gens within group context
 */
    int get_parity(
        const tc::Group &context,
        const Symbol &g_gens,
        const Symbol &sg_gens
    ) {
        if (g_gens.size() != sg_gens.size() + 1) return 0;

        const auto proper_sg_gens = recontext_gens(context.rank(), g_gens, sg_gens);

        int i = 0;
        for (; i < sg_gens.size(); ++i) {
            if (proper_sg_gens[i] != i) {
                break;
            }
        }

        return i & 1;
    }

/**
 * Reverse the orientation of all primitives in this mesh.
 */
    void flip(std::vector<Symbol> &prims) {
        for (auto &prim: prims) {
            flip(prim);
        }
    }

/**
 * Convert the indexes of this mesh to those of a different context, using g_gens to build the parent context and sg_gens to build this context.
 */
    void recontext(
        std::vector<Symbol> &prims,
        const tc::Group &context,
        const Symbol &g_gens,
        const Symbol &sg_gens
    ) {
        const auto proper_sg_gens = recontext_gens(context.rank(), g_gens, sg_gens);

        const auto table = solve(context, g_gens, Symbol(0));
        const auto path = solve(context, sg_gens, Symbol(0)).path();

        auto map = path.walk(0U, proper_sg_gens, [&table](auto coset, auto gen) {
            return table.get(coset, gen);
        });

        for (Symbol &prim: prims) {
            for (auto &ind: prim) {
                ind = map[ind];
            }
        }

        if (get_parity(context, g_gens, sg_gens) == 1)
            flip(prims);
    }

/**
 * Union several meshes of the same dimension
 */
    std::vector<Symbol> merge(const std::vector<std::vector<Symbol>> &meshes) {
        size_t size = 0;
        for (const auto &mesh: meshes) {
            size += mesh.size();
        }

        std::vector<Symbol> res;
        res.reserve(size);
        for (const auto &mesh: meshes) {
            res.insert(res.end(), mesh.begin(), mesh.end());
        }

        return res;
    }

    std::vector<std::vector<Symbol>> each_tile(
        std::vector<Symbol> base,
        const tc::Group &context,
        const Symbol &g_gens,
        const Symbol &sg_gens
    ) {
        recontext(base, context, g_gens, sg_gens);

        const auto table = solve(context, g_gens, Symbol(0));
        const auto path = solve(context, g_gens, sg_gens).path();

        auto _gens = generators(context);

        auto res = path.walk(base, _gens, [&table](auto from, auto &gen) {
            apply(table, gen, from);
            return from;
        });

        return res;
    }

    [[nodiscard]]
    std::vector<Symbol> tile(
        std::vector<Symbol> base,
        const tc::Group &context,
        const Symbol &g_gens,
        const Symbol &sg_gens
    ) {
        auto res = each_tile(std::move(base), context, g_gens, sg_gens);

        return merge(res);
    }

/**
 * Produce a mesh of primitives that fill out the volume of the subgroup generated by generators g_gens within the group context
 */
    std::vector<Symbol> triangulate(
        const tc::Group &context,
        const Symbol &g_gens
    ) {
        if (g_gens.size() == 0) {
            std::vector<Symbol> res;
            Symbol prims(1);
            prims.setZero();
            res.push_back(prims);
            return res;
        }

        const auto &combos = combinations(g_gens, g_gens.size() - 1);

        std::vector<std::vector<Symbol>> meshes;

        for (const auto &sg_gens: combos) {
            auto base = triangulate(context, sg_gens);
            auto raised = tile(base, context, g_gens, sg_gens);
            raised.erase(raised.begin(), raised.begin() + base.size());
            auto fanned = fan(raised, 0);
            meshes.push_back(fanned);
        }

        const std::vector<Symbol> &result = merge(meshes);
        return result;
    }

    template<class T>
    auto hull(const tc::Group &group, T all_sg_gens, const std::vector<Symbol> &exclude) {
        std::vector<std::vector<Symbol>> parts;
        auto g_gens = group.gens;
        for (const Symbol &sg_gens: all_sg_gens) {
            bool excluded = false;
            for (const auto &test: exclude) {
                if (sg_gens == test) {
                    excluded = true;
                    break;
                }
            }
            if (excluded) continue;

            const auto &base = triangulate(group, sg_gens);
            const auto &tiles = each_tile(base, group, g_gens, sg_gens);
            for (const auto &tile: tiles) {
                parts.push_back(tile);
            }
        }
        return parts;
    }
}
