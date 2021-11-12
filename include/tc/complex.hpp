#pragma once

#include <tc/group.hpp>
#include <tc/solver.hpp>
#include <cmath>
#include <optional>
#include <numeric>
#include <iostream>

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

/**
 * An primitive stage N indices.
 * @tparam N
 */
    template<unsigned N>
    struct Primitive {
        static_assert(N > 0, "Primitives must contain at least one point. Primitive<0> or lower is impossible.");

        std::array<unsigned, N> inds;

        Primitive() = default;

        Primitive(const Primitive<N> &) = default;

        Primitive(const Primitive<N - 1> &sub, unsigned root) {
            std::copy(sub.inds.begin(), sub.inds.end(), inds.begin());
            inds[N - 1] = root;
        }

        ~Primitive() = default;

        inline void flip() {
            if (N > 1) std::swap(inds[0], inds[1]);
        }

        void apply(const tc::Cosets &table, unsigned int gen) {
            for (auto &ind: inds) {
                ind = table.get(ind, gen);
            }
            flip();
        }
    };

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
 * Apply some context transformation to all primitives of this mesh.
 */
    template<unsigned N>
    std::vector<Primitive<N>> apply(std::vector<Primitive<N>> prims, const tc::Cosets &table, unsigned int gen) {
        for (auto &prim: prims) {
            prim.apply(table, gen);
        }
        return prims;
    }

/**
 * Reverse the orientation of all primitives in this mesh.
 */
    template<unsigned N>
    void flip(std::vector<Primitive<N>> prims) {
        for (auto &prim: prims) {
            prim.flip();
        }
    }

/**
 * Convert the indexes of this mesh to those of a different context, using g_gens to build the parent context and sg_gens to build this context.
 */
    template<unsigned N>
    [[nodiscard]]
    std::vector<Primitive<N>> recontext(
        std::vector<Primitive<N>> prims,
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

        std::vector<Primitive<N>> res(prims);
        for (Primitive<N> &prim: res) {
            for (auto &ind: prim.inds) {
                ind = map[ind];
            }
        }

        if (get_parity(context, g_gens, sg_gens) == 1)
            flip(res);

        return res;
    }

/**
 * Union several meshes of the same dimension
 */
    template<unsigned N>
    std::vector<Primitive<N>> merge(const std::vector<std::vector<Primitive<N>>> &meshes) {
        size_t size = 0;
        for (const auto &mesh: meshes) {
            size += mesh.size();
        }

        std::vector<Primitive<N>> res;
        res.reserve(size);
        for (const auto &mesh: meshes) {
            res.insert(res.end(), mesh.begin(), mesh.end());
        }

        return res;
    }

    template<unsigned N>
    [[nodiscard]]
    std::vector<std::vector<Primitive<N>>> each_tile(
        std::vector<Primitive<N>> prims,
        const tc::Group &context,
        const Symbol &g_gens,
        const Symbol &sg_gens
    ) {
        std::vector<Primitive<N>> base = recontext(prims, context, g_gens, sg_gens);

        const auto table = solve(context, g_gens, Symbol(0));
        const auto path = solve(context, g_gens, sg_gens).path();

        auto _gens = generators(context);

        auto res = path.walk(base, _gens, [&table](auto from, auto gen){
            return apply(from, table, gen);
        });

        return res;
    }

    template<unsigned N>
    [[nodiscard]]
    std::vector<Primitive<N>> tile(
        std::vector<Primitive<N>> prims,
        const tc::Group &context,
        const Symbol &g_gens,
        const Symbol &sg_gens
    ) {
        auto res = each_tile<N>(prims, context, g_gens, sg_gens);

        return merge(res);
    }

/**
 * Produce a mesh of higher dimension by fanning a single point to all primitives in this mesh.
 */
    template<unsigned N>
    [[nodiscard]]
    std::vector<Primitive<N + 1>> fan(std::vector<Primitive<N>> prims, int root) {
        std::vector<Primitive<N + 1>> res(prims.size());
        std::transform(prims.begin(), prims.end(), res.begin(),
            [root](const Primitive<N> &prim) {
                return Primitive<N + 1>(prim, root);
            }
        );
        return res;
    }

/**
 * Produce a mesh of primitives that fill out the volume of the subgroup generated by generators g_gens within the group context
 */
    template<unsigned N>
    std::vector<Primitive<N>> triangulate(
        const tc::Group &context,
        const Symbol &g_gens
    ) {
        if (g_gens.size() + 1 != N) // todo make static assert
            throw std::logic_error("g_gens size must be one less than N");

        const auto &combos = combinations(g_gens, g_gens.size() - 1);

        std::vector<std::vector<Primitive<N>>> meshes;

        for (const auto &sg_gens: combos) {
            auto base = triangulate<N - 1>(context, sg_gens);
            auto raised = tile(base, context, g_gens, sg_gens);
            raised.erase(raised.begin(), raised.begin() + base.size());
            meshes.push_back(fan(raised, 0));
        }

        return merge(meshes);
    }

/**
 * Single-index primitives should not be further triangulated.
 */
    template<>
    std::vector<Primitive<1>> triangulate(
        const tc::Group &context,
        const Symbol &g_gens
    ) {
        if (g_gens.size() != 0) // todo make static assert
            throw std::logic_error("g_gens must be empty for a trivial Mesh");

        std::vector<Primitive<1>> res;
        res.emplace_back();
        return res;
    }

    template<unsigned N, class T>
    auto hull(const tc::Group &group, T all_sg_gens, const std::vector<Symbol> &exclude) {
        std::vector<std::vector<Primitive<N>>> parts;
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

            const auto &base = triangulate<N>(group, sg_gens);
            const auto &tiles = each_tile(base, group, g_gens, sg_gens);
            for (const auto &tile: tiles) {
                parts.push_back(tile);
            }
        }
        return parts;
    }
}
