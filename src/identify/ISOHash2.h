/*************************************************************************************************
CNFTools -- Copyright (c) 2026, Ashlin Iser, Frederick Gehm, KIT - Karlsruhe Institute of Technology

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 **************************************************************************************************/

#ifndef ISOHASH2_H_
#define ISOHASH2_H_

#ifndef ISOHASH2_STABILIZATION_VARIANT
#define ISOHASH2_STABILIZATION_VARIANT 1
#endif

#ifndef ISOHASH2_DEDUP_CLAUSES
#define ISOHASH2_DEDUP_CLAUSES 1
#endif

#include <vector>
#include <algorithm>
#include <iostream>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <sstream>

#if ISOHASH2_STABILIZATION_VARIANT == 0
#include <unordered_map>
#endif

#define XXH_INLINE_ALL
#include "src/external/xxhash/xxhash.h"
#include "src/util/CNFFormula.h"

namespace CNF {

struct IsoHash2Settings {
    int max_iterations = 31; // 0 = until stabilized
    bool print_stats = false;
};

class IsoHash2 {
public:
    using Hash = uint64_t;
    using Clause = Cl;
    using Literal = Lit;

    struct Stats {
        Hash hash = 0;
        int round = 0;
        bool stabilized = false;
    };

private:
    struct LitColors {
        Hash val[2] = {1, 1};
    };

    struct ColorFunction {
        std::vector<LitColors> colors_by_var;
        explicit ColorFunction(size_t n_vars) : colors_by_var(n_vars + 1) {}

        inline Hash& operator()(Lit lit) {
            return colors_by_var[lit.var()].val[(int)lit.sign()];
        }
        inline const Hash& operator()(Lit lit) const {
            return colors_by_var[lit.var()].val[(int)lit.sign()];
        }
    };

    // mix64variant13 [Steele et al. 2014]
    static inline Hash fast_mix(Hash k) {
        k ^= k >> 30; k *= 0xbf58476d1ce4e5b9ULL;
        k ^= k >> 27; k *= 0x94d049bb133111ebULL;
        k ^= k >> 31;
        return k;
    }

    static inline Hash rotl64(Hash x, int r) {
        return (x << r) | (x >> (64 - r));
    }

#if ISOHASH2_STABILIZATION_VARIANT == 0
    struct StateKey {
        Hash p;
        Hash n;

        bool operator==(const StateKey& other) const {
            return p == other.p && n == other.n;
        }
    };

    struct StateKeyHasher {
        size_t operator()(const StateKey& key) const {
            Hash x = key.p ^ rotl64(key.n, 1);
            return (size_t)fast_mix(x + 0x9e3779b97f4a7c15ULL);
        }
    };
#endif

    struct VarState {
        Hash p;
        Hash n;
        size_t var_idx;

        bool operator<(const VarState& other) const {
            if (p != other.p) return p < other.p;
            return n < other.n;
        }

        bool operator==(const VarState& other) const {
            return p == other.p && n == other.n;
        }
    };

    static inline Hash literal_key(Literal lit) {
        return (((Hash)lit.var()) << 1) ^ (Hash)lit.sign();
    }

    static Hash clause_fingerprint(const Clause& clause) {
        Hash h = fast_mix((Hash)clause.size() + 0x9e3779b97f4a7c15ULL);

        for (const Literal lit : clause) {
            Hash x = literal_key(lit);
            h ^= fast_mix(x + 0xbf58476d1ce4e5b9ULL + rotl64(h, 7));
        }

        return fast_mix(h);
    }

    static bool same_clause(const Clause& a, const Clause& b) {
        return a.size() == b.size() && std::equal(a.begin(), a.end(), b.begin());
    }

    IsoHash2Settings settings;
    const CNFFormula& cnf;
    ColorFunction color_functions[2];
    Stats stats;

    std::vector<Hash> partition_buffer;
    std::vector<size_t> prev_partition_ids;
    std::vector<size_t> curr_partition_ids;
    std::vector<VarState> rank_buffer;
    size_t prev_partition_count = 0;

    std::vector<const Clause*> active_clauses;

    inline ColorFunction& old_color() {
        return color_functions[stats.round % 2];
    }

    inline const ColorFunction& old_color() const {
        return color_functions[stats.round % 2];
    }

    inline ColorFunction& new_color() {
        return color_functions[(stats.round + 1) % 2];
    }

    inline const ColorFunction& new_color() const {
        return color_functions[(stats.round + 1) % 2];
    }

    void build_view() {
        active_clauses.clear();
        active_clauses.reserve(cnf.nClauses());

#if ISOHASH2_DEDUP_CLAUSES == 0
        for (const auto& clause_ptr : cnf) {
            active_clauses.push_back(clause_ptr);
        }
        return;
#else
        struct ClauseData {
            Hash hash;
            const Clause* ptr;
        };

        auto clause_less = [](const Clause& a, const Clause& b) {
            if (a.size() != b.size()) return a.size() < b.size();

            for (size_t i = 0; i < a.size(); ++i) {
                Hash ka = literal_key(a[i]);
                Hash kb = literal_key(b[i]);

                if (ka != kb) return ka < kb;
            }

            return false;
        };

        std::vector<ClauseData> flat_list;
        flat_list.reserve(cnf.nClauses());

        for (const auto& clause_ptr : cnf) {
            flat_list.push_back({clause_fingerprint(*clause_ptr), clause_ptr});
        }

        std::sort(flat_list.begin(), flat_list.end(),
            [&](const ClauseData& a, const ClauseData& b) {
                if (a.hash != b.hash) return a.hash < b.hash;
                return clause_less(*a.ptr, *b.ptr);
            }
        );

        const Clause* last_kept = nullptr;

        for (const ClauseData& item : flat_list) {
            if (last_kept != nullptr && same_clause(*last_kept, *item.ptr)) {
                continue;
            }

            active_clauses.push_back(item.ptr);
            last_kept = item.ptr;
        }
#endif
    }

    inline Hash state_hash_oriented(const LitColors& lc) const {
        Hash p = lc.val[0];
        Hash n = lc.val[1];
        Hash x = p ^ rotl64(n, 1);
        return fast_mix(x + 0x9e3779b97f4a7c15ULL);
    }

    inline Hash state_hash_canonical(const LitColors& lc) const {
        Hash p = lc.val[0];
        Hash n = lc.val[1];
        if (p > n) std::swap(p, n);
        Hash x = p ^ rotl64(n, 1);
        return fast_mix(x + 0x9e3779b97f4a7c15ULL);
    }

    Hash clause_hash(const Clause& clause) const {
        Hash a = 0, b = 0;
        const auto& c_func = old_color();

        for (const Literal lit : clause) {
            Hash y = fast_mix(c_func(lit) + 0x9e3779b97f4a7c15ULL);
            a += y;
            b ^= rotl64(y, 23);
        }

        Hash combined = a ^ fast_mix(b + 0xbf58476d1ce4e5b9ULL) ^ (Hash)clause.size();
        return fast_mix(combined);
    }

    void finalize_literal_colors() {
        const size_t num_vars = cnf.nVars();
        auto* agg_vec = &new_color().colors_by_var;
        const auto* old_vec = &old_color().colors_by_var;

        for (size_t i = 1; i <= num_vars; ++i) {
            Hash old_p = (*old_vec)[i].val[0];
            Hash old_n = (*old_vec)[i].val[1];
            Hash agg_p = (*agg_vec)[i].val[0];
            Hash agg_n = (*agg_vec)[i].val[1];

            Hash x_p = old_p + fast_mix(agg_p) + rotl64(old_n, 1);
            Hash x_n = old_n + fast_mix(agg_n) + rotl64(old_p, 1);

            (*agg_vec)[i].val[0] = fast_mix(x_p);
            (*agg_vec)[i].val[1] = fast_mix(x_n);
        }
    }

    void iteration_step() {
        auto& nc_vec = new_color().colors_by_var;
        std::memset(nc_vec.data(), 0, nc_vec.size() * sizeof(LitColors));

        for (const Clause* clause_ptr : active_clauses) {
            Hash ch = clause_hash(*clause_ptr);
            for (const Literal lit : *clause_ptr) {
                new_color()(lit) += ch;
            }
        }
        finalize_literal_colors();
    }

    bool check_stabilization() {
        const size_t n = cnf.nVars();

        if (n == 0) {
            prev_partition_count = 0;
            return prev_partition_ids.empty();
        }

        const auto& current_colors = old_color().colors_by_var;

#if ISOHASH2_STABILIZATION_VARIANT == 0
        curr_partition_ids.assign(n, 0);

        std::unordered_map<StateKey, size_t, StateKeyHasher> ids;
        ids.reserve(n * 2 + 1);

        size_t next_id = 0;

        for (size_t i = 1; i <= n; ++i) {
            StateKey key { current_colors[i].val[0], current_colors[i].val[1] };
            auto result = ids.emplace(key, next_id);

            if (result.second) {
                ++next_id;
            }

            curr_partition_ids[i - 1] = result.first->second;
        }

        bool stable = prev_partition_ids.size() == curr_partition_ids.size() && prev_partition_ids == curr_partition_ids;
        prev_partition_ids.swap(curr_partition_ids);
        prev_partition_count = next_id;
        return stable;
#else
        for (size_t i = 1; i <= n; ++i) {
            rank_buffer[i - 1] = {
                current_colors[i].val[0],
                current_colors[i].val[1],
                i - 1
            };
        }

        std::sort(rank_buffer.begin(), rank_buffer.end());

        size_t partition_count = 0;
        size_t begin = 0;

        while (begin < n) {
            size_t end = begin + 1;
            size_t min_var_idx = rank_buffer[begin].var_idx;

            while (end < n && rank_buffer[end] == rank_buffer[begin]) {
                min_var_idx = std::min(min_var_idx, rank_buffer[end].var_idx);
                ++end;
            }

            for (size_t i = begin; i < end; ++i) {
                curr_partition_ids[rank_buffer[i].var_idx] = min_var_idx;
            }

            ++partition_count;
            begin = end;
        }

        bool stable = prev_partition_ids.size() == curr_partition_ids.size()
                   && prev_partition_ids == curr_partition_ids;

        prev_partition_ids.swap(curr_partition_ids);
        prev_partition_count = partition_count;
        return stable;
#endif
    }

public:
    IsoHash2(const CNFFormula& formula, const IsoHash2Settings& s) :
        settings(s),
        cnf(formula),
        color_functions{ColorFunction(cnf.nVars()), ColorFunction(cnf.nVars())},
        partition_buffer(cnf.nVars()),
        prev_partition_ids(cnf.nVars(), 0),
        curr_partition_ids(cnf.nVars(), 0),
        rank_buffer(cnf.nVars())
    {
        build_view();
    }

    Stats run() {
        stats = Stats{};
        prev_partition_count = cnf.nVars() > 0 ? 1 : 0;

        std::fill(prev_partition_ids.begin(), prev_partition_ids.end(), 0);
        std::fill(curr_partition_ids.begin(), curr_partition_ids.end(), 0);
        std::fill(color_functions[0].colors_by_var.begin(), color_functions[0].colors_by_var.end(), LitColors{});
        std::fill(color_functions[1].colors_by_var.begin(), color_functions[1].colors_by_var.end(), LitColors{});

        while (stats.round < settings.max_iterations || settings.max_iterations == 0) {
            iteration_step();
            stats.round++;

            bool stable = check_stabilization();

            if (settings.print_stats) std::cerr << "c Round " << stats.round << " partitions: " << prev_partition_count << "\n";

            if (stable) {
                stats.stabilized = true;
                // if (settings.print_stats) { 
                    std::cerr << "c Stabilized after " << stats.round << " rounds.\n";
                // }
                break;
            }
        }

        if (!stats.stabilized && settings.print_stats) {
            std::cerr << "c Reached max iterations (" << settings.max_iterations << ").\n";
        }

        const size_t n = cnf.nVars();
        if (partition_buffer.size() != n) partition_buffer.resize(n);

        const auto& final_colors = old_color().colors_by_var;
        for (size_t i = 1; i <= n; ++i) {
            partition_buffer[i - 1] = state_hash_canonical(final_colors[i]);
        }
        std::sort(partition_buffer.begin(), partition_buffer.end());
        stats.hash = XXH3_64bits(partition_buffer.data(), partition_buffer.size() * sizeof(Hash));
        return stats;
    }
};

inline IsoHash2::Stats isohash2_stats(const char* filename, const IsoHash2Settings& s = {}) {
    CNFFormula cnf(filename);
    cnf.normalizeVariableNames();

    IsoHash2 hasher(cnf, s);
    return hasher.run();
}

inline std::string isohash2(const char* filename, const IsoHash2Settings& s = {}) {
    const auto stats = isohash2_stats(filename, s);

    std::ostringstream oss;
    oss << std::hex << std::setw(16) << std::setfill('0') << std::nouppercase << stats.hash;

    return oss.str();
}

} // namespace CNF

#endif // ISOHASH2_H_