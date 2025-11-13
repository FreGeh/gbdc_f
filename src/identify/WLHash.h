#ifndef WLF_WLHASHER_H
#define WLF_WLHASHER_H

#include <vector>
#include <numeric>
#include <algorithm>
#include <type_traits>
#include <iostream>
#include <unordered_set>

#define XXH_INLINE_ALL
#include "src/external/xxhash/xxhash.h"
#include "src/util/CNFFormula.h"

namespace WLF {

struct WLSettings {
    uint64_t max_iterations = 100;
    bool print_stats = false;
};

template <
    bool one_pass,
    bool use_xor,
    bool use_sort,
    bool mate_ref,
    bool half_bit_hash
>
class WLHasher {
public:
    using Hash = std::conditional_t<half_bit_hash, uint32_t, uint64_t>;
    using Clause = Cl;

    struct Stats {
        Hash hash = 0;
        int round = 0;
        bool stabilized = false;
    };

private:
    struct LitColors {
        Hash p = 1;
        Hash n = 1;
    };
    struct ColorFunction {
        std::vector<LitColors> colors_by_var;

        explicit ColorFunction(size_t n_vars) : colors_by_var(n_vars + 1) {}

        inline Hash& operator()(Lit lit) {
            return lit.sign() ? colors_by_var[lit.var()].n : colors_by_var[lit.var()].p;
        }
        inline const Hash& operator()(Lit lit) const {
            return lit.sign() ? colors_by_var[lit.var()].n : colors_by_var[lit.var()].p;
        }
    };

    struct Add { inline void operator()(Hash& h, Hash v) const { h += v; } };
    struct Xor  { inline void operator()(Hash& h, Hash v) const { h ^= v; } };
    using Combiner = std::conditional_t<use_xor, Xor, Add>;

    const WLSettings& settings;
    const CNFFormula& cnf;
    ColorFunction color_functions[2];
    Stats stats;
    size_t prev_partition_count = 0;

    inline ColorFunction& old_color() { return color_functions[stats.round % 2]; }
    inline const ColorFunction& old_color() const { return color_functions[stats.round % 2]; }

    inline ColorFunction& new_color() { return color_functions[(stats.round + 1) % 2]; }
    inline const ColorFunction& new_color() const { return color_functions[(stats.round + 1) % 2]; }

    Hash variable_hash(const LitColors& lc) const {
        if (lc.p > lc.n) {
            const LitColors canonical_pair = {lc.n, lc.p};
            return XXH3_64bits(&canonical_pair, sizeof(canonical_pair));
        }
        return XXH3_64bits(&lc, sizeof(lc));
    }

    Hash clause_hash(const Clause& clause) const {
        if constexpr (use_sort) {
            std::vector<Hash> sorted_colors;
            sorted_colors.reserve(clause.size());
            for (const Lit lit : clause) {
                sorted_colors.push_back(old_color()(lit));
            }
            std::sort(sorted_colors.begin(), sorted_colors.end());
            return XXH3_64bits(sorted_colors.data(), sorted_colors.size() * sizeof(Hash));
        } else {
            Hash combined_hash = 0;
            Combiner combine;
            for (const Lit lit : clause) {
                combine(combined_hash, old_color()(lit));
            }
            return XXH3_64bits(&combined_hash, sizeof(combined_hash));
        }
    }

    void finalize_literal_colors() {
        auto* agg_colors_vec = &new_color().colors_by_var;
        const auto* old_colors_vec = &old_color().colors_by_var;

        for (size_t var_id = 1; var_id <= cnf.nVars(); ++var_id) {
            auto& agg_lc = (*agg_colors_vec)[var_id];
            const auto& old_lc = (*old_colors_vec)[var_id];

            const Hash old_p_color = old_lc.p;
            const Hash old_n_color = old_lc.n;
            const Hash agg_p_color = agg_lc.p;
            const Hash agg_n_color = agg_lc.n;

            {
                Hash features[3] = { old_p_color, agg_p_color, 0 };
                if constexpr (mate_ref) { features[2] = old_n_color; }
                agg_lc.p = XXH3_64bits(features, sizeof(features));
            }
            {
                Hash features[3] = { old_n_color, agg_n_color, 0 };
                if constexpr (mate_ref) { features[2] = old_p_color; }
                agg_lc.n = XXH3_64bits(features, sizeof(features));
            }
        }
    }


    void iteration_step_one_pass() {
        for(auto& lc : new_color().colors_by_var) { lc.p = 0; lc.n = 0; }
        Combiner combine;

        for (const auto& clause_ptr : cnf) {
            Hash ch = clause_hash(*clause_ptr);
            for (const Lit lit : *clause_ptr) {
                combine(new_color()(lit), ch);
            }
        }
        finalize_literal_colors();
    }

    void iteration_step_two_pass() {
        std::vector<Hash> clause_hashes(cnf.nClauses());
        for (size_t i = 0; i < cnf.nClauses(); ++i) {
            clause_hashes[i] = clause_hash(*cnf[i]);
        }
        
        for(auto& lc : new_color().colors_by_var) { lc.p = 0; lc.n = 0; }
        Combiner combine;

        for (size_t i = 0; i < cnf.nClauses(); ++i) {
            for (const Lit lit : *cnf[i]) {
                combine(new_color()(lit), clause_hashes[i]);
            }
        }
        finalize_literal_colors();
    }
    
    bool check_stabilization() {
        std::unordered_set<Hash> partitions;
        partitions.reserve(cnf.nVars());

        for (size_t var_id = 1; var_id <= cnf.nVars(); ++var_id) {
            partitions.insert(variable_hash(new_color().colors_by_var[var_id]));
        }

        if (partitions.size() == prev_partition_count) {
            return true;
        }

        prev_partition_count = partitions.size();
        return false;
    }

public:
    WLHasher(const CNFFormula& formula, const WLSettings& s) :
        settings(s),
        cnf(formula),
        color_functions{ColorFunction(cnf.nVars()), ColorFunction(cnf.nVars())}
    {}

    Stats run() {
        while (stats.round < settings.max_iterations) {
            if constexpr (one_pass) {
                iteration_step_one_pass();
            } else {
                iteration_step_two_pass();
            }
            stats.round++;

            if (settings.print_stats) {
                std::cerr << "c Round " << stats.round << " completed. Partitions: " << prev_partition_count << std::endl;
            }
            
            if (check_stabilization()) {
                stats.stabilized = true;
                if (settings.print_stats) {
                    std::cerr << "c Stabilized after " << stats.round << " rounds." << std::endl;
                }
                break;
            }
        }

        if (!stats.stabilized && settings.print_stats) {
            std::cerr << "c Reached max iterations (" << settings.max_iterations << ")." << std::endl;
        }

        Hash final_hash_accumulator = 0;
        Combiner combine;
        for (size_t var_id = 1; var_id <= cnf.nVars(); ++var_id) {
            combine(final_hash_accumulator, variable_hash(new_color().colors_by_var[var_id]));
        }
        stats.hash = XXH3_64bits(&final_hash_accumulator, sizeof(final_hash_accumulator));
        
        return stats;
    }
};

} // namespace WLF

#endif // WLF_WLHASHER_H