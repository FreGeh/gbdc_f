#include "bits/stdc++.h"
#define XXH_INLINE_ALL
#include "src/external/xxhash/xxhash.h"
#include "src/util/CNFFormula.h"

using namespace std;

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
    using Hash = conditional_t<half_bit_hash, uint32_t, uint64_t>;
    using Clause = typename CNFFormula::Cl;
    using Literal = typename CNFFormula::Lit;

    struct Stats {
        Hash hash;
        int round = 0;
        int vertex_color_classes = 0;
        int edge_color_classes = 0;
        bool stabilized = false;
    };

private:
    struct ColorFunction {
        vector<LitColors> colors;
        explicit ColorFunction(size_t n) : colors(2*n + 2, 0) {}
        inline Hash &operator ()(Literal lit) { 
            return literal_colors[static_cast<size_t>(lit)];
        }
    };

ColorFunction color_functions[2];
const CNFFormula &cnf;
const WLSettings settings;
Stats stats;

inline ColorFunction &old_color() { return color_functions[stats.round & 1U]; }
inline ColorFunction &new_color() { return color_functions[(stats.round + 1U) & 1U]; }

template <typename D>
static inline Hash hash(const D data) {
    return XXH3_64bits(&data, sizeof(data));
}

struct Add {
    inline void operator()(Hash &h, Hash v) const { h += v; }
};

struct Xor {
    inline void operator()(Hash &h, Hash v) const { h ^= v; }
};

template <typename D, typename C, typename M, typename F>
static inline Hash hash_container(const C &container, M map, F combine) {
    Hash hash = 0;
    for (const D &element : container)
        combine(hash, map(element));
    return hash;
}

WLHasher(const char* filename, const WLSettings settings) : 
    settings(settings),
    cnf(filename),
    color_functions {ColorFunction(cnf.nVars()), ColorFunction(cnf.nVars())}
{}

void mate_reference() {
    for (LitColors& lit_colors : old_color().colors) {
        lit_colors.mate_reference();
    }
}

using Comb = conditional_t<use_xor, Xor, Add>;
Hash clause_hash(const Clause cl) {
    if (!settings.use_sort) {
        Hash clause_hash = hash_container<Literal>(cl, [&](Literal lit){ return old_color()(lit); },Comb{});
        // hash again to preserve clause structure
        if (settings.rehash) clause_hash = hash(clause_hash);
        return clause_hash;
    } else {
        vector<Hash> sorted;
        sorted.reserve(cl.size());
        for (const Lit lit : cl) {
            sorted.push_back(old_color()(lit));
        }
        std::sort(sorted.begin(), sorted.end());
        return XXH3_64bits(sorted.data(), sorted.size() * sizeof(Hash));
    }
};

Hash variable_hash() {
    return hash_container<LitColors>(old_color().colors, [&](LitColors lc){ return lc.variable_hash(); });
}

void iteration_step() {
    mate_reference();
    for (const Clause cl : cnf.clauses()) {
        const Hash clh = clause_hash(cl) : cfg.rehash ? hash(cl.size()) : cl.size();
        for (const Lit lit : cl) {
            hash_container<LitColors>(&new_color()(lit), clh);
        }
    }
    stats.round++;
}

inline bool same_partition(const vector<Hash> &old_colors, const vector<Hash> &new_colors) {
    if (old_colors.size() != new_colors.size()) return false;
    unordered_map<Hash, Hash> mapping;
    unordered_set<Hash> used_new;
    mapping.reserve(old_colors.size());
    used_new.reserve(old_colors.size());

    for (int i=0;i<old_colors.size();i++) {
        Hash x = old_colors[i];
        Hash y = new_colors[i];
        auto it = mapping.find(x);
        
        if (it == mapping.end()) {
            if (!used_new.insert(y).second) return false;
            else mapping.emplace(x,y);
        } else {
            if (it->second != y) return false;
        }
    }
    return true;
}

inline int get_color_classes(const vector<Hash> &colors) {
    unordered_set<Hash> unique;
    for (auto c : colors) {
        unique.insert(c);
    }
    return unique.size();
}

// iteration step of assigning new colors
inline void update_colors(const CNFFormula &formula, const WLSettings &settings, vector<Hash> &old_V, vector<Hash> &new_V, vector<Hash> &new_E) {
    const int num_E = formula.nClauses();
    const int num_V = formula.nVars()*2;

    XXH3StatePtr state(XXH3_createState());

    vector<Accumulator> acc(num_V);
    vector<vector<uint64_t>> in_clauses;
    if (!settings.sum_encoding && !settings.stat_encoding) {
        in_clauses.assign(num_V, {});   // allocate per literal
    }

    // update clause (edge) colors
    for (int e = 0; e < num_E; e++) {
        const auto &clause = *formula[e];
        const uint64_t clause_size = clause.size();

        XXH3_64bits_reset(state.get());
        XXH3_64bits_update(state.get(), &clause_size, sizeof(clause_size));

        if (settings.stat_encoding) {
            Accumulator stats;
            for (const Lit &lit : clause) {
                uint64_t c = old_V[literal_to_vertex_idx(lit)];
                stats.sum += c;
                stats.sum2 += c * c;
                stats.xor_sum ^= c;
                if (c < stats.mn) stats.mn = c;
                if (c > stats.mx) stats.mx = c;
            }
            if (clause.empty()) { stats.mn = 0; }

            // Hash stats directly from a stack-allocated array. NO push_backs.
            uint64_t features[5] = {stats.sum, stats.sum2, stats.xor_sum, stats.mn, stats.mx};
            XXH3_64bits_update(state.get(), features, sizeof(features));
            new_E[e] = XXH3_64bits_digest(state.get());
        } 
        else if (settings.sum_encoding) {
            uint64_t sum = 0;
            for (const Lit &lit : clause) {
                sum += old_V[literal_to_vertex_idx(lit)];
            }
            XXH3_64bits_update(state.get(), &sum, sizeof(sum));
            new_E[e] = XXH3_64bits_digest(state.get());
        }
        else {
            vector<Hash> neigh;
            neigh.reserve(clause.size());
            for (const Lit &lit : clause) {
                neigh.push_back(old_V[literal_to_vertex_idx(lit)]);
            }
            sort(neigh.begin(), neigh.end());
            XXH3_64bits_update(state.get(), neigh.data(), neigh.size() * sizeof(Hash));
            new_E[e] = XXH3_64bits_digest(state.get());
        }

        for (const Lit& lit : *formula[e]) {
            unsigned v = literal_to_vertex_idx(lit);

            if (settings.stat_encoding) {
                acc[v].deg += 1;
                acc[v].sum  += new_E[e];
                acc[v].sum2 += new_E[e] * new_E[e];
                acc[v].xor_sum ^= new_E[e];
                if (new_E[e] < acc[v].mn) acc[v].mn = new_E[e];
                if (new_E[e] > acc[v].mx) acc[v].mx = new_E[e];
            } 
            else if (settings.sum_encoding) {
                acc[v].deg += 1;
                acc[v].sum += new_E[e];
            }
            else {
                in_clauses[v].push_back(new_E[e]);
            }
        }
    }

    // update literal (vertex) colors
    for (int v = 0; v < num_V; v++) {
        const uint64_t self  = old_V[v];
        const uint64_t mate  = old_V[mate_of_vertex_idx(v)];

        if (settings.stat_encoding) {
            const auto &A = acc[v];
            uint64_t feat[8] = {
                self, mate, (uint64_t)A.deg, A.sum, A.sum2, A.xor_sum,
                (A.mn == UINT64_MAX ? 0 : A.mn), A.mx
            };
            new_V[v] = hash_words(feat, 8);
        }
        else if (settings.sum_encoding) {
            const auto &A = acc[v];
            uint64_t feat[4] = {self, mate, (uint64_t)A.deg, A.sum};
            new_V[v] = hash_words(feat, 4);
        }
        else {
            XXH3_64bits_reset(state.get());

            auto &lits = in_clauses[v];
            sort(lits.begin(), lits.end());

            XXH3_64bits_update(state.get(), &self, sizeof(self));
            XXH3_64bits_update(state.get(), &mate, sizeof(mate));
            XXH3_64bits_update(state.get(), lits.data(), lits.size() * sizeof(Hash));
            new_V[v] = XXH3_64bits_digest(state.get());
        }
    }
}

inline void update_colors_one_pass(const CNFFormula &formula, vector<Hash> &old_V, vector<Hash> &new_V, vector<Hash> &new_E) {
    const int num_E = formula.nClauses();
    const int num_V = formula.nVars() * 2;
    new_V.assign(num_V, 0);

    for (int e = 0; e < num_E; e++) {
        const auto &clause = *formula[e];

        uint64_t sum = 0;
        for (const Lit &lit : clause) {
            sum += old_V[literal_to_vertex_idx(lit)];
        }
        new_E[e] = XXH3_64bits(&sum, sizeof(sum)); 

        for (const Lit& lit : clause) {
            unsigned v_idx = literal_to_vertex_idx(lit);
            new_V[v_idx] += new_E[e];
        }
    }
}

// calculating final hash
inline uint64_t standard_digest(const vector<Hash> &new_V, const vector<Hash> &new_E) {
    map<Hash, uint64_t> cntV, cntE;
    for (Hash c : new_V) ++cntV[c];
    for (Hash c : new_E) ++cntE[c];

    vector<uint64_t> buf;
    buf.reserve(2 * (cntV.size() + cntE.size()) + 2);

    buf.push_back((uint64_t)cntV.size());
    for (const auto &p : cntV) { 
        buf.push_back(p.first); 
        buf.push_back(p.second); 
    }

    buf.push_back((uint64_t)cntE.size());
    for (const auto &p : cntE) { 
        buf.push_back(p.first); 
        buf.push_back(p.second); 
    }

    return XXH3_64bits(buf.data(), buf.size()*sizeof(uint64_t));
}

inline uint64_t quick_digest(const vector<Hash> &new_V, const vector<Hash> &new_E) {
    unordered_map<Hash, uint64_t> cntV, cntE;
    for (Hash c : new_V) ++cntV[c];
    for (Hash c : new_E) ++cntE[c];

    uint64_t final_digest = 0;

    for (const auto &p : cntV) { 
        uint64_t pair[2] = {p.first, p.second};
        uint64_t pair_hash = XXH3_64bits(pair, sizeof(pair));
        final_digest ^= pair_hash;
    }

    for (const auto &p : cntE) { 
        uint64_t pair[2] = {p.first, p.second};
        uint64_t pair_hash = XXH3_64bits(pair, sizeof(pair));
        final_digest ^= pair_hash;
    }

    uint64_t class_counts[2] = {cntV.size(), cntE.size()};
    final_digest ^= XXH3_64bits(class_counts, sizeof(class_counts));

    return final_digest;
}

inline uint64_t super_quick_digest(const vector<Hash> &V, const vector<Hash> &E) {
    uint64_t xor_sum = 0;
    for (Hash c : V) {
        xor_sum ^= c;
    }

    uint64_t e_xor_sum = 0;
    for (Hash c : E) {
        e_xor_sum ^= c;
    }
    
    uint64_t final_sums[2] = {xor_sum, e_xor_sum};
    return XXH3_64bits(final_sums, sizeof(final_sums));
}

// main function
Stats run(const CNFFormula &formula, const WLSettings &settings) {
    Stats stats;
    const size_t num_V = formula.nVars() * 2;
    const size_t num_E = formula.nClauses();
    vector<Hash> old_V(num_V), new_V(num_V);
    vector<Hash> old_E(num_E, 0), new_E(num_E, 0);
    int prev_V_classes = 0, prev_E_classes = 0;

    vector<unsigned> degrees(num_V, 0);
    for (unsigned e = 0; e < num_E; ++e) {
        for (const Lit &lit : *formula[e]) {
            degrees[literal_to_vertex_idx(lit)]++;
        }
    }

    for (int e = 0; e < num_E; e++) old_E[e] = formula[e]->size();
    for (int v = 0; v < num_V; v++) old_V[v] = degrees[v];

    while (stats.round < settings.max_iterations && !stats.stabilized) {
        // refinement
        if (settings.one_pass) {
            update_colors_one_pass(formula, old_V, new_V, new_E);
        } else {
            update_colors(formula, settings, old_V, new_V, new_E);
        }
        stats.round++;

        // per iteration stats
        stats.vertex_color_classes = get_color_classes(new_V);
        stats.edge_color_classes = get_color_classes(new_E);
        if (settings.print_stats) {
            cerr << "round " << (stats.round)
                 << "  V_classes=" << stats.vertex_color_classes
                 << "  E_classes=" << stats.edge_color_classes << "\n";
        }

        // stabilization
        bool same_V=false, same_E=false;
        if (settings.simple_stab_check) {
            same_V = (prev_V_classes == stats.vertex_color_classes);
            same_E = (prev_E_classes == stats.edge_color_classes);
        } else {
            same_V = same_partition(old_V, new_V); 
            same_E = same_partition(old_E, new_E);
        }
        if (same_V && same_E) {
            if (settings.quick_digest) {
                stats.hash = quick_digest(new_V, new_E);
            } else {
                stats.hash = standard_digest(new_V, new_E);
            }
            stats.stabilized=true;
            cerr << "c stabilized: " << stats.round << "\n";
            return stats;
        }

        prev_V_classes = stats.vertex_color_classes;
        prev_E_classes = stats.edge_color_classes;

        old_V.swap(new_V); old_E.swap(new_E);
    }
    if (settings.quick_digest) {
        stats.hash = quick_digest(old_V, old_E);
    } 
    else if (settings.super_quick_digest) {
        stats.hash = super_quick_digest(old_V, old_E);
    }
    else {
        stats.hash = standard_digest(old_V, old_E);
    }
    return stats;
}


inline string to_hex64(uint64_t x) {
    ostringstream oss;
    oss << hex << nouppercase << setfill('0') << setw(16)
        << static_cast<unsigned long long>(x);
    return oss.str();
}

inline string wlhash(const char* filename, WLSettings settings) {
    CNFFormula f(filename);
    cerr << "c parsed vars=" << f.nVars() << " clauses=" << f.nClauses() << "\n";

    f.normalizeVariableNames();

    Stats res = run(f, settings);
    return to_hex64(res.hash);
}
};
}
