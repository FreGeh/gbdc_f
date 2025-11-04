#include "bits/stdc++.h"
#define XXH_INLINE_ALL
#include "src/external/xxhash/xxhash.h"
#include "src/util/CNFFormula.h"

using namespace std;

namespace WLF {

struct WLSettings {
    uint64_t max_iterations = 100;
    bool print_stats = false;
    bool sum_encoding = false;
    bool stat_encoding = false;
    bool quick_digest = false;
};

struct WLResult {
    uint64_t hash;
    int round = 0;
    int vertex_color_classes = 0;
    int edge_color_classes = 0;
    bool stabilized = false;
};

struct Accumulator {
    uint64_t sum = 0, sum2 = 0, xor_sum = 0, mn = UINT64_MAX, mx = 0;
    uint32_t deg = 0;
};

inline unsigned literal_to_vertex_idx(const Lit& lit) {
    return ((lit.var() - 1) << 1) | (lit.sign() ? 1u : 0u);
}

inline unsigned mate_of_vertex_idx(unsigned v_idx) {
    return v_idx ^ 1u;
}

using Color = uint64_t;

inline uint64_t hash_words(const uint64_t* data, size_t n_words) {
    return XXH3_64bits(static_cast<const void*>(data), n_words * sizeof(uint64_t));
}
inline uint64_t hash_words(const vector<uint64_t>& v) {
    return hash_words(v.data(), v.size());
}

inline bool same_partition(const vector<Color> &old_colors, const vector<Color> &new_colors) {
    if (old_colors.size() != new_colors.size()) return false;
    unordered_map<Color, Color> mapping;
    unordered_set<Color> used_new;
    mapping.reserve(old_colors.size());
    used_new.reserve(old_colors.size());

    for (int i=0;i<old_colors.size();i++) {
        Color x = old_colors[i];
        Color y = new_colors[i];
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

inline int get_color_classes(const vector<Color> &colors) {
    unordered_set<Color> unique;
    for (auto c : colors) {
        unique.insert(c);
    }
    return unique.size();
}

inline void standard_encoding(vector<uint64_t> &neigh, vector<uint64_t> &features) {
    sort(neigh.begin(), neigh.end());

    features.insert(features.end(), neigh.begin(), neigh.end());
}

inline void sum_encoding(vector<uint64_t> &neigh, vector<uint64_t> &features) {
    uint64_t sum = 0;

    for (uint64_t c : neigh) {
        sum += c; 
    }

    features.push_back(sum);
}

inline void stat_encoding(vector<uint64_t> &neigh, vector<uint64_t> &features) {
    Accumulator stats;
    for (uint64_t c : neigh) {
        stats.sum += c; 
        stats.sum2 += c * c; 
        stats.xor_sum ^= c;
        if (c < stats.mn) stats.mn = c; 
        else if (c > stats.mx) stats.mx = c;
    }
    if (neigh.empty()) { stats.mn = 0; stats.mx = 0; }

    features.push_back(stats.sum); 
    features.push_back(stats.sum2); 
    features.push_back(stats.xor_sum); 
    features.push_back(stats.mn); 
    features.push_back(stats.mx);
}


// iteration step of assigning new colors
inline void update_colors(const CNFFormula &formula, const WLSettings &settings, vector<Color> &old_V, vector<Color> &new_V, vector<Color> &new_E) {
    const int num_E = formula.nClauses();
    const int num_V = formula.nVars()*2;
    vector<Color> neigh;
    vector<uint64_t> clause_features;
    vector<Accumulator> acc(num_V);
    vector<vector<uint64_t>> in_clauses;
    if (!settings.sum_encoding && !settings.stat_encoding) {
        in_clauses.assign(num_V, {});   // allocate per literal
    }

    // clauses (edges)
    for (int e = 0; e < num_E; e++) {
        neigh.clear();
        clause_features.clear();

        for (const Lit &lit : *formula[e]) {
            neigh.push_back(old_V[literal_to_vertex_idx(lit)]);
        }
        clause_features.push_back(formula[e]->size()); 

        if (settings.stat_encoding) {
            stat_encoding(neigh, clause_features);
        } 
        else if (settings.sum_encoding) {
            sum_encoding(neigh, clause_features);
        }
        else {
            standard_encoding(neigh, clause_features);
        }

        new_E[e] = hash_words(clause_features);

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

    // literals (vertices)
    vector<uint64_t> literal_features;
    for (int v = 0; v < num_V; v++) {
        const uint64_t self  = old_V[v];
        const uint64_t mate  = old_V[mate_of_vertex_idx(v)];

        if (settings.stat_encoding) {
            const auto &A = acc[v];
            uint64_t feat[2 + 6];
            int k = 0;

            feat[k++] = self;
            feat[k++] = mate;
            feat[k++] = A.deg;
            feat[k++] = A.sum;
            feat[k++] = A.sum2;
            feat[k++] = A.xor_sum;
            feat[k++] = (A.mn == UINT64_MAX ? 0 : A.mn);
            feat[k++] = A.mx;

            new_V[v] = hash_words(feat, k);
        }
        else if (settings.sum_encoding) {
            const auto &A = acc[v];
            uint64_t feat[2 + 2];

            feat[0] = self;
            feat[1] = mate;
            feat[2] = A.deg;
            feat[3] = A.sum;

            new_V[v] = hash_words(feat, 4);
        }
        else {
            auto &lits = in_clauses[v];
            sort(lits.begin(), lits.end());

            literal_features.insert(literal_features.end(), lits.begin(), lits.end());

            new_V[v] = hash_words(literal_features);
        }
    }
}

// calculating final hash
inline uint64_t standard_digest(const vector<Color> &new_V, const vector<Color> &new_E) {
    map<Color, uint64_t> cntV, cntE;
    for (Color c : new_V) ++cntV[c];
    for (Color c : new_E) ++cntE[c];

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

inline uint64_t quick_digest(const vector<Color> &new_V, const vector<Color> &new_E) {
    unordered_map<Color, uint64_t> cntV, cntE;
    for (Color c : new_V) ++cntV[c];
    for (Color c : new_E) ++cntE[c];

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

// main function
WLResult run(const CNFFormula &formula, const WLSettings &settings) {
    WLResult stats;
    const size_t num_V = formula.nVars() * 2;
    const size_t num_E = formula.nClauses();
    vector<Color> old_V(num_V), new_V(num_V);
    vector<Color> old_E(num_E, 0), new_E(num_E, 0);

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
        update_colors(formula, settings, old_V, new_V, new_E);
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
        same_V = same_partition(old_V, new_V); 
        same_E = same_partition(old_E, new_E);
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

        old_V.swap(new_V); old_E.swap(new_E);
    }
    if (settings.quick_digest) {
        stats.hash = quick_digest(old_V, old_E);
    } else {
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

    WLResult res = run(f, settings);
    return to_hex64(res.hash);
}

}
