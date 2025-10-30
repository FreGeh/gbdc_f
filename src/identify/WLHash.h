#include "bits/stdc++.h"
#define XXH_INLINE_ALL
#include "src/external/xxhash/xxhash.h"
#include "src/util/CNFFormula.h"
#include "src/transform/cnf2hypergraph.h"

using namespace std;

namespace WLF {

struct WLSettings {
    uint64_t max_iterations = 100;
    bool print_stats = false;
    bool sum_sort = false;
    bool quick_digest = false;
};

struct WLResult {
    uint64_t hash;
    int round = 0;
    int vertex_color_classes = 0;
    int edge_color_classes = 0;
    bool stabilized = false;
};

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
    uint64_t sum = 0, sum2 = 0, xor_sum = 0, mn = UINT64_MAX, mx = 0;
    for (uint64_t c : neigh) {
        sum += c; 
        sum2 += c * c; 
        xor_sum ^= c;
        if (c < mn) mn = c; 
        if (c > mx) mx = c;
    }
    if (neigh.empty()) { mn = 0; mx = 0; }

    features.push_back(sum); 
    features.push_back(sum2); 
    features.push_back(xor_sum); 
    features.push_back(mn); 
    features.push_back(mx);
}


inline void initialize_colors(const CNF::IncidenceHypergraph &graph, vector<Color> &old_V, vector<Color> &old_E) {
    // edges
    for (int e = 0; e < graph.nClauses(); e++) {
        old_E[e]=graph.clauseSize(e);
    }

    // vertices
    for (int v = 0; v < graph.nVertices(); v++) {
        old_V[v]=graph.degree(v);
    }
}

// iteration step of assigning new colors
inline void update_colors(const CNF::IncidenceHypergraph &graph, const WLSettings &settings, vector<Color> &old_V, vector<Color> &old_E, vector<Color> &new_V, vector<Color> &new_E) {
    vector<Color> neigh;
    const int num_E = graph.nClauses();
    const int num_V = graph.nVertices();
    vector<vector<uint64_t>> edge_features(num_E), vertex_features(num_V);

    // edges (clauses)
    for (int e = 0; e < num_E; e++) {
        neigh.clear();

        auto verts = graph.literalsOfClause(e);
        for (int v_id : verts) {
            neigh.push_back(old_V[v_id]);
        }

        edge_features[e].push_back(old_E[e]);
        edge_features[e].push_back(neigh.size());

        if (settings.sum_sort) {
            sum_encoding(neigh, edge_features[e]);
        } 
        else {
            standard_encoding(neigh, edge_features[e]);
        }
    }

    // vertices (literals)
    for (int v = 0; v < num_V; v++) {
        neigh.clear();

        auto edges = graph.clausesOfLiteral(v);
        for (auto e_id : edges) {
            neigh.push_back(old_E[e_id]);
        }

        vertex_features[v].push_back(old_V[v]);
        vertex_features[v].push_back(old_V[graph.mateOf(v)]);
        vertex_features[v].push_back(neigh.size());
        if (settings.sum_sort) {
            sum_encoding(neigh, vertex_features[v]);
        } 
        else {
            standard_encoding(neigh, vertex_features[v]);
        }
    }

    new_E.resize(num_E);
    new_V.resize(num_V);

    //raw
    for (int e = 0; e < num_E; ++e) new_E[e] = hash_words(edge_features[e]);
    for (int v = 0; v < num_V; ++v) new_V[v] = hash_words(vertex_features[v]);
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
WLResult run(const CNF::IncidenceHypergraph &graph, const WLSettings &settings) {
    WLResult stats;
    vector<Color> old_V(graph.nVertices()), old_E(graph.nClauses()), new_V(graph.nVertices()), new_E(graph.nClauses());

    initialize_colors(graph, old_V, old_E);

    while (stats.round < settings.max_iterations && !stats.stabilized) {
        // refinement
        update_colors(graph, settings, old_V, old_E, new_V, new_E);
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
    CNF::IncidenceHypergraph g(f);

    WLResult res = run(g, settings);
    return to_hex64(res.hash);
}

}
