#include "bits/stdc++.h"
#define XXH_INLINE_ALL
#include "src/external/xxhash/xxhash.h"
#include "src/util/CNFFormula.h"
#include "src/transform/cnf2hypergraph.h"

using namespace std;

namespace WLF {

struct WLSettings {
    uint64_t max_iterations = 10;
    bool print_stats = false;
    bool canonicalize_color_classes = false;
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
inline uint64_t hash_words(const std::vector<uint64_t>& v) {
    return hash_words(v.data(), v.size());
}

inline void encode(const std::vector<Color>& sorted, std::vector<uint64_t>& pairs) {
    pairs.clear();
    if (sorted.empty()) return;
    Color cur = sorted[0];
    uint64_t cnt = 1;
    for (size_t i = 1; i < sorted.size(); i++) {
        if (sorted[i] == cur) { 
            ++cnt; 
        }
        else { 
            pairs.push_back(cur); 
            pairs.push_back(cnt); 
            cur = sorted[i]; cnt = 1; 
        }
    }
    pairs.push_back(cur); pairs.push_back(cnt);
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

inline uint64_t digest(const vector<Color> &new_V, const vector<Color> &new_E) {
    unordered_map<Color, uint64_t> cntV, cntE;
    for (Color c : new_V) ++cntV[c];
    for (Color c : new_E) ++cntE[c];

    vector<pair<Color,uint64_t>> histV, histE;
    for (auto &kv : cntV) histV.push_back(kv);
    for (auto &kv : cntE) histE.push_back(kv);
    sort(histV.begin(), histV.end(),[](auto &a, auto &b){ return a.first < b.first; });
    sort(histE.begin(), histE.end(),[](auto &a, auto &b){ return a.first < b.first; });

    vector<uint64_t> buf;
    buf.push_back((uint64_t)histV.size());
    for (auto &p : histV) { buf.push_back(p.first); buf.push_back(p.second); }
    buf.push_back((uint64_t)histE.size());
    for (auto &p : histE) { buf.push_back(p.first); buf.push_back(p.second); }

    return XXH3_64bits(buf.data(), buf.size()*sizeof(uint64_t));
}

inline void update_colors(const CNF::IncidenceHypergraph &graph, const WLSettings &settings, vector<Color> &old_V, vector<Color> &old_E, vector<Color> &new_V, vector<Color> &new_E) {
    std::vector<Color> neigh;
    std::vector<uint64_t> pairs;
    const int num_E = graph.nClauses();
    const int num_V = graph.nVertices();
    vector<vector<uint64_t>> sign_E(num_E), sign_V(num_V);

    // edges
    for (int e = 0; e < num_E; e++) {
        neigh.clear();
        auto verts = graph.literalsOfClause(e);
        for (int v_id : verts) {
            neigh.push_back(old_V[v_id]);
        }
        sort(neigh.begin(), neigh.end());
        encode(neigh, pairs);

        // layout: old_E[e], degree, pairs of (color,count)
        auto &buf = sign_E[e];
        buf.clear();
        buf.push_back(old_E[e]);
        buf.push_back(neigh.size());
        buf.insert(buf.end(), pairs.begin(), pairs.end());
    }

    // vertices
    for (int v = 0; v < num_V; v++) {
        neigh.clear();
        auto edges = graph.clausesOfLiteral(v);

        for (auto e_id : edges) {
            neigh.push_back(old_E[e_id]);
        }
        sort(neigh.begin(), neigh.end());
        encode(neigh, pairs);

        // layout: old_V[v], old_V[mate(v)], degree, pairs of (color,count)
        auto &buf = sign_V[v];
        buf.clear();
        buf.push_back(old_V[v]);
        buf.push_back(old_V[graph.mateOf(v)]);
        buf.push_back(neigh.size());
        buf.insert(buf.end(), pairs.begin(), pairs.end());
    }

    new_E.resize(num_E);
    new_V.resize(num_V);

    //raw
    if (!settings.canonicalize_color_classes) {
        for (int e = 0; e < num_E; ++e) new_E[e] = hash_words(sign_E[e]);
        for (int v = 0; v < num_V; ++v) new_V[v] = hash_words(sign_V[v]);
        return;
    }

    //canonical
    // edges
    if (num_E > 0) {
        vector<int> ordE(num_E); iota(ordE.begin(), ordE.end(), 0);
        sort(ordE.begin(), ordE.end(), [&](int a, int b){ return sign_E[a] < sign_E[b]; });
        Color color_id_e = 0;
        new_E[ordE[0]] = color_id_e;
        for (int i = 1; i < num_E; ++i) {
            if (sign_E[ordE[i]] != sign_E[ordE[i-1]]) color_id_e++;
            new_E[ordE[i]] = color_id_e;
        }
    }

    if (num_V > 0) {
        // vertices
        vector<int> ordV(num_V); iota(ordV.begin(), ordV.end(), 0);
        sort(ordV.begin(), ordV.end(), [&](int a, int b){ return sign_V[a] < sign_V[b]; });
        Color color_id_v = 0;
        new_V[ordV[0]] = color_id_v;
        for (int i = 1; i < num_V; ++i) {
            if (sign_V[ordV[i]] != sign_V[ordV[i-1]]) color_id_v++;
            new_V[ordV[i]] = color_id_v;
        }
    }
}

inline void initialize(const CNF::IncidenceHypergraph &graph, vector<Color> &old_V, vector<Color> &old_E) {
    // edges
    for (int e = 0; e < graph.nClauses(); e++) {
        old_E[e]=graph.clauseSize(e);
    }

    // vertices
    for (int v = 0; v < graph.nVertices(); v++) {
        old_V[v]=graph.degree(v);
    }
}

// main function
WLResult run(const CNF::IncidenceHypergraph &graph, const WLSettings &settings) {
    WLResult out;
    std::vector<Color> old_V(graph.nVertices()), old_E(graph.nClauses()), new_V(graph.nVertices()), new_E(graph.nClauses());

    initialize(graph, old_V, old_E);

    while (out.round < settings.max_iterations && !out.stabilized) {
        // refinement
        update_colors(graph, settings, old_V, old_E, new_V, new_E);
        out.round++;

        // per iteration stats
        out.vertex_color_classes = get_color_classes(new_V);
        out.edge_color_classes = get_color_classes(new_E);
        if (settings.print_stats) {
            std::cerr << "round " << (out.round)
                      << "  V_classes=" << out.vertex_color_classes
                      << "  E_classes=" << out.edge_color_classes << "\n";
        }

        // stabilization
        bool same_V=false, same_E=false;
        if (!settings.canonicalize_color_classes) {
            same_V = same_partition(old_V, new_V); 
            same_E = same_partition(old_E, new_E);
        } else {
            same_V = (old_V == new_V);
            same_E = (old_E == new_E);
        }
        if (same_V && same_E) { 
            out.hash = digest(new_V, new_E);
            out.stabilized=true;
            return out;
        }

        old_V.swap(new_V); old_E.swap(new_E);
    }
    out.hash = digest(old_V, old_E);
    return out;
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
