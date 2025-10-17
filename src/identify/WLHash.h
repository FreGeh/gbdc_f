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
    bool neigh_hop = false;
    bool advanced_sort = false;
    bool sum_sort = false;
};

struct WLResult {
    uint64_t hash;
    int round = 0;
    int vertex_color_classes = 0;
    int edge_color_classes = 0;
    bool stabilized = false;
};

// helper functions

inline void radix_sort(vector<uint64_t>& a) {
    const size_t n = a.size();
    if (n <= 1) return;
    if (n < 64) { sort(a.begin(), a.end()); return; }

    static vector<uint64_t> buf;
    if (buf.size() < n) buf.resize(n);

    for (int pass = 0; pass < 8; ++pass) {
        size_t cnt[256] = {0};
        const unsigned shift = static_cast<unsigned>(pass * 8);

        // count
        for (size_t i = 0; i < n; ++i) {
            const unsigned byte = static_cast<unsigned>((a[i] >> shift) & 0xFFull);
            ++cnt[byte];
        }

        // prefix sums -> positions
        size_t pos[256];
        pos[0] = 0;
        for (int i = 1; i < 256; ++i) pos[i] = pos[i - 1] + cnt[i - 1];

        // stable scatter
        for (size_t i = 0; i < n; ++i) {
            const unsigned byte = static_cast<unsigned>((a[i] >> shift) & 0xFFull);
            buf[pos[byte]++] = a[i];
        }
        a.swap(buf);
    }
}

using Color = uint64_t;

inline uint64_t hash_words(const uint64_t* data, size_t n_words) {
    return XXH3_64bits(static_cast<const void*>(data), n_words * sizeof(uint64_t));
}
inline uint64_t hash_words(const vector<uint64_t>& v) {
    return hash_words(v.data(), v.size());
}

inline void encode(const vector<Color>& sorted, vector<uint64_t>& pairs) {
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

inline void sort_neigh(vector<uint64_t> &neigh, const WLSettings& settings) {
    if (settings.advanced_sort) radix_sort(neigh);
    else sort(neigh.begin(), neigh.end());
}

inline void append_neigh(vector<uint64_t> &neigh, const WLSettings &settings, vector<uint64_t> &buf, vector<uint64_t> &pairs) {
    if (settings.sum_sort) {
        uint64_t sum = 0, sum2 = 0, x = 0, mn = UINT64_MAX, mx = 0;
        for (uint64_t c : neigh) {
            sum += c; 
            sum2 += c * c; 
            x ^= c;
            if (c < mn) mn = c; 
            if (c > mx) mx = c;
        }
        if (neigh.empty()) { mn = 0; mx = 0; }
        buf.push_back(sum); 
        buf.push_back(sum2); 
        buf.push_back(x); 
        buf.push_back(mn); 
        buf.push_back(mx);
    } else {
        sort_neigh(neigh, settings);
        encode(neigh, pairs);
        buf.insert(buf.end(), pairs.begin(), pairs.end());
    }
}

inline void append_twohop_vertex(const CNF::IncidenceHypergraph& g, int v, const vector<Color>& old_V, vector<uint64_t>& buf) {
    unordered_map<Color, uint64_t> freq;
    auto edges = g.clausesOfLiteral(v);
    for (int e_id : edges) {
        auto lits = g.literalsOfClause(e_id);
        for (int u : lits) {
            if (u != v) ++freq[ old_V[u] ];
        }
    }
    if (freq.empty()) { 
        buf.push_back(0); 
        return; 
    }
    vector<pair<Color,uint64_t>> tmp; tmp.reserve(freq.size());
    for (auto& kv : freq) tmp.emplace_back(kv.first, kv.second);
    sort(tmp.begin(), tmp.end(), [](auto& a, auto& b){ return a.first < b.first; });
    buf.push_back(static_cast<uint64_t>(tmp.size()));
    for (auto& p : tmp) { 
        buf.push_back(p.first); 
        buf.push_back(p.second); 
    }
}

inline void append_twohop_edge(const CNF::IncidenceHypergraph& g, int e, const vector<Color>& old_E, vector<uint64_t>& buf) {
    unordered_map<Color, uint64_t> freq;
    auto verts = g.literalsOfClause(e);
    for (int v_id : verts) {
        auto edges2 = g.clausesOfLiteral(v_id);
        for (int e2 : edges2) {
            if (e2 != e) ++freq[ old_E[e2] ];
        }
    }
    if (freq.empty()) { 
        buf.push_back(0); 
        return; 
    }
    vector<pair<Color,uint64_t>> tmp; tmp.reserve(freq.size());
    for (auto& kv : freq) tmp.emplace_back(kv.first, kv.second);
    sort(tmp.begin(), tmp.end(), [](auto& a, auto& b){ return a.first < b.first; });
    buf.push_back(static_cast<uint64_t>(tmp.size()));
    for (auto& p : tmp) { 
        buf.push_back(p.first); 
        buf.push_back(p.second); 
    }
}

// important main methods

// initializing colors
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

// iteration step of assigning new colors
inline void update_colors(const CNF::IncidenceHypergraph &graph, const WLSettings &settings, vector<Color> &old_V, vector<Color> &old_E, vector<Color> &new_V, vector<Color> &new_E) {
    vector<Color> neigh;
    vector<uint64_t> pairs;
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

        // layout: old_E[e], degree, pairs of (color,count)
        auto &buf = sign_E[e];
        buf.clear();
        buf.push_back(old_E[e]);
        buf.push_back(neigh.size());
        append_neigh(neigh, settings, buf, pairs);
        if (settings.neigh_hop) {
            append_twohop_edge(graph, e, old_E, buf);
        }
    }

    // vertices
    for (int v = 0; v < num_V; v++) {
        neigh.clear();
        auto edges = graph.clausesOfLiteral(v);

        for (auto e_id : edges) {
            neigh.push_back(old_E[e_id]);
        }

        // layout: old_V[v], old_V[mate(v)], degree, pairs of (color,count)
        auto &buf = sign_V[v];
        buf.clear();
        buf.push_back(old_V[v]);
        buf.push_back(old_V[graph.mateOf(v)]);
        buf.push_back(neigh.size());
        append_neigh(neigh, settings, buf, pairs);
        if (settings.neigh_hop) {
            append_twohop_vertex(graph, v, old_V, buf);
        }
    }

    new_E.resize(num_E);
    new_V.resize(num_V);

    //raw
    for (int e = 0; e < num_E; ++e) new_E[e] = hash_words(sign_E[e]);
    for (int v = 0; v < num_V; ++v) new_V[v] = hash_words(sign_V[v]);
}

// calculating final hash
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

// main function
WLResult run(const CNF::IncidenceHypergraph &graph, const WLSettings &settings) {
    WLResult out;
    vector<Color> old_V(graph.nVertices()), old_E(graph.nClauses()), new_V(graph.nVertices()), new_E(graph.nClauses());

    initialize(graph, old_V, old_E);

    while (out.round < settings.max_iterations && !out.stabilized) {
        // refinement
        update_colors(graph, settings, old_V, old_E, new_V, new_E);
        out.round++;

        // per iteration stats
        out.vertex_color_classes = get_color_classes(new_V);
        out.edge_color_classes = get_color_classes(new_E);
        if (settings.print_stats) {
            cerr << "round " << (out.round)
                 << "  V_classes=" << out.vertex_color_classes
                 << "  E_classes=" << out.edge_color_classes << "\n";
        }

        // stabilization
        bool same_V=false, same_E=false;
        same_V = same_partition(old_V, new_V); 
        same_E = same_partition(old_E, new_E);
        if (same_V && same_E) { 
            out.hash = digest(new_V, new_E);
            out.stabilized=true;
            cerr << "c stabilized: " << out.round << "\n";
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
