#include "bits/stdc++.h"
#include "src/external/xxhash/xxhash.h"
#include "src/util/CNFFormula.h"
#include "src/transform/cnf2hypergraph.h"

using namespace std;

namespace WLF {

struct WLSettings {
    uint64_t max_iterations = 10;
    bool print_stats = false;
};

struct WLResult {
    uint64_t hash;
    int round = 0;
    int node_color_classes = 0;
    int edge_color_classes = 0;
    bool stabilized = false;
};

// helper functions
using Color = uint64_t;

inline uint64_t hash_words(const uint64_t* data, size_t n_words) {
    return XXH3_64bits(static_cast<const void*>(data), n_words * sizeof(uint64_t));
}
inline uint64_t hash_words(const std::vector<uint64_t>& v) {
    return hash_words(v.data(), v.size());
}

inline void encode(const std::vector<Color>& sorted, std::vector<uint64_t>& out_pairs) {
    out_pairs.clear();
    if (sorted.empty()) return;
    Color cur = sorted[0];
    uint64_t cnt = 1;
    for (size_t i = 1; i < sorted.size(); i++) {
        if (sorted[i] == cur) { 
            ++cnt; 
        }
        else { 
            out_pairs.push_back(cur); 
            out_pairs.push_back(cnt); 
            cur = sorted[i]; cnt = 1; 
        }
    }
    out_pairs.push_back(cur); out_pairs.push_back(cnt);
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

//hashed histogram of the resulting labels
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

inline void update_colors(const CNF::IncidenceHypergraph &graph, vector<Color> &old_V, vector<Color> &old_E, vector<Color> &new_V, vector<Color> &new_E) {
    std::vector<Color> neigh;
    std::vector<uint64_t> pairs, buf;

    // edges
    for (int e = 0; e < graph.numClauses(); e++) {
        neigh.clear();
        auto lits = graph.literalsOfClause(e);

        for (int lit : lits) {
            neigh.push_back(old_V[lit]);
        }
        sort(neigh.begin(), neigh.end());
        encode(neigh, pairs);

        buf.clear();

        // layout: [old_E[e], degree, (color,count)*]
        buf.push_back(old_E[e]);
        buf.push_back(neigh.size());
        buf.insert(buf.end(), pairs.begin(), pairs.end());

        new_E[e] = hash_words(buf);
    }

    //nodes
    for (int v = 0; v < graph.numVertices(); v++) {
        neigh.clear();
        auto clauses = graph.clausesOfLiteral(v);

        for (auto cl : clauses) {
            neigh.push_back(old_E[cl]);
        }
        sort(neigh.begin(), neigh.end());
        encode(neigh, pairs);

        buf.clear();
        // layout: [old_V[v], old_V[mate(v)], degree, (color,count)*]
        buf.push_back(old_V[v]);
        buf.push_back(old_V[graph.mateOf(v)]);
        buf.push_back(neigh.size());
        buf.insert(buf.end(), pairs.begin(), pairs.end());

        new_V[v] = hash_words(buf);
    }
}

// initial coloring
inline void initialize(const CNF::IncidenceHypergraph &graph, vector<Color> &old_V, vector<Color> &old_E) {
    // edges
    for (int e = 0; e < graph.numClauses(); e++) {
        old_E[e]=graph.clauseSize(e);
    }

    //nodes
    for (int v = 0; v < graph.numVertices(); v++) {
        old_V[v]=graph.degree(v);
    }
}

// iteration steps
void refine_round(const CNF::IncidenceHypergraph &graph, vector<Color> &old_V, vector<Color> &old_E, vector<Color> &new_V, vector<Color> &new_E) {
    update_colors(graph, old_V, old_E, new_V, new_E);
}

// main function
WLResult run(const CNF::IncidenceHypergraph &graph, const WLSettings &settings) {
    WLResult out;
    std::vector<Color> old_V(graph.numVertices()), old_E(graph.numClauses()), new_V(graph.numVertices()), new_E(graph.numClauses());

    initialize(graph, old_V, old_E);

    while (out.round < settings.max_iterations && !out.stabilized) {
        refine_round(graph, old_V, old_E, new_V, new_E);
        out.round++;

        out.node_color_classes = get_color_classes(new_V);
        out.edge_color_classes = get_color_classes(new_E);
        if (settings.print_stats) {
            std::cerr << "round " << (out.round)
                      << "  V_classes=" << out.node_color_classes
                      << "  E_classes=" << out.edge_color_classes << "\n";
        }

        bool same_V = same_partition(old_V, new_V); 
        bool same_E = same_partition(old_E, new_E);
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