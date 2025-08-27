#pragma once
#include "bits/stdc++.h"
using namespace std;

#include "src/external/xxhash/xxhash.h"
#include "src/util/CNFFormula.h"
#include "src/transform/cnf2hypergraph.h"

namespace WL {

struct WLResult {
    vector<uint64_t> vertex_labels;
    vector<uint64_t> edge_labels;
    int rounds = 0;
    bool stable = false;
    uint64_t hash = 0;
};

struct WLSettings {
    uint64_t max_iterations=10;
    bool remove_duplicate_literals = true;
    bool early_stopping = false;
    bool split_polarity = true;
    bool edge_labels = true;
    bool node_labels = true;
    bool print_stats = false;
};

inline static vector<uint64_t> canonicalize(const vector<vector<uint64_t>>& sigs) {
    const int n = (int)sigs.size();
    if (n == 0) return {};
    vector<int> perm(n);
    iota(perm.begin(), perm.end(), 0);

    // Deterministic order of signatures
    sort(perm.begin(), perm.end(),
         [&](int a, int b){ return sigs[a] < sigs[b]; });

    vector<uint64_t> ids(n);
    uint64_t cur = 0;
    ids[perm[0]] = 0;
    for (int i = 1; i < n; ++i) {
        if (sigs[perm[i]] != sigs[perm[i - 1]]) ++cur;
        ids[perm[i]] = cur;
    }
    return ids;
}

// counting
inline static uint64_t nonzero_count(const vector<uint64_t>& v) {
    return uint64_t(std::count_if(v.begin(), v.end(),
                                  [](uint64_t c){ return c != 0; }));
}


// hashing array of 64 bit ints with xxh3
inline uint64_t hashing(const vector<uint64_t>& data) {
    return data.empty()
         ? XXH3_64bits(nullptr, 0)
         : XXH3_64bits(data.data(), data.size() * sizeof(uint64_t));
}

inline WLResult refine_1WL(const CNF::cnf2hypergraph& G, WLSettings settings) {
    WLResult out;
    out.vertex_labels.resize(G.nverts);
    out.edge_labels.resize(G.nedges);

    // 0) Initial Coloring
    vector<uint64_t> new_edge_labels(G.nedges);
    vector<uint64_t> new_vertex_labels(G.nverts);

    // vertices: [polarity, degree]
    {
        vector<vector<uint64_t>> v_sig(G.nverts);
        for (int v = 0; v < G.nverts; ++v) {
            const uint64_t degree   = uint64_t(G.inc_ofs[v + 1] - G.inc_ofs[v]);
            const uint64_t polarity = settings.split_polarity ? uint64_t(G.vertex_polarity[v]) : 0;
            v_sig[v] = { polarity, degree };
        }
        new_vertex_labels = canonicalize(v_sig);
    }

    // edges: [kind, size]
    {
        vector<vector<uint64_t>> e_sig(G.nedges);
        for (int e = 0; e < G.nedges; ++e) {
            const uint64_t size = uint64_t(G.edge_verts_ofs[e + 1] - G.edge_verts_ofs[e]);
            const uint64_t kind = uint64_t(int(G.edge_kind[e]));
            e_sig[e] = { kind, size };
        }
        new_edge_labels = canonicalize(e_sig);
    }

    // assigning new labels
    out.vertex_labels = new_vertex_labels;
    out.edge_labels = new_edge_labels;
    vector<array<uint64_t,3>> round_fingerprint;
    round_fingerprint.reserve(settings.max_iterations);
    vector<uint64_t> final_vertex_count, final_clause_count, final_varpair_count;

    // 1) Refinement loop
    for (int round = 1; round <= settings.max_iterations; ++round) {
        // a) Recompute each edge label from incident vertex labels
        vector<vector<uint64_t>> e_sig(G.nedges);
        for (int e = 0; e < G.nedges; ++e) {
            const int begin = G.edge_verts_ofs[e];
            const int end   = G.edge_verts_ofs[e + 1];

            // Collect and sort incident vertex labels
            vector<uint64_t> neighbor_labels;
            neighbor_labels.reserve(size_t(end - begin));
            for (int i = begin; i < end; ++i) {
                const int vtx = G.edge_verts[i];
                neighbor_labels.push_back(out.vertex_labels[vtx]);
            }
            sort(neighbor_labels.begin(), neighbor_labels.end());

            // Build signature = [kind, size, sorted neighbor labels]
            e_sig[e].reserve(2 + neighbor_labels.size());
            e_sig[e].push_back(uint64_t(int(G.edge_kind[e])));
            e_sig[e].push_back(uint64_t(end - begin));
            e_sig[e].insert(e_sig[e].end(), neighbor_labels.begin(), neighbor_labels.end());
        }
        // assign new canonicalized edge
        new_edge_labels = canonicalize(e_sig);

        // b) Recompute each vertex label from incident edge labels
        vector<vector<uint64_t>> v_sig(G.nverts);
        for (int v = 0; v < G.nverts; ++v) {
            const int begin = G.inc_ofs[v];
            const int end   = G.inc_ofs[v + 1];

            // Collect and sort incident edge labels
            vector<uint64_t> neighbor_labels;
            neighbor_labels.reserve(size_t(end - begin));
            for (int i = begin; i < end; ++i) {
                const int edge = G.inc_edges[i];
                neighbor_labels.push_back(out.edge_labels[edge]);
            }
            sort(neighbor_labels.begin(), neighbor_labels.end());

            const uint64_t degree   = uint64_t(end - begin);
            const uint64_t polarity = settings.split_polarity ? uint64_t(G.vertex_polarity[v]) : 0;


            // Build signature = [polarity, degree, sorted neighbor labels]
            v_sig[v].reserve(2 + neighbor_labels.size());
            v_sig[v].push_back(polarity);
            v_sig[v].push_back(degree);
            v_sig[v].insert(v_sig[v].end(), neighbor_labels.begin(), neighbor_labels.end());
        }
        new_vertex_labels = canonicalize(v_sig);

        // creating histogramm
        // Count vertex labels
        uint64_t vertex_classes=0;
        for (auto x : new_vertex_labels) {
            if (x + 1 > vertex_classes) {
                vertex_classes = x + 1;
            }
        }
        vector<uint64_t> vertex_count(vertex_classes, 0);
        for (auto x : new_vertex_labels) ++vertex_count[x];

        // Count edge labels (seperately for both kinds)
        uint64_t edge_classes = 0;
        for (auto x : new_edge_labels) {
            if (x + 1 > edge_classes) { 
                edge_classes = x + 1;
            }
        }
        vector<uint64_t> clause_count(edge_classes, 0), varpair_count(edge_classes, 0);                   // C=Clause, P=VarPair
        for (int e = 0; e < G.nedges; ++e) {
            const uint64_t id = new_edge_labels[e];
            if (G.edge_kind[e] == CNF::cnf2hypergraph::Clause) {
                ++clause_count[id];
            }
            else {
                ++varpair_count[id];
            }
        }

        final_vertex_count=vertex_count;
        final_clause_count=clause_count;
        final_varpair_count=varpair_count;

        // put the counts per color class into this rounds fingerprint
        uint64_t V_colors = nonzero_count(vertex_count);
        uint64_t C_colors = nonzero_count(clause_count);    
        uint64_t P_colors = nonzero_count(varpair_count);
        round_fingerprint.push_back({V_colors, C_colors, P_colors});

        // --print-stats
        if (settings.print_stats) {
            cerr << "c iteration " << round
                << "  V_colors="  << V_colors
                << "  E_colors="  << (C_colors + P_colors)
                << "  (Clause="   << C_colors
                << ", VarPair="   << P_colors << ")\n";
        }

        // stabilization check
        const bool edges_unchanged   = (new_edge_labels   == out.edge_labels);
        const bool vertices_unchanged= (new_vertex_labels == out.vertex_labels);

        // Update state and check convergence
        out.rounds = round;
        out.edge_labels.swap(new_edge_labels);
        out.vertex_labels.swap(new_vertex_labels);

        if (edges_unchanged && vertices_unchanged) {
            out.stable = true;
            std::cerr << "stabilized after " << out.rounds << " iterations\n";
            std::cout << "stabilized\n";
            break;
        }
    }

    // 2) Final hash
    vector<uint64_t> buffer;
    buffer.reserve(16 + 2*final_vertex_count.size() + 4*max(final_clause_count.size(),final_varpair_count.size()) + 3*round_fingerprint.size());

    // helper to append nonzero (id,count) pairs
    auto push_pairs = [&](const vector<uint64_t>& cnt) {
        uint64_t nz = 0; for (auto c : cnt) if (c) ++nz;
        buffer.push_back(nz);
        for (uint64_t id = 0; id < cnt.size(); ++id) if (cnt[id]) {
            buffer.push_back(id);
            buffer.push_back(cnt[id]);
        }
    };

    // vertices, clause-edges, varpair-edges
    push_pairs(final_vertex_count);
    push_pairs(final_clause_count);
    push_pairs(final_varpair_count);

    // add all round fingerprints
    buffer.push_back(uint64_t(round_fingerprint.size()));
    for (const auto& fp : round_fingerprint) {
        buffer.push_back(fp[0]); buffer.push_back(fp[1]); buffer.push_back(fp[2]);
    }
    out.hash = hashing(buffer);

    return out;
}

inline std::string to_hex64(uint64_t x) {
    std::ostringstream oss;
    oss << std::hex << std::nouppercase << std::setfill('0') << std::setw(16)
        << static_cast<unsigned long long>(x);
    return oss.str();
}

inline std::string wlhyphash(const char* filename, WLSettings settings) {
    CNFFormula f(filename);
    std::cerr << "c parsed vars=" << f.nVars() << " clauses=" << f.nClauses() << "\n";
    if (f.nClauses() == 0 && f.nVars() == 0) {
        std::cerr << "c ERROR: nothing parsed\n";
    }
    f.normalizeVariableNames();
    CNF::cnf2hypergraph g(f);
    WLResult res = refine_1WL(g, settings);
    return to_hex64(res.hash);
}



} // namespace WL
