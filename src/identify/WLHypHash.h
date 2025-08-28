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
    bool split_polarity = true;
    bool print_stats = false;
    bool final_multiset = true; // sorted multisets of final vertex labels and final edge labels
    bool color_inc_matrix = true; // incidence counts between vertex color classes and edge color classes, M[i, j] = number of incidences where a vertex of class i touches an edge of class j
    bool detailed_histogram = false; // for each vertex color i, histogram the degrees of vertices in that class. For edges, histogram the sizes of edges in that class
    bool per_round_multiset = false; // at each round, hash the sorted multiset of vertex labels and of edge labels, append these per-round digests
    bool whole_final_hash = false; // serialize the colored quotient as rows of canonicalized neighbor multisets per class
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
    return uint64_t(std::count_if(v.begin(), v.end(),[](uint64_t c){ return c != 0; }));
}

// hashing array of 64 bit ints with xxh3
inline uint64_t hashing(const vector<uint64_t>& data) {
    return data.empty()
         ? XXH3_64bits(nullptr, 0)
         : XXH3_64bits(data.data(), data.size() * sizeof(uint64_t));
}

// final hash buffer functions
// Push nonzero (id, count) pairs of a histogram, deterministically ordered by id.
inline static void push_pairs(vector<uint64_t>& buffer, const vector<uint64_t>& cnt) {
    uint64_t nz = 0; for (auto c : cnt) if (c) ++nz;
    buffer.push_back(nz);
    for (uint64_t id = 0; id < cnt.size(); ++id) if (cnt[id]) {
        buffer.push_back(id);
        buffer.push_back(cnt[id]);
    }
}

// Push a sorted multiset of labels: size followed by elements in ascending order.
inline static void push_sorted_multiset(vector<uint64_t>& buffer, const vector<uint64_t>& labels) {
    vector<uint64_t> tmp = labels;
    sort(tmp.begin(), tmp.end());
    buffer.push_back(uint64_t(tmp.size()));
    buffer.insert(buffer.end(), tmp.begin(), tmp.end());
}

// Push a sparse matrix encoded as (rows, cols, nnz, r1,c1,val1, r2,c2,val2, ...),
// iterating rows then cols in increasing order to ensure determinism.
inline static void push_sparse_matrix(
    vector<uint64_t>& buffer, const vector<uint64_t>& mat, uint64_t rows, uint64_t cols)
{
    uint64_t nnz = 0; for (auto x : mat) if (x) ++nnz;
    buffer.push_back(rows);
    buffer.push_back(cols);
    buffer.push_back(nnz);
    for (uint64_t r = 0; r < rows; ++r) {
        const uint64_t base = r * cols;
        for (uint64_t c = 0; c < cols; ++c) {
            uint64_t v = mat[base + c];
            if (v) { buffer.push_back(r); buffer.push_back(c); buffer.push_back(v); }
        }
    }
}

// Append the two sorted multisets of final labels (vertices, edges).
inline static void append_final_label_multisets(
    vector<uint64_t>& buffer, const vector<uint64_t>& vertex_labels, const vector<uint64_t>& edge_labels)
{
    push_sorted_multiset(buffer, vertex_labels);
    push_sorted_multiset(buffer, edge_labels);
}

// Build and append the colored quotient incidence matrix M[i,j] over final colors.
inline static void append_color_incidence_matrix(vector<uint64_t>& buffer, 
    const CNF::cnf2hypergraph& G, const vector<uint64_t>& vertex_labels, const vector<uint64_t>& edge_labels)
{
    uint64_t VC = 0, EC = 0;
    if (!vertex_labels.empty())
        VC = 1 + *max_element(vertex_labels.begin(), vertex_labels.end());
    if (!edge_labels.empty())
        EC = 1 + *max_element(edge_labels.begin(), edge_labels.end());

    // Edge case: empty graph
    if (VC == 0 || EC == 0) {
        push_sparse_matrix(buffer, /*mat*/{}, /*rows*/VC, /*cols*/EC);
        return;
    }

    vector<uint64_t> M(VC * EC, 0);

    // Count incidences between classes
    for (int e = 0; e < G.nedges; ++e) {
        const uint64_t j = edge_labels[e];
        const int b = G.edge_verts_ofs[e], t = G.edge_verts_ofs[e + 1];
        for (int k = b; k < t; ++k) {
            const int v = G.edge_verts[k];
            const uint64_t i = vertex_labels[v];
            M[i * EC + j] += 1;
        }
    }

    push_sparse_matrix(buffer, M, VC, EC);
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
            break;
        }
    }

    // 2) Final hash
    vector<uint64_t> buffer;

    // vertices, clause-edges, varpair-edges
    push_pairs(buffer, final_vertex_count);
    push_pairs(buffer, final_clause_count);
    push_pairs(buffer, final_varpair_count);

    if (settings.final_multiset) {
        append_final_label_multisets(buffer, out.vertex_labels, out.edge_labels);
    }
    if (settings.color_inc_matrix) {
        append_color_incidence_matrix(buffer, G, out.vertex_labels, out.edge_labels);
    }

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
