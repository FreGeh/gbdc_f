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

    // buffers reused each iteration
    vector<uint64_t> neighbor_labels;   neighbor_labels.reserve(32); // multiset to sort
    vector<uint64_t> signature_parts;   signature_parts.reserve(32); // features + multiset

    // 0) Initial Coloring
    // vertices (polarity and degree)
    for (int v = 0; v < G.nverts; ++v) {
        const uint64_t degree   = uint64_t(G.inc_ofs[v + 1] - G.inc_ofs[v]);
        const uint64_t polarity = uint64_t(v & 1); // 0 = positive literal, 1 = negative literal

        signature_parts.clear();
        signature_parts.push_back(polarity);
        signature_parts.push_back(degree);

        out.vertex_labels[v] = hashing(signature_parts);
    }

    // edges (kind and size)
    for (int e = 0; e < G.nedges; ++e) {
        const uint64_t size = uint64_t(G.edge_verts_ofs[e + 1] - G.edge_verts_ofs[e]);
        const uint64_t kind = uint64_t(int(G.edge_kind[e])); // 0 = clause, 1 = variable-pair

        signature_parts.clear();
        signature_parts.push_back(kind);
        signature_parts.push_back(size);

        out.edge_labels[e] = hashing(signature_parts);
    }

    // 1) Refinement loop
    vector<uint64_t> new_edge_labels(G.nedges);
    vector<uint64_t> new_vertex_labels(G.nverts);

    for (int round = 1; round <= settings.max_iterations; ++round) {
        bool edges_unchanged   = true;
        bool vertices_unchanged = true;

        // a) Recompute each edge label from incident vertex labels
        for (int e = 0; e < G.nedges; ++e) {
            const int begin = G.edge_verts_ofs[e];
            const int end   = G.edge_verts_ofs[e + 1];

            // Collect and sort incident vertex labels
            neighbor_labels.clear();
            neighbor_labels.reserve(size_t(end - begin));
            for (int i = begin; i < end; ++i) {
                const int vtx = G.edge_verts[i];
                neighbor_labels.push_back(out.vertex_labels[vtx]);
            }
            sort(neighbor_labels.begin(), neighbor_labels.end());

            // Build signature = [kind, size, sorted neighbor labels]
            signature_parts.clear();
            signature_parts.push_back(uint64_t(int(G.edge_kind[e])));
            signature_parts.push_back(uint64_t(end - begin));
            signature_parts.insert(signature_parts.end(), neighbor_labels.begin(), neighbor_labels.end());

            const uint64_t label = hashing(signature_parts);
            edges_unchanged &= (label == out.edge_labels[e]);
            new_edge_labels[e] = label;
        }

        // b) Recompute each vertex label from incident edge labels
        for (int v = 0; v < G.nverts; ++v) {
            const int begin = G.inc_ofs[v];
            const int end   = G.inc_ofs[v + 1];

            // Collect and sort incident edge labels
            neighbor_labels.clear();
            neighbor_labels.reserve(size_t(end - begin));
            for (int i = begin; i < end; ++i) {
                const int edge = G.inc_edges[i];
                neighbor_labels.push_back(new_edge_labels[edge]);
            }
            sort(neighbor_labels.begin(), neighbor_labels.end());

            // Build signature = [polarity, degree, sorted neighbor labels]
            signature_parts.clear();
            signature_parts.push_back(uint64_t(v & 1));
            signature_parts.push_back(uint64_t(end - begin));
            signature_parts.insert(signature_parts.end(), neighbor_labels.begin(), neighbor_labels.end());

            const uint64_t label = hashing(signature_parts);
            vertices_unchanged &= (label == out.vertex_labels[v]);
            new_vertex_labels[v] = label;
        }

        // --print-stats
        if (settings.print_stats) {
            cerr << "c iteration: " << out.round;
            cerr << ", edge labels: " << new_edge_labels.unique();
            cerr << ", clause labels: " << new_vertex_labels.unique() << "\n";
        }

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

    // 2) Final hash: hash sorted stabilized labels plus sizes
    vector<uint64_t> edges = out.edge_labels;
    vector<uint64_t> verts = out.vertex_labels;
    sort(edges.begin(), edges.end());
    sort(verts.begin(), verts.end());

    vector<uint64_t> buffer;
    buffer.reserve(2 + edges.size() + verts.size());
    buffer.push_back(uint64_t(G.nedges));
    buffer.push_back(uint64_t(G.nverts));
    buffer.insert(buffer.end(), edges.begin(), edges.end());
    buffer.insert(buffer.end(), verts.begin(), verts.end());
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
