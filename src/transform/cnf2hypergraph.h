#pragma once
#include "bits/stdc++.h"
using namespace std;
#include "src/util/CNFFormula.h"

namespace CNF {

class cnf2hypergraph {
public:
    enum EdgeKind : uint8_t { Clause = 0, VariablePair = 1 };

    // sizes
    int nvars = 0;         // number of CNF variables
    int nverts = 0;        // = 2 * nvars, literal vertices
    int nedges = 0;        // = nvars + nClauses, var-pair edges first, then clauses

    // edge metadata, size = nedges
    vector<EdgeKind> edge_kind;      // VariablePair or Clause
    vector<int>      edge_clause_id; // -1 for var-pair edges, clause id for clause edges

    // CSR: for each edge, the list of its vertices
    vector<int> edge_verts_ofs; // size = nedges + 1
    vector<int> edge_verts;     // flat storage of vertex ids

    // CSR: for each vertex, the list of incident edges
    vector<int> inc_ofs;        // size = nverts + 1
    vector<int> inc_edges;      // flat storage of edge ids

    cnf2hypergraph() = default;
    explicit cnf2hypergraph(const CNFFormula& F) { build(F); }

    static inline int literalVertexId(int var0, bool positive) {
        // var0 is 0-based
        return (var0 << 1) + (positive ? 0 : 1);
    }
    static inline pair<int,bool> decodeLiteralVertexId(int vid) {
        return { vid >> 1, (vid & 1) == 0 };
    }

    inline int numVertices()  const { return nverts; }
    inline int numHyperedges() const { return nedges; }

    // access helpers, return [begin, end) pointer ranges
    inline pair<const int*, const int*> edgeVertices(int e) const {
        int s = edge_verts_ofs[e], t = edge_verts_ofs[e + 1];
        return { &edge_verts[s], &edge_verts[t] };
    }
    inline pair<const int*, const int*> vertexIncidentEdges(int v) const {
        int s = inc_ofs[v], t = inc_ofs[v + 1];
        return { &inc_edges[s], &inc_edges[t] };
    }

    void clear() {
        nvars = nverts = nedges = 0;
        edge_kind.clear();
        edge_clause_id.clear();
        edge_verts_ofs.clear();
        edge_verts.clear();
        inc_ofs.clear();
        inc_edges.clear();
    }

    void build(const CNFFormula& F) {
        clear();
        nvars  = static_cast<int>(F.nVars());
        nverts = 2 * nvars;
        const int m = static_cast<int>(F.nClauses());
        nedges = nvars + m; // var-pair edges [0..nvars-1], clause edges [nvars..nvars+m-1]

        edge_kind.resize(nedges);
        edge_clause_id.assign(nedges, -1);
        edge_verts_ofs.resize(nedges + 1);

        // First pass, sizes and degrees
        vector<int> edge_sizes(nedges, 0);
        vector<int> deg(nverts, 1); // each literal participates in exactly one var-pair edge

        // var-pair edges have size 2
        for (int v0 = 0; v0 < nvars; ++v0) {
            edge_sizes[v0] = 2;
        }

        // clause edges sizes, and count literal occurrences for vertex degrees
        for (int cid = 0; cid < m; ++cid) {
            const Cl* clause = F[cid];
            const int eid = nvars + cid;
            edge_sizes[eid] = static_cast<int>(clause->size());
            for (const Lit& lit : *clause) {
                const int var0 = static_cast<int>(lit.var()); // CNFFormula is 0-based
                const bool pos = !lit.sign();
                const int vid  = literalVertexId(var0, pos);
                ++deg[vid];
            }
        }

        // CSR offsets for edges -> vertices
        edge_verts_ofs[0] = 0;
        for (int e = 0; e < nedges; ++e) edge_verts_ofs[e + 1] = edge_verts_ofs[e] + edge_sizes[e];
        edge_verts.resize(edge_verts_ofs[nedges]);

        // CSR offsets for vertices -> incident edges
        inc_ofs.resize(nverts + 1);
        inc_ofs[0] = 0;
        for (int v = 0; v < nverts; ++v) inc_ofs[v + 1] = inc_ofs[v] + deg[v];
        inc_edges.resize(inc_ofs[nverts]);

        // write cursors
        vector<int> cur_edge_pos = edge_verts_ofs; // where to write next vertex for each edge
        vector<int> cur_inc_pos  = inc_ofs;        // where to write next incident edge for each vertex

        // Fill var-pair edges
        for (int var0 = 0; var0 < nvars; ++var0) {
            const int eid = var0;
            edge_kind[eid] = VariablePair;

            const int vp = literalVertexId(var0, true);
            const int vn = literalVertexId(var0, false);

            edge_verts[cur_edge_pos[eid]++] = vp;
            edge_verts[cur_edge_pos[eid]++] = vn;

            inc_edges[cur_inc_pos[vp]++] = eid;
            inc_edges[cur_inc_pos[vn]++] = eid;
        }

        // Fill clause edges
        for (int cid = 0; cid < m; ++cid) {
            const int eid = nvars + cid;
            edge_kind[eid] = Clause;
            edge_clause_id[eid] = cid;

            const Cl* clause = F[cid];
            for (const Lit& lit : *clause) {
                const int var0 = static_cast<int>(lit.var()); // 1-based to 0-based
                const bool pos = !lit.sign();
                const int vid  = literalVertexId(var0, pos);

                edge_verts[cur_edge_pos[eid]++] = vid;
                inc_edges[cur_inc_pos[vid]++]    = eid;
            }
        }
    }
};

} // namespace CNF
