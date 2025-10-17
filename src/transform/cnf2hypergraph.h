#ifndef SRC_TRANSFORM_CNF2HYPERGRAPH_H_
#define SRC_TRANSFORM_CNF2HYPERGRAPH_H_

#include <vector>
#include <algorithm>
#include "src/util/CNFFormula.h"

using namespace std;

namespace CNF {

class IncidenceHypergraph {
    unsigned vars_ = 0;
    unsigned verts_ = 0;
    unsigned edges_ = 0;
    vector<unsigned> edge_offsets_;
    vector<unsigned> edge_vertices_;
    vector<unsigned> inc_offsets_;
    vector<unsigned> inc_edges_;

public:
    IncidenceHypergraph() = default;

    explicit IncidenceHypergraph(const CNFFormula& f) {
        build(f);
    }

    static inline unsigned literalIndex(unsigned var0, bool positive) {
        return (var0 << 1) | (positive ? 0u : 1u);
    }

    static inline unsigned varOf(unsigned lit_idx) {
        return lit_idx >> 1;
    }

    static inline bool isPositive(unsigned lit_idx) {
        return (lit_idx & 1u) == 0u;
    }

    static inline unsigned mateOf(unsigned lit_idx) {
        return lit_idx ^ 1u;
    }

    inline unsigned nVars() const {
        return vars_;
    }

    inline unsigned nVertices() const {
        return verts_;
    }

    inline unsigned nClauses() const {
        return edges_;
    }

    struct Span {
        const unsigned* b;
        const unsigned* e;

        inline const unsigned* begin() const {
            return b;
        }

        inline const unsigned* end() const {
            return e;
        }

        inline unsigned size() const {
            return (unsigned)(e - b);
        }

        inline bool empty() const {
            return b == e;
        }

        inline unsigned operator[](unsigned i) const {
            return b[i];
        }
    };

    inline Span literalsOfClause(unsigned e) const {
        return { &edge_vertices_[edge_offsets_[e]], &edge_vertices_[edge_offsets_[e + 1]] };
    }

    inline Span clausesOfLiteral(unsigned v) const {
        return { &inc_edges_[inc_offsets_[v]], &inc_edges_[inc_offsets_[v + 1]] };
    }

    inline unsigned degree(unsigned v) const {
        return inc_offsets_[v + 1] - inc_offsets_[v];
    }

    inline unsigned clauseSize(unsigned e) const {
        return edge_offsets_[e + 1] - edge_offsets_[e];
    }

    void build(const CNFFormula& f) {
        clear();

        vars_ = (unsigned)f.nVars();
        verts_ = 2u * vars_;
        edges_ = (unsigned)f.nClauses();

        edge_offsets_.assign(edges_ + 1, 0u);
        inc_offsets_.assign(verts_ + 1, 0u);

        if (edges_ == 0u) {
            return;
        }

        vector<unsigned> edge_size(edges_, 0u);
        vector<unsigned> deg(verts_, 0u);

        for (unsigned cid = 0; cid < edges_; ++cid) {
            const Cl* c = f[cid];
            edge_size[cid] = (unsigned)c->size();
            for (const Lit& lit : *c) {
                const unsigned v0 = (unsigned)lit.var() - 1u;
                const bool pos = !lit.sign();
                ++deg[literalIndex(v0, pos)];
            }
        }

        for (unsigned e = 0; e < edges_; ++e) {
            edge_offsets_[e + 1] = edge_offsets_[e] + edge_size[e];
        }

        edge_vertices_.resize(edge_offsets_[edges_]);

        for (unsigned v = 0; v < verts_; ++v) {
            inc_offsets_[v + 1] = inc_offsets_[v] + deg[v];
        }

        inc_edges_.resize(inc_offsets_[verts_]);

        vector<unsigned> cur_edge = edge_offsets_;
        vector<unsigned> cur_inc = inc_offsets_;
        vector<unsigned> lits;

        for (unsigned cid = 0; cid < edges_; ++cid) {
            const Cl* c = f[cid];

            lits.clear();
            lits.reserve(c->size());

            for (const Lit& lit : *c) {
                const unsigned v0 = (unsigned)lit.var() - 1u;
                const bool pos = !lit.sign();
                lits.push_back(literalIndex(v0, pos));
            }

            sort(lits.begin(), lits.end());

            for (unsigned li : lits) {
                edge_vertices_[cur_edge[cid]++] = li;
                inc_edges_[cur_inc[li]++] = cid;
            }
        }

        for (unsigned v = 0; v < verts_; ++v) {
            sort(&inc_edges_[inc_offsets_[v]], &inc_edges_[inc_offsets_[v + 1]]);
        }
    }

    void clear() {
        vars_ = 0u;
        verts_ = 0u;
        edges_ = 0u;
        edge_offsets_.clear();
        edge_vertices_.clear();
        inc_offsets_.clear();
        inc_edges_.clear();
    }
};

}  // namespace CNF

#endif  // SRC_TRANSFORM_CNF2HYPERGRAPH_H_
