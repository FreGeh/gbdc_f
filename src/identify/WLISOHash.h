/*************************************************************************************************
CNFTools -- Copyright (c) 2021, Markus Iser, KIT - Karlsruhe Institute of Technology
WLISOHash -- Copyright (c) 2025, Frederick Gehm, KIT - Karlsruhe Institute of Technology

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 **************************************************************************************************/

#ifndef WLISOHash_H_
#define WLISOHash_H_

#include <vector>
#include <algorithm>
#include <stdio.h>

#include "src/external/md5/md5.h"
#include "src/external/xxhash/xxhash.h"


#include "src/util/StreamBuffer.h"
#include "src/util/SolverTypes.h"
#include "src/util/CNFFormula.h"


namespace CNF {
    // one node for each literal and each clause
    struct Node {
        enum class Type { Literal, Clause } type;
        std::size_t index;
        bool polarity; // true if positive (only relevant for Literals, false for Clauses)
        
        // --- optional node features for round 0 ---
        std::size_t sizeMetric = 0; // degree for literals, length for clauses

    };

    // bipartite graph
    struct Graph {
        std::vector<Node> nodes;
        std::vector<std::vector<std::size_t>> adj;  // adj[u] = vector with neighbors of u

        // turn CNFFormula into bipartite Graph
        void toGraph(const CNFFormula& cnf) {
            nodes.clear();
            adj.clear();

            // --- NODES ---
            std::size_t literalCount = 2 * cnf.nVars();
            nodes.reserve(literalCount + cnf.nClauses());

            // literals
            for (std::size_t var = 1; var <= cnf.nVars(); ++var) {
                nodes.push_back({ Node::Type::Literal, var, false, 0});
                nodes.push_back({ Node::Type::Literal, var, true, 0});
            }

            // clauses
            std::size_t nodeOffset = nodes.size();
            for (std::size_t i = 0; i < cnf.nClauses(); ++i) {
                nodes.push_back({ Node::Type::Clause, i, false, 0});
            }

            // --- ADJ LIST ---
            adj.resize(nodes.size());

            for (std::size_t i = 0; i < cnf.nClauses(); ++i) {
                const Cl* clause = cnf[i];
                std::size_t clauseNode = nodeOffset + i;

                nodes[clauseNode].sizeMetric = clause->size();

                for (const Lit& lit : *clause) {
                    std::size_t litNode = 2 * (lit.var() - 1) + (lit.sign() ? 1 : 0);
                    adj[litNode].push_back(clauseNode);
                    adj[clauseNode].push_back(litNode);

                    ++nodes[litNode].sizeMetric;
                }
            }
        };
    };

    bool loopCondition(std::vector<size_t> &color_old, std::vector<size_t> &color_new, std::size_t node_size) {
        for (std::size_t u = 0; u < node_size; ++u) {
            if (color_old[u] != color_new[u]) {
                return false;
            }
            return true;
        }
    }

    std::size_t hashNewColor(std::size_t &old_color, std::vector<std::size_t> &neighbours_colors) {
        size_t new_color;
        for (std::size_t n : neighbours_colors) {
            new_color = new_color + n;
        }
        return new_color;
    }

    std::size_t hashOutput(std::vector<std::size_t> &colored_nodes) {
        size_t outputHash;
        for (std::size_t n : colored_nodes) {
            outputHash = outputHash + n;
        }
        return outputHash;
    }

    std::size_t wlhash(const char* filename) {
        CNFFormula cnf = CNFFormula(filename);
        cnf.normalizeVariableNames();

        Graph g;
        g.toGraph(cnf);

        std::vector<size_t> color_old(g.nodes.size(), 0);
        std::vector<size_t> color_new(g.nodes.size(), 0);

        // round 0
        for (std::size_t u = 0; u < g.nodes.size(); ++u) {
            const Node& n = g.nodes[u];
            color_new[u] = n.sizeMetric; // initial coloring
        }

        std::vector<size_t> buffer;
        // loop
        while (true) {
            // 1. compute raw hashes
            std::vector<size_t> raw(g.nodes.size());
            for (std::size_t u = 0; u < g.nodes.size(); ++u) {
                buffer.clear();
                for (size_t w : g.adj[u]) {
                    buffer.push_back(color_old[w]);
                }
                std::sort(buffer.begin(), buffer.end());
                std::size_t h = hashNewColor(color_old[u], buffer);
                raw[u] = h;
            }

            // 2. canonise hashes into good color ids
            std::unordered_map<size_t,size_t> canon;
            size_t next = 0;
            for (size_t u = 0; u < g.nodes.size(); ++u) {
                auto [it,ins] = canon.emplace(raw[u], next);
                if (ins) ++next;
                color_new[u] = it->second;
            }

            // 3. loop condition
            if (color_new == color_old) break;
            color_old.swap(color_new);
        }

        // output
        return hashOutput(color_new);
    }
} // namespace CNF

#endif  // WLISOHash_H_
