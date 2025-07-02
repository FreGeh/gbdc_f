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

#define XXH_INLINE_ALL
#include "src/external/xxhash/xxhash.h"

#include "src/util/StreamBuffer.h"
#include "src/util/SolverTypes.h"
#include "src/util/CNFFormula.h"
#include <map>
#include <cstdint>


namespace CNF {
    // one node for each literal and each clause
    struct Node {
        enum class Type { Literal, Clause } type;
        std::uint64_t index;
        bool polarity; // true if positive (only relevant for Literals, false for Clauses)
        
        // --- optional node features for round 0 ---
        std::uint64_t sizeMetric = 0; // degree for literals, length for clauses

    };

    // bipartite graph
    struct Graph {
        std::vector<Node> nodes;
        std::vector<std::vector<std::uint64_t>> adj;  // adj[u] = vector with neighbors of u

        // turn CNFFormula into bipartite Graph
        void toGraph(const CNFFormula& cnf) {
            nodes.clear();
            adj.clear();

            // --- NODES ---
            const std::uint64_t literalCount = 2 * cnf.nVars();
            nodes.reserve(literalCount + cnf.nClauses());

            // literals
            for (std::uint64_t var = 0; var < cnf.nVars(); ++var) {
                nodes.push_back({ Node::Type::Literal, var, false, 0});
                nodes.push_back({ Node::Type::Literal, var, true, 0});
            }

            // clauses
            for (std::uint64_t i = 0; i < cnf.nClauses(); ++i) {
                nodes.push_back({ Node::Type::Clause, i, false, 0});
            }

            // --- ADJ LIST ---
            adj.resize(nodes.size());

            for (std::uint64_t i = 0; i < cnf.nClauses(); ++i) {
                const Cl* clause = cnf[i];
                std::uint64_t clauseNode = literalCount + i;
                nodes[clauseNode].sizeMetric = clause->size();

                // REMOVE DUPLICATE LITERALS (make optional later on for testing)
                std::vector<Lit> tmp(clause->begin(), clause->end());
                std::sort(tmp.begin(), tmp.end(),[](Lit a, Lit b){ return a.x < b.x; });
                tmp.erase(std::unique(tmp.begin(), tmp.end(),[](Lit a, Lit b){ return a.x == b.x; }), tmp.end());

                for (const Lit& lit : tmp) {
                    std::uint64_t litNode = 2 * lit.var() + (lit.sign() ? 1 : 0);
                    if (litNode >= adj.size() || clauseNode >= adj.size()) {
                        std::cerr << "ERROR: litNode or clauseNode out of bounds: litNode=" << litNode << ", clauseNode=" << clauseNode << ", adj.size()=" << adj.size() << std::endl;
                    }
                    adj[litNode].push_back(clauseNode);
                    adj[clauseNode].push_back(litNode);

                    ++nodes[litNode].sizeMetric;
                }
            }
        };
    };

    enum class HashMode {
        XXH3_64,
        XXH3_128
    };

    static constexpr HashMode DefaultHashMode = HashMode::XXH3_64;

    struct Hash {
        static constexpr HashMode mode = DefaultHashMode;
        
        static std::uint64_t hash(const void* data, std::uint64_t len) {
            if constexpr (mode == HashMode::XXH3_64) {
                return static_cast<std::uint64_t>( XXH3_64bits(data,len) );
            }
        }

        static std::string toHex(std::uint64_t v) {
            char buf[17];
            std::snprintf(buf, sizeof(buf), "%016llx",
                            static_cast<unsigned long long>(v));
            return std::string(buf);
        }
    };

    inline std::uint64_t hashNewColor(std::uint64_t old_color, const std::vector<std::uint64_t> &neighbours) {
        std::vector<std::uint64_t> buf;
        buf.reserve(neighbours.size()+1);
        buf.push_back(old_color);
        buf.insert(buf.end(), neighbours.begin(), neighbours.end());
        return Hash::hash(buf.data(), buf.size()*sizeof(buf[0]));
    }

    inline std::uint64_t hashOutput(const std::vector<std::uint64_t> &colored_nodes) {
        return Hash::hash(colored_nodes.data(),colored_nodes.size()*sizeof(colored_nodes[0]));
    }

    inline std::uint64_t wlhash(const char* filename) {
        CNFFormula cnf = CNFFormula(filename);
        cnf.normalizeVariableNames();

        Graph g;
        g.toGraph(cnf);

        std::vector<uint64_t> color_old(g.nodes.size(), 0);
        std::vector<uint64_t> color_new(g.nodes.size(), 0);
        std::vector<uint64_t> raw(g.nodes.size());

        // round 0
        for (std::uint64_t u = 0; u < g.nodes.size(); ++u) {
            color_new[u] = g.nodes[u].sizeMetric; // initial coloring
        }

        std::vector<uint64_t> buffer;
        // loop
        int iter = 0;
        while (true) {
            // 1. loop condition
            if (color_new == color_old) {
                break;
            }
            color_old.swap(color_new);

            // 2. compute raw hashes
            std::vector<uint64_t> raw(g.nodes.size());
            for (std::uint64_t u = 0; u < g.nodes.size(); ++u) {
                buffer.clear();
                for (uint64_t w : g.adj[u]) {
                    buffer.push_back(color_old[w]);
                }
                std::sort(buffer.begin(), buffer.end());
                raw[u] = hashNewColor(color_old[u], buffer);
            }

            // 3. canonise hashes into good color ids
            std::map<uint64_t,uint64_t> canon;
            uint64_t next = 0;
            for (uint64_t u = 0; u < g.nodes.size(); ++u) {
                auto [it,ins] = canon.emplace(raw[u], next);
                if (ins) ++next;
                color_new[u] = it->second;
            }
        }

        // output
        return hashOutput(color_new);
    }
} // namespace CNF

#endif // WLISOHash_H_