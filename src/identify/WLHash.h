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
#include <cstdio>
#include <limits>
#include <map>
#include <cstdint>

#include "src/external/md5/md5.h"

#define XXH_INLINE_ALL
#include "src/external/xxhash/xxhash.h"

#include "src/util/StreamBuffer.h"
#include "src/util/SolverTypes.h"
#include "src/util/CNFFormula.h"

namespace CNF {

struct Config {
    std::size_t k_dimension = 1;
    std::size_t max_iterations = 100;
    bool print_stats = false;
    bool remove_duplicate_literals = true;
    bool early_stopping = true;
    bool split_nodes = true;
    bool edge_labels = true;
    bool use128 = false;
};

struct Stats {
    std::uint64_t iterations = 0;
    std::uint64_t num_literal_nodes = 0;
    std::uint64_t num_clause_nodes = 0;
    std::uint64_t num_tuples = 0;
    std::vector<std::uint64_t> num_colors_per_iteration;
};

struct Graph {
    struct Node {
        enum class Type { Literal, Clause } type;
        std::uint64_t idx = 0;
        bool pol = false;
        std::uint64_t size = 0;
    };

    enum class Sign : std::uint8_t { Pos = 0, Neg = 1 };

    struct Edge {
        std::uint64_t dst = 0;
        Sign sign = Sign::Pos;
    };

    std::vector<Node> nodes;
    std::vector<std::vector<Edge>> adj;

    void build(const CNFFormula& cnf, const Config& cfg, Stats& stats) {
        nodes.clear();
        adj.clear();

        const std::uint64_t lit_count = cfg.split_nodes ? 2 * cnf.nVars() : cnf.nVars();
        const std::uint64_t cl_count = cnf.nClauses();

        stats.num_literal_nodes = lit_count;
        stats.num_clause_nodes = cl_count;

        nodes.reserve(lit_count + cl_count);

        for (std::uint64_t var = 0; var < cnf.nVars(); ++var) {
            nodes.push_back({Node::Type::Literal, var, false, 0});
            if (cfg.split_nodes) {
                nodes.push_back({Node::Type::Literal, var, true, 0});
            }
        }

        for (std::uint64_t i = 0; i < cl_count; ++i) {
            nodes.push_back({Node::Type::Clause, i, false, 0});
        }

        adj.resize(nodes.size());

        for (std::uint64_t i = 0; i < cl_count; ++i) {
            const Cl* clause = cnf[i];
            std::uint64_t clause_idx = lit_count + i;
            nodes[clause_idx].size = clause->size();

            std::vector<Lit> lits(clause->begin(), clause->end());
            std::sort(lits.begin(), lits.end(),[](Lit a, Lit b) { 
                return a.x < b.x; 
            });
            if (cfg.remove_duplicate_literals) {
                lits.erase(std::unique(lits.begin(), lits.end(), [](Lit a, Lit b) { 
                    return a.x == b.x; 
                }), lits.end());
            }

            for (const Lit& lit : lits) {
                std::uint64_t lit_idx;
                if (cfg.split_nodes) {
                    lit_idx = 2 * lit.var() + (lit.sign() ? 1 : 0);
                } else {
                    lit_idx = lit.var();
                }

                Sign edge_sign = cfg.edge_labels
                    ? (lit.sign() ? Sign::Neg : Sign::Pos)
                    : Sign::Pos;

                adj[lit_idx].push_back({clause_idx, edge_sign});
                adj[clause_idx].push_back({lit_idx, edge_sign});
                ++nodes[lit_idx].size;
            }
        }
    }
};

struct Color {
    std::uint64_t hi = 0;
    std::uint64_t lo = 0;

    constexpr Color() = default;
    constexpr Color(std::uint64_t v) : hi(0), lo(v) {}
    constexpr Color(std::uint64_t hi_, std::uint64_t lo_) : hi(hi_), lo(lo_) {}

    friend constexpr bool operator==(const Color& a, const Color& b) noexcept {
        return a.hi == b.hi && a.lo == b.lo;
    }
    friend constexpr bool operator<(const Color& a, const Color& b) noexcept {
        if (a.hi != b.hi) { return a.hi < b.hi; }
        return a.lo < b.lo;
    }

    constexpr explicit operator std::uint64_t() const noexcept {
        return lo;
    }
};

struct Tools {
    // overflow safe pow
    static std::uint64_t ipow(std::uint64_t base, std::size_t exp) {
        constexpr std::uint64_t MAX = std::numeric_limits<std::uint64_t>::max();
        std::uint64_t result = 1;
        while (exp--) {
            if (base != 0 && result > MAX / base) {
                throw std::overflow_error("ipow overflow");
            }
            result *= base;
        }
        return result;
    }

    // Tuple Helper Functions: tuple  <--id_to_tuple--  id  <--tuple_to_id--  tuple
    static std::uint64_t tuple_to_id(const std::vector<std::uint64_t>& tuple, std::uint64_t nodes) {
        std::uint64_t id = 0;
        std::uint64_t base = 1;
        for (std::size_t i = 0; i < tuple.size(); ++i) {
            id += base * tuple[i];
            base *= nodes;
        }
        return id;
    }

    static void id_to_tuple(std::uint64_t id, std::size_t k, std::uint64_t nodes, std::vector<std::uint64_t>& tuple) {
        tuple.resize(k);
        for (std::size_t i = 0; i < k; ++i) {
            tuple[i] = id % nodes;
            id /= nodes;
        }
    }

    static Color hash(const void* data, std::size_t len, bool use128) {
        if (use128) {
            auto h = XXH3_128bits(data, len);
            return Color{h.high64, h.low64};
        }
        return Color{XXH3_64bits(data, len)};
    }

    static Color recolor(Color old_own_color, const std::vector<Color>& neighbours, bool use128) {
        std::vector<Color> buf;
        buf.reserve(neighbours.size() + 1);
        buf.push_back(old_own_color);
        buf.insert(buf.end(), neighbours.begin(), neighbours.end());
        return hash(buf.data(), buf.size() * sizeof(Color), use128);
    }

    static Color combine_sign(Color color, std::uint8_t sign) {
        Color out;
        out.hi = (color.hi << 1) | (color.lo >> 63);
        out.lo = (color.lo << 1) | static_cast<std::uint64_t>(sign);
        return out;
    }

    static std::string to_hex(Color color) {
        char buf[33];
        if (color.hi == 0) {
            std::snprintf(buf, sizeof(buf), "%016llx",
                          static_cast<unsigned long long>(color.lo));
        } else {
            std::snprintf(buf, sizeof(buf), "%016llx%016llx",
                          static_cast<unsigned long long>(color.hi),
                          static_cast<unsigned long long>(color.lo));
        }
        return std::string(buf);
    }
};

class Hasher {
public:
    explicit Hasher(Config cfg = {}) : cfg_(cfg) {}

    std::string operator()(const CNFFormula& cnf) {
        graph_.build(cnf, cfg_, stats_);
        return run_refinement();
    }

    const Stats& stats() const noexcept { 
        return stats_; 
    }

private:
    Config cfg_;
    Stats stats_;
    Graph graph_;

    std::string run_refinement() {
        const std::uint64_t nodes = graph_.nodes.size();
        const std::size_t k = cfg_.k_dimension;
        const std::uint64_t num_tuple = Tools::ipow(nodes, k);

        stats_.num_tuples = num_tuple;

        std::vector<Color> old_colors(num_tuple, Color{0});
        std::vector<Color> new_colors(num_tuple, Color{0});
        std::vector<Color> raw(num_tuple);

        std::vector<std::uint64_t> tuple;

        // initialize colors
        for (std::uint64_t id = 0; id < num_tuple; ++id) {
            Tools::id_to_tuple(id, k, nodes, tuple);

            Color color = Color{0};
            for (auto v : tuple) {
                Color initial_color = Color{graph_.nodes[v].size};
                color = Tools::recolor(color, {initial_color}, cfg_.use128);
            }
            new_colors[id] = color;
        }

        // refinement loop
        for (stats_.iterations = 0; stats_.iterations < cfg_.max_iterations; ++stats_.iterations) {
            // loop conditions
            if (new_colors == old_colors) {
                break;
            }
            old_colors.swap(new_colors);


            // compute raw hashes
            std::vector<std::uint64_t> child(k);
            std::vector<Color> neighbours;

            for (std::uint64_t id = 0; id < num_tuple; ++id) {
                Tools::id_to_tuple(id, k, nodes, tuple);
                neighbours.clear();

                for (std::size_t dim = 0; dim < k; ++dim) {
                    child = tuple;
                    for (const auto& e : graph_.adj[tuple[dim]]) {
                        child[dim] = e.dst;
                        std::uint64_t child_id = Tools::tuple_to_id(child, nodes);
                        neighbours.push_back(Tools::combine_sign(old_colors[child_id], static_cast<std::uint8_t>(e.sign)));
                    }
                }

                std::sort(neighbours.begin(), neighbours.end());
                raw[id] = Tools::recolor(old_colors[id], neighbours, cfg_.use128);
            }

            // canonise hashes into good color ids
            std::map<Color, std::uint64_t> canon;
            std::uint64_t next = 0;

            for (std::uint64_t id = 0; id < num_tuple; ++id) {
                auto [it, ins] = canon.emplace(raw[id], next);
                if (ins) {
                    ++next;
                }
                new_colors[id] = Color{it->second};
            }

            stats_.num_colors_per_iteration.push_back(next);
        }

        // --print-stats
        if (cfg_.print_stats) {
            std::cout << stats_.iterations << "\n";
            std::cout << stats_.num_literal_nodes << "," << stats_.num_clause_nodes << " => " << stats_.num_tuples << "\n";
            for (auto it : stats_.num_colors_per_iteration) {
                std::cout << it << "\n";
            }
        }

        // Hash in Hex as output
        return Tools::to_hex(Tools::hash(new_colors.data(), new_colors.size() * sizeof(Color), cfg_.use128));
    }
};

inline std::string wlhash(const char* filename, const Config& cfg = {}) {
    CNFFormula f(filename);
    f.normalizeVariableNames();
    return Hasher(cfg)(f);
}

}

#endif // WLISOHash_H_