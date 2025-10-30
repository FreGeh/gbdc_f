/*************************************************************************************************
CNFTools -- Copyright (c) 2021, Markus Iser, KIT - Karlsruhe Institute of Technology

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

#ifndef SRC_UTIL_INTERVALCNFFORMULA
#define SRC_UTIL_INTERVALCNFFORMULA

#include <vector>
#include <algorithm>
#include <memory>
#include <string>

#include "src/util/StreamBuffer.h"
#include "src/util/SolverTypes.h"

class IntervalCNFFormula {
    std::vector<Lit> literals;
    unsigned variables = 0;
    unsigned n_clauses = 0;
    unsigned n_literals = 0;

 public:
    explicit inline IntervalCNFFormula(const char* filename, const bool shrink_to_fit) {
        readDimacsFromFile(filename, shrink_to_fit);
    }

    inline size_t nVars() const {
        return variables;
    }
    inline size_t nClauses() const {
        return n_clauses;
    }
    inline size_t nLiterals() const {
        return n_literals;
    }
    inline size_t maxClauseLength() const {
        return 0; // dummy, only interesting for SizeGroupedCNFFormula
    }

    struct Clause {
        using It = std::vector<Lit>::const_iterator;
        const It begin_, end_;
        inline It begin() const {
            return begin_;
        }
        inline It end() const {
            return end_;
        }
        inline unsigned size() const {
            return end_ - begin_;
        }
    };
    struct ClauseIt {
        std::vector<Lit>::const_iterator length_slot;
        inline ClauseIt& operator ++ () {
            length_slot += length_slot->x;
            return *this;
        }
        inline Clause operator * () {
            return Clause {length_slot + 1, length_slot + length_slot->x};
        }
        inline bool operator != (const ClauseIt o) {
            return length_slot != o.length_slot;
        }
    };
    struct Clauses {
        const ClauseIt begin_, end_;
        inline ClauseIt begin() const {
            return begin_;
        }
        inline ClauseIt end() const {
            return end_;
        }
    };
    Clauses clauses() const {
        return Clauses {{literals.begin()}, {literals.end()}};
    }

 private:
 struct MutClause {
        using It = std::vector<Lit>::iterator;
        const It begin_, end_;
        inline It begin() const {
            return begin_;
        }
        inline It end() const {
            return end_;
        }
        inline unsigned size() const {
            return end_ - begin_;
        }
    };
    struct MutClauseIt {
        std::vector<Lit>::iterator length_slot;
        inline MutClauseIt& operator ++ () {
            length_slot += length_slot->x;
            return *this;
        }
        inline MutClause operator * () {
            return MutClause {length_slot + 1, length_slot + length_slot->x};
        }
        inline bool operator != (const MutClauseIt o) {
            return length_slot != o.length_slot;
        }
    };
    struct MutClauses {
        const MutClauseIt begin_, end_;
        inline MutClauseIt begin() const {
            return begin_;
        }
        inline MutClauseIt end() const {
            return end_;
        }
    };
    MutClauses mut_clauses() {
        return MutClauses {{literals.begin()}, {literals.end()}};
    }
    // create gapless representation of variables
    void normalizeVariableNames() {
        std::vector<unsigned> name;
        constexpr unsigned empty = ~0U;
        name.resize(variables+1, empty);
        unsigned max = 0;
        for (MutClause cl : mut_clauses()) {
            for (Lit& lit : cl) {
                if (name[lit.var()] == empty) name[lit.var()] = max++;
                lit = Lit(name[lit.var()], lit.sign());
            }
        }
        variables = max;
    }

    void readDimacsFromFile(const char* filename, const bool shrink_to_fit) {
        StreamBuffer in(filename);
        while (in.skipWhitespace()) {
            if (*in == 'p' || *in == 'c') {
                if (!in.skipLine()) break;
            } else {
                literals.push_back({});
                unsigned length = 1; // including length slot
                int plit;
                while (in.readInteger(&plit)) {
                    if (plit == 0) break;
                    const unsigned var = abs(plit);
                    literals.push_back(Lit(var, plit < 0));
                    ++length;
                    if (var > variables) variables = var;
                }
                if (length != 1) {
                    literals[literals.size() - length].x = length;
                    ++n_clauses;
                    n_literals += length - 1;
                } else
                    literals.pop_back();
            }
        }
        if (shrink_to_fit) literals.shrink_to_fit();
        normalizeVariableNames();
    }
};

#endif  // SRC_UTIL_INTERVALCNFFORMULA
