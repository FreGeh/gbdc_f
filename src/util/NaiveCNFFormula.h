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

#ifndef SRC_UTIL_NAIVECNFFORMULA_H_
#define SRC_UTIL_NAIVECNFFORMULA_H_

#include <vector>
#include <algorithm>
#include <memory>
#include <string>

#include "src/util/StreamBuffer.h"
#include "src/util/SolverTypes.h"

class NaiveCNFFormula {
    For formula;
    unsigned variables = 0;
    unsigned literals = 0;

 public:
    explicit inline NaiveCNFFormula(const char* filename, const bool shrink_to_fit) {
        readDimacsFromFile(filename, shrink_to_fit);
    }
    NaiveCNFFormula(const NaiveCNFFormula&) = delete;
    NaiveCNFFormula& operator=(const NaiveCNFFormula&) = delete;
    ~NaiveCNFFormula() {
        for (Cl* clause : formula) {
            delete clause;
        }
    }

    inline size_t nVars() const {
        return variables;
    }
    inline size_t nClauses() const {
        return formula.size();
    }
    inline size_t nLiterals() const {
        return literals;
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
        For::const_iterator it;
        inline ClauseIt& operator ++ () {
            ++it;
            return *this;
        }
        inline Clause operator * () {
            return Clause {(*it)->begin(), (*it)->end()};
        }
        inline bool operator != (const ClauseIt o) {
            return it != o.it;
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
        return Clauses {
            {formula.begin()},
            {formula.end()}
        };
    }

 private:
    // create gapless representation of variables
    void normalizeVariableNames() {
        std::vector<unsigned> name;
        constexpr unsigned empty = ~0U;
        name.resize(variables+1, empty);
        unsigned int max = 0;
        for (Cl* clause : formula) {
            for (Lit& lit : *clause) {
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
                Cl* clause = new Cl;
                int plit;
                while (in.readInteger(&plit)) {
                    if (plit == 0) break;
                    const unsigned var = abs(plit);
                    clause->push_back(Lit(var, plit < 0));
                    if (var > variables) variables = var;
                }
                if (clause->size() == 0) {
                    delete clause;
                    continue;
                }
                if (shrink_to_fit) clause->shrink_to_fit();
                formula.push_back(clause);
                literals += clause->size();
            }
        }
        normalizeVariableNames();
    }
};

#endif  // SRC_UTIL_CNFFORMULA_H_
