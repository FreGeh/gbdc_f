/*************************************************************************************************
CNFTools -- Copyright (c) 2021, Ashlin Iser, KIT - Karlsruhe Institute of Technology

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

#ifndef SRC_UTIL_CNFFORMULA_H_
#define SRC_UTIL_CNFFORMULA_H_

#include <vector>
#include <algorithm>
#include <memory>
#include <string>

#include "src/util/StreamBuffer.h"
#include "src/util/SolverTypes.h"

class CNFFormula {
    For formula;
    unsigned variables;
    unsigned total_literals;

 public:
    CNFFormula() : formula(), variables(0), total_literals(0) { }

    explicit CNFFormula(const char* filename) : CNFFormula() {
        readDimacsFromFile(filename);
    }

    ~CNFFormula() {
        for (Cl* clause : formula) {
            delete clause;
        }
    }

    typedef For::const_iterator const_iterator;

    inline const_iterator begin() const {
        return formula.begin();
    }

    inline const_iterator end() const {
        return formula.end();
    }

    inline const Cl* operator[] (int i) const {
        return formula[i];
    }

    inline size_t nVars() const {
        return variables;
    }

    inline size_t nLits() const {
        return total_literals;
    }

    inline size_t nClauses() const {
        return formula.size();
    }

    inline int newVar() {
        return ++variables;
    }

    inline void clear() {
        for (Cl* clause : formula) delete clause;
        formula.clear();
        variables = 0;
        total_literals = 0;
    }

    void normalizeVariableNames() {
        std::vector<unsigned> seen(variables + 1, 0);
        for (Cl* clause : formula) {
            for (Lit& lit : *clause) {
                seen[lit.var()] = 1;
            }
        }

        std::vector<unsigned> map(variables + 1, 0);
        unsigned next = 1;
        for (unsigned v = 1; v <= variables; ++v) {
            if (seen[v]) map[v] = next++;
        }

        for (Cl* clause : formula) {
            for (Lit& lit : *clause) {
                lit = Lit(map[lit.var()], lit.sign());
            }
        }

        variables = next - 1;
    }

    bool hasEmptyClause() const {
        for (const Cl* clause : formula) {
            if (clause->empty()) return true;
        }

        return false;
    }

    void removeDuplicateClauses() {
        std::sort(formula.begin(), formula.end(), [](const Cl* a, const Cl* b) {
            return std::lexicographical_compare(a->begin(), a->end(), b->begin(), b->end());
        });

        size_t write = 0;
        total_literals = 0;
        variables = 0;
        for (Cl* clause : formula) {
            bool duplicate = write > 0 && clause->size() == formula[write-1]->size() && std::equal(clause->begin(), clause->end(), formula[write-1]->begin());
            if (duplicate) {
                delete clause;
            } else {
                formula[write++] = clause;
                total_literals += clause->size();
                for (const Lit& lit : *clause) {
                    variables = std::max(variables, (unsigned)lit.var());
                }
            }
        }

        formula.resize(write);
        formula.shrink_to_fit();
    }

    void normalizeLogicalCNF() {
        removeDuplicateClauses();

        // needed if there is empty clause -> whole formula becomes false
        if (hasEmptyClause()) {
            clear();
            Cl* empty = new Cl();
            formula.push_back(empty);
            variables = 0;
            total_literals = 0;
            return;
        }

        normalizeVariableNames();
    }

    void readDimacsFromFile(const char* filename) {
        StreamBuffer in(filename);
        Cl clause;
        while (in.skipWhitespace()) {
            if (*in == 'p' || *in == 'c') {
                if (!in.skipLine()) break;
            } else {
                int plit;
                while (in.readInteger(&plit)) {
                    if (plit == 0) break;
                    clause.push_back(Lit(abs(plit), plit < 0));
                }
                readClause(clause.begin(), clause.end());
                clause.clear();
            }
        }
    }

    void readClause(std::initializer_list<Lit> list) {
        readClause(list.begin(), list.end());
    }

    void readClauses(const For& formula) {
        for (Cl* clause : formula) {
            readClause(clause->begin(), clause->end());
        }
    }

    template <typename Iterator>
    void readClause(Iterator begin, Iterator end) {
        Cl* clause = new Cl { begin, end };
        if (clause->size() > 0) {
            // remove redundant literals
            std::sort(clause->begin(), clause->end());
            unsigned dup = 0;
            for (auto it = clause->begin(), jt = clause->begin()+1; jt != clause->end(); ++jt) {
                if (*it != *jt) {  // unique
                    if (it->var() == jt->var()) {
                        delete clause;
                        return;  // no tautologies
                    }
                    ++it;
                    *it = *jt;
                } else {
                    ++dup;
                }
            }
            clause->resize(clause->size() - dup);
            clause->shrink_to_fit();
            variables = std::max(variables, (unsigned int)clause->back().var());
            total_literals += clause->size();
        }
        formula.push_back(clause);
    }
};

#endif  // SRC_UTIL_CNFFORMULA_H_

