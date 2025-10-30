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

#ifndef SRC_UTIL_SIZEGROUPEDCNFFORMULA
#define SRC_UTIL_SIZEGROUPEDCNFFORMULA

#include <vector>
#include <algorithm>
#include <memory>
#include <string>

#include "src/util/StreamBuffer.h"
#include "src/util/SolverTypes.h"

class SizeGroupedCNFFormula {
    std::vector<std::vector<Lit>*> clause_length_literals;
    unsigned variables = 0;
    unsigned n_clauses = 0;
    unsigned literals = 0;

 public:
    explicit inline SizeGroupedCNFFormula(const char* filename, const bool shrink_to_fit) {
        readDimacsFromFile(filename, shrink_to_fit);
    }
    ~SizeGroupedCNFFormula() {
        for (const std::vector<Lit>* clause_length : clause_length_literals)
            delete clause_length;
    }

    inline size_t nVars() const {
        return variables;
    }
    inline size_t nClauses() const {
        return n_clauses;
    }
    inline size_t nLiterals() const {
        return literals;
    }
    inline size_t maxClauseLength() const {
        return clause_length_literals.size();
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
        const std::vector<std::vector<Lit>*>& clause_length_literals;
        std::size_t length;
        std::vector<Lit>::const_iterator clause_begin;
        std::vector<Lit>::const_iterator subvector_end;
        inline ClauseIt& operator ++ () {
            clause_begin += length;
            if (clause_begin == subvector_end) {
                while (++length != clause_length_literals.size() && clause_length_literals[length]->size() == 0)
                    ;
                if (length != clause_length_literals.size()) {
                    clause_begin = clause_length_literals[length]->begin();
                    subvector_end = clause_length_literals[length]->end();
                }
            }
            return *this;
        }
        inline Clause operator * () {
            return Clause {clause_begin, clause_begin + length};
        }
        inline bool operator != (const ClauseIt o) {
            return length != o.length || clause_begin != o.clause_begin;
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
            ++ClauseIt{clause_length_literals, 0, clause_length_literals[0]->begin(), clause_length_literals[0]->end()},
            {clause_length_literals, clause_length_literals.size(), (*(clause_length_literals.end() - 1))->end()}
        };
    }

 private:
    // create gapless representation of variables
    void normalizeVariableNames() {
        std::vector<unsigned> name;
        constexpr unsigned empty = ~0U;
        name.resize(variables+1, empty);
        unsigned max = 0;
        for (std::vector<Lit>* clause_length : clause_length_literals) {
            for (Lit& lit : *clause_length) {
                if (name[lit.var()] == empty) name[lit.var()] = max++;
                lit = Lit(name[lit.var()], lit.sign());
            }
        }
        variables = max;
    }

    // https://stackoverflow.com/a/1322548/27720282
    static inline unsigned next_power_of_2(unsigned n) {
        n--;
        n |= n >> 0b00001;
        n |= n >> 0b00010;
        n |= n >> 0b00100;
        n |= n >> 0b01000;
        n |= n >> 0b10000;
        n++;
        return n;
    }
    void readDimacsFromFile(const char* filename, const bool shrink_to_fit) {
        StreamBuffer in(filename);
        std::vector<Lit> clause;
        while (in.skipWhitespace()) {
            if (*in == 'p' || *in == 'c') {
                if (!in.skipLine()) break;
            } else {
                int plit;
                while (in.readInteger(&plit)) {
                    if (plit == 0) break;
                    const unsigned var = abs(plit);
                    clause.push_back(Lit(var, plit < 0));
                    if (var > variables) variables = var;
                }
                if (clause.size() >= clause_length_literals.size()) {
                    const unsigned old_size = clause_length_literals.size();
                    const unsigned new_size = clause.size() + 1;
                    clause_length_literals.reserve(next_power_of_2(new_size));
                    clause_length_literals.resize(new_size);
                    for (auto it = clause_length_literals.begin() + old_size; it != clause_length_literals.end(); ++it)
                        *it = new std::vector<Lit>;
                }
                std::vector<Lit>& insert_here = *clause_length_literals[clause.size()];
                insert_here.reserve(next_power_of_2(insert_here.size() + clause.size()));
                insert_here.insert(insert_here.end(), clause.begin(), clause.end());
                ++n_clauses;
                literals += clause.size();
                clause.clear();
            }
        }
        if (shrink_to_fit)
            for (std::vector<Lit>* clause_length : clause_length_literals)
                clause_length->shrink_to_fit();
        normalizeVariableNames();
    }
};

#endif  // SRC_UTIL_SIZEGROUPEDCNFFORMULA
