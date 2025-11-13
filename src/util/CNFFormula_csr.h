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
#ifndef SRC_UTIL_CSR_H_
#define SRC_UTIL_CSR_H_

#include <vector>
#include <algorithm>
#include <memory>
#include <string>
#include <cstddef> // For std::ptrdiff_t

#include "src/util/StreamBuffer.h"
#include "src/util/SolverTypes.h"

class CNFFormula_csr {
private:
    std::vector<Lit> _all_literals;
    std::vector<size_t> _clause_offsets;
    unsigned variables;

public:
    class ClauseView {
    private:
        const Lit* _begin;
        const Lit* _end;
    public:
        ClauseView(const Lit* b, const Lit* e) : _begin(b), _end(e) {}
        const Lit* begin() const { return _begin; }
        const Lit* end() const { return _end; }
        size_t size() const { return _end - _begin; }
    };

public: // Constructor and main accessors
    CNFFormula_csr() : variables(0) {
        _clause_offsets.push_back(0);
    }

    explicit CNFFormula_csr(const char* filename) : CNFFormula_csr() {
        readDimacsFromFile(filename);
    }

    ~CNFFormula_csr() = default;

    inline size_t nVars() const { return variables; }
    inline size_t nLits() const { return _all_literals.size(); }
    inline size_t nClauses() const { return _clause_offsets.size() - 1; }

    inline ClauseView operator[](size_t i) const {

        size_t start = _clause_offsets[i];
        size_t end = _clause_offsets[i+1];
        return ClauseView(&_all_literals[start], &_all_literals[end]);
    }

public: // Iterators for range-based for loops (for (auto clause : cnf))
    class ClauseIterator {
    public:
        using iterator_category = std::forward_iterator_tag;
        using value_type = ClauseView;
        using difference_type = std::ptrdiff_t;
        using pointer = const ClauseView*;
        using reference = ClauseView; // Return by value is efficient for this small view

    private:
        const CNFFormula_csr* _formula;
        size_t _clause_index;

    public:
        ClauseIterator(const CNFFormula_csr* formula, size_t index)
            : _formula(formula), _clause_index(index) {}

        reference operator*() const { return (*_formula)[_clause_index]; }

        ClauseIterator& operator++() { ++_clause_index; return *this; }
        ClauseIterator operator++(int) { ClauseIterator tmp = *this; ++(*this); return tmp; }
        bool operator==(const ClauseIterator& other) const { return _clause_index == other._clause_index; }
        bool operator!=(const ClauseIterator& other) const { return !(*this == other); }
    };

    ClauseIterator begin() const { return ClauseIterator(this, 0); }
    ClauseIterator end() const { return ClauseIterator(this, nClauses()); }

public: // Modification and parsing methods
    inline int newVar() { return ++variables; }

    inline void clear() {
        _all_literals.clear();
        _clause_offsets.clear();
        variables = 0;
        _clause_offsets.push_back(0);
    }

    inline void normalizeVariableNames() {
        if (variables == 0 || _all_literals.empty()) return;

        std::vector<char> seen(variables + 1, 0);
        for (const Lit& lit : _all_literals) {
            seen[lit.var().id] = 1;
        }

        std::vector<unsigned> map(variables + 1, 0);
        unsigned next = 1;
        for (unsigned v = 1; v <= variables; ++v) {
            if (seen[v]) map[v] = next++;
        }

        for (Lit& lit : _all_literals) {
            lit = Lit(map[lit.var().id], lit.sign());
        }
        variables = next - 1;
    }

    template <typename Iterator>
    void readClause(Iterator begin, Iterator end) {
        Cl temp_clause(begin, end);

        if (temp_clause.empty()) return;

        std::sort(temp_clause.begin(), temp_clause.end());
        auto last = std::unique(temp_clause.begin(), temp_clause.end());
        temp_clause.erase(last, temp_clause.end());

        for (size_t i = 0; i + 1 < temp_clause.size(); ++i) {
            if (temp_clause[i].var() == temp_clause[i+1].var()) {
                return; // Tautology, ignore clause.
            }
        }

        variables = std::max(variables, (unsigned)temp_clause.back().var().id);

        _all_literals.insert(_all_literals.end(), temp_clause.begin(), temp_clause.end());
        _clause_offsets.push_back(_all_literals.size());
    }

    void readDimacsFromFile(const char* filename) {
        StreamBuffer in(filename);
        Cl clause_buffer;
        
        unsigned parsed_vars = 0;
        unsigned parsed_clauses = 0;

        while (in.skipWhitespace()) {
            if (*in == 'p') {
                if (!in.skipString("p cnf")) {
                    return;
                }
                if (!in.readInteger((int*)&parsed_vars) || !in.readInteger((int*)&parsed_clauses)) {
                    return;
                }
                
                _clause_offsets.reserve(parsed_clauses + 1);
                size_t estimated_literals = static_cast<size_t>(parsed_clauses) * 3;
                _all_literals.reserve(estimated_literals);
                
                if (!in.skipLine()) break;

            } else if (*in == 'c') {
                if (!in.skipLine()) break;
            } else {
                int plit;
                while (in.readInteger(&plit)) {
                    if (plit == 0) break;
                    clause_buffer.push_back(Lit(abs(plit), plit < 0));
                }
                readClause(clause_buffer.begin(), clause_buffer.end());
                clause_buffer.clear();
            }
        }
    }

    void readClause(std::initializer_list<Lit> list) { readClause(list.begin(), list.end()); }
    void readClauses(const For& formula) {
        for (Cl* clause : formula) { readClause(clause->begin(), clause->end()); }
    }
};

#endif // SRC_UTIL_CNFFORMULA_CSR_H_