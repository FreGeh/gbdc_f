#pragma once

#include <vector>

#include "src/extract/IExtractor.h"
#include "src/util/CNFFormula.h"

namespace CNF {

// SAT Isomorphism Graph (SIG)
class cnf2sig : public IExtractor {
 private:
    CNFFormula F;
    const char* filename_;
    const char* output_;

 public:
    cnf2sig(const char* filename, const char* output = nullptr);
    virtual ~cnf2sig();
    virtual void run();
};

}  // namespace CNF