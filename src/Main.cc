/*************************************************************************************************
CNFTools -- Copyright (c) 2020, Ashlin Iser, KIT - Karlsruhe Institute of Technology

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

#include <iostream>
#include <array>
#include <cstdio>
#include <filesystem>
#include <iomanip>

#include "src/external/argparse/argparse.h"
#include "src/external/ipasir.h"

#include "src/identify/GBDHash.h"
#include "src/identify/ISOHash.h"
#include "src/identify/WLHash.h"
#include "src/identify/TIMONHash.h"

#include "src/util/CNFFormula.h"
#include "src/util/SolverTypes.h"

#include "src/util/ResourceLimits.h"
#include "src/transform/cnf2bip.h"
#include "src/transform/cnf2kis.h"
#include "src/transform/cnf2cnf.h"

#include "src/extract/CNFSaniCheck.h"
#include "src/extract/CNFBaseFeatures.h"
#include "src/extract/CNFGateFeatures.h"
#include "src/extract/WCNFBaseFeatures.h"
#include "src/extract/OPBBaseFeatures.h"

#include "src/util/StreamCompressor.h"

template <bool one_pass, bool use_xor, bool use_sort, bool mate_ref, bool half_bit_hash>
void run_wlhash(const std::string& filename, const WLF::WLSettings& settings) {
    CNFFormula cnf(filename.c_str());
    WLF::WLHasher<one_pass, use_xor, use_sort, mate_ref, half_bit_hash> hasher(cnf, settings);
    
    auto stats = hasher.run();

    std::cerr << "c Hash: " << std::hex << std::setw(half_bit_hash ? 8 : 16) << std::setfill('0') 
            << stats.hash 
            << std::dec << std::setfill(' ') << std::endl;
    if (settings.print_stats) {
        std::cerr << "c Stabilized: " << (stats.stabilized ? "yes" : "no") << std::endl;
        std::cerr << "c Rounds: " << stats.round << std::endl;
    }
}

using WLHashFuncPtr = void (*)(const std::string&, const WLF::WLSettings&);

// [one_pass | use_xor | use_sort | mate_ref | 32bit_hash]
constexpr WLHashFuncPtr g_dispatch_table[32] = {
    &run_wlhash<0,0,0,0,0>, // 00000 | Index 0
    &run_wlhash<0,0,0,0,1>, // 00001 | Index 1
    &run_wlhash<0,0,0,1,0>, // 00010 | Index 2
    &run_wlhash<0,0,0,1,1>, // 00011 | Index 3

    &run_wlhash<0,0,1,0,0>, // 00100 | Index 4
    &run_wlhash<0,0,1,0,1>, // 00101 | Index 5
    &run_wlhash<0,0,1,1,0>, // 00110 | Index 6
    &run_wlhash<0,0,1,1,1>, // 00111 | Index 7

    &run_wlhash<0,1,0,0,0>, // 01000 | Index 8
    &run_wlhash<0,1,0,0,1>, // 01001 | Index 9
    &run_wlhash<0,1,0,1,0>, // 01010 | Index 10
    &run_wlhash<0,1,0,1,1>, // 01011 | Index 11

    &run_wlhash<0,1,1,0,0>, // 01100 | Index 12
    &run_wlhash<0,1,1,0,1>, // 01101 | Index 13
    &run_wlhash<0,1,1,1,0>, // 01110 | Index 14
    &run_wlhash<0,1,1,1,1>, // 01111 | Index 15

    &run_wlhash<1,0,0,0,0>, // 10000 | Index 16
    &run_wlhash<1,0,0,0,1>, // 10001 | Index 17
    &run_wlhash<1,0,0,1,0>, // 10010 | Index 18
    &run_wlhash<1,0,0,1,1>, // 10011 | Index 19

    &run_wlhash<1,0,1,0,0>, // 10100 | Index 20
    &run_wlhash<1,0,1,0,1>, // 10101 | Index 21
    &run_wlhash<1,0,1,1,0>, // 10110 | Index 22
    &run_wlhash<1,0,1,1,1>, // 10111 | Index 23

    &run_wlhash<1,1,0,0,0>, // 11000 | Index 24
    &run_wlhash<1,1,0,0,1>, // 11001 | Index 25
    &run_wlhash<1,1,0,1,0>, // 11010 | Index 26
    &run_wlhash<1,1,0,1,1>, // 11011 | Index 27

    &run_wlhash<1,1,1,0,0>, // 11100 | Index 28
    &run_wlhash<1,1,1,0,1>, // 11101 | Index 29
    &run_wlhash<1,1,1,1,0>, // 11110 | Index 30
    &run_wlhash<1,1,1,1,1>, // 11111 | Index 31
};

int main(int argc, char** argv) {
    argparse::ArgumentParser argparse("CNF Tools");

    argparse.add_argument("tool").help("Select Tool: id, isohash, wlhash, timonhash, normalize, sanitize, checksani, cnf2kis, cnf2bip, extract, gates")
        .default_value("identify")
        .action([](const std::string& value) {
            static const std::vector<std::string> choices = { "id", "isohash", "wlhash", "timonhash", "normalize", "sanitize", "checksani", "cnf2kis", "cnf2bip", "extract", "gates" };
            if (std::find(choices.begin(), choices.end(), value) != choices.end()) {
                return value;
            }
            return std::string{ "identify" };
        });

    argparse.add_argument("file").help("Path to Input File");
    argparse.add_argument("-o", "--output").default_value(std::string("-")).help("Path to Output File (used by cnf2* transformers, default is stdout)");
    argparse.add_argument("-t", "--timeout").default_value(0).scan<'i', int>().help("Time limit in seconds");
    argparse.add_argument("-m", "--memout").default_value(0).scan<'i', int>().help("Memory limit in MB");
    argparse.add_argument("-f", "--fileout").default_value(0).scan<'i', int>().help("File size limit in MB"); 
    
    // flags for WL Hash
    argparse.add_argument("--max-iters").default_value(100u).scan<'i', unsigned>().help("Maximum WL iterations before stopping");
    argparse.add_argument("--print-stats").default_value(false).implicit_value(true).help("Print useful stats");
    argparse.add_argument("--one-pass").default_value(false).implicit_value(true).help("Use only one pass");
    argparse.add_argument("--use-xor").default_value(false).implicit_value(true).help("use xor instead of normal sums");
    argparse.add_argument("--use-sort").default_value(false).implicit_value(true).help("sort instead of adding or xoring");
    argparse.add_argument("--mate-ref").default_value(false).implicit_value(true).help("always build a reference hash of a literals other polarity hash");
    argparse.add_argument("--half-bit").default_value(false).implicit_value(true).help("use only 32 bit");
    
    // flags for timon hash
    argparse.add_argument("--no-sort").default_value(false).implicit_value(true).help("disable sorting for timon isohash");
    argparse.add_argument("--no-rehash").default_value(false).implicit_value(true).help("disable rehashing for timon isohash");
    
    try {
        argparse.parse_args(argc, argv);
    }
    catch (const std::runtime_error& err) {
        std::cout << err.what() << std::endl;
        std::cout << argparse;
        exit(0);
    }

    std::string filename = argparse.get("file");
    std::string toolname = argparse.get("tool");
    std::string output = argparse.get("output");

    ResourceLimits limits(argparse.get<int>("timeout"), argparse.get<int>("memout"), argparse.get<int>("fileout"));
    limits.set_rlimits();
    std::cerr << "c Running: " << toolname << " " << filename << std::endl;

    std::string ext = std::filesystem::path(filename).extension();
    if (ext == ".xz" || ext == ".lzma" || ext == ".bz2" || ext == ".gz") {
        ext = std::filesystem::path(filename).stem().extension();
    }

    if (ext == ".cnf" || ext == ".wecnf") {
        std::cerr << "Detected CNF" << std::endl;
    } else if (ext == ".opb") {
        std::cerr << "Detected OPB" << std::endl;
    } else if (ext == ".qcnf" || ext == ".qdimacs") {
        std::cerr << "Detected QBF" << std::endl;
    } else if (ext == ".wcnf") {
        std::cerr << "Detected WCNF" << std::endl;
    }

    try {
        if (toolname == "id") {
            if (ext == ".cnf" || ext == ".wecnf") {
                std::cout << CNF::gbdhash(filename.c_str()) << std::endl;
            }
            else if (ext == ".opb") {
                std::cout << OPB::gbdhash(filename.c_str()) << std::endl;
            }
            else if (ext == ".qcnf" || ext == ".qdimacs") {
                std::cout << PQBF::gbdhash(filename.c_str()) << std::endl;
            }
            else if (ext == ".wcnf") {
                std::cout << WCNF::gbdhash(filename.c_str()) << std::endl;
            }
        }
        else if (toolname == "wlhash") {
            if (ext == ".cnf") {
                WLF::WLSettings config;
                config.max_iterations = argparse.get<unsigned>("--max-iters");
                config.print_stats = argparse.get<bool>("--print-stats");
                size_t index = 0;
                index |= (1 << 4); // one_pass default true
                index |= (1 << 1); // mate_ref default true

                if (argparse.get<bool>("--one-pass"))   index ^= (1 << 4); // turns it false
                if (argparse.get<bool>("--use-xor"))    index ^= (1 << 3);
                if (argparse.get<bool>("--use-sort"))   index ^= (1 << 2);
                if (argparse.get<bool>("--mate-ref"))   index ^= (1 << 1); // turns it false
                if (argparse.get<bool>("--half-bit"))   index ^= (1 << 0);
                
                g_dispatch_table[index](filename, config);
            }
        }
        else if (toolname == "timonhash") {
            if (ext == ".cnf") {
                unsigned depth = 200;
                bool sort = false;
                bool rehash = false;
                depth = 2 * argparse.get<unsigned>("--max-iters");
                sort = !argparse.get<bool>("--no-sort");
                rehash = !argparse.get<bool>("--no-rehash");

                std::cout << "c timonhash: max_iters=" << depth << " sort=" << (sort ? 1 : 0) << " rehash=" << (rehash ? 1 : 0) << std::endl;
                std::cerr << CNF::weisfeiler_leman_hash(filename.c_str(), 0, true, true, false, depth, true, rehash, true, 6u, false, false, sort) << std::endl;
            }
        }
        else if (toolname == "isohash") {
            if (ext == ".cnf") {
                std::cerr << "c Hash: " << CNF::isohash(filename.c_str()) << std::endl;
            }
            else if (ext == ".wcnf") {
                std::cout << "c Hash: " << WCNF::isohash(filename.c_str()) << std::endl;
            }
        } 
        else if (toolname == "normalize") {
            std::cerr << "Normalizing " << filename << std::endl;
            CNF::Normaliser norm(filename.c_str(), output == "-" ? nullptr : output.c_str());
            norm.run();
        } 
        else if (toolname == "sanitize") {
            CNF::Sanitiser sani(filename.c_str(), output == "-" ? nullptr : output.c_str());
            sani.run();
        }
        else if (toolname == "checksani") {
            CNF::SaniCheck ana(filename.c_str(), true);
            ana.run();
            bool sani = true;
            std::cout << "hash" << " " << CNF::gbdhash(filename.c_str()) << std::endl;
            std::cout << "filename " << filename << std::endl;
            sani = ana.getFeature("head_vars") == ana.getFeature("norm_vars") && ana.getFeature("head_clauses") == ana.getFeature("norm_clauses");
            // std::cout << ana.getFeature("head_vars") << " " << ana.getFeature("norm_vars") << std::endl;
            // std::cout << ana.getFeature("head_clauses") << " " << ana.getFeature("norm_clauses") << std::endl;
            std::cout << "header_consistent " << (sani ? "yes" : "no") << std::endl;
            sani = ana.getFeature("whitespace_normalised") == 1.0;
            std::cout << "whitespace_normalised " << (sani ? "yes" : "no") << std::endl;
            sani = ana.getFeature("has_comment") == 0.0;
            std::cout << "no_comment " << (sani ? "yes" : "no") << std::endl;
            sani = ana.getFeature("has_tautological_clause") == 0.0;
            std::cout << "no_tautological_clause " << (sani ? "yes" : "no") << std::endl;
            sani = ana.getFeature("has_duplicate_literals") == 0.0;
            std::cout << "no_duplicate_literals " << (sani ? "yes" : "no") << std::endl;
            sani = ana.getFeature("has_empty_clause") == 0.0;
            std::cout << "no_empty_clause " << (sani ? "yes" : "no") << std::endl;
        }
        else if (toolname == "cnf2kis") {
            std::cerr << "Generating Independent Set Problem " << filename << std::endl;
            IndependentSetFromCNF gen(filename.c_str());
            gen.generate_independent_set_problem(output == "-" ? nullptr : output.c_str());
        }
        else if (toolname == "cnf2bip") {
            std::cerr << "Generating Bipartite Graph " << filename << std::endl;
            CNF::cnf2bip gen(filename.c_str(), output == "-" ? nullptr : output.c_str());
            gen.run();
        }
        else if (toolname == "extract" || toolname == "gates") {
            IExtractor *stats;
            if (toolname == "extract") {
                if (ext == ".cnf") {
                    stats = new CNF::BaseFeatures(filename.c_str());
                }
                else if (ext == ".wcnf") {
                    stats = new WCNF::BaseFeatures(filename.c_str());
                }
                else if (ext == ".opb") {
                    stats = new OPB::BaseFeatures(filename.c_str());
                }
                else {
                    std::cerr << "Format " << ext << " not supported by extract" << std::endl;
                    return 1;
                }
            }
            else if (toolname == "gates") {
                stats = new CNF::GateFeatures(filename.c_str());
                if (ext != ".cnf") {
                    std::cerr << "Format " << ext << " not supported by extract" << std::endl;
                    return 1;
                }
            }
            stats->run();
            for (std::string name : stats->getNames()) {
                std::cout << name << " ";
            }
            std::cout << std::endl;
            for (double feature : stats->getFeatures()) {
                std::cout << feature << " ";
            }
            std::cout << std::endl;
        }
    }
    catch (std::bad_alloc& e) {
        std::cerr << "Memory Limit Exceeded" << std::endl;
        return 1;
    }
    catch (MemoryLimitExceeded& e) {
        std::cerr << "Memory Limit Exceeded" << std::endl;
        return 1;
    }
    catch (TimeLimitExceeded& e) {
        std::cerr << "Time Limit Exceeded" << std::endl;
        return 1;
    }
    catch (FileSizeLimitExceeded& e) {
        std::remove(output.c_str());
        std::cerr << "File Size Limit Exceeded" << std::endl;
        return 1;
    }
    return 0;
}
