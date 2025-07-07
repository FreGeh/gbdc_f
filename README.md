# GBDC: Global Benchmark Database, C++ Extension Module

[![Build-Test](https://github.com/Udopia/gbdc/actions/workflows/build_test.yml/badge.svg?branch=master)](https://github.com/Udopia/gbdc/actions/workflows/build_test.yml)

[GBDC](https://github.com/Udopia/gbdc) provides efficient implementations of functions for benchmark instance identification, instance feature extraction and instance transformation.
GBDC provides a command-line tool as well as the Python package `gbdc`.
The Python package `gbdc` is used by [Global Benchmark Database](https://github.com/Udopia/gbd).

## [Documentation](https://udopia.github.io/gbdc/doc/Index.html)

GBDC provides benchmark instance identifiers, feature extractors, and instance transformers for several problem domains, including propositional satisfiability (SAT) and optimization (MaxSAT), as well as Pseudo-Boolean Optimization (PBO).
A description of the supported domains, feature extractors, and instance transformers can be found in the [documentation](https://udopia.github.io/gbdc/doc/Index.html).

## Installation from PyPI
* Pre-built distributions for Linux and MacOS.
* Requires at least Python 3.8.0 (3.10.0 for Apple Silicon).
* Installation via `pip install gbdc`

## Installation from Source

* GBDC uses `libarchive` for reading from a large variety of compressed formats (in some systems provided by the package `libarchive-dev`).
* Some GBDC functions use an [IPASIR](https://github.com/biotomas/ipasir) SAT Solver. GBDC's build-system pulls the external SAT Solver [CaDiCaL](http://fmv.jku.at/cadical/) by A. Biere (MIT licensed).

### Steps:
1. **Install Dependencies** (libarchive, pybind, ninja)
   - For Ubuntu: `apt install libarchive-dev pybind11-dev ninja-build`
   - For macOS: `brew install libarchive pybind11 ninja`

2. Run `pip install . --user` in the repository directory.


### Build & Use Executeable
1. `mkdir build && cd build`
2. `cmake ..`
3. `make`
4. `./gbdc wlhash ../test/resources/test_files/0a4ed112f2cdc0a524976a15d1821097-cliquecoloring_n12_k9_c8.cnf.xz --print-stats`

<!-- #### Shipped Dependencies

* A copy of the command-line argument parser by P. S. Kumar [`argparse.h`](https://github.com/p-ranav/argparse) (MIT licensed) resides in the `lib` folder.

* A copy of the [MD5 hash](https://github.com/CommanderBubble/MD5) implementation by M. Lloyd (MIT licensed) resides in the `lib` folder. -->

<!-- ## Publications

* Gate feature extraction uses our gate recognition algorithm which is described in the following publications:

    * [*Recognition of Nested Gates in CNF Formulas* (SAT 2015, Iser et al.)](https://rdcu.be/czCr1)

    * [*Recognition and Exploitation of Gate Structure in SAT Solving* (2020, Iser)](https://d-nb.info/1209199122/34)

* The Python module `gbdc` is used in our project [GBD Benchmark Database](https://github.com/Udopia/gbd)

    * [*Collaborative Management of Benchmark Instances and their Attributes* (2020, Iser et al.)](https://arxiv.org/pdf/2009.02995.pdf) -->


Requirements: CMake, gcc, py-pybind11, libarchive
1. spack cmake packages installieren
2. gbdc bauen etc
2. gbd mit pip installieren
3. `gbd -d /nfs/share/instances/gbd/cnf_local.db get -r local -c min > gbd.ae.txt`
4. mit sbatch und slurm server 
5. analysieren mit pandas, etc

Ausgangslage: Datei mit jeweils [GBDHash] [datei] space seperated

- Starten: `sbatch ~/gbdc_f/gbdc_wlhash.slurm`
- Status: `scontrol show job` und `ls ~/gbdc_f/out | wc -l` und 

`f=$(ls -t ~/gbdc_f/out/*.out | head -1); id=${f##*/}; id=${id%.out}; wl=$(grep -Eo '^[0-9a-f]{16}$' "$f"); orig=$(awk -v h="$id" '$1==h{ $1=""; sub(/^ /,""); print; exit }' ~/gbdc_f/instances.lst); printf 'GBD: %s\nWL : %s\nCNF: %s\n' "$id" "$wl" "$orig"; ~/gbdc_f/build/gbdc wlhash "$orig"` um zu überprüfen ob letztgeschriebenes richtig oder falsch ist

Um jetzt wlhash wirklich zu extrahieren
```
OUTDIR=$HOME/gbdc_f/out
LIST=$HOME/gbdc_f/instances.lst
RESULT=$HOME/gbdc_f/wlhash_results.tsv

paste \
  <(awk '{print $1}' "$LIST") \
  <(for f in "$OUTDIR"/*.out; do
        grep -Eo '^[0-9a-f]{16}$' "$f"
    done) \
  > "$RESULT"

echo "→ Ergebnis in $RESULT"
head "$RESULT"
```
