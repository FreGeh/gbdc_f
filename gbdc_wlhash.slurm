#!/bin/bash -l
#SBATCH --job-name=gbdc_wlhash
#SBATCH --nodelist=hoare
#SBATCH --exclusive
#SBATCH --cpus-per-task=20
#SBATCH --mem=0
#SBATCH --time=72:00:00
#SBATCH --chdir=/nfs/home/fgehm/gbdc_f/build
#SBATCH --output=/nfs/home/fgehm/gbdc_f/logs/slurm-%j.out
#SBATCH --error=/nfs/home/fgehm/gbdc_f/logs/slurm-%j.err

# variablen
PARALLEL_JOBS=${PARALLEL_JOBS:-160}
TIME_LIMIT=${TIME_LIMIT:-120}
LIST=${LIST:-/nfs/home/fgehm/gbdc_f/instances.lst}
RESULT=/nfs/home/fgehm/gbdc_f/wlhash_results_variable.tsv

BIN=/nfs/home/fgehm/gbdc_f/build/gbdc
OUTDIR=/nfs/home/fgehm/gbdc_f/out_variable
ERRDIR=/nfs/home/fgehm/gbdc_f/err_variable

export PATH=$HOME/local/bin:$PATH
which runsolver  || { echo "runsolver not in PATH";  exit 1; }
which parallel   || { echo "parallel  not in PATH";  exit 1; }

mkdir -p "$OUTDIR" "$ERRDIR"

while read -r ID INST; do
    OUT="$OUTDIR/${ID}.out"
    ERR="$ERRDIR/${ID}.err"
    if [[ ! -s "$OUT" ]]; then
        echo "runsolver -C $TIME_LIMIT \"$BIN\" wlhash --k 2 \"$INST\" > \"$OUT\" 2> \"$ERR\""
    fi
done < "$LIST" | parallel -j "$PARALLEL_JOBS"

# Merge
echo "Merging WL-Hashes into $RESULT …"
paste \
  <(awk '{print $1}' "$LIST") \
  <(for f in "$ERRDIR"/*.out; do
        grep -m1 -Eo '^[0-9a-f]{16}$' "$f"
    done) \
  > "$RESULT"
echo "Done — see $RESULT"
