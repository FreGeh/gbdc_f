#!/bin/bash -l
#SBATCH --job-name=gbdc_wlhyphash_iterations
#SBATCH --nodelist=hoare
#SBATCH --exclusive
#SBATCH --cpus-per-task=80
#SBATCH --ntasks-per-node=1
#SBATCH --mem=0
#SBATCH --time=24:00:00
#SBATCH --chdir=/nfs/home/fgehm/gbdc_f/build
#SBATCH --output=/nfs/home/fgehm/gbdc_f/logs/slurm-%j.out
#SBATCH --error=/nfs/home/fgehm/gbdc_f/logs/slurm-%j.err

set -euo pipefail

# run vars
max_iterations=(4 12)

# folders
BIN=/nfs/home/fgehm/gbdc_f/build/gbdc
RUN_DIR="/nfs/home/fgehm/gbdc_f/run/${SLURM_JOB_ID}"
OUTDIR="${RUN_DIR}/out"
ERRDIR="${RUN_DIR}/err"

# vars for execution
TOTAL_CPUS=${SLURM_CPUS_PER_TASK:-${SLURM_CPUS_ON_NODE:-$(nproc)}}
PARALLEL_JOBS=${PARALLEL_JOBS:-$TOTAL_CPUS}
TIME_LIMIT=${TIME_LIMIT:-6000}
LIST=${LIST:-/nfs/home/fgehm/gbdc_f/instances.lst}
RESULT="${RUN_DIR}/results_${SLURM_JOB_ID}.tsv"
SUMMARY="${RUN_DIR}/hash_summary_${SLURM_JOB_ID}.txt"

export PATH="$HOME/local/bin:$PATH"

which runsolver  >/dev/null || { echo "runsolver not in PATH";  exit 1; }
which parallel   >/dev/null || { echo "parallel not in PATH";   exit 1; }
[[ -f "$LIST" ]] || { echo "List file not found at $LIST"; exit 1; }

mkdir -p "$OUTDIR" "$ERRDIR"

echo "Starting job generation and execution for ${SLURM_JOB_ID}"

# one run per instance per max_iterations
while read -r ID INST; do
  for MI in "${max_iterations[@]}"; do
    OUT="$OUTDIR/${ID}_mi${MI}.out"
    ERR="$ERRDIR/${ID}_mi${MI}.err"
    [[ -s "$OUT" || -s "$ERR" ]] && continue
    printf 'runsolver -C %d "%s" wlhyphash --max-iters %d "%s" > "%s" 2> "%s"\n' \
           "$TIME_LIMIT" "$BIN" "$MI" "$INST" "$OUT" "$ERR"
  done
done < "$LIST" | parallel -j "$PARALLEL_JOBS"

echo "All jobs completed"

# parse results
echo "Starting to parse results"
echo -e "ID\tmax_iters\thash" > "$RESULT"

while read -r ID INST; do
  for MI in "${max_iterations[@]}"; do
    ERR="$ERRDIR/${ID}_mi${MI}.err"
    if [[ -s "$ERR" ]]; then
      # extract the hash line like: c Hash: <value>
      HASH=$(awk '/^c Hash:/{print $3; exit}' "$ERR")
      if [[ -n "${HASH:-}" ]]; then
        echo -e "${ID}\t${MI}\t${HASH}" >> "$RESULT"
      else
        echo -e "${ID}\t${MI}\tNA" >> "$RESULT"
      fi
    else
      echo -e "${ID}\t${MI}\tNA" >> "$RESULT"
    fi
  done
done < "$LIST"

echo "Result parsing finished. Final results are in $RESULT"

# summarize (across all max_iters)
echo "Computing unique hash count"

TOTAL=$(awk 'NR>1{c++} END{print c+0}' "$RESULT")
NA=$(awk -F'\t' 'NR>1 && $3=="NA"{c++} END{print c+0}' "$RESULT")
UNIQUE=$(awk -F'\t' 'NR>1 && $3!="NA"{print $3}' "$RESULT" | LC_ALL=C sort -u | wc -l | awk '{print $1}')
SAME=$(( TOTAL - NA - UNIQUE ))

{
  echo "TOTAL: $TOTAL"
  echo "UNIQUE: $UNIQUE"
  echo "NA: $NA"
  echo "SAME: $SAME"
} > "$SUMMARY"

echo "Summary written to: $SUMMARY"
