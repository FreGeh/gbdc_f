#!/bin/bash -l
#SBATCH --job-name=gbdc_wlhash
#SBATCH --nodelist=hoare
#SBATCH --exclusive
#SBATCH --cpus-per-task=160          # all cores on hoare
#SBATCH --ntasks=1                   # we drive parallel ourselves
#SBATCH --mem=0                      # full node memory
#SBATCH --time=30:00:00
#SBATCH --chdir=/nfs/home/fgehm/gbdc_f/build
#SBATCH --output=/nfs/home/fgehm/gbdc_f/logs/slurm-%j.out
#SBATCH --error=/nfs/home/fgehm/gbdc_f/logs/slurm-%j.err

set -euo pipefail

############################
# Adjustable parameters
############################
K_LIST=(1)
MAX_ITERS_LIST=(1 2 3 4 5 6 7 8 9 10)       # 1000 = “unlimited” for your study
PARALLEL_JOBS=${PARALLEL_JOBS:-140}          # leave some headroom on a 160-core node
TIME_LIMIT=${TIME_LIMIT:-9600}               # CPU seconds (soft)
WALL_LIMIT=${WALL_LIMIT:-18000}               # Wall seconds (hard)
LIST=${LIST:-/nfs/home/fgehm/gbdc_f/instances.lst}

BIN=/nfs/home/fgehm/gbdc_f/build/gbdc
OUTDIR=/nfs/home/fgehm/gbdc_f/out_variable_new3
ERRDIR=/nfs/home/fgehm/gbdc_f/err_variable_new3
LOGDIR=/nfs/home/fgehm/gbdc_f/logs

HASH_REGEX='[0-9a-fA-F]{16}'
JOBLOG="$LOGDIR/parallel_joblog-${SLURM_JOB_ID}.tsv"
FAILLOG="$LOGDIR/failures-${SLURM_JOB_ID}.log"
CMDFILE="$LOGDIR/cmds-${SLURM_JOB_ID}.txt"

mkdir -p "$OUTDIR" "$ERRDIR" "$LOGDIR"

# Check tools
export PATH="$HOME/local/bin:$PATH"
which runsolver >/dev/null || { echo "runsolver not in PATH"; exit 1; }
which parallel  >/dev/null || { echo "parallel not in PATH";  exit 1; }

############################
# 1) Generate command list
############################
: > "$CMDFILE"
: > "$FAILLOG"

while read -r ID INST; do
  [[ -z "$ID" ]] && continue
  for k in "${K_LIST[@]}"; do
    for it in "${MAX_ITERS_LIST[@]}"; do
      OUT="$OUTDIR/${ID}_k${k}_i${it}.out"
      ERR="$ERRDIR/${ID}_k${k}_i${it}.err"

      # Skip if OUT already contains a 16-hex WLHash
      if [[ -s "$OUT" ]] && grep -qE "$HASH_REGEX" "$OUT"; then
        continue
      fi

      # One-line command for GNU parallel
      printf 'runsolver -C %d -W %d "%s" wlhash --k %d --max-iters %d "%s" > "%s" 2> "%s" || echo -e "%s\t%d\t%d\t$?" >> "%s"\n' \
        "$TIME_LIMIT" "$WALL_LIMIT" "$BIN" "$k" "$it" "$INST" "$OUT" "$ERR" \
        "$ID" "$k" "$it" "$FAILLOG" \
        >> "$CMDFILE"
    done
  done
done < "$LIST"

echo "Commands to run: $(wc -l < "$CMDFILE")"

############################
# 2) Run commands
############################
parallel --joblog "$JOBLOG" -j "$PARALLEL_JOBS" < "$CMDFILE"

echo "All commands finished."
echo "Joblog:     $JOBLOG"
echo "Failures:   $FAILLOG (empty if none)"
echo "Outputs in: $OUTDIR / $ERRDIR"
