#!/bin/bash -l
#SBATCH --job-name=gbdc_wlhash_run
#SBATCH --nodelist=hoare
#SBATCH --exclusive
#SBATCH --cpus-per-task=80
#SBATCH --ntasks-per-node=1
#SBATCH --mem=0
#SBATCH --time=58:00:00
#SBATCH --chdir=/nfs/home/fgehm/gbdc_f/build
#SBATCH --output=/nfs/home/fgehm/gbdc_f/logs/slurm-%j.out
#SBATCH --error=/nfs/home/fgehm/gbdc_f/logs/slurm-%j.err

set -euo pipefail

# binaries and paths
BIN=/nfs/home/fgehm/gbdc_f/build/gbdc
RUN_DIR="/nfs/home/fgehm/gbdc_f/run/${SLURM_JOB_ID}"
OUTDIR="${RUN_DIR}/out"
ERRDIR="${RUN_DIR}/err"
mkdir -p "$OUTDIR" "$ERRDIR"

# exec vars
TOTAL_CPUS=${SLURM_CPUS_PER_TASK:-${SLURM_CPUS_ON_NODE:-$(nproc)}}
PARALLEL_JOBS=${PARALLEL_JOBS:-$TOTAL_CPUS}
TIME_LIMIT=${TIME_LIMIT:-12000}

# NEW: headered TSV with many columns (col1=local, col2=hash/ID)
LIST=${LIST:-/nfs/home/fgehm/gbdc_f/instances.lst}

# env / deps
export PATH="$HOME/local/bin:$PATH"
command -v runsolver >/dev/null || { echo "runsolver not in PATH"; exit 1; }
command -v parallel  >/dev/null || { echo "parallel not in PATH";  exit 1; }
[[ -f "$LIST" ]] || { echo "List file not found at $LIST"; exit 1; }
ulimit -n 4096

PAR_FLAGS="--will-cite --jobs $PARALLEL_JOBS --delay 0 --no-notice --joblog ${RUN_DIR}/joblog.tsv"
CMD_LOG="${RUN_DIR}/cmds.txt"; : > "$CMD_LOG"

echo "Starting job generation and execution for ${SLURM_JOB_ID}"

# Generate commands per row: none | adv | sum
# IMPORTANT: skip header; read TABs; LOCAL=$1, ID=$2
{
  tail -n +2 "$LIST" | while IFS=$'\t' read -r LOCAL ID _; do
    [[ -z "${LOCAL:-}" || -z "${ID:-}" ]] && continue

    # 1) none
    OUT="${OUTDIR}/${ID}_none.out"; ERR="${ERRDIR}/${ID}_none.err"
    if [[ ! -s "$OUT" && ! -s "$ERR" ]]; then
      printf 'runsolver -C %d "%s" wlhash "%s" > "%s" 2> "%s"\n' \
             "$TIME_LIMIT" "$BIN" "$LOCAL" "$OUT" "$ERR"
    fi

    # 2) advanced-sort
    OUT="${OUTDIR}/${ID}_adv.out"; ERR="${ERRDIR}/${ID}_adv.err"
    if [[ ! -s "$OUT" && ! -s "$ERR" ]]; then
      printf 'runsolver -C %d "%s" wlhash --advanced-sort "%s" > "%s" 2> "%s"\n' \
             "$TIME_LIMIT" "$BIN" "$LOCAL" "$OUT" "$ERR"
    fi

    # 3) sum-sort
    OUT="${OUTDIR}/${ID}_sum.out"; ERR="${ERRDIR}/${ID}_sum.err"
    if [[ ! -s "$OUT" && ! -s "$ERR" ]]; then
      printf 'runsolver -C %d "%s" wlhash --sum-sort "%s" > "%s" 2> "%s"\n' \
             "$TIME_LIMIT" "$BIN" "$LOCAL" "$OUT" "$ERR"
    fi
  done
} | tee -a "$CMD_LOG" | parallel $PAR_FLAGS

echo "All jobs completed"
echo "Commands generated: $(wc -l < "$CMD_LOG")"
echo "Artifacts in: ${RUN_DIR}"
echo "Joblog: ${RUN_DIR}/joblog.tsv"
