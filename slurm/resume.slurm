#!/bin/bash -l
#SBATCH --job-name=gbdc_resume_maxiters
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

# --- Fixed run ID so we keep resuming in the same folder ---
RUN_ID="49379"

# --- Binaries and paths ---
BIN=/nfs/home/fgehm/gbdc_f/build/gbdc
RUN_DIR="/nfs/home/fgehm/gbdc_f/run/${RUN_ID}"
OUTDIR="${RUN_DIR}/out"
ERRDIR="${RUN_DIR}/err"
mkdir -p "$OUTDIR" "$ERRDIR"

# --- Exec vars ---
TOTAL_CPUS=${SLURM_CPUS_PER_TASK:-${SLURM_CPUS_ON_NODE:-$(nproc)}}
PARALLEL_JOBS=${PARALLEL_JOBS:-$TOTAL_CPUS}
TIME_LIMIT=${TIME_LIMIT:-12000}

# --- Max-iteration set to continue/run ---
# Adjust as needed; add new values to extend the run.
MAX_ITER=(1 2 3 4 5 6 7 8 9 10 11 12)

# --- Input list (TSV: col1=LOCAL path, col2=ID) ---
LIST=${LIST:-/nfs/home/fgehm/gbdc_f/instances.lst}

# --- Env / deps ---
export PATH="$HOME/local/bin:$PATH"
command -v runsolver >/dev/null || { echo "runsolver not in PATH"; exit 1; }
command -v parallel  >/dev/null || { echo "parallel not in PATH";  exit 1; }
[[ -f "$LIST" ]] || { echo "List file not found at $LIST"; exit 1; }
ulimit -n 4096

CMD_LOG="${RUN_DIR}/cmds.txt"
JOBLOG="${RUN_DIR}/joblog.tsv"

# Do not truncate cmds on resume; create if missing
[[ -e "$CMD_LOG" ]] || : > "$CMD_LOG"

# GNU parallel flags with resume support
PAR_FLAGS="--will-cite --jobs $PARALLEL_JOBS --delay 0 --no-notice --joblog $JOBLOG"
if [[ -s "$JOBLOG" ]]; then
  PAR_FLAGS="$PAR_FLAGS --resume-failed"
  echo "Resuming run in ${RUN_DIR} with --resume-failed"
else
  echo "Starting new run in ${RUN_DIR}"
fi

echo "Generating commands for RUN_ID=${RUN_ID}"
echo "Max-iters set(s): ${MAX_ITER[*]}"

# Generate commands per row and per max-iters: modes {none, sum}
# IMPORTANT: skip header; read TABs; LOCAL=$1, ID=$2
{
  tail -n +2 "$LIST" | while IFS=$'\t' read -r LOCAL ID _; do
    [[ -z "${LOCAL:-}" || -z "${ID:-}" ]] && continue

    for d in "${MAX_ITER[@]}"; do
      suffix="mi${d}"

      # none
      printf 'runsolver -C %d "%s" wlhash --max-iters %d "%s" > "%s/%s_none_%s.out" 2> "%s/%s_none_%s.err"\n' \
             "$TIME_LIMIT" "$BIN" "$d" "$LOCAL" "$OUTDIR" "$ID" "$suffix" "$ERRDIR" "$ID" "$suffix"

      # sum-sort
      printf 'runsolver -C %d "%s" wlhash --sum-sort --max-iters %d "%s" > "%s/%s_sum_%s.out" 2> "%s/%s_sum_%s.err"\n' \
             "$TIME_LIMIT" "$BIN" "$d" "$LOCAL" "$OUTDIR" "$ID" "$suffix" "$ERRDIR" "$ID" "$suffix"
    done
  done
} | tee -a "$CMD_LOG" | parallel $PAR_FLAGS

echo "All jobs completed"
echo "Commands logged in: ${CMD_LOG}"
echo "Artifacts in: ${RUN_DIR}"
echo "Joblog: ${JOBLOG}"
