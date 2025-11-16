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
TIME_LIMIT=${TIME_LIMIT:-6000}

MAX_ITER=(1 2 4 6 8)

LIST=${LIST:-/nfs/home/fgehm/gbdc_f/instances_new.lst}

# env / deps
export PATH="$HOME/local/bin:$PATH"
command -v runsolver >/dev/null || { echo "runsolver not in PATH"; exit 1; }
command -v parallel  >/dev/null || { echo "parallel not in PATH";  exit 1; }
[[ -f "$LIST" ]] || { echo "List file not found at $LIST"; exit 1; }
ulimit -n 4096

PAR_FLAGS="--will-cite --jobs $PARALLEL_JOBS --delay 0 --no-notice --joblog ${RUN_DIR}/joblog.tsv"
CMD_LOG="${RUN_DIR}/cmds.txt"; : > "$CMD_LOG"

echo "Starting job generation and execution for ${SLURM_JOB_ID}"
echo "Max-iters set(s): ${MAX_ITER[*]}"

{
  tail -n +2 "$LIST" | while IFS=$'\t' read -r LOCAL ID _; do
    [[ -z "${LOCAL:-}" || -z "${ID:-}" ]] && continue

    for d in "${MAX_ITER[@]}"; do
      suffix="mi${d}"

      # 1) timon
      OUT="${OUTDIR}/${ID}_timon_${suffix}.out"; ERR="${ERRDIR}/${ID}_timon_${suffix}.err"
      if [[ ! -s "$OUT" && ! -s "$ERR" ]]; then
        printf 'runsolver -C %d "%s" timonhash --max-iters %d --no-sort "%s" > "%s" 2> "%s"\n' \
              "$TIME_LIMIT" "$BIN" "$d" "$LOCAL" "$OUT" "$ERR"
      fi

      OUT="${OUTDIR}/${ID}_xor_${suffix}.out"; ERR="${ERRDIR}/${ID}_xor_${suffix}.err"
      if [[ ! -s "$OUT" && ! -s "$ERR" ]]; then
        printf 'runsolver -C %d "%s" wlhash --use-xor --max-iters %d "%s" > "%s" 2> "%s"\n' \
              "$TIME_LIMIT" "$BIN" "$d" "$LOCAL" "$OUT" "$ERR"
      fi

      OUT="${OUTDIR}/${ID}_sort_${suffix}.out"; ERR="${ERRDIR}/${ID}_sort_${suffix}.err"
      if [[ ! -s "$OUT" && ! -s "$ERR" ]]; then
        printf 'runsolver -C %d "%s" wlhash --use-sort --max-iters %d "%s" > "%s" 2> "%s"\n' \
              "$TIME_LIMIT" "$BIN" "$d" "$LOCAL" "$OUT" "$ERR"
      fi

      OUT="${OUTDIR}/${ID}_half_${suffix}.out"; ERR="${ERRDIR}/${ID}_half_${suffix}.err"
      if [[ ! -s "$OUT" && ! -s "$ERR" ]]; then
        printf 'runsolver -C %d "%s" wlhash --half-bit --max-iters %d "%s" > "%s" 2> "%s"\n' \
              "$TIME_LIMIT" "$BIN" "$d" "$LOCAL" "$OUT" "$ERR"
      fi
    done
  done
} | tee -a "$CMD_LOG" | parallel $PAR_FLAGS

echo "All jobs completed"
echo "Commands generated: $(wc -l < "$CMD_LOG")"
echo "Artifacts in: ${RUN_DIR}"
echo "Joblog: ${RUN_DIR}/joblog.tsv"
