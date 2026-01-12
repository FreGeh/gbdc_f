#!/bin/bash -l
#SBATCH --job-name=gbdc_wlhash_run
#SBATCH --nodelist=liskov
#SBATCH --exclusive
#SBATCH --cpus-per-task=128
#SBATCH --ntasks-per-node=1
#SBATCH --mem=0
#SBATCH --time=80:00:00
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

MAX_ITER=(1 2 3 4 5 6 7 8 9 10)

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

    # isohash
    OUT="${OUTDIR}/${ID}_isohash.out"; ERR="${ERRDIR}/${ID}_isohash.err"
    if [[ ! -s "$OUT" && ! -s "$ERR" ]]; then
      printf 'runsolver -C %d "%s" isohash "%s" > "%s" 2> "%s"\n' \
            "$TIME_LIMIT" "$BIN" "$LOCAL" "$OUT" "$ERR"
    fi

    for d in "${MAX_ITER[@]}"; do
      suffix="mi${d}"

      # default: mate ref, use sum, one pass
      OUT="${OUTDIR}/${ID}_onepass_${suffix}.out"; ERR="${ERRDIR}/${ID}_onepass_${suffix}.err"
      if [[ ! -s "$OUT" && ! -s "$ERR" ]]; then
        printf 'runsolver -C %d "%s" wlhash --mate-ref --one-pass --max-iters %d "%s" > "%s" 2> "%s"\n' \
              "$TIME_LIMIT" "$BIN" "$d" "$LOCAL" "$OUT" "$ERR"
      fi

      # wlhash xor - default
      OUT="${OUTDIR}/${ID}_xor_${suffix}.out"; ERR="${ERRDIR}/${ID}_xor_${suffix}.err"
      if [[ ! -s "$OUT" && ! -s "$ERR" ]]; then
        printf 'runsolver -C %d "%s" wlhash --mate-ref --one-pass --use-xor --max-iters %d "%s" > "%s" 2> "%s"\n' \
              "$TIME_LIMIT" "$BIN" "$d" "$LOCAL" "$OUT" "$ERR"
      fi

      # wlhash sort - default
      OUT="${OUTDIR}/${ID}_sort_${suffix}.out"; ERR="${ERRDIR}/${ID}_sort_${suffix}.err"
      if [[ ! -s "$OUT" && ! -s "$ERR" ]]; then
        printf 'runsolver -C %d "%s" wlhash --mate-ref --one-pass --use-sort --max-iters %d "%s" > "%s" 2> "%s"\n' \
              "$TIME_LIMIT" "$BIN" "$d" "$LOCAL" "$OUT" "$ERR"
      fi
    done
  done
} | tee -a "$CMD_LOG" | parallel $PAR_FLAGS

echo "All jobs completed"
echo "Commands generated: $(wc -l < "$CMD_LOG")"
echo "Artifacts in: ${RUN_DIR}"
echo "Joblog: ${RUN_DIR}/joblog.tsv"
