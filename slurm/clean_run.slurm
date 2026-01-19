#!/bin/bash -l
#SBATCH --job-name=gbdc_isohash2
#SBATCH --nodelist=liskov
#SBATCH --exclusive
#SBATCH --cpus-per-task=128
#SBATCH --hint=nomultithread
#SBATCH --ntasks-per-node=1
#SBATCH --mem=0
#SBATCH --time=140:00:00
#SBATCH --chdir=/nfs/home/fgehm/isohash2_implementation/build
#SBATCH --output=/nfs/home/fgehm/isohash2_implementation/logs/slurm-%j.out
#SBATCH --error=/nfs/home/fgehm/isohash2_implementation/logs/slurm-%j.err

set -euo pipefail

# binaries and paths
BIN=/nfs/home/fgehm/isohash2_implementation/build/gbdc
RUN_DIR="/nfs/home/fgehm/isohash2_implementation/run/${SLURM_JOB_ID}"
OUTDIR="${RUN_DIR}/out"
ERRDIR="${RUN_DIR}/err"
mkdir -p "$OUTDIR" "$ERRDIR"

# exec vars
TOTAL_CPUS=${SLURM_CPUS_PER_TASK:-${SLURM_CPUS_ON_NODE:-$(nproc)}}
PARALLEL_JOBS=${PARALLEL_JOBS:-$TOTAL_CPUS}
TIME_LIMIT=${TIME_LIMIT:-18000}

MAX_ITER=(1 8)

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

      # default: mate ref, use sum, one pass
      OUT="${OUTDIR}/${ID}_${suffix}.out"; ERR="${ERRDIR}/${ID}_${suffix}.err"
      if [[ ! -s "$OUT" && ! -s "$ERR" ]]; then
        printf 'runsolver -C %d "%s" isohash2 --max-iters %d "%s" > "%s" 2> "%s"\n' \
              "$TIME_LIMIT" "$BIN" "$d" "$LOCAL" "$OUT" "$ERR"
      fi

    done
  done
} | tee -a "$CMD_LOG" | parallel $PAR_FLAGS

echo "All jobs completed"
echo "Commands generated: $(wc -l < "$CMD_LOG")"
echo "Artifacts in: ${RUN_DIR}"
echo "Joblog: ${RUN_DIR}/joblog.tsv"
