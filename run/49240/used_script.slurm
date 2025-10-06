#!/bin/bash -l
#SBATCH --job-name=gbdc_wlhash_run
#SBATCH --nodelist=hoare
#SBATCH --exclusive
#SBATCH --cpus-per-task=80
#SBATCH --ntasks-per-node=1
#SBATCH --mem=0
#SBATCH --time=60:00:00
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

# execution vars
TOTAL_CPUS=${SLURM_CPUS_PER_TASK:-${SLURM_CPUS_ON_NODE:-$(nproc)}}
PARALLEL_JOBS=${PARALLEL_JOBS:-$TOTAL_CPUS}
TIME_LIMIT=${TIME_LIMIT:-6000}

# NOTE: keep the same two-column LIST format as your old script: "ID INST"
# If you want the first 1000, point LIST to that file; otherwise to the full list.
LIST=${LIST:-/nfs/home/fgehm/gbdc_f/instances_first_1000.lst}

# ensure runsolver/parallel as before
export PATH="$HOME/local/bin:$PATH"
which runsolver  >/dev/null || { echo "runsolver not in PATH";  exit 1; }
which parallel   >/dev/null || { echo "parallel not in PATH";   exit 1; }
[[ -f "$LIST" ]] || { echo "List file not found at $LIST"; exit 1; }

echo "Starting job generation and execution for ${SLURM_JOB_ID}"

# Generate 4 modes per (ID, INST): none | neigh | adv | sum
# Pipe the commands directly to parallel (exactly like your old script)
while read -r ID INST; do
  # 1) none
  OUT_NONE="$OUTDIR/${ID}_none.out"
  ERR_NONE="$ERRDIR/${ID}_none.err"
  if [[ ! -s "$OUT_NONE" && ! -s "$ERR_NONE" ]]; then
    printf 'runsolver -C %d "%s" wlhash "%s" > "%s" 2> "%s"\n' \
           "$TIME_LIMIT" "$BIN" "$INST" "$OUT_NONE" "$ERR_NONE"
  fi

  # 2) neigh-hop
  OUT_NEIGH="$OUTDIR/${ID}_neigh.out"
  ERR_NEIGH="$ERRDIR/${ID}_neigh.err"
  if [[ ! -s "$OUT_NEIGH" && ! -s "$ERR_NEIGH" ]]; then
    printf 'runsolver -C %d "%s" wlhash --neigh-hop "%s" > "%s" 2> "%s"\n' \
           "$TIME_LIMIT" "$BIN" "$INST" "$OUT_NEIGH" "$ERR_NEIGH"
  fi

  # 3) advanced-sort
  OUT_ADV="$OUTDIR/${ID}_adv.out"
  ERR_ADV="$ERRDIR/${ID}_adv.err"
  if [[ ! -s "$OUT_ADV" && ! -s "$ERR_ADV" ]]; then
    printf 'runsolver -C %d "%s" wlhash --advanced-sort "%s" > "%s" 2> "%s"\n' \
           "$TIME_LIMIT" "$BIN" "$INST" "$OUT_ADV" "$ERR_ADV"
  fi

  # 4) sum-sort
  OUT_SUM="$OUTDIR/${ID}_sum.out"
  ERR_SUM="$ERRDIR/${ID}_sum.err"
  if [[ ! -s "$OUT_SUM" && ! -s "$ERR_SUM" ]]; then
    printf 'runsolver -C %d "%s" wlhash --sum-sort "%s" > "%s" 2> "%s"\n' \
           "$TIME_LIMIT" "$BIN" "$INST" "$OUT_SUM" "$ERR_SUM"
  fi
done < "$LIST" | parallel -j "$PARALLEL_JOBS"

echo "All jobs completed"
echo "Artifacts in: ${RUN_DIR}"
