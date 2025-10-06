#!/bin/bash -l
#SBATCH --job-name=gbdc_wlhash_iterations
#SBATCH --nodelist=hoare
#SBATCH --exclusive
#SBATCH --cpus-per-task=80
#SBATCH --ntasks-per-node=1
#SBATCH --mem=0
#SBATCH --time=40:00:00
#SBATCH --chdir=/nfs/home/fgehm/gbdc_f/build
#SBATCH --output=/nfs/home/fgehm/gbdc_f/logs/slurm-%j.out
#SBATCH --error=/nfs/home/fgehm/gbdc_f/logs/slurm-%j.err

set -euo pipefail

# iteration caps
max_iterations=(0 1 2 4 6 8 10)

# binaries and paths
BIN=/nfs/home/fgehm/gbdc_f/build/gbdc
RUN_DIR="/nfs/home/fgehm/gbdc_f/run/${SLURM_JOB_ID}"
OUTDIR="${RUN_DIR}/out"
ERRDIR="${RUN_DIR}/err"

# execution vars
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

# Generate one raw and one canonical job per (ID, MI)
while read -r ID INST; do
  for MI in "${max_iterations[@]}"; do
    # RAW
    OUT_RAW="$OUTDIR/${ID}_mi${MI}_raw.out"
    ERR_RAW="$ERRDIR/${ID}_mi${MI}_raw.err"
    if [[ ! -s "$OUT_RAW" && ! -s "$ERR_RAW" ]]; then
      printf 'runsolver -C %d "%s" wlhash --max-iters %d "%s" > "%s" 2> "%s"\n' \
             "$TIME_LIMIT" "$BIN" "$MI" "$INST" "$OUT_RAW" "$ERR_RAW"
    fi

    # CANONICAL
    OUT_CAN="$OUTDIR/${ID}_mi${MI}_canon.out"
    ERR_CAN="$ERRDIR/${ID}_mi${MI}_canon.err"
    if [[ ! -s "$OUT_CAN" && ! -s "$ERR_CAN" ]]; then
      printf 'runsolver -C %d "%s" wlhash --max-iters %d --canonicalize "%s" > "%s" 2> "%s"\n' \
             "$TIME_LIMIT" "$BIN" "$MI" "$INST" "$OUT_CAN" "$ERR_CAN"
    fi
  done
done < "$LIST" | parallel -j "$PARALLEL_JOBS"

echo "All jobs completed"

# Parse results
echo "Starting to parse results"
echo -e "ID\tmax_iters\tmode\thash" > "$RESULT"

extract_hash() {
  # prefer stderr (stats printed to stderr), fallback to stdout
  local err_file="$1"; local out_file="$2"
  local h
  h=$(awk '/^c Hash:/{print $3; exit}' "$err_file" 2>/dev/null || true)
  if [[ -z "${h:-}" ]]; then
    h=$(awk '/^c Hash:/{print $3; exit}' "$out_file" 2>/dev/null || true)
  fi
  if [[ -z "${h:-}" ]]; then echo "NA"; else echo "$h"; fi
}

while read -r ID INST; do
  for MI in "${max_iterations[@]}"; do
    # RAW
    OUT_RAW="$OUTDIR/${ID}_mi${MI}_raw.out"
    ERR_RAW="$ERRDIR/${ID}_mi${MI}_raw.err"
    HASH=$(extract_hash "$ERR_RAW" "$OUT_RAW")
    echo -e "${ID}\t${MI}\traw\t${HASH}" >> "$RESULT"

    # CANON
    OUT_CAN="$OUTDIR/${ID}_mi${MI}_canon.out"
    ERR_CAN="$ERRDIR/${ID}_mi${MI}_canon.err"
    HASH=$(extract_hash "$ERR_CAN" "$OUT_CAN")
    echo -e "${ID}\t${MI}\tcanon\t${HASH}" >> "$RESULT"
  done
done < "$LIST"

echo "Result parsing finished. Final results are in $RESULT"

# Summaries (combined and per mode)
TOTAL=$(awk 'NR>1{c++} END{print c+0}' "$RESULT")
NA=$(awk -F'\t' 'NR>1 && $4=="NA"{c++} END{print c+0}' "$RESULT")
UNIQUE_ALL=$(awk -F'\t' 'NR>1 && $4!="NA"{print $4}' "$RESULT" | LC_ALL=C sort -u | wc -l | awk '{print $1}')

UNIQUE_RAW=$(awk -F'\t' 'NR>1 && $3=="raw"   && $4!="NA"{print $4}' "$RESULT" | LC_ALL=C sort -u | wc -l | awk '{print $1}')
UNIQUE_CAN=$(awk -F'\t' 'NR>1 && $3=="canon" && $4!="NA"{print $4}' "$RESULT" | LC_ALL=C sort -u | wc -l | awk '{print $1}')

MATCH_RAW_CAN=$(awk -F'\t' 'NR>1 {key=$1 FS $2; mode=$3; hash=$4; a[key,mode]=hash}
  END{
    m=0;
    for (k in a) {
      split(k, parts, SUBSEP);
      idmi=parts[1];
      # if both modes exist and both not NA and equal
      if ((idmi, "raw") in a && (idmi, "canon") in a) {
        if (a[idmi, "raw"] != "NA" && a[idmi, "canon"] != "NA" && a[idmi, "raw"]==a[idmi, "canon"]) m++
      }
    }
    print m+0
  }' "$RESULT")

{
  echo "TOTAL rows: $TOTAL"
  echo "NA rows: $NA"
  echo "UNIQUE hashes (all modes): $UNIQUE_ALL"
  echo "UNIQUE hashes (raw):   $UNIQUE_RAW"
  echo "UNIQUE hashes (canon): $UNIQUE_CAN"
  echo "RAW==CANON matching pairs (same ID & max_iters): $MATCH_RAW_CAN"
} > "$SUMMARY"

echo "Summary written to: $SUMMARY"
