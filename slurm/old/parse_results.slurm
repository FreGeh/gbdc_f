#!/bin/bash -l
#SBATCH --job-name=gbdc_parsing
#SBATCH --nodelist=hoare
#SBATCH --exclusive
#SBATCH --cpus-per-task=80
#SBATCH --ntasks-per-node=1
#SBATCH --mem=0
#SBATCH --time=2:00:00

set -euo pipefail

############################
# Configure here
JOB_ID="49122"
BASE_RUN="/nfs/home/fgehm/gbdc_f/run"
OUT_SUB="out"
ERR_SUB="err"
RESULT_NAME="results.tsv"
############################

RUN_DIR="${BASE_RUN}/${JOB_ID}"
OUTDIR="${RUN_DIR}/${OUT_SUB}"
ERRDIR="${RUN_DIR}/${ERR_SUB}"
RESULT="${RUN_DIR}/${RESULT_NAME}"

[[ -d "$OUTDIR" && -d "$ERRDIR" ]] || { echo "Missing out or err dir in ${RUN_DIR}" >&2; exit 3; }

# header
printf "ID\tmax_iters\tmode\thash\tcpu_time_s\tmax_mem_kib\n" > "$RESULT"

# prefer stderr for the hash, fallback to stdout
extract_hash() {
  local err_file="$1" out_file="$2" h=""
  if [[ -s "$err_file" ]]; then
    h=$(awk '/^c[[:space:]]+Hash:/{print $3; exit}' "$err_file" || true)
  fi
  if [[ -z "${h:-}" && -s "$out_file" ]]; then
    h=$(awk '/^c[[:space:]]+Hash:/{print $3; exit}' "$out_file" || true)
  fi
  [[ -z "${h:-}" ]] && echo "NA" || echo "$h"
}

shopt -s nullglob
for out_file in "${OUTDIR}"/*_mi*_*.out; do
  base=$(basename "$out_file")
  # expected: <ID>_mi<iters>_(raw|canon).out
  if [[ "$base" =~ ^(.+)_mi([0-9]+)_(raw|canon)\.out$ ]]; then
    id="${BASH_REMATCH[1]}"
    mi="${BASH_REMATCH[2]}"
    mode="${BASH_REMATCH[3]}"
    err_file="${ERRDIR}/${id}_mi${mi}_${mode}.err"

    hash_val=$(extract_hash "$err_file" "$out_file")

    # from runsolver stdout
    cpu_time=$(awk -F': ' '/^CPU time \(s\):/{print $2; exit}' "$out_file" || true)
    [[ -z "${cpu_time:-}" ]] && cpu_time="NA"

    max_mem=$(awk -F': ' '/^Max\.[[:space:]]+memory.*\(KiB\):/{print $2; exit}' "$out_file" || true)
    [[ -z "${max_mem:-}" ]] && max_mem="NA"

    printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$id" "$mi" "$mode" "$hash_val" "$cpu_time" "$max_mem" >> "$RESULT"
  fi
done

echo "Wrote ${RESULT}"
