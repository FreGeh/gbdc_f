#!/bin/bash -l
#SBATCH --job-name=gbdc_iso_wl_wlhyp
#SBATCH --nodelist=hoare
#SBATCH --exclusive
#SBATCH --cpus-per-task=80
#SBATCH --ntasks-per-node=1
#SBATCH --mem=0
#SBATCH --time=24:00:00
#SBATCH --chdir=/nfs/home/fgehm/gbdc_f/build
#SBATCH --output=/nfs/home/fgehm/gbdc_f/logs/slurm-%j.out
#SBATCH --error=/nfs/home/fgehm/gbdc_f/logs/slurm-%j.err

# folders
BIN=/nfs/home/fgehm/gbdc_f/build/gbdc
OUTDIR="/nfs/home/fgehm/gbdc_f/run/${SLURM_JOB_ID}/out"
ERRDIR="/nfs/home/fgehm/gbdc_f/run/${SLURM_JOB_ID}/err"

# vars for execution
TOTAL_CPUS=${SLURM_CPUS_PER_TASK:-${SLURM_CPUS_ON_NODE:-$(nproc)}}
PARALLEL_JOBS=${PARALLEL_JOBS:-${SLURM_CPUS_PER_TASK}}
TIME_LIMIT=${TIME_LIMIT:-600}
# use the first-1000 list
LIST=${LIST:-/nfs/home/fgehm/gbdc_f/instances_first_1000.lst}

# per-tool results + summaries
RESULT_ISO="/nfs/home/fgehm/gbdc_f/run/${SLURM_JOB_ID}/results_isohash_${SLURM_JOB_ID}.tsv"
RESULT_WL="/nfs/home/fgehm/gbdc_f/run/${SLURM_JOB_ID}/results_wlhash_${SLURM_JOB_ID}.tsv"
RESULT_WLHYP="/nfs/home/fgehm/gbdc_f/run/${SLURM_JOB_ID}/results_wlhyphash_${SLURM_JOB_ID}.tsv"

SUMMARY_ISO="/nfs/home/fgehm/gbdc_f/run/${SLURM_JOB_ID}/hash_summary_isohash_${SLURM_JOB_ID}.txt"
SUMMARY_WL="/nfs/home/fgehm/gbdc_f/run/${SLURM_JOB_ID}/hash_summary_wlhash_${SLURM_JOB_ID}.txt"
SUMMARY_WLHYP="/nfs/home/fgehm/gbdc_f/run/${SLURM_JOB_ID}/hash_summary_wlhyphash_${SLURM_JOB_ID}.txt"

# wl parameters (kept minimal)
MAX_ITERS_LIST=(10)
K_LIST=(1)

export PATH=$HOME/local/bin:$PATH
which runsolver  || { echo "runsolver not in PATH";  exit 1; }
which parallel   || { echo "parallel  not in PATH";  exit 1; }

mkdir -p "$OUTDIR" "$ERRDIR"

echo "Starting job generation and execution for ${SLURM_JOB_ID}"
echo "Using list: $LIST"

# ---- generate commands for all three tools; pipe to GNU parallel ----
while read -r ID INST; do
  # isohash
  OUT_ISO="$OUTDIR/${ID}.iso.out"
  ERR_ISO="$ERRDIR/${ID}.iso.err"
  [[ -s "$OUT_ISO" ]] || printf 'runsolver -C %d "%s" isohash "%s" > "%s" 2> "%s"\n' \
                         "$TIME_LIMIT" "$BIN" "$INST" "$OUT_ISO" "$ERR_ISO"

  # wlhash (k / max-iters grid kept for consistency)
  for k in "${K_LIST[@]}"; do
    for it in "${MAX_ITERS_LIST[@]}"; do
      OUT_WL="$OUTDIR/${ID}_k${k}_i${it}.wl.out"
      ERR_WL="$ERRDIR/${ID}_k${k}_i${it}.wl.err"
      [[ -s "$OUT_WL" ]] || printf 'runsolver -C %d "%s" wlhash --k %d --max-iters %d "%s" > "%s" 2> "%s"\n' \
                           "$TIME_LIMIT" "$BIN" "$k" "$it" "$INST" "$OUT_WL" "$ERR_WL"
    done
  done

  # wlhyphash (no k/iters needed, but keep naming for alignment)
  for k in "${K_LIST[@]}"; do
    for it in "${MAX_ITERS_LIST[@]}"; do
      OUT_HYP="$OUTDIR/${ID}_k${k}_i${it}.wlhyp.out"
      ERR_HYP="$ERRDIR/${ID}_k${k}_i${it}.wlhyp.err"
      [[ -s "$OUT_HYP" ]] || printf 'runsolver -C %d "%s" wlhyphash "%s" > "%s" 2> "%s"\n' \
                             "$TIME_LIMIT" "$BIN" "$INST" "$OUT_HYP" "$ERR_HYP"
    done
  done
done < "$LIST" | parallel -j "${PARALLEL_JOBS}"

echo "All jobs completed"

# =========================
# parse: ISOHASH
# =========================
echo "Parsing isohash results"
echo -e "ID\tisohash" > "$RESULT_ISO"

while read -r ID INST; do
  ERR_ISO="$ERRDIR/${ID}.iso.err"
  if [[ -s "$ERR_ISO" ]]; then
    HASH=$(grep "c Hash:" "$ERR_ISO" | awk '{print $3}')
    [[ -n "$HASH" ]] || HASH=NA
    echo -e "${ID}\t${HASH}" >> "$RESULT_ISO"
  else
    echo -e "${ID}\tNA" >> "$RESULT_ISO"
  fi
done < "$LIST"

# summary: ISOHASH
TOTAL=$(awk 'NR>1{c++} END{print c+0}' "$RESULT_ISO")
UNIQUE=$(awk -F'\t' 'NR>1 && $2!="NA"{print $2}' "$RESULT_ISO" | LC_ALL=C sort -u | wc -l | awk '{print $1}')
NA=$(awk -F'\t' 'NR>1 && $2=="NA"{c++} END{print c+0}' "$RESULT_ISO")
SAME=$(( TOTAL - NA - UNIQUE ))
{
  echo "TOTAL: $TOTAL"
  echo "UNIQUE: $UNIQUE"
  echo "NA: $NA"
  echo "SAME: $SAME"
} > "$SUMMARY_ISO"
echo "Summary (isohash): $SUMMARY_ISO"

# =========================
# parse: WLHASH
# =========================
echo "Parsing wlhash results"
echo -e "ID\tk\tmax_iters\twlhash" > "$RESULT_WL"

while read -r ID INST; do
  for k in "${K_LIST[@]}"; do
    for it in "${MAX_ITERS_LIST[@]}"; do
      ERR_WL="$ERRDIR/${ID}_k${k}_i${it}.wl.err"
      if [[ -s "$ERR_WL" ]]; then
        HASH=$(grep "c Hash:" "$ERR_WL" | awk '{print $3}')
        [[ -n "$HASH" ]] || HASH=NA
        echo -e "${ID}\t${k}\t${it}\t${HASH}" >> "$RESULT_WL"
      else
        echo -e "${ID}\t${k}\t${it}\tNA" >> "$RESULT_WL"
      fi
    done
  done
done < "$LIST"

# summary: WLHASH
TOTAL=$(awk 'NR>1{c++} END{print c+0}' "$RESULT_WL")
UNIQUE=$(awk -F'\t' 'NR>1 && $4!="NA"{print $4}' "$RESULT_WL" | LC_ALL=C sort -u | wc -l | awk '{print $1}')
NA=$(awk -F'\t' 'NR>1 && $4=="NA"{c++} END{print c+0}' "$RESULT_WL")
SAME=$(( TOTAL - NA - UNIQUE ))
{
  echo "TOTAL: $TOTAL"
  echo "UNIQUE: $UNIQUE"
  echo "NA: $NA"
  echo "SAME: $SAME"
} > "$SUMMARY_WL"
echo "Summary (wlhash): $SUMMARY_WL"

# =========================
# parse: WLHYPHASH
# =========================
echo "Parsing wlhyphash results"
echo -e "ID\tk\tmax_iters\twlhyphash" > "$RESULT_WLHYP"

while read -r ID INST; do
  for k in "${K_LIST[@]}"; do
    for it in "${MAX_ITERS_LIST[@]}"; do
      ERR_HYP="$ERRDIR/${ID}_k${k}_i${it}.wlhyp.err"
      if [[ -s "$ERR_HYP" ]]; then
        HASH=$(grep "c Hash:" "$ERR_HYP" | awk '{print $3}')
        [[ -n "$HASH" ]] || HASH=NA
        echo -e "${ID}\t${k}\t${it}\t${HASH}" >> "$RESULT_WLHYP"
      else
        echo -e "${ID}\t${k}\t${it}\tNA" >> "$RESULT_WLHYP"
      fi
    done
  done
done < "$LIST"

# summary: WLHYPHASH
TOTAL=$(awk 'NR>1{c++} END{print c+0}' "$RESULT_WLHYP")
UNIQUE=$(awk -F'\t' 'NR>1 && $4!="NA"{print $4}' "$RESULT_WLHYP" | LC_ALL=C sort -u | wc -l | awk '{print $1}')
NA=$(awk -F'\t' 'NR>1 && $4=="NA"{c++} END{print c+0}' "$RESULT_WLHYP")
SAME=$(( TOTAL - NA - UNIQUE ))
{
  echo "TOTAL: $TOTAL"
  echo "UNIQUE: $UNIQUE"
  echo "NA: $NA"
  echo "SAME: $SAME"
} > "$SUMMARY_WLHYP"
echo "Summary (wlhyphash): $SUMMARY_WLHYP"

echo "All summaries written. Per-tool results:"
echo "  $RESULT_ISO"
echo "  $RESULT_WL"
echo "  $RESULT_WLHYP"
