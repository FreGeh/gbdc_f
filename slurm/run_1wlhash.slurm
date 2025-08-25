#!/bin/bash -l
#SBATCH --job-name=gbdc_wlhash_k1
#SBATCH --nodelist=hoare
#SBATCH --exclusive
#SBATCH --cpus-per-task=80
#SBATCH --ntasks-per-node=1
#SBATCH --mem=0
#SBATCH --time=72:00:00
#SBATCH --chdir=/nfs/home/fgehm/gbdc_f/build
#SBATCH --output=/nfs/home/fgehm/gbdc_f/logs/slurm-%j.out
#SBATCH --error=/nfs/home/fgehm/gbdc_f/logs/slurm-%j.err

# folders
BIN=/nfs/home/fgehm/gbdc_f/build/gbdc
OUTDIR="/nfs/home/fgehm/gbdc_f/run/${SLURM_JOB_ID}/out"
ERRDIR="/nfs/home/fgehm/gbdc_f/run/${SLURM_JOB_ID}/err"

# vars for execution
TOTAL_CPUS=${SLURM_CPUS_PER_TASK:-${SLURM_CPUS_ON_NODE:-$(nproc)}}
PARALLEL_JOBS=${PARALLEL_JOBS:-$TOTAL_CPUS}
TIME_LIMIT=${TIME_LIMIT:-600}
LIST=${LIST:-/nfs/home/fgehm/gbdc_f/instances.lst}
RESULT="/nfs/home/fgehm/gbdc_f/run/${SLURM_JOB_ID}/results_${SLURM_JOB_ID}.tsv"
UNIQUE_SUMMARY_FILE="/nfs/home/fgehm/gbdc_f/run/${SLURM_JOB_ID}/hash_summary_${SLURM_JOB_ID}.txt"

# run vars
MAX_ITERS_LIST=(1 2 3 4 5 6 7 8 9 10 15 25 50 100)
K_LIST=(1)

export PATH=$HOME/local/bin:$PATH
which runsolver  || { echo "runsolver not in PATH";  exit 1; }
which parallel   || { echo "parallel  not in PATH";  exit 1; }

mkdir -p "$OUTDIR" "$ERRDIR"

# run jobs

echo "Starting job generation and execution for ${SLURM_JOB_ID}"

## for each instance, for each k, for each iteration run the command
while read -r ID INST; do
  for k in "${K_LIST[@]}"; do
    for it in "${MAX_ITERS_LIST[@]}"; do
      OUT="$OUTDIR/${ID}_k${k}_i${it}.out"
      ERR="$ERRDIR/${ID}_k${k}_i${it}.err"
      [[ -s "$OUT" ]] && continue
      printf 'runsolver -C %d "%s" wlhash --k %d --max-iters %d "%s" > "%s" 2> "%s"\n' \
             "$TIME_LIMIT" "$BIN" "$k" "$it" "$INST" "$OUT" "$ERR"
    done
  done
done < "$LIST" | parallel -j "$PARALLEL_JOBS"

echo "All jobs completed"

# -- parse results --

echo "Starting to parse results"

## nice header
echo -e "ID\tk\tmax_iters\twlhash" > "$RESULT"

## reading all the hashes and parsing into result file
while read -r ID INST; do
  for k in "${K_LIST[@]}"; do
    for it in "${MAX_ITERS_LIST[@]}"; do
      ERR="$ERRDIR/${ID}_k${k}_i${it}.err"
      # Check if the output file exists and is not empty
      if [[ -s "$ERR" ]]; then
        # Read hash
        HASH=$(grep "c Hash:" "$ERR" | awk '{print $3}')
        # Write hash into correct result file, if not found write NA
        if [[ -n "$HASH" ]]; then
            echo -e "${ID}\t${k}\t${it}\t${HASH}" >> "$RESULT"
        else
            echo -e "${ID}\t${k}\t${it}\tNA" >> "$RESULT"
        fi
      fi
    done
  done
done < "$LIST"

echo "Result parsing finished. Final results are in $RESULT"

## analyzing hashes

echo "Computing unique hash count"

# overall stats: use column 4 for the hash
TOTAL=$(awk 'NR>1{c++} END{print c+0}' "$RESULT")
NA=$(awk -F'\t' 'NR>1 && $4=="NA"{c++} END{print c+0}' "$RESULT")
UNIQUE=$(awk -F'\t' 'NR>1 && $4!="NA"{print $4}' "$RESULT" | LC_ALL=C sort -u | wc -l | awk '{print $1}')
SAME=$(( TOTAL - NA - UNIQUE ))

{
  echo "TOTAL: $TOTAL"
  echo "UNIQUE: $UNIQUE"
  echo "NA: $NA"
  echo "SAME: $SAME"
  echo ""
  echo "# UNIQUE hashes by k"
} > "$UNIQUE_SUMMARY_FILE"

# UNIQUE per k
awk -F'\t' 'NR>1 && $4!="NA"{print $2, $4}' "$RESULT" \
| awk '{seen[$1"|" $2]=1} END{for (k in seen){split(k,a,"|"); cnt[a[1]]++} for (k in cnt) print k, cnt[k]}' \
| LC_ALL=C sort -n -k1,1 \
| awk '{printf "k=%s\tUNIQUE=%s\n",$1,$2}' >> "$UNIQUE_SUMMARY_FILE"

echo "" >> "$UNIQUE_SUMMARY_FILE"
echo "# UNIQUE hashes by k and max_iters" >> "$UNIQUE_SUMMARY_FILE"

# UNIQUE per (k, max_iters)
awk -F'\t' 'NR>1 && $4!="NA"{print $2, $3, $4}' "$RESULT" \
| awk '{seen[$1"|" $2"|" $3]=1} END{for (key in seen){split(key,a,"|"); combo=a[1] "\t" a[2]; cnt[combo]++} for (c in cnt) print c "\t" cnt[c]}' \
| LC_ALL=C sort -n -k1,1 -k2,2 \
| awk -F'\t' '{printf "k=%s\tmax_iters=%s\tUNIQUE=%s\n",$1,$2,$3}' >> "$UNIQUE_SUMMARY_FILE"

echo "" >> "$UNIQUE_SUMMARY_FILE"
echo "# NA counts by k and max_iters" >> "$UNIQUE_SUMMARY_FILE"

# NA per (k, max_iters)
awk -F'\t' 'NR>1 && $4=="NA"{na[$2"\t"$3]++} END{for (c in na) print c "\t" na[c]}' "$RESULT" \
| LC_ALL=C sort -n -k1,1 -k2,2 \
| awk -F'\t' '{printf "k=%s\tmax_iters=%s\tNA=%s\n",$1,$2,$3}' >> "$UNIQUE_SUMMARY_FILE"

echo "Summary written to: $UNIQUE_SUMMARY_FILE"