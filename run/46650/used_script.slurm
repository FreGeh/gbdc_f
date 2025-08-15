#!/bin/bash -l
#SBATCH --job-name=gbdc_wlhash_k1
#SBATCH --nodelist=hoare
#SBATCH --exclusive
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=20
#SBATCH --mem=0
#SBATCH --time=72:00:00
#SBATCH --chdir=/nfs/home/fgehm/gbdc_f/build
#SBATCH --output=/nfs/home/fgehm/gbdc_f/logs/slurm-%j.out
#SBATCH --error=/nfs/home/fgehm/gbdc_f/logs/slurm-%j.err

# folders
BIN=/nfs/home/fgehm/gbdc_f/build/gbdc
OUTDIR="/nfs/home/fgehm/gbdc_f/out_${SLURM_JOB_ID}"
ERRDIR="/nfs/home/fgehm/gbdc_f/err_${SLURM_JOB_ID}"

# vars for execution
PARALLEL_JOBS=${SLURM_NTASKS_PER_NODE:-20} 
TIME_LIMIT=${TIME_LIMIT:-600}
LIST=${LIST:-/nfs/home/fgehm/gbdc_f/instances.lst}
RESULT="/nfs/home/fgehm/gbdc_f/results_${SLURM_JOB_ID}.tsv"

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
    # fehlerhaft hier müsste errdir und .err sein
      OUT="$OUTDIR/${ID}_k${k}_i${it}.out"
      # Check if the output file exists and is not empty
      if [[ -s "$OUT" ]]; then
        # Read hash
        HASH=$(grep "c Hash:" "$OUT" | awk '{print $3}')
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
