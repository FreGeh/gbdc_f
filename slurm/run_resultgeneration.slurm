#!/bin/bash -l
#SBATCH --job-name=gbdc_wlhash_summary
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
OUTDIR="/nfs/home/fgehm/gbdc_f/out_46650"
ERRDIR="/nfs/home/fgehm/gbdc_f/err_46650"

# vars for execution
TOTAL_CPUS=${SLURM_CPUS_PER_TASK:-${SLURM_CPUS_ON_NODE:-$(nproc)}}
PARALLEL_JOBS=${PARALLEL_JOBS:-$TOTAL_CPUS}
TIME_LIMIT=${TIME_LIMIT:-600}
LIST=${LIST:-/nfs/home/fgehm/gbdc_f/instances.lst}
RESULT="/nfs/home/fgehm/gbdc_f/results_46650.tsv"
UNIQUE_SUMMARY_FILE="/nfs/home/fgehm/gbdc_f/hash_summary_46650.txt"

# run vars
MAX_ITERS_LIST=(1 2 3 4 5 6 7 8 9 10 15 25 50 100)
K_LIST=(1)


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

## analyzing hashes => do with pandas