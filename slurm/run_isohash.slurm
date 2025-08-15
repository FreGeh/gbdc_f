#!/bin/bash -l
#SBATCH --job-name=gbdc_isohash
#SBATCH --nodelist=liskov
#SBATCH --exclusive
#SBATCH --cpus-per-task=256
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
PARALLEL_JOBS=${PARALLEL_JOBS:-$TOTAL_CPUS}
TIME_LIMIT=${TIME_LIMIT:-60}
LIST=${LIST:-/nfs/home/fgehm/gbdc_f/instances.lst}
RESULT="/nfs/home/fgehm/gbdc_f/run/${SLURM_JOB_ID}/results_${SLURM_JOB_ID}.tsv"
UNIQUE_SUMMARY_FILE="/nfs/home/fgehm/gbdc_f/run/${SLURM_JOB_ID}/hash_summary_${SLURM_JOB_ID}.txt"


export PATH=$HOME/local/bin:$PATH
which runsolver  || { echo "runsolver not in PATH";  exit 1; }
which parallel   || { echo "parallel  not in PATH";  exit 1; }

mkdir -p "$OUTDIR" "$ERRDIR"

# run jobs

echo "Starting job generation and execution for ${SLURM_JOB_ID}"

## for each instance, for each k, for each iteration run the command
while read -r ID INST; do
      OUT="$OUTDIR/${ID}.out"
      ERR="$ERRDIR/${ID}.err"
      [[ -s "$OUT" ]] && continue
      printf 'runsolver -C %d "%s" isohash "%s" > "%s" 2> "%s"\n' \
             "$TIME_LIMIT" "$BIN" "$INST" "$OUT" "$ERR"
done < "$LIST" | parallel -j "$PARALLEL_JOBS"

echo "All jobs completed"

# -- parse results --

echo "Starting to parse results"

## nice header
echo -e "ID\tisohash" > "$RESULT"

## reading all the hashes and parsing into result file
while read -r ID INST; do
      ERR="$ERRDIR/${ID}.err"
      # Check if the output file exists and is not empty
      if [[ -s "$ERR" ]]; then
        # Read hash
        HASH=$(grep "c Hash:" "$ERR" | awk '{print $3}')
        # Write hash into correct result file, if not found write NA
        if [[ -n "$HASH" ]]; then
            echo -e "${ID}\t${HASH}" >> "$RESULT"
        else
            echo -e "${ID}\tNA" >> "$RESULT"
        fi
      else
        echo -e "${ID}\tNA" >> "$RESULT"
      fi
done < "$LIST"

echo "Result parsing finished. Final results are in $RESULT"

## get amount of unique hashes

echo "Computing unique hash count (excluding NA)"

# extract hash column (skip header) and print amount of entries, amount of unique hashes, amount of NA into SUMMARY FILE
TOTAL=$(tail -n +2 "$RESULT" | wc -l | awk '{print $1}')
UNIQUE=$(tail -n +2 "$RESULT" | awk -F'\t' '{print $2}' | grep -v '^NA$' | LC_ALL=C sort -u | wc -l | awk '{print $1}')
NA=$(tail -n +2 "$RESULT" | awk -F'\t' '{print $2}' | grep -c '^NA$' | awk '{print $1}')
{
  echo "TOTAL: $TOTAL"
  echo "UNIQUE: $UNIQUE"
  echo "NA: $NA"
} > "$UNIQUE_SUMMARY_FILE"

echo "Unique hash count written to: $UNIQUE_SUMMARY_FILE"