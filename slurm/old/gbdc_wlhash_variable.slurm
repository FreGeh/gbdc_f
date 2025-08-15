#!/bin/bash -l
#SBATCH --job-name=gbdc_wlhash
#SBATCH --nodelist=hoare
#SBATCH --exclusive
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=20
#SBATCH --mem=0
#SBATCH --time=72:00:00
#SBATCH --chdir=/nfs/home/fgehm/gbdc_f/build
#SBATCH --output=/nfs/home/fgehm/gbdc_f/logs/slurm-%j.out
#SBATCH --error=/nfs/home/fgehm/gbdc_f/logs/slurm-%j.err

set -euo pipefail

# Konfigurierbare Parameter
K_LIST=(1 2 3 4)
MAX_ITERS_LIST=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 25 50 100)
PARALLEL_JOBS=${PARALLEL_JOBS:-20}
TIME_LIMIT=${TIME_LIMIT:-600}
LIST=${LIST:-/nfs/home/fgehm/gbdc_f/instances.lst}

BIN=/nfs/home/fgehm/gbdc_f/build/gbdc
OUTDIR=/nfs/home/fgehm/gbdc_f/out_variable
ERRDIR=/nfs/home/fgehm/gbdc_f/err_variable
RESULT=/nfs/home/fgehm/gbdc_f/wlhash_results_variable.tsv
HASH_REGEX='[0-9a-fA-F]{16}'

# Vorbereitungen
export PATH=$HOME/local/bin:$PATH
which runsolver  >/dev/null || { echo "runsolver not in PATH";  exit 1; }
which parallel   >/dev/null || { echo "parallel  not in PATH";  exit 1; }

mkdir -p "$OUTDIR" "$ERRDIR"

############################
# 1. Jobs generieren
############################
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

############################
# 2. Ergebnisse zusammenführen
############################
{
  printf "id\tk\titers\thash\n"
  while read -r ID _; do
    for k in "${K_LIST[@]}"; do
      for it in "${MAX_ITERS_LIST[@]}"; do
        OUT="$OUTDIR/${ID}_k${k}_i${it}.out"
        HASH=$(grep -m1 -Eo "$HASH_REGEX" "$OUT" || true)
        printf "%s\t%d\t%d\t%s\n" "$ID" "$k" "$it" "$HASH"
      done
    done
  done < "$LIST"
} > "$RESULT"

echo "Fertig – Resultate in $RESULT"
