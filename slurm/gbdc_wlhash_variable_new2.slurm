#!/bin/bash -l
#SBATCH --job-name=gbdc_wlhash
#SBATCH --nodelist=hoare
#SBATCH --exclusive
#SBATCH --cpus-per-task=160
#SBATCH --ntasks=1
#SBATCH --mem=0
#SBATCH --time=30:00:00
#SBATCH --chdir=/nfs/home/fgehm/gbdc_f/build2
#SBATCH --output=/nfs/home/fgehm/gbdc_f/logs/slurm-%j.out
#SBATCH --error=/nfs/home/fgehm/gbdc_f/logs/slurm-%j.err

set -euo pipefail

############################
# 0) Parameter
############################
K_LIST=(1)                                       # später beliebig erweitern
MAX_ITERS_LIST=(1 2 4 6 8 10)            # aufsteigend!
PARALLEL_JOBS=${PARALLEL_JOBS:-140}
TIME_LIMIT=${TIME_LIMIT:-600}                   # CPU-Sekunden
WALL_LIMIT=${WALL_LIMIT:-900}                  # Wand-Sekunden
LIST=${LIST:-/nfs/home/fgehm/gbdc_f/instances.lst}

BIN=/nfs/home/fgehm/gbdc_f/build2/gbdc
OUTDIR=/nfs/home/fgehm/gbdc_f/out_variable_new5
ERRDIR=/nfs/home/fgehm/gbdc_f/err_variable_new5
LOGDIR=/nfs/home/fgehm/gbdc_f/logs
RESULT=/nfs/home/fgehm/gbdc_f/wlhash_results_variable5.tsv

HASH_RE='[0-9a-fA-F]{16}'                        # 64-Bit-Hash
CMDFILE="$LOGDIR/cmds-${SLURM_JOB_ID}.txt"
JOBLOG="$LOGDIR/parallel_joblog-${SLURM_JOB_ID}.tsv"

mkdir -p "$OUTDIR" "$ERRDIR" "$LOGDIR"

export PATH="$HOME/local/bin:$PATH"
which runsolver  >/dev/null || { echo "runsolver not in PATH"; exit 1; }
which parallel   >/dev/null || { echo "parallel  not in PATH";  exit 1; }

############################
# 1) Command-Liste erzeugen
############################
: > "$CMDFILE"

while read -r ID INST; do
  [[ -z $ID ]] && continue
  for k in "${K_LIST[@]}"; do
    # Prüfen, ob schon eine Break-Iteration vorliegt
    for it in "${MAX_ITERS_LIST[@]}"; do
      ERR="$ERRDIR/${ID}_k${k}_i${it}.err"
      if [[ -s $ERR ]]; then
        # Hash + Break schon vorhanden?
        break_it=$(grep -m1 '^c Iteration Loop Break:' "$ERR" 2>/dev/null || true)
        [[ -n $break_it && $it -gt $break_it ]] && continue   # Datei wurde schon kopiert
      fi

      # Kommando nur anlegen, wenn wir die Iteration wirklich noch brauchen
      printf 'runsolver -C %d -W %d "%s" wlhash --k %d --max-iters %d "%s" > "%s" 2> "%s"\n' \
             "$TIME_LIMIT" "$WALL_LIMIT" "$BIN" "$k" "$it" "$INST" \
             "$OUTDIR/${ID}_k${k}_i${it}.out" "$ERR" \
             >> "$CMDFILE"
      # Sobald das erste noch offene it gesammelt ist, brechen wir ab;
      # höhere it-Werte werden ggf. später per Copy-Trick erzeugt
      break
    done
  done
done < "$LIST"

############################
# 2) Solver laufen lassen
############################
parallel --joblog "$JOBLOG" -j "$PARALLEL_JOBS" < "$CMDFILE"

############################
# 3) Break-Iteration auswerten
############################
echo "Post-processing Copy-Logs …"
while read -r ID INST; do
  [[ -z $ID ]] && continue
  for k in "${K_LIST[@]}"; do
    # Kleinste vorhandene Iteration mit gültigem Hash suchen
    for it in "${MAX_ITERS_LIST[@]}"; do
      ERR="$ERRDIR/${ID}_k${k}_i${it}.err"
      [[ -s $ERR ]] || continue
      HASH=$(grep -m1 -E "^c Hash: $HASH_RE" "$ERR" | awk '{print $3}' || grep -m1 -E "$HASH_RE" "$OUT" || true)
      BREAK=$(grep -m1 '^c Iteration Loop Break:' "$ERR" | awk '{print $5}')
      [[ -z $HASH || -z $BREAK ]] && continue   # unvollständig → nächste it testen

      # 3.1 – fehlende höheren it-Files anlegen
      for it2 in "${MAX_ITERS_LIST[@]}"; do
        [[ $it2 -le $BREAK ]] && continue       # noch „echte“ Runs
        ERR2="$ERRDIR/${ID}_k${k}_i${it2}.err"
        OUT2="$OUTDIR/${ID}_k${k}_i${it2}.out"
        [[ -s $ERR2 ]] && continue              # existiert bereits

        printf 'Detected CNF\nc Iteration Loop Break: %d\nc Hash: %s\n' \
               "$BREAK" "$HASH" >  "$ERR2"
        : > "$OUT2"                             # leere OUT-Datei
      done
      break   # ID,k abgeschlossen
    done
  done
done < "$LIST"

############################
# 4) TSV-Ergebnis bauen
############################
{
  printf "id\tk\titers\thash\n"
  while read -r ID _; do
    for k in "${K_LIST[@]}"; do
      for it in "${MAX_ITERS_LIST[@]}"; do
        ERR="$ERRDIR/${ID}_k${k}_i${it}.err"
        HASH=$(grep -m1 -E "^c Hash: $HASH_RE" "$ERR" | awk '{print $3}' || grep -m1 -E "$HASH_RE" "$OUT" || true)
        printf "%s\t%d\t%d\t%s\n" "$ID" "$k" "$it" "$HASH"
      done
    done
  done < "$LIST"
} > "$RESULT"

echo "Fertig – Resultate in $RESULT"
echo "Joblog:   $JOBLOG"
