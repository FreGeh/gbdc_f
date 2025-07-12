# Basis-Dateien
RESULT=$HOME/gbdc_f/wlhash_results.tsv      # 3 Spalten: GBD  CNF  WL
NEWLIST=$HOME/gbdc_f/instances_dupes.lst    # Ziel: GBD  CNF

# 1) Hash-Häufigkeiten ermitteln  ➜  wlhash_freq.txt
cut -f3 "$RESULT" | sort | uniq -c | sort -nr \
  > ~/gbdc_f/wlhash_freq.txt
echo "→ Hash-Häufigkeiten in wlhash_freq.txt (Spalte1 = Count, Spalte2 = WL)"

# 2) nur Hashes mit Count ≥ 2 herausfiltern  ➜  wlhash_dupes.txt
awk '$1 >= 2 {print $2}' ~/gbdc_f/wlhash_freq.txt \
  > ~/gbdc_f/wlhash_dupes.txt
echo "→ Duplikat-Hashes in wlhash_dupes.txt"

# 3) Zu jeder Dupe-WL den passenden GBD-Hash & CNF-Pfad heraussuchen
awk -v dupfile="$HOME/gbdc_f/wlhash_dupes.txt" '
    BEGIN {
        while ((getline h < dupfile) > 0) dup[h]=1
        close(dupfile)
    }
    dup[$3] { print $1, $2 }        # GBDHash <TAB> CNF-Pfad
' "$RESULT" > "$NEWLIST"

echo "→ Neue Instanzliste für Re-Run: $NEWLIST"
wc -l "$NEWLIST"
head "$NEWLIST"
