#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 JOBID [RUN_BASE_DIR]"
  echo "Example: $0 46654 /nfs/home/fgehm/gbdc_f/run"
  exit 1
fi

JOBID="$1"
RUN_BASE="${2:-/nfs/home/fgehm/gbdc_f/run}"
RUN="${RUN_BASE%/}/${JOBID}"
RESULT="$RUN/results_${JOBID}.tsv"
SUMMARY="$RUN/hash_summary_${JOBID}.txt"

# sanity checks
if [[ ! -d "$RUN/err" || ! -d "$RUN/out" ]]; then
  echo "Expected $RUN/err and $RUN/out to exist."
  exit 2
fi

# if no .err files, avoid iterating over a literal pattern
shopt -s nullglob

echo -e "ID\tisohash" > "$RESULT"

# parse all .err files (hash is printed to stderr). We also pass the .out as fallback.
count=0
for f in "$RUN/err"/*.err; do
  id="$(basename "$f" .err)"
  out="$RUN/out/${id}.out"
  hash="$(awk '/^c[[:space:]]+Hash:/{ if ($3 ~ /^[0-9a-fA-F]{32}$/) h=$3 } END{ print h }' \
            "$f" "$out" 2>/dev/null || true)"
  [[ -z "$hash" ]] && hash=NA
  printf "%s\t%s\n" "$id" "$hash" >> "$RESULT"
  ((count++)) || true
done

if [[ $count -eq 0 ]]; then
  echo "No .err files found in $RUN/err. Nothing to do."
  exit 0
fi

# summary
TOTAL="$(tail -n +2 "$RESULT" | wc -l | awk '{print $1}')"
UNIQUE="$(tail -n +2 "$RESULT" | awk -F'\t' '{print $2}' | grep -v '^NA$' | LC_ALL=C sort -u | wc -l | awk '{print $1}')"
NA="$(tail -n +2 "$RESULT" | awk -F'\t' '{print $2}' | grep -c '^NA$' | awk '{print $1}')"

{
  echo "TOTAL: $TOTAL"
  echo "UNIQUE: $UNIQUE"
  echo "NA: $NA"
} > "$SUMMARY"

echo "Wrote $RESULT and $SUMMARY"
