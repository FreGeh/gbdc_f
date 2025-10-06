#!/bin/bash -l
#SBATCH --job-name=join
#SBATCH --nodelist=iverson
#SBATCH --exclusive
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=0
#SBATCH --time=58:00:00
#SBATCH --chdir=/nfs/home/fgehm/gbdc_f/build
#SBATCH --output=/nfs/home/fgehm/gbdc_f/logs/slurm-%j.out
#SBATCH --error=/nfs/home/fgehm/gbdc_f/logs/slurm-%j.err


set -euo pipefail

# >>> Set your run id here <<<
RUN_ID="49265"

RUN_DIR="/nfs/home/fgehm/gbdc_f/run/${RUN_ID}"
RESULTS="${RUN_DIR}/results.tsv"                    # /nfs/home/fgehm/gbdc_f/run/<RUN_ID>/results.tsv
INSTANCES="/nfs/home/fgehm/gbdc_f/instances_new.lst"    # project root
OUT="${RUN_DIR}/results_joined.tsv"                 # new file

# Resolve absolute paths (avoids awk seeing a bare 'results.tsv')
RESULTS_ABS=$(readlink -f "$RESULTS"    || true)
INSTANCES_ABS=$(readlink -f "$INSTANCES"|| true)
OUT_ABS="${OUT}"  # parent dir exists already as part of the run dir

echo "Resolved paths:"
echo "  RESULTS   = $RESULTS_ABS"
echo "  INSTANCES = $INSTANCES_ABS"
echo "  OUT       = $OUT_ABS"

[[ -f "$RESULTS"   ]] || { echo "Missing $RESULTS";   exit 1; }
[[ -f "$INSTANCES" ]] || { echo "Missing $INSTANCES"; exit 1; }

awk -v out="$OUT_ABS" '
BEGIN {
  FS = OFS = "\t"
  ireq_str = "hash local isohash family result clauses variables bytes cls1 cls2 cls3 cls4 cls5 cls6 cls7 cls8 cls9 cls10p horn invhorn positive negative"
  rreq_str = "hash mode wlhash stabilized_at cpu_time_s max_mem_kib"
  ni = split(ireq_str, ireq, " ")
  nr = split(rreq_str, rreq, " ")
}

# Pass 1: read instances.lst, build table keyed by hash
FNR==1 && NR==FNR {
  for (i=1; i<=NF; i++) IIdx[$i] = i
  for (i=1; i<=ni; i++) if (!(ireq[i] in IIdx)) { print "Missing column in instances.lst:", ireq[i] > "/dev/stderr"; exit 1 }
  next
}
NR==FNR {
  key = $(IIdx["hash"])
  v = $(IIdx["local"])
  v = v OFS $(IIdx["isohash"])
  v = v OFS $(IIdx["family"])
  v = v OFS $(IIdx["result"])
  v = v OFS $(IIdx["clauses"])
  v = v OFS $(IIdx["variables"])
  v = v OFS $(IIdx["bytes"])
  v = v OFS $(IIdx["cls1"])
  v = v OFS $(IIdx["cls2"])
  v = v OFS $(IIdx["cls3"])
  v = v OFS $(IIdx["cls4"])
  v = v OFS $(IIdx["cls5"])
  v = v OFS $(IIdx["cls6"])
  v = v OFS $(IIdx["cls7"])
  v = v OFS $(IIdx["cls8"])
  v = v OFS $(IIdx["cls9"])
  v = v OFS $(IIdx["cls10p"])
  v = v OFS $(IIdx["horn"])
  v = v OFS $(IIdx["invhorn"])
  v = v OFS $(IIdx["positive"])
  v = v OFS $(IIdx["negative"])
  inst_map[key] = v
  next
}

# Pass 2: read results.tsv, write joined rows in results order
FNR==1 {
  for (i=1; i<=NF; i++) RIdx[$i] = i
  if (!("hash" in RIdx) && ("id" in RIdx)) RIdx["hash"] = RIdx["id"]  # tolerate id/hash swap

  for (i=1; i<=nr; i++) if (!(rreq[i] in RIdx)) { print "Missing column in results.tsv:", rreq[i] > "/dev/stderr"; exit 1 }

  print "hash","mode","wlhash","stabilized_at","cpu_time_s","max_mem_kib",
        "local","isohash","family","result","clauses","variables","bytes",
        "cls1","cls2","cls3","cls4","cls5","cls6","cls7","cls8","cls9",
        "cls10p","horn","invhorn","positive","negative" > out
  next
}

{
  key = $(RIdx["hash"])
  if (key in inst_map) {
    print $(RIdx["hash"]),
          $(RIdx["mode"]),
          $(RIdx["wlhash"]),
          $(RIdx["stabilized_at"]),
          $(RIdx["cpu_time_s"]),
          $(RIdx["max_mem_kib"]),
          inst_map[key] >> out
  }
}
' "$INSTANCES_ABS" "$RESULTS_ABS"

echo "Wrote: $OUT_ABS"
