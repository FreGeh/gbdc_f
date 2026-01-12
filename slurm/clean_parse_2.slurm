#!/bin/bash -l
#SBATCH --job-name=gbdc_wlhash_parse
#SBATCH --partition=hoare
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=80
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=04:00:00
#SBATCH --chdir=/nfs/home/fgehm/gbdc_f/build
#SBATCH --output=/nfs/home/fgehm/gbdc_f/logs/parse-%j.out
#SBATCH --error=/nfs/home/fgehm/gbdc_f/logs/parse-%j.err

set -euo pipefail

# Usage: sbatch parse.sh [RUN_ID]
RUN_ID="50758"

RUN_DIR="/nfs/home/fgehm/gbdc_f/run/${RUN_ID}"
OUTDIR="${RUN_DIR}/out"
ERRDIR="${RUN_DIR}/err"
RESULTS="${RUN_DIR}/results_isohash.tsv"

echo "Parsing run directory: ${RUN_DIR}"
[[ -d "$OUTDIR" && -d "$ERRDIR" ]] || { echo "Missing out/ or err/ in ${RUN_DIR}"; exit 1; }

export RUN_DIR OUTDIR ERRDIR RESULTS

command -v python3 >/dev/null || { echo "python3 not found"; exit 1; }

python3 - << 'PY'
import os, re, sys

outdir = os.environ["OUTDIR"]
errdir = os.environ["ERRDIR"]
results_path = os.environ["RESULTS"]

# Only accept isohash files, with optional _mi<iters>
fname_re = re.compile(
    r"^(?P<id>[0-9a-fA-F]{8,})_isohash(?:_mi(?P<iters>\d+))?\.(?P<ext>err|out)$"
)

hash_re        = re.compile(r"(?i)\b(Hash|Digest):\s*([0-9a-f]{8,64})")
plain_hash_re  = re.compile(r"^\s*([0-9a-fA-F]{8,}|[0-9]{8,})\s*$")
stab_re        = re.compile(r"(?i)stabil\w*\D+(\d+)")
cpu_time_1     = re.compile(r"(?im)^CPU time \(s\):\s*([0-9.]+)\s*$")
cpu_time_2     = re.compile(r"(?im)^total CPU time \(s\):\s*([0-9.]+)\s*$")
max_mem_kib    = re.compile(r"(?im)^Max\. memory .*?\(KiB\):\s*([0-9]+)\s*$")
max_rss_fall   = re.compile(r"(?i)maximum resident set size=\s*([0-9]+)\b")

def build_name(id_, iters, ext):
    if iters == "NA":
        return f"{id_}_isohash.{ext}"
    return f"{id_}_isohash_mi{iters}.{ext}"

def parse_err(path):
    wlhash = "NA"
    stab   = "NA"
    try:
        with open(path, "r", errors="ignore") as f:
            for line in f:
                if wlhash == "NA":
                    m = hash_re.search(line)
                    if m:
                        wlhash = m.group(2).lower()
                    else:
                        pm = plain_hash_re.match(line)
                        if pm:
                            wlhash = pm.group(1).lower()
                if stab == "NA":
                    s = stab_re.search(line)
                    if s:
                        stab = s.group(1)
                if wlhash != "NA" and stab != "NA":
                    break
    except FileNotFoundError:
        pass
    except Exception as e:
        sys.stderr.write(f"[warn] err parse failed {path}: {e}\n")
    return wlhash, stab

def parse_out(path):
    cpu = "NA"
    mem = "NA"
    try:
        with open(path, "r", errors="ignore") as f:
            for line in f:
                if cpu == "NA":
                    m = cpu_time_1.search(line) or cpu_time_2.search(line)
                    if m:
                        cpu = m.group(1)
                if mem == "NA":
                    mm = max_mem_kib.search(line)
                    if mm:
                        mem = mm.group(1)
                if cpu != "NA" and mem != "NA":
                    break
            if mem == "NA":
                f.seek(0)
                for line in f:
                    mm2 = max_rss_fall.search(line)
                    if mm2:
                        mem = mm2.group(1)
                        break
    except FileNotFoundError:
        pass
    except Exception as e:
        sys.stderr.write(f"[warn] out parse failed {path}: {e}\n")
    return cpu, mem

def write_header(fp):
    fp.write("id\tmode\tmax_iters\twlhash\tstabilized_at\tcpu_time_s\tmax_mem_kib\n")

# Collect keys from both dirs, but only isohash matches
keys = set()

for d in (errdir, outdir):
    try:
        with os.scandir(d) as it:
            for entry in it:
                if not entry.is_file():
                    continue
                m = fname_re.match(entry.name)
                if not m:
                    continue
                id_ = m.group("id")
                iters = m.group("iters") or "NA"
                keys.add((id_, iters))
    except FileNotFoundError:
        pass

rows = 0
with open(results_path, "w") as outfp:
    write_header(outfp)
    for id_, iters in sorted(keys):
        err_path = os.path.join(errdir, build_name(id_, iters, "err"))
        out_path = os.path.join(outdir, build_name(id_, iters, "out"))

        wlhash, stab = parse_err(err_path)
        cpu, mem     = parse_out(out_path)

        outfp.write(f"{id_}\tisohash\t{iters}\t{wlhash}\t{stab}\t{cpu}\t{mem}\n")
        rows += 1
        if rows % 1000 == 0:
            outfp.flush()

print(f"Wrote {rows} isohash rows to: {results_path}")
PY


echo "Done. Results at: ${RESULTS}"
