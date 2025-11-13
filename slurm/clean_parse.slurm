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

# >>> Set your run id here <<<
RUN_ID="50359"

RUN_DIR="/nfs/home/fgehm/gbdc_f/run/${RUN_ID}"
OUTDIR="${RUN_DIR}/out"
ERRDIR="${RUN_DIR}/err"
RESULTS="${RUN_DIR}/results.tsv"

echo "Parsing run directory: ${RUN_DIR}"
[[ -d "$OUTDIR" && -d "$ERRDIR" ]] || { echo "Missing out/ or err/ in ${RUN_DIR}"; exit 1; }

export RUN_DIR OUTDIR ERRDIR RESULTS

command -v python3 >/dev/null || { echo "python3 not found"; exit 1; }

python3 - << 'PY'
import os, re, sys

run_dir  = os.environ['RUN_DIR']
outdir   = os.environ['OUTDIR']
errdir   = os.environ['ERRDIR']
results_path = os.environ['RESULTS']

# ---------- configurable mode mapping (edit if your naming changes) ----------
OUT2ERR = {'timon':'timon', 'xor':'xor', 'sort':'sort', 'half':'half'}
ERR2OUT = {v:k for k,v in OUT2ERR.items()}
# -----------------------------------------------------------------------------

# <id>_<mode>_mi<iters>.(err|out)  (allow underscores in mode)
fname_re = re.compile(r'^(?P<id>[0-9a-fA-F]{8,})_(?P<mode>[A-Za-z0-9_\-]+)_mi(?P<iters>\d+)\.(?P<ext>err|out)$')

# Patterns inside files
hash_re        = re.compile(r'(?i)\b(Hash|Digest):\s*([0-9a-f]{8,64})')
plain_hash_re  = re.compile(r'^\s*([0-9a-fA-F]{8,}|[0-9]{8,})\s*$')  # bare hex or decimal
stab_re        = re.compile(r'(?i)stabil\w*\D+(\d+)')                 # "stabilized" info
cpu_time_1     = re.compile(r'(?im)^CPU time \(s\):\s*([0-9.]+)\s*$')
cpu_time_2     = re.compile(r'(?im)^total CPU time \(s\):\s*([0-9.]+)\s*$')
max_mem_kib    = re.compile(r'(?im)^Max\. memory .*?\(KiB\):\s*([0-9]+)\s*$')
max_rss_fall   = re.compile(r'(?i)maximum resident set size=\s*([0-9]+)\b')

def parse_err(path):
    wlhash = "NA"
    stab   = "NA"
    try:
        with open(path, 'r', errors='ignore') as f:
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
        with open(path, 'r', errors='ignore') as f:
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

processed = set()
rows = 0

with open(results_path, 'w') as outfp:
    write_header(outfp)

# Pass 1: every .err file (pair to its .out using ERR2OUT mapping)
try:
    with os.scandir(errdir) as it:
        for entry in it:
            if not entry.is_file(): continue
            name = entry.name
            m = fname_re.match(name)
            if not m or m.group('ext') != 'err': continue

            id_   = m.group('id')
            mode_err = m.group('mode')
            iters = m.group('iters')
            stem_err = name[:-4]

            out_mode = ERR2OUT.get(mode_err, mode_err)
            out_name = f"{id_}_{out_mode}_mi{iters}.out"

            err_path = os.path.join(errdir, name)
            out_path = os.path.join(outdir, out_name)

            wlhash, stab = parse_err(err_path)
            cpu, mem     = parse_out(out_path)

            with open(results_path, 'a') as outfp:
                outfp.write(f"{id_}\t{out_mode}\t{iters}\t{wlhash}\t{stab}\t{cpu}\t{mem}\n")
                rows += 1
                if rows % 1000 == 0:
                    outfp.flush()

            processed.add(f"{id_}_{out_mode}_mi{iters}")
except FileNotFoundError:
    pass

# Pass 2: .out files without .err counterpart (pair using OUT2ERR mapping)
try:
    with os.scandir(outdir) as it:
        for entry in it:
            if not entry.is_file(): continue
            name = entry.name
            m = fname_re.match(name)
            if not m or m.group('ext') != 'out': continue

            id_   = m.group('id')
            mode_out = m.group('mode')
            iters = m.group('iters')
            key = f"{id_}_{mode_out}_mi{iters}"
            if key in processed:
                continue

            err_mode = OUT2ERR.get(mode_out, mode_out)
            err_name = f"{id_}_{err_mode}_mi{iters}.err"

            out_path = os.path.join(outdir, name)
            err_path = os.path.join(errdir, err_name)

            wlhash, stab = parse_err(err_path)
            cpu, mem     = parse_out(out_path)

            with open(results_path, 'a') as outfp:
                outfp.write(f"{id_}\t{mode_out}\t{iters}\t{wlhash}\t{stab}\t{cpu}\t{mem}\n")
                rows += 1
                if rows % 1000 == 0:
                    outfp.flush()
except FileNotFoundError:
    pass

print(f"Wrote {rows} rows to: {results_path}")
PY

echo "Done. Results at: ${RESULTS}"
