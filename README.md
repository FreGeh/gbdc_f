# Server
## Setup
1. Spack env profil erstellen
2. Gewollte Spack module installieren
3. gbdc clonen

## Build instance list
1. go into venv
2. set db var `export GBD_DB="/nfs/share/instances/gbd/cnf_local.db:/nfs/home/fgehm/gbd/data/base.db:/nfs/home/fgehm/gbd/data/meta.db"` 
3. get all wanted features in one instance file
```bash
gbd get -g cnf_local:local -r meta:hash meta:isohash meta:family meta:result meta:track base:clauses base:variables base:bytes base:cls1 base:cls2 base:cls3 base:cls4 base:cls5 base:cls6 base:cls7 base:cls8 base:cls9 base:cls10p base:horn base:invhorn base:positive base:negative -c min -H -d $'\t' > instances_new.lst
```


## Scramble 
`./scranfilize (-p/-P/-r/-R) -f 0.5 source_instances/source_cnf.xz scrambled_instances/r_noindividualflips.cnf`
also `./scranfilize -r -f 0 -s 0 source_instances/source_cnf.xz scrambled_instances/r_noindividualflips.cnf` zb

## Build & Use Executeable
1. `mkdir build && cd build`
2. `cmake ..`
3. `make`
4. `./gbdc wlhash ../test/resources/test_files/0a4ed112f2cdc0a524976a15d1821097-cliquecoloring_n12_k9_c8.cnf.xz --print-stats` schneller test ob funktioniert (`1>/dev/null` für stdout ausblenden und `2>/dev/null` für stderr ausblenden) `./gbdc wlhash ../test/scrambled/halfflips/p_half_flips.cnf --print-stats`

## testing
`perf record -F 999 --call-graph dwarf -g -- ./gbdc wlhash ../test/resources/test_files/0a4ed112f2cdc0a524976a15d1821097-cliquecoloring_n12_k9_c8.cnf.xz --stat-encoding`

## Server Useage
Ausgangslage: Datei mit jeweils [GBDHash] [datei] space seperated in instances.lst

0. Immer `spack env activate gbdc` bei Neuanmeldung für builds

1. Server Status Checken
- für hoare - 10/20 parallel jobs
- schauen wer idle wer verfügbar `sinfo -p all -o "%N %T %C"`
- resourcen anzeigen lassen `sinfo -o '| %.10P || %.4c || %.10z || %.8m || %.8f' | sed 'a |-'`
- architektur anzeigen lassen `sinfo -h -N -p all -t idle -o "%N" | xargs -I{} bash -lc 'printf "%-12s " "{}"; scontrol show node {} | sed -n -e "s/.*Arch=\\([^ ]*\\).*/\\1/p" -e "s/.*Features=\\([^ ]*\\).*/ (features: \\1)/p"'`

2. DIRs und so im Slurm Script anpassen

3. Slurm Script mit Sbatch starten
- `sbatch script.slurm`

4. Status prüfen (Server)
- `squeue -u $USER` alle
- `squeue -j <jobid>` bestimmten
- `scontrol show job <jobid>`
- `sacct -j <jobid> --format=JobID,State,ExitCode,Elapsed`
- `sstat -j <jobid>.0 --format=AveCPU,AveRSS,MaxRSS`

4. Status prüfen (Ergebnisse)
- `tail -f /nfs/home/fgehm/gbdc_f/logs/slurm-<jobid>.out` um logs letzte ausgaben zu sehen
- `ls /nfs/home/fgehm/gbdc_f/run/<jobid>/err | wc -l` Anzahl der schon ausgeführten Commands
- `ls run/49265/err/ | wc -l`

5. Bei Problem stoppen
- `scancel <jobid>` bestimmten
- `scancel -u $USER` alle

## Slurm usage

* `cpus-per-task` = threads per SLURM task.
* `ntasks-per-node` = how many SLURM tasks on a node.
* `PARALLEL_JOBS` = how many shell jobs GNU Parallel runs at once.