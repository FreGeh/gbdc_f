# Server
## Setup
1. Spack env profil erstellen
2. Gewollte Spack module installieren
3. gbdc clonen

## Build & Use Executeable
1. `mkdir build && cd build`
2. `cmake ..`
3. `make`
4. `./gbdc wlhash ../test/resources/test_files/0a4ed112f2cdc0a524976a15d1821097-cliquecoloring_n12_k9_c8.cnf.xz --print-stats` schneller test ob funktioniert (`1>/dev/null` für stdout ausblenden und `2>/dev/null` für stderr ausblenden)

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

5. Bei Problem stoppen
- `scancel <jobid>` bestimmten
- `scancel -u $USER` alle

## Slurm usage

* `cpus-per-task` = threads per SLURM task.
* `ntasks-per-node` = how many SLURM tasks on a node.
* `PARALLEL_JOBS` = how many shell jobs GNU Parallel runs at once.