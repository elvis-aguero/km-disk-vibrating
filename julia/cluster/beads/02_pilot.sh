#!/bin/bash
#SBATCH --job-name=bd-pilot
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=output/beads/02_pilot/%x-%j.out
#SBATCH --error=output/beads/02_pilot/%x-%j.err
#
# Stage 02: 3-case pilot on the REAL domain, to measure time and memory per case before
# sizing the big jobs. PTFE at radius R, BDF2, Gamma = 0.1, f = 20, 60 and 90 Hz, so it also
# exercises the fit path (`good` / `unconverged` / `wave_cut` labels) for real.
# Needs the D25Q200 cache (stage 00), unless DOMAIN is overridden.
#
# Resources: one case needs about 6.2 GB for the LU cache (README, nr=2500, 60 steps per
# period). The sweep runs min(cpus, floor(0.7 * mem / 6.2 GB)) cases at once:
# 3 cores and 32G -> 3 concurrent cases, one per core, so all 3 cases run together.
#
#     cd /oscar/data/dharri15/eaguerov/Github/km-disk-vibrating
#     sbatch --dependency=afterok:8262774 --kill-on-invalid-dep=yes julia/cluster/beads/02_pilot.sh
#   (the dependency makes it start only if the cache job ended with exit code 0)
#
# After it ends:  sacct -j <jobid> --format=JobID,Elapsed,MaxRSS,State
#                 and the `elapsed_s` column of <run dir>/ptfe/summary.csv

source "${SLURM_SUBMIT_DIR:?run with sbatch}/julia/cluster/beads/_common.sh"
make_run_dir 02_pilot

julia -t "$THREADS" "$BEADS_DIR/run_bead_sweep.jl" stage=02_pilot outdir="$RUN_DIR" \
    beads=ptfe basis=R scheme="${SCHEME:-bdf2}" gammas=0.1 freqs=20,60,90 \
    domain="${DOMAIN:-25,200}" tres="${TRES:-60}" require_finite=true
