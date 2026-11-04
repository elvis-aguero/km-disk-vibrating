#!/bin/bash
#SBATCH --job-name=bd-sweep-R
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=80G
#SBATCH --time=04:00:00
#SBATCH --output=output/beads/04_sweep_R/%x-%j.out
#SBATCH --error=output/beads/04_sweep_R/%x-%j.err
#
# Stage 04: four beads at their ACTUAL radius R (PDF table 2), 5-90 Hz, BDF2.
# Disk mass = sphere mass. Output: gain and phase lag, lab and shaker frame.
#
# GAMMAS is REQUIRED. Choose it from the stage-03 acceleration test; there is no default.
#     GAMMAS=0.5 sbatch julia/cluster/beads/04_sweep_R.sh
#
# Resources (ESTIMATE; refine from the pilot): 4 beads x 34 frequencies x (number of Gamma)
# cases. For one Gamma that is 136 cases x ~150 s / 6 concurrent = ~1 h. The 4 h limit
# covers up to about 3 Gamma values; for more, raise it:  sbatch --time=08:00:00 ...
# Memory: measured peak 57-58 GB for 6 concurrent cases (~9.5 GB per case; hit 99% of 56G in the first campaign), so 80G.
#
# One bead failing, or any case failing, does not stop the others: see gain_phase.csv
# (column `quality`) and the end-of-job "Quality counts" line.
# Resume: RUN_DIR=<same run dir> GAMMAS=... sbatch ...
# Optional overrides: BEADS (default all four), SCHEME, TRES, DOMAIN, FREQS.

: "${GAMMAS:?set GAMMAS first, e.g.  GAMMAS=0.5 sbatch julia/cluster/beads/04_sweep_R.sh}"
source "${SLURM_SUBMIT_DIR:?run with sbatch}/julia/cluster/beads/_common.sh"
make_run_dir 04_sweep_R

julia -t "$THREADS" "$BEADS_DIR/run_bead_sweep.jl" stage=04_sweep_R outdir="$RUN_DIR" \
    beads="${BEADS:-ptfe,al32,al16,nylon}" basis=R scheme="${SCHEME:-bdf2}" gammas="$GAMMAS" \
    freqs="${FREQS:-grid}" domain="${DOMAIN:-25,200}" tres="${TRES:-60}"
