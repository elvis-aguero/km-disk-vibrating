#!/bin/bash
#SBATCH --job-name=bd-sweep-rc
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=56G
#SBATCH --time=04:00:00
#SBATCH --output=output/beads/05_sweep_rc/%x-%j.out
#SBATCH --error=output/beads/05_sweep_rc/%x-%j.err
#
# Stage 05: same as stage 04, but diskRadius = the CONTACT-LINE radius rc (PDF table 3).
# Disk mass is still the sphere mass.
#
# GAMMAS is REQUIRED (use the same value as stage 04, so the two sets can be compared).
#     GAMMAS=0.5 sbatch julia/cluster/beads/05_sweep_rc.sh
#
# Resources, failure handling, resume and overrides: identical to 04_sweep_R.sh.

: "${GAMMAS:?set GAMMAS first, e.g.  GAMMAS=0.5 sbatch julia/cluster/beads/05_sweep_rc.sh}"
source "${SLURM_SUBMIT_DIR:?run with sbatch}/julia/cluster/beads/_common.sh"
make_run_dir 05_sweep_rc

julia -t "$THREADS" "$BEADS_DIR/run_bead_sweep.jl" stage=05_sweep_rc outdir="$RUN_DIR" \
    beads="${BEADS:-ptfe,al32,al16,nylon}" basis=rc scheme="${SCHEME:-bdf2}" gammas="$GAMMAS" \
    freqs="${FREQS:-grid}" domain="${DOMAIN:-25,200}" tres="${TRES:-60}"
