#!/bin/bash
#SBATCH --job-name=bd-smoke
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=output/beads/01_smoke/%x-%j.out
#SBATCH --error=output/beads/01_smoke/%x-%j.err
#
# Stage 01: plumbing test of the driver and the post-processor on the SMALL cache
# (D5Quant20, nr=50). It does not need the D25Q200 cache, so it can run while stage 00 runs.
# It checks that: the driver starts, both time schemes run, per-bead failure handling and
# the gain/phase CSVs work. It says nothing about the physics: the 20-radius bath is tiny
# (waves reach the wall), so `wave_cut` rows are expected here.
#
# Resources: nr=50 matrices are tiny (a few MB); 2 cores, 4G and 20 min are a generous
# cap for 6 short cases plus Julia start-up.
#
#     cd /oscar/data/dharri15/eaguerov/Github/km-disk-vibrating
#     sbatch julia/cluster/beads/01_smoke.sh

source "${SLURM_SUBMIT_DIR:?run with sbatch}/julia/cluster/beads/_common.sh"
make_run_dir 01_smoke

# BDF2, two beads, two frequencies (4 cases).
julia -t "$THREADS" "$BEADS_DIR/run_bead_sweep.jl" stage=01_smoke outdir="$RUN_DIR/bdf2" \
    beads=ptfe,al32 basis=R scheme=bdf2 gammas=0.1 freqs=50,60 domain=5,20 tres=60 require_finite=true

# Euler, one bead, one frequency (1 case): the other time-stepping branch.
julia -t "$THREADS" "$BEADS_DIR/run_bead_sweep.jl" stage=01_smoke outdir="$RUN_DIR/euler" \
    beads=ptfe basis=rc scheme=euler gammas=0.1 freqs=60 domain=5,20 tres=60 require_finite=true

echo
echo "Smoke test finished. Look at:"
echo "  $RUN_DIR/bdf2/gain_phase_all.csv"
echo "  $RUN_DIR/euler/gain_phase_all.csv"
