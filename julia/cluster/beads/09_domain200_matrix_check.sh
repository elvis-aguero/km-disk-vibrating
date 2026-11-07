#!/bin/bash
#SBATCH --job-name=dom200-matrix
#SBATCH --account=mbessa-condo
#SBATCH --qos=mbessa-condo
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=01:00:00
#SBATCH --output=output/beads/09_schur_check/%x-%j.out
#SBATCH --error=output/beads/09_schur_check/%x-%j.err
#
# Echo-sensitive rows in the small domain (25 x 200, wall echo prediction 0.07-0.15): al32 rc at 52, 60, 70 Hz, uniform nr = 2500,
# once with the paper's polar DtN matrix (should reproduce the old CSV rows) and once with the single-layer matrix (the current default).
# Separates "matrix" from "graded grid" for the C7b difference at R = 100 a. Solver: :auto (Schur). Six cases, a few minutes each
# (ESTIMATE; nr = 2500 Schur case about 1-2 min). Memory: ~3 GB per case.
source "${SLURM_SUBMIT_DIR:?run with sbatch}/julia/cluster/beads/_common.sh"
make_run_dir 09_schur_check
for METHOD in polar single_layer; do
  for FREQ in 52 60 70; do
    echo "=== al32 rc ${FREQ} Hz, DtN method ${METHOD}"
    KMDISK_DTN_METHOD="$METHOD" julia -t "$THREADS" "$BEADS_DIR/run_bead_sweep.jl" stage="09_dom200_${METHOD}_${FREQ}" \
        outdir="$RUN_DIR/${METHOD}_${FREQ}" beads="al32" basis="rc" scheme=bdf2 gammas=0.5 freqs="$FREQ" domain=25,200 tres=60 \
        || echo "case ${METHOD} ${FREQ} exited with code $? (continuing)"
  done
done
