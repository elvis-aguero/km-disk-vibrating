#!/bin/bash
#SBATCH --job-name=bd-converge
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=80G
#SBATCH --time=03:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=elvis_vera@brown.edu
#SBATCH --output=output/beads/07_convergence/%x-%j.out
#SBATCH --error=output/beads/07_convergence/%x-%j.err
#
# Time-step convergence study. PTFE, radius R, 2 frequencies x 2 accelerations = 4 cases per level.
# Levels (steps per bath period): BDF2 at 30, 60, 120, 240; Euler at 60, 120, 240.
# Three or more levels give the observed order of accuracy and a Richardson error estimate.
# Analysis: python3 -I julia/cluster/beads/analyze_convergence.py <RUN_DIR>
#
# Resources (ESTIMATE, not measured): stage 06 needed ~150 s per case at 60 steps/period and
# 47 GB peak for 6 concurrent cases. Cost scales with steps, so the 240-step cases may need
# ~600 s each. Total ~2200 s of work per (f, Gamma) -> ~25 min on 6 threads if balanced.
# Memory is not known for 240 steps: 80 GB is a margin over 47 GB; check MaxRSS afterwards.
# Cases run on 6 threads inside one Julia process, one level at a time. Finished cases are
# skipped on resume (set RUN_DIR).

source "${SLURM_SUBMIT_DIR:?run with sbatch}/julia/cluster/beads/_common.sh"
make_run_dir 07_convergence

GAMMAS="${GAMMAS:-0.5,1.0}"
FREQS_CONV="${FREQS:-40,90}"     # both in the 5 Hz grid; both 'good' in the main campaign
DOM="${DOMAIN:-25,200}"

run() {   # scheme tres
    julia -t "$THREADS" "$BEADS_DIR/run_bead_sweep.jl" stage="07_$1_tres$2" outdir="$RUN_DIR/$1_tres$2" \
        beads=ptfe basis=R scheme="$1" gammas="$GAMMAS" freqs="$FREQS_CONV" domain="$DOM" tres="$2" \
        || echo "level $1/$2 exited with code $? (continuing)"
}

for t in 30 60 120 240; do run bdf2 "$t"; done
for t in 60 120 240;    do run euler "$t"; done
