#!/bin/bash
#SBATCH --job-name=bd-accel
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=80G
#SBATCH --time=04:00:00
#SBATCH --output=output/beads/03_accel/%x-%j.out
#SBATCH --error=output/beads/03_accel/%x-%j.err
#
# Stage 03: does the gain / phase-lag curve depend on the acceleration Gamma?
# One bead (PTFE, the reference), diskRadius = R, BDF2, full 34-frequency grid.
#
# Default Gamma list: 0.1, 0.25, 0.5, 1.0   (4 x 34 = 136 cases).
# Gamma > 1 runs as a SECOND job, so a failure there cannot touch these results:
#     GAMMAS=1.5,2.0 sbatch julia/cluster/beads/03_accel.sh        (2 x 34 = 68 cases)
# Rows with Gamma > 1 carry gamma_gt_1 = true: the effective gravity changes sign during
# the cycle and the model may not be valid there.
#
# Resources (an ESTIMATE from the committed reference run, 116-203 s per case at nr=2500;
# refine it from the stage-02 pilot):
#   - Memory: measured peak 55 GB (campaign 20261102_170942) for 6 concurrent cases = ~9 GB per case; it hit 99% of 56G, so 80G.
#   - CPUs: 6, one per concurrent case.
#   - Time: 136 cases x ~150 s / 6 = ~1.0 h. 4 h is a ~4x margin. For the Gamma > 1 job
#     (68 cases) the same limit is more than enough.
#
#     cd /oscar/data/dharri15/eaguerov/Github/km-disk-vibrating
#     sbatch --dependency=afterok:<pilot jobid> julia/cluster/beads/03_accel.sh
# Resume after a time-out or a crash: RUN_DIR=<the same run dir> sbatch ... (finished cases are skipped).
#
# Optional overrides (environment): GAMMAS, BEADS, SCHEME, TRES, DOMAIN, FREQS (default grid).

source "${SLURM_SUBMIT_DIR:?run with sbatch}/julia/cluster/beads/_common.sh"
make_run_dir 03_accel

julia -t "$THREADS" "$BEADS_DIR/run_bead_sweep.jl" stage=03_accel outdir="$RUN_DIR" \
    beads="${BEADS:-ptfe}" basis=R scheme="${SCHEME:-bdf2}" gammas="${GAMMAS:-0.1,0.25,0.5,1.0}" \
    freqs="${FREQS:-grid}" domain="${DOMAIN:-25,200}" tres="${TRES:-60}"
