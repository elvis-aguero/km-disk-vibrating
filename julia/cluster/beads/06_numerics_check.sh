#!/bin/bash
#SBATCH --job-name=bd-numerics
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=56G
#SBATCH --time=03:00:00
#SBATCH --output=output/beads/06_euler_check/%x-%j.out
#SBATCH --error=output/beads/06_euler_check/%x-%j.err
#
# Stage 06: how much do the numerics matter? PTFE at radius R, 7 frequencies (10,30,50,58,60,62,90) including
# 55-65 Hz, two variants that differ from the BDF2 / 60-step reference of stage 03:
#   (a) Euler, 60 steps per period     -> time-scheme error
#   (b) BDF2,  30 steps per period     -> time-step error
# The reference (BDF2, 60 steps) is the stage-03 run at the same Gamma; compare against it.
# These two numbers are the "noise floor" for the acceleration test: the simulation is
# deterministic, so a spread across accelerations only matters if it exceeds them.
#
# GAMMAS is REQUIRED and must be one of the Gamma values used in stage 03.
#     GAMMAS=0.5 sbatch julia/cluster/beads/06_numerics_check.sh
#
# Resources (ESTIMATE): 2 variants x 7 freqs x (number of Gamma) = 14 cases per Gamma,
# ~150 s each / 6 concurrent = ~10 min. 3 h is a very wide margin. The 30-step variant
# needs about half the LU memory per case, so the same 56G is safe.

: "${GAMMAS:?set GAMMAS first, e.g.  GAMMAS=0.5 sbatch julia/cluster/beads/06_numerics_check.sh}"
source "${SLURM_SUBMIT_DIR:?run with sbatch}/julia/cluster/beads/_common.sh"
make_run_dir 06_euler_check

FREQS_CHECK="${FREQS:-10,30,50,58,60,62,90}"   # 7 points: low, mid, three around the peak, high. All are in the 34-point grid.
DOM="${DOMAIN:-25,200}"

julia -t "$THREADS" "$BEADS_DIR/run_bead_sweep.jl" stage=06_euler_tres60 outdir="$RUN_DIR/euler_tres60" \
    beads=ptfe basis=R scheme=euler gammas="$GAMMAS" freqs="$FREQS_CHECK" domain="$DOM" tres=60

julia -t "$THREADS" "$BEADS_DIR/run_bead_sweep.jl" stage=06_bdf2_tres30 outdir="$RUN_DIR/bdf2_tres30" \
    beads=ptfe basis=R scheme=bdf2 gammas="$GAMMAS" freqs="$FREQS_CHECK" domain="$DOM" tres=30
