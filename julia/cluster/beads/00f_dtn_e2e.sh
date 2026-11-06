#!/bin/bash
#SBATCH --job-name=bd-dtn-e2e
#SBATCH --account=mbessa-condo
#SBATCH --qos=mbessa-condo
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=32G
#SBATCH --time=01:30:00
#SBATCH --output=output/beads/00_dtn/%x-%j.out
#SBATCH --error=output/beads/00_dtn/%x-%j.err
#
# End-to-end test of the single-layer DtN (cubic, 6th-order FD) at the PRODUCTION domain 25 x 200:
# PTFE, radius R, Gamma 0.5, 40/60/90 Hz, BDF2, 60 steps/period, run through the real solver with the
# alternative matrix (KMDISK_DTN_FILE). Compare gain_phase.csv with the production campaign
# output/beads/campaign_20261104_174348/04_sweep_R/ptfe/gain_phase.csv (same code commit, same settings).
# Resources: ~24 s to build the matrix (measured). Cases: stage 06 measured ~150 s and ~8 GB per case at
# 60 steps/period (ESTIMATE for these 3 cases in parallel: ~5 min, ~25 GB). 3 threads, 32G, 1.5 h cap.
# From the repo root:  sbatch julia/cluster/beads/00f_dtn_e2e.sh
source "${SLURM_SUBMIT_DIR:?run with sbatch}/julia/cluster/beads/_common.sh"
make_run_dir 00_dtn_e2e
julia --startup-file=no julia/scripts/build_dtn_singlelayer.jl 25 200 6 "$RUN_DIR/D25_Quant200_singlelayer_ord6.jld2"
KMDISK_DTN_FILE="$RUN_DIR/D25_Quant200_singlelayer_ord6.jld2" julia -t "$THREADS" "$BEADS_DIR/run_bead_sweep.jl" \
    stage=00_dtn_e2e outdir="$RUN_DIR/sweep" beads=ptfe basis=R scheme=bdf2 gammas=0.5 freqs=40,60,90 domain=25,200 tres=60
