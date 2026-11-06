#!/bin/bash
#SBATCH --job-name=bd-dtn-ctrl-al32rc
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
# CONTROL for the single-layer DtN: Al 1/32" at rc (the case that failed), production domain 25 x 200, Gamma 0.5,
# 40/60/70 Hz, BDF2/60, run with the single-layer matrix (cubic, 6th-order FD). Compare with the production campaign
# output/beads/campaign_20261104_174348/05_sweep_rc/al32/gain_phase.csv. Same resources as 00f (measured ~150 s and ~8 GB per
# case at nr = 2500; ESTIMATE ~5 min for 3 cases in parallel). From the repo root:  sbatch julia/cluster/beads/00g_dtn_control.sh
source "${SLURM_SUBMIT_DIR:?run with sbatch}/julia/cluster/beads/_common.sh"
make_run_dir 00_dtn_ctrl
julia --startup-file=no julia/scripts/build_dtn_singlelayer.jl 25 200 6 "$RUN_DIR/D25_Quant200_singlelayer_ord6.jld2"
KMDISK_DTN_FILE="$RUN_DIR/D25_Quant200_singlelayer_ord6.jld2" julia -t "$THREADS" "$BEADS_DIR/run_bead_sweep.jl" \
    stage=00_dtn_ctrl outdir="$RUN_DIR/sweep" beads=al32 basis=rc scheme=bdf2 gammas=0.5 freqs=40,60,70 domain=25,200 tres=60
