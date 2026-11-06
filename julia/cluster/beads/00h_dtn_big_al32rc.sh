#!/bin/bash
#SBATCH --job-name=bd-dtn-big-al32rc
#SBATCH --account=mbessa-condo
#SBATCH --qos=mbessa-condo
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=360G
#SBATCH --time=04:00:00
#SBATCH --output=output/beads/00_dtn/%x-%j.out
#SBATCH --error=output/beads/00_dtn/%x-%j.err
#
# Al 1/32" at rc in the LARGE domain 25 x 700 (nr = 8750, 6.0 cm radius instead of 1.71 cm), Gamma 0.5, 40/60/70 Hz,
# BDF2/60, with the single-layer DtN (cubic, 6th-order FD; the control job 00g checks this matrix against production at
# 25 x 200). Compare gain/phase with 05_sweep_rc/al32 of the production campaign: if the wall echo caused the
# oscillation, the large-domain values should be smooth and match the periods-2-4 estimate (1.458 lab gain at 60 Hz).
# ESTIMATES (nothing measured at nr = 8750): matrix build ~5-10 min (24 s at nr = 2500, cost ~ nr^2, 1 thread);
# per case ~8 GB x (8750/2500)^2 = ~100 GB if memory ~ nr^2, so 3 cases ~300 GB -> 360G; per case ~30 min if time ~ nr^2.
# 4 h cap. From the repo root:  sbatch julia/cluster/beads/00h_dtn_big_al32rc.sh
source "${SLURM_SUBMIT_DIR:?run with sbatch}/julia/cluster/beads/_common.sh"
make_run_dir 00_dtn_big
julia --startup-file=no julia/scripts/build_dtn_singlelayer.jl 25 700 6 "$RUN_DIR/D25_Quant200_singlelayer_ord6_D700.jld2"
KMDISK_DTN_FILE="$RUN_DIR/D25_Quant200_singlelayer_ord6_D700.jld2" julia -t "$THREADS" "$BEADS_DIR/run_bead_sweep.jl" \
    stage=00_dtn_big outdir="$RUN_DIR/sweep" beads=al32 basis=rc scheme=bdf2 gammas=0.5 freqs=40,60,70 domain=25,700 tres=60
