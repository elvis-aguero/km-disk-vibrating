#!/bin/bash
#SBATCH --job-name=km-validate
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=02:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=elvis_vera@brown.edu
#SBATCH --output=output/validation/logs/%x-%j.out
#SBATCH --error=output/validation/logs/%x-%j.err
#
# One 30-case validation sweep (Gamma 0.05, 0.2, 0.5 x 10..100 Hz, disk R=0.2 cm, D50Quant100,
# nr 2500) against the digitized experiment (matlab/0_data/external/experimental_measurements).
# Used for the 2x2 comparison: force balance (pre-fix / post-fix of PR 7) x time scheme (euler / bdf2).
#
#   CODE_DIR=<checkout to run> SCHEME=euler|bdf2 LABEL=<name> sbatch julia/cluster/validate_vs_experiment.sh
#
# DtN matrix: the registry now resolves the single-layer cache by default. Add KMDISK_DTN_METHOD=polar to the environment
# to reproduce the earlier runs with the polar (paper) matrix; DTN_BUILD=single_layer builds the matrix inside the job.
#
# Output: <submit dir>/output/validation/<LABEL>/sweep (summary.csv, validation_overlay.csv, plots).
# Resources (ESTIMATE from the bead campaigns, same nr=2500 LU size): ~9.5 GB and ~160 s per case;
# 8 cases run at a time -> ~76 GB and ~10-15 min. 96G and 2 h leave a margin.
set -euo pipefail
MAIN="${SLURM_SUBMIT_DIR:?run with sbatch}"
: "${CODE_DIR:?}" "${SCHEME:?}" "${LABEL:?}"
module load julia/1.12.4s-adme
OUT="$MAIN/output/validation/$LABEL"
mkdir -p "$OUT/sweep" "$OUT/logs"
export KMDISK_LOG_DIR="$OUT/logs" KMDISK_TIME_SCHEME="$SCHEME"
cd "$CODE_DIR"
# A separate checkout has no (gitignored) DTN cache: link the one of the main checkout.
if [[ "$CODE_DIR" != "$MAIN" && ! -e julia/data/dtn_cache ]]; then
    mkdir -p julia/data && ln -s "$MAIN/julia/data/dtn_cache" julia/data/dtn_cache
fi
echo "code: $(git rev-parse --short HEAD) in $CODE_DIR; scheme $SCHEME; label $LABEL; threads $SLURM_CPUS_PER_TASK"
# DTN_BUILD=single_layer: build the 50 x 100 matrix with the single-layer construction (docs/dtn_single_layer.md,
# ~24 s at nr = 2500, measured) and use it instead of the registered cache.
if [[ "${DTN_BUILD:-}" == "single_layer" ]]; then
    julia --project=julia/scripts julia/scripts/build_dtn_singlelayer.jl 50 100 6 "$OUT/D50_Quant100_singlelayer_ord6.jld2"
    export KMDISK_DTN_FILE="$OUT/D50_Quant100_singlelayer_ord6.jld2"
fi
julia -t "$SLURM_CPUS_PER_TASK" --project=julia/scripts julia/scripts/run_matlab_baseline_check.jl "$OUT/sweep"
