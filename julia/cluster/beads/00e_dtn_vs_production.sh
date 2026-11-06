#!/bin/bash
#SBATCH --job-name=bd-dtn-vs-prod
#SBATCH --account=mbessa-condo
#SBATCH --qos=mbessa-condo
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=output/beads/00_dtn/%x-%j.out
#SBATCH --error=output/beads/00_dtn/%x-%j.err
#
# Single-layer DtN (cubic) vs the production cache D25_Quant200_nr2500. ESTIMATE: ~1.5e8 kernel evaluations
# per matrix (~10 s; measured 0.16 s at nr = 400, cost ~ nr^2); 16G for a few 2500^2 double matrices (50 MB each).
# From the repo root:  sbatch julia/cluster/beads/00e_dtn_vs_production.sh

set -euo pipefail
module load julia/1.12.4s-adme
cd "${SLURM_SUBMIT_DIR:?run this script with sbatch}"
echo "Host: $(hostname)   commit: $(git rev-parse --short HEAD)"
julia --startup-file=no julia/scripts/pilot_dtn_vs_production.jl julia/data/dtn_cache/D25_Quant200_nr2500.jld2 2500 200
