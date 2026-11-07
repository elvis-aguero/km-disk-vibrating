#!/bin/bash
#SBATCH --job-name=dtn-pointsource
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:20:00
#SBATCH --output=output/tests/%x-%j.out
#SBATCH --error=output/tests/%x-%j.err
# Point-source (algebraic tail) test of the DtN matrix; scripts/dtn_pointsource_check.jl. ESTIMATE: under 5 minutes (the polar nr = 200 matrix takes about a minute).
set -euo pipefail
MAIN="${SLURM_SUBMIT_DIR:?run with sbatch}"
module load julia/1.12.4s-adme
mkdir -p "$MAIN/output/tests"; cd "$MAIN"
julia -t "$SLURM_CPUS_PER_TASK" --project=julia/scripts julia/scripts/dtn_pointsource_check.jl
