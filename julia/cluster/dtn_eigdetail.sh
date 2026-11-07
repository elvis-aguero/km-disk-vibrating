#!/bin/bash
#SBATCH --job-name=dtn-eigdetail
#SBATCH --account=mbessa-condo
#SBATCH --qos=mbessa-condo
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=output/tests/%x-%j.out
#SBATCH --error=output/tests/%x-%j.err
# Eigenvalues with eigenvectors at nr = 1000 and 2500 (eigvals took 1 and 3 s in job 8489793; eigen is a few times that).
set -euo pipefail
MAIN="${SLURM_SUBMIT_DIR:?run with sbatch}"
module load julia/1.12.4s-adme
mkdir -p "$MAIN/output/tests"; cd "$MAIN"
for nr in 1000 2500; do julia -t "$SLURM_CPUS_PER_TASK" --project=julia/scripts julia/scripts/dtn_eigdetail.jl $nr 6; done
