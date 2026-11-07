#!/bin/bash
#SBATCH --job-name=dtn-eigcheck
#SBATCH --account=mbessa-condo
#SBATCH --qos=mbessa-condo
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --output=output/tests/%x-%j.out
#SBATCH --error=output/tests/%x-%j.err
#
# Eigenvalues of the single-layer DtN matrix at nr = 1000, 2500, 5000 (and 8750), looking for Re < 0 or complex
# eigenvalues (growing modes; task 40 of array 8482751 diverged after 29 periods). Dense nonsymmetric eigvals is O(n^3):
# ESTIMATE, not measured: n = 2500 a few minutes, 8750 about 1-3 h with LAPACK's mostly serial Hessenberg QR. Memory:
# a few n^2 doubles, < 10 GB at 8750. Output prints after each size, so a time-out still keeps the smaller sizes.
set -euo pipefail
MAIN="${SLURM_SUBMIT_DIR:?run with sbatch}"
module load julia/1.12.4s-adme
mkdir -p "$MAIN/output/tests"
cd "$MAIN"
julia -t "$SLURM_CPUS_PER_TASK" --project=julia/scripts julia/scripts/dtn_eigcheck.jl 1000 2500 5000 8750
