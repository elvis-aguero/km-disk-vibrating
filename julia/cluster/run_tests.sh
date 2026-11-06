#!/bin/bash
#SBATCH --job-name=km-tests
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=elvis_vera@brown.edu
#SBATCH --output=output/tests/%x-%j.out
#SBATCH --error=output/tests/%x-%j.err
#
# The fast test tier (unit + integration, julia/test/runtests.jl) on a compute node, for when GitHub Actions is not
# available. KMDISK_RUN_SLOW_TESTS=1 adds the slow tier (nr = 2500 comparisons and the 30-case validation sweep).
#   CODE_DIR=<checkout> sbatch julia/cluster/run_tests.sh      (default: the submit directory)
# Resources (ESTIMATE, not measured): the fast tier has nr <= 150 solves plus two ~40 s spot rows of the nr = 2500 cache;
# a few tens of minutes on 8 threads. 32G covers nr = 2500 LU caches (about 6 GB per case) with margin; 2 h cap.
set -euo pipefail
MAIN="${SLURM_SUBMIT_DIR:?run with sbatch}"
CODE_DIR="${CODE_DIR:-$MAIN}"
module load julia/1.12.4s-adme
mkdir -p "$MAIN/output/tests"
cd "$CODE_DIR"
echo "code: $(git rev-parse --short HEAD)$(git status --porcelain | grep -q . && echo ' (dirty tree)') in $CODE_DIR; threads $SLURM_CPUS_PER_TASK"
julia -t "$SLURM_CPUS_PER_TASK" --project=julia/scripts -e 'include("julia/scripts/_bootstrap.jl"); include("julia/test/runtests.jl")'
