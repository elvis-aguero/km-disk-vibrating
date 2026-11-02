#!/bin/bash
#SBATCH --job-name=bd-dtn-selftest
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=6G
#SBATCH --time=00:20:00
#SBATCH --output=output/beads/00_dtn/%x-%j.out
#SBATCH --error=output/beads/00_dtn/%x-%j.err
#
# Stage 00a: self-test of the multi-process DTN generator on the SMALL domain (5, 20): nr=50.
# It runs `regenerate_dtn.jl` in test mode (DTN_CACHE_DIR set): the matrix is written to this
# job's own folder, the manifest is not touched, and the result must equal the registered
# serial library generator to a relative difference below 1e-12 (the registered D5Quant20 cache is stale, see regenerate_dtn.jl). The generator was checked against the
# MATLAB reference to machine precision (see dtn_generator.jl), so this tests the worker
# assembly, the write/reload/verify path and the progress logging, in about a minute.
#
# Resources: nr=50 is tiny. 4 workers, 6G (master ~3.5 GB measured for the big job's master,
# plus 4 small workers), 20 min cap for package loading and worker start-up.

set -euo pipefail

module load julia/1.12.4s-adme

cd "${SLURM_SUBMIT_DIR:?run this script with sbatch}"
TEST_DIR="output/beads/00_dtn/selftest_${SLURM_JOB_ID}"
mkdir -p "$TEST_DIR"

echo "Host: $(hostname)  git: $(git rev-parse --short HEAD)  workers: $SLURM_CPUS_PER_TASK"
julia --version

DTN_WORKERS="$SLURM_CPUS_PER_TASK" DTN_CACHE_DIR="$TEST_DIR" julia julia/cluster/beads/regenerate_dtn.jl 5 20

echo "Self-test passed."
