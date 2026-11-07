#!/bin/bash
#SBATCH --job-name=km-graded-first
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00:40:00
#SBATCH --output=output/graded/%x-%j.out
#SBATCH --error=output/graded/%x-%j.err
#
# First, cheapest job: the two new test files (map, uniform bit-identity, Laplacian order, DtN vs exact, small end-to-end case)
# and C1 + C2 tables (graded_pilot.jl parts A and B). From the worktree root:
#   mkdir -p output/graded && sbatch julia/cluster/graded_first.sh
# Resources (ESTIMATES, not measured): minutes; the end-to-end test builds one uniform nr = 600 DtN matrix.
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:?run with sbatch from the worktree root}"
module load julia/1.12.4s-adme
MAN=julia/scripts/Manifest.toml
if [[ ! -f "$MAN" ]]; then
    sed "s|^path = \".*julia/scripts/\.\.\"|path = \"$PWD/julia/scripts/..\"|" \
        /oscar/data/dharri15/eaguerov/Github/km-disk-vibrating/julia/scripts/Manifest.toml > "$MAN"
fi
OUT="output/graded/first_${SLURM_JOB_ID}"
mkdir -p "$OUT/logs"
export KMDISK_LOG_DIR="$OUT/logs"
echo "code: $(git rev-parse --short HEAD) in $PWD"
echo "=== new tests"
julia -t "$SLURM_CPUS_PER_TASK" --project=julia/scripts -e '
include("julia/scripts/_bootstrap.jl")
using Test, LinearAlgebra
@testset "graded" begin
    include("julia/test/unit/test_graded_grid.jl")
    include("julia/test/integration/test_graded_grid_dtn.jl")
end'
echo "tests exit status: $?"
echo "=== C1, C2"
julia -t "$SLURM_CPUS_PER_TASK" --project=julia/scripts julia/scripts/graded_pilot.jl "$OUT" parts=AB
