#!/bin/bash
#SBATCH --job-name=km-graded-c7b
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=output/graded/%x-%j.out
#SBATCH --error=output/graded/%x-%j.err
# C7b (informational): al32 rc at 52, 60, 70 Hz, graded (order 2, ppw 32, q 1.04) at R = 100, 200, 350, 1000, 3000 disk radii, printed next to the
# uniform 25x200 and 25x700 rows of the campaign CSV. From the worktree root:   mkdir -p output/graded && sbatch julia/cluster/graded_c7b.sh
# Optional: BUDGET_MIN=45. Resources (ESTIMATES): 15 cases, nr 100-1000, a few minutes each; 32G is ample.

set -euo pipefail
MAIN="${SLURM_SUBMIT_DIR:?run with sbatch from the repo root}"
cd "$MAIN"
module load julia/1.12.4s-adme
# A fresh worktree has no (gitignored) Manifest: reuse the main checkout's, pointing the FaradayDisk path at this tree.
MAN=julia/scripts/Manifest.toml
if [[ ! -f "$MAN" ]]; then
    sed "s|^path = \".*julia/scripts/\.\.\"|path = \"$PWD/julia/scripts/..\"|" \
        /oscar/data/dharri15/eaguerov/Github/km-disk-vibrating/julia/scripts/Manifest.toml > "$MAN"
    grep -n '^path = ' "$MAN" | head -3
fi
OUT="output/graded/c7b_${SLURM_JOB_ID}"
mkdir -p "$OUT/logs"
export KMDISK_LOG_DIR="$OUT/logs"
echo "code: $(git rev-parse --short HEAD)$(git status --porcelain | grep -q . && echo ' (dirty tree)') in $PWD; threads $SLURM_CPUS_PER_TASK; out $OUT"
julia -t "$SLURM_CPUS_PER_TASK" --project=julia/scripts julia/scripts/graded_pilot.jl "$OUT" parts="${PARTS:-G}" budget_min="${BUDGET_MIN:-45}"
