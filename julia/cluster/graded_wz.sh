#!/bin/bash
#SBATCH --job-name=km-graded-wz
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:30:00
#SBATCH --output=output/graded/%x-%j.out
#SBATCH --error=output/graded/%x-%j.err
# C3 (28 rows) and C3b (5 Hz rows) with the fine wall zone, configs (order 2, ppw 32 and 24, q 1.04, wall_zone 2 disk radii); same criteria as before.
# PARTS=I (C3) and J (C3b). Resources (ESTIMATES): about 60 cases, minutes each; 32G ample.
# Optional: PARTS=EF (default), BUDGET_MIN=75. Resources (ESTIMATES): C3b minutes; C7 at R = 8000 has nr about 3400 (n = 6800, LU cache
# about 12 GB), a few minutes per case; 32G is enough, 1.5 h cap.

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
OUT="output/graded/wz_${SLURM_JOB_ID}"
mkdir -p "$OUT/logs"
export KMDISK_LOG_DIR="$OUT/logs"
echo "code: $(git rev-parse --short HEAD)$(git status --porcelain | grep -q . && echo ' (dirty tree)') in $PWD; threads $SLURM_CPUS_PER_TASK; out $OUT"
julia -t "$SLURM_CPUS_PER_TASK" --project=julia/scripts julia/scripts/graded_pilot.jl "$OUT" parts="${PARTS:-IJ}" budget_min="${BUDGET_MIN:-75}"
