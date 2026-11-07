#!/bin/bash
#SBATCH --job-name=km-graded-5hz
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:30:00
#SBATCH --output=output/graded/%x-%j.out
#SBATCH --error=output/graded/%x-%j.err
# C3b (all quality=good 5 Hz rows of 25x700, configs order 2 / ppw 24, 32 / q 1.04, same criterion as C3) and C7 (informational:
# al16 R, 5 Hz, R = 350, 1000, 3000, 8000 disk radii, ppw 32). From the worktree root:   mkdir -p output/graded && sbatch julia/cluster/graded_5hz.sh
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
OUT="output/graded/fivehz_${SLURM_JOB_ID}"
mkdir -p "$OUT/logs"
export KMDISK_LOG_DIR="$OUT/logs"
echo "code: $(git rev-parse --short HEAD)$(git status --porcelain | grep -q . && echo ' (dirty tree)') in $PWD; threads $SLURM_CPUS_PER_TASK; out $OUT"
julia -t "$SLURM_CPUS_PER_TASK" --project=julia/scripts julia/scripts/graded_pilot.jl "$OUT" parts="${PARTS:-EF}" budget_min="${BUDGET_MIN:-75}"
