#!/bin/bash
#SBATCH --job-name=km-graded-pilot
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=output/graded/%x-%j.out
#SBATCH --error=output/graded/%x-%j.err
#
# Graded radial grid pilot (docs/graded_grid.md): C1 (DtN vs exact packets), C2 (Laplacian order), C3 (28 bead cases x 9
# configurations vs the uniform rows of the campaign), C6 (cost; one uniform nr = 2500 case timed in the same job).
# From the repo (worktree) root:   mkdir -p output/graded && sbatch julia/cluster/graded_pilot.sh
# Optional: PARTS=ABCD (default), BUDGET_MIN=95 (no new C3 configuration after this many minutes).
# Resources (ESTIMATES, not measured): C1/C2 a few minutes; each C3 case 5-20 s at nr 100-1000 (one LU cache of 31 factorizations
# of n = 2 nr + 2 <= 2000: < 1 GB); 28 cases x 9 configurations = 252 cases, about 30-60 min; the uniform nr = 2500 case about
# 3-5 min and ~7 GB of LU cache. 32G leaves a large margin; the 2 h cap and BUDGET_MIN stop the job before the cap.
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
OUT="output/graded/pilot_${SLURM_JOB_ID}"
mkdir -p "$OUT/logs"
export KMDISK_LOG_DIR="$OUT/logs"
echo "code: $(git rev-parse --short HEAD)$(git status --porcelain | grep -q . && echo ' (dirty tree)') in $PWD; threads $SLURM_CPUS_PER_TASK; out $OUT"
julia -t "$SLURM_CPUS_PER_TASK" --project=julia/scripts julia/scripts/graded_pilot.jl "$OUT" parts="${PARTS:-ABCD}" budget_min="${BUDGET_MIN:-95}"
