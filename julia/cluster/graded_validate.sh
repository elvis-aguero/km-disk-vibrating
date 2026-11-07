#!/bin/bash
#SBATCH --job-name=km-graded-valid
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=output/graded/%x-%j.out
#SBATCH --error=output/graded/%x-%j.err
#
# (1) The fast test tier (362 existing tests + the new graded-grid tests). (2) C4 of docs/graded_grid.md: the 30-case experimental
# validation (50 x 100 domain, disk R = 0.2 cm, BDF2) on graded grids, one run per configuration "<ppw>,<q>,<laplacian_order>".
# From the repo (worktree) root:   mkdir -p output/graded && sbatch julia/cluster/graded_validate.sh
# Optional: C4_CONFIGS="12,1.04,4 16,1.04,2" (default below), SKIP_TESTS=1, BUDGET_MIN=100 (no new C4 configuration after this).
# Criterion (fixed before the run): amplitude RMSE <= 0.022578 and wrapped-phase RMSE <= 3.1557 deg (uniform single-layer: 0.022135 / 3.0557).
# Resources (ESTIMATES, not measured): the test tier takes tens of minutes (run_tests.sh); one C4 configuration, 30 cases at nr 300-700,
# about 5-10 min including the plots; 32G is far more than needed.
set -uo pipefail
MAIN="${SLURM_SUBMIT_DIR:?run with sbatch from the repo root}"
cd "$MAIN"
module load julia/1.12.4s-adme
MAN=julia/scripts/Manifest.toml
if [[ ! -f "$MAN" ]]; then
    sed "s|^path = \".*julia/scripts/\.\.\"|path = \"$PWD/julia/scripts/..\"|" \
        /oscar/data/dharri15/eaguerov/Github/km-disk-vibrating/julia/scripts/Manifest.toml > "$MAN"
fi
OUT="output/graded/validate_${SLURM_JOB_ID}"
mkdir -p "$OUT/logs"
export KMDISK_LOG_DIR="$OUT/logs"
T0=$(date +%s)
echo "code: $(git rev-parse --short HEAD)$(git status --porcelain | grep -q . && echo ' (dirty tree)') in $PWD; threads $SLURM_CPUS_PER_TASK; out $OUT"
if [[ "${SKIP_TESTS:-0}" != "1" ]]; then
    echo "=== fast test tier"
    julia -t "$SLURM_CPUS_PER_TASK" --project=julia/scripts -e 'include("julia/scripts/_bootstrap.jl"); include("julia/test/runtests.jl")' 2>&1 | tee "$OUT/tests.log" | tail -n 60
    echo "test tier exit status: ${PIPESTATUS[0]}"
fi
for cfg in ${C4_CONFIGS:-"12,1.04,4 16,1.04,4 16,1.04,2 24,1.04,2"}; do
    if (( ($(date +%s) - T0) / 60 > ${BUDGET_MIN:-100} )); then echo "BUDGET REACHED: C4 config $cfg not started"; continue; fi
    echo "=== C4 config (ppw,q,order) = $cfg"
    KMDISK_TIME_SCHEME=bdf2 KMDISK_GRADED="$cfg" julia -t "$SLURM_CPUS_PER_TASK" --project=julia/scripts \
        julia/scripts/run_matlab_baseline_check.jl "$OUT/c4_${cfg//,/_}" 2>&1 | tail -n 8
    echo "C4 criterion: amplitude RMSE <= 0.022578 and phase RMSE <= 3.1557 deg (config $cfg)"
done
echo "done after $(( ($(date +%s) - T0) / 60 )) min"
