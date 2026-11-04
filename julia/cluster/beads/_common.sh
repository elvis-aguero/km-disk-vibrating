# Shared setup for the bead-campaign stage scripts. Source it; do not run it.
#   source "${SLURM_SUBMIT_DIR:?run with sbatch}/julia/cluster/beads/_common.sh"
# Provides: THREADS, BEADS_DIR, make_run_dir <stage>, and a pre-flight syntax check.

set -euo pipefail

module load julia/1.12.4s-adme

cd "${SLURM_SUBMIT_DIR:?run this script with sbatch from the repo root}"
BEADS_DIR=julia/cluster/beads
[[ -f "$BEADS_DIR/run_bead_sweep.jl" ]] || {
    echo "Error: run sbatch from the repo root ($BEADS_DIR/run_bead_sweep.jl not found in $PWD)." >&2
    exit 1
}
THREADS="${SLURM_CPUS_PER_TASK:-1}"

echo "Host:          $(hostname)"
echo "Repo root:     $PWD"
echo "Git commit:    $(git rev-parse --short HEAD)$(git status --porcelain | grep -q . && echo ' (dirty tree)' || true)"
echo "Julia threads: $THREADS"
echo "Mem (MB):      ${SLURM_MEM_PER_NODE:-unset}"
julia --version

# Pre-flight: fail in a second, not after package loading, on a syntax error.
julia --startup-file=no -e '
    for f in ARGS
        ex = Meta.parseall(read(f, String))
        any(a -> a isa Expr && a.head in (:error, :incomplete), ex.args) && error("syntax error in $f")
    end
    println("pre-flight parse OK: ", length(ARGS), " files")
' "$BEADS_DIR"/*.jl

# RUN_DIR may be set by the caller to resume an earlier run (finished cases are skipped).
make_run_dir() {
    local stage="$1"
    RUN_DIR="${RUN_DIR:-output/beads/$stage/${SLURM_JOB_NAME}_${SLURM_JOB_ID}_$(date +%Y%m%d_%H%M%S)}"
    mkdir -p "$RUN_DIR/logs"
    export KMDISK_LOG_DIR="$RUN_DIR/logs"
    echo "Run directory: $RUN_DIR"
}
