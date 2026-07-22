#!/bin/bash
#SBATCH --job-name=km-disk-baseline-sweep
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#
# Submits the full 30-case (gamma x frequency) validation sweep at the D50Quant100 domain
# (nr=2500 — the slow, real-accuracy domain, not a test fixture) and compares the result
# against the digitized experimental data, producing the amplitude/phase RMSE numbers and
# overlay plots. This is the Julia-port equivalent of running MATLAB's sweeper.m +
# overlay_validation.py end to end.
#
# Single node, single task, many CPU threads: this port's parallelism is Threads.jl
# (Threads.@threads over the sweep's Cartesian grid), not Distributed.jl/MPI, so there is
# nothing to gain from --nodes>1 or --ntasks>1 — see julia/README.md's design notes.
#
# Adjust for your cluster before submitting:
#   - --partition / --account: uncomment and fill in below; every cluster's accounting
#     setup is different and this can't be guessed.
#   - --cpus-per-task / --mem / --time: 16 cores / 32GB / 24h is a starting guess. The
#     nr=2500 domain's system matrices are ~2500x2500 dense LU factors (~50MB each,
#     cached per cycle-phase); MATLAB historically took ~200-1500s per case serially, so
#     30 cases threaded across N cores should land well under a day, but tune from your
#     own first run's actual wall-clock (see the "elapsed_s" field in summary.csv).
#   - `module load julia`: this loads Julia via environment modules, the common setup on
#     clusters like Brown's Oscar. Run `module avail julia` on your cluster and adjust the
#     name/version, or delete this line entirely if Julia is already on PATH (e.g. via
#     juliaup in your home directory).
#
# #SBATCH --partition=batch
# #SBATCH --account=your_allocation_here

set -euo pipefail

module load julia 2>/dev/null || true

# --- locate the repo root, independent of the submission directory ---
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$REPO_ROOT"

THREADS="${SLURM_CPUS_PER_TASK:-4}"
echo "Repo root: $REPO_ROOT"
echo "Julia threads: $THREADS"
julia --version

# --- output layout: everything from this run under one timestamped/job-tagged folder ---
#
# output/<job_id>_<timestamp>/
#   sweep/                        -- run_matlab_baseline_check.jl's own output
#     cases/gamma..._f...Hz.csv     -- per-case (time_s, CoM_cm, eta_boundary_cm) traces
#     summary.csv                    -- one row per case: gamma, freq, amplitude, status, ...
#     val_amp_<name>.{png,pdf}       -- amplitude-ratio overlay vs. experiment
#     val_phase_<name>.{png,pdf}     -- phase-difference overlay vs. experiment
#   logs/                          -- structured .log (human-readable) + .jsonl
#                                     (machine-parseable) logs, one pair per run
#   <job-name>-<job_id>.out/.err   -- this SLURM job's own stdout/stderr (written directly
#                                     by SLURM next to this script, per the #SBATCH
#                                     --output/--error directives above; not moved into
#                                     output/ so you can tail them while the job is queued
#                                     or still starting up)
RUN_ID="${SLURM_JOB_ID:-local}_$(date +%Y%m%d_%H%M%S)"
OUTDIR="$REPO_ROOT/output/$RUN_ID"
mkdir -p "$OUTDIR/sweep" "$OUTDIR/logs"
export KMDISK_LOG_DIR="$OUTDIR/logs"
echo "Output directory: $OUTDIR"

# --- DTN operator cache: ported once from the legacy MATLAB .mat caches, then reused
# from the JLD2 cache on every run (never regenerated) ---
# See julia/scripts/migrate_dtn_caches.jl's docstring: this reads the existing
# matlab/1_code/D*Quant*/DTNnew345*.mat files (already committed in the repo) and writes
# julia/data/dtn_cache/manifest.toml + *.jld2 once. Every later run, including this one,
# resolves the cached matrix via the registry (src/dtn/dtn_registry.jl) — the sweep
# below never touches MATLAB or regenerates anything.
MANIFEST="$REPO_ROOT/julia/data/dtn_cache/manifest.toml"
if [[ -f "$MANIFEST" ]]; then
    echo "DTN cache manifest already present at $MANIFEST — skipping migration."
else
    echo "No DTN cache manifest found — migrating legacy MATLAB .mat caches (one-off)..."
    julia -t "$THREADS" --project=julia/scripts julia/scripts/migrate_dtn_caches.jl
fi

# --- the actual run: full sweep + RMSE-vs-experiment + overlay plots ---
julia -t "$THREADS" --project=julia/scripts julia/scripts/run_matlab_baseline_check.jl "$OUTDIR/sweep"

echo ""
echo "Done."
echo "Results:  $OUTDIR/sweep/summary.csv, $OUTDIR/sweep/cases/*.csv"
echo "Plots:    $OUTDIR/sweep/val_amp_*.{png,pdf}, $OUTDIR/sweep/val_phase_*.{png,pdf}"
echo "Logs:     $OUTDIR/logs/"
