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
#   - --cpus-per-task / --mem / --time: 16 cores / 32GB / 24h is a starting guess, but read
#     this before assuming more cores helps. The nr=2500 domain's LU factors are dense
#     n=5002 matrices -- a FULL LU cache (one factor per cycle-phase, 40-60 of them) is
#     ~10-12GB. GMRES (the fallback when a cache doesn't fit) is NOT a mildly-slower
#     alternative here -- measured multiple *minutes per step* at this scale vs. a
#     fraction of a second for the LU-cached path -- so run_sweep caps how many cases run
#     *concurrently* to however many full LU caches the budget affords, rather than ever
#     dividing memory evenly across every thread and accepting GMRES for everyone (see
#     run_sweep's docstring in src/sweep.jl). Concretely, at --mem=32G only ONE case runs
#     at a time here (32GB * 0.7 safety margin / ~12GB per cache = 1) regardless of
#     --cpus-per-task -- extra threads just queue, they don't help. If you want N cases
#     running with :lu in parallel, request roughly `N * 12GB / 0.7` of memory (e.g.
#     --mem=140G for 4-way) with --cpus-per-task >= N; if you can't get that much memory,
#     lower --cpus-per-task to match what --mem actually supports instead of requesting
#     idle cores. MATLAB historically took ~200-1500s per case (that baseline used LU
#     caching too), so even fully serial (1-way) at --mem=32G should land in the 2-12 hour
#     range for the full 30-case sweep; tune --time from your own first run's actual
#     wall-clock (see the "elapsed_s" field in summary.csv).
#   - `module load julia`: this loads Julia via environment modules, the common setup on
#     clusters like Brown's Oscar. Run `module avail julia` on your cluster and adjust the
#     name/version, or delete this line entirely if Julia is already on PATH (e.g. via
#     juliaup in your home directory).
#
# Julia's Sys.free_memory() is NOT cgroup-aware — under SLURM it reports the whole node's
# free memory, not this job's actual --mem allocation, which caused a real OOM kill before
# available_memory_bytes() started preferring $SLURM_MEM_PER_NODE (which SLURM sets
# automatically from --mem above) over it. If you ever need to override the detected
# budget entirely (e.g. debugging on a shared node where even that isn't reliable), export
# KMDISK_RAM_BUDGET_BYTES before this script's `julia` invocations.
#
# #SBATCH --partition=batch
# #SBATCH --account=your_allocation_here

set -euo pipefail

module load julia 2>/dev/null || true

# --- locate the repo root ---
# NOT via ${BASH_SOURCE[0]}/dirname alone: under `sbatch`, Slurm copies this script into a
# spool directory on the compute node (e.g. /var/spool/slurmd/...) and executes it from
# there, so BASH_SOURCE[0] resolves to that spooled copy, not the file in your checkout —
# deriving the repo root from it lands on something like /var/spool, where mkdir below
# fails with "Permission denied" (confirmed against a real cluster run). Slurm sets
# SLURM_SUBMIT_DIR to wherever `sbatch` was actually invoked from, which is the right
# anchor for a batch job; only fall back to BASH_SOURCE when running this script directly
# with bash (no Slurm involved, so no spooling to work around).
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
    START_DIR="$SLURM_SUBMIT_DIR"
else
    START_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

# Walk up from START_DIR to find the repo root, rather than assuming a fixed number of
# ".." levels — this way it's correct whether sbatch was run from the repo root itself or
# from a subdirectory of it (e.g. `cd julia/cluster && sbatch run_baseline_sweep.sh`).
REPO_ROOT=""
CANDIDATE="$START_DIR"
while true; do
    if [[ -f "$CANDIDATE/julia/cluster/run_baseline_sweep.sh" ]]; then
        REPO_ROOT="$CANDIDATE"
        break
    fi
    [[ "$CANDIDATE" == "/" ]] && break
    CANDIDATE="$(dirname "$CANDIDATE")"
done

if [[ -z "$REPO_ROOT" ]]; then
    echo "Error: could not locate the km-disk-vibrating repo root starting from '$START_DIR'." >&2
    echo "Submit this job with 'sbatch' from inside your repo checkout (a subdirectory is" >&2
    echo "fine too), or run this script directly with bash from its own location." >&2
    exit 1
fi

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
#
# This job's own stdout/stderr (<job-name>-<job_id>.out/.err, per the #SBATCH
# --output/--error directives above) land in whatever directory you ran `sbatch` from —
# Slurm resolves relative --output/--error paths against the submission directory, not
# this script's own location — so they won't be inside output/<job_id>_<timestamp>/
# itself; that's fine, and lets you tail them while the job is still queued or starting,
# before this script has even created that directory.
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
