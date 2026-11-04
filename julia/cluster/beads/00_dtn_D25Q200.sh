#!/bin/bash
#SBATCH --job-name=bd-dtn-D25Q200
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=output/beads/00_dtn/%x-%j.out
#SBATCH --error=output/beads/00_dtn/%x-%j.err
#
# Stage 00: regenerate the DTN operator for spatialResolution=25, bathDiameter=200 (nr=2500).
# The old cache was incomplete (rows 564:2500 were zero) and was deleted. From the repo root:
#
#     sbatch julia/cluster/beads/00_dtn_D25Q200.sh
#
# History: the first attempt (job 8262774, 32 THREADS, 8G) was cancelled after 39 min. Measured
# on its node: 1.00 busy core of 32 over 20 s, memory growing 3.3 -> 3.5 GB in ~2 min. Probable
# cause, NOT proven (see regenerate_dtn.jl): Julia's GC needs all threads at a safepoint, and
# each row ends in a long allocation-free loop. What IS measured (benchmark 8264099): 8 threads
# are 2.66x faster than serial on the heaviest rows, 8 separate processes 7.49x, with identical
# results. So this version uses 24 worker PROCESSES.
#
# Resources:
#   - CPUs: 24 worker processes, one per core; the master process mostly waits.
#   - Memory: MEASURED in job 8264099 (heaviest rows of this exact problem): every worker peaked
#     at 1.4 GiB. 24 workers = 33.6 GiB, plus the master (1.4 GiB after the benchmark's thread
#     run, 3.5 GiB in the self-test) = about 35-37 GiB. 48G leaves about 30% margin; the first
#     request, 40G, had under 10%. A worker's heap can still grow over ~100 rows instead of the
#     benchmark's 8; the campaign's final report records the real peak (sacct MaxRSS).
#   - Time: MEASURED in job 8264099 (serial, one thread): row 4 = 41.7 s, row 1250 = 52.0 s,
#     row 2500 = 58.1 s, so about 35.6 core-hours for all 2,500 rows (my first estimate, 14
#     core-hours, was 2.5x too low). / 24 cores = about 1.5 h if the scaling is ideal; 6 h is
#     the cap (4x margin).
#   - No checkpoint: a time-limit kill loses everything, so the job stops with a clear
#     "progress" log (every 5% of rows, in the .err file) that shows the ETA early.

set -euo pipefail

module load julia/1.12.4s-adme

cd "${SLURM_SUBMIT_DIR:?run this script with sbatch}"
[[ -f julia/cluster/beads/regenerate_dtn.jl ]] || {
    echo "Error: run sbatch from the repo root (julia/cluster/beads/regenerate_dtn.jl not found in $PWD)." >&2
    exit 1
}

echo "Host:          $(hostname)"
echo "Repo root:     $PWD"
echo "Git commit:    $(git rev-parse --short HEAD)"
echo "Workers:       $SLURM_CPUS_PER_TASK"
julia --version

julia --startup-file=no -e '
    ex = Meta.parseall(read(ARGS[1], String))
    any(a -> a isa Expr && a.head in (:error, :incomplete), ex.args) && error("syntax error in $(ARGS[1])")
    println("pre-flight parse OK: ", ARGS[1])
' julia/cluster/beads/regenerate_dtn.jl

DTN_WORKERS="$SLURM_CPUS_PER_TASK" julia julia/cluster/beads/regenerate_dtn.jl 25 200

echo
echo "Done. Cache:    julia/data/dtn_cache/D25_Quant200_nr2500.jld2"
echo "      Manifest: julia/data/dtn_cache/manifest.toml"
