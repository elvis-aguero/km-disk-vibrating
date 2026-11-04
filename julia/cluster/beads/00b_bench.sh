#!/bin/bash
#SBATCH --job-name=bd-bench
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00:20:00
#SBATCH --output=output/beads/00_dtn/%x-%j.out
#SBATCH --error=output/beads/00_dtn/%x-%j.err
#
# Stage 00b: falsification test, threads vs processes on the real DTN rows (see
# bench_threads_vs_procs.jl). Expected time: serial rows ~1-2 min, thread run ~0.5-3 min,
# process run ~0.5-1 min plus worker start-up. 8 cores = N; 16G: 8 threads' row temporaries
# (a few hundred MB each) + 8 workers (~1 GB each) + master (3.5 GB measured) is about 14 GB
# at the very most; the job prints the measured peaks, and `sacct -j <id> -o MaxRSS` too.

set -euo pipefail
module load julia/1.12.4s-adme
cd "${SLURM_SUBMIT_DIR:?run with sbatch}"
echo "Host: $(hostname)  git: $(git rev-parse --short HEAD)  N=$SLURM_CPUS_PER_TASK"
BENCH_N="$SLURM_CPUS_PER_TASK" julia -t "$SLURM_CPUS_PER_TASK" julia/cluster/beads/bench_threads_vs_procs.jl
