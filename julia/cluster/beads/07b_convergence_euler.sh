#!/bin/bash
#SBATCH --job-name=bd-conv-euler
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=70G
#SBATCH --time=02:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=elvis_vera@brown.edu
#SBATCH --output=output/beads/07_convergence/%x-%j.out
#SBATCH --error=output/beads/07_convergence/%x-%j.err
#
# Euler levels of the convergence study, run in parallel with 07_convergence.sh (job 8296798)
# on a separate node. 4 cases per level fill 4 threads. Memory: 07_convergence.sh peaked at
# 57 GB with 4 concurrent cases (MaxRSS, measured); 70 GB is a margin over that.
# Same cases and code as 07_convergence.sh; only the levels differ.
source "${SLURM_SUBMIT_DIR:?run with sbatch}/julia/cluster/beads/_common.sh"
make_run_dir 07_convergence
for t in 60 120 240; do
    julia -t "$THREADS" "$BEADS_DIR/run_bead_sweep.jl" stage="07_euler_tres$t" outdir="$RUN_DIR/euler_tres$t" \
        beads=ptfe basis=R scheme=euler gammas=0.5,1.0 freqs=40,90 domain=25,200 tres="$t" \
        || echo "level euler/$t exited with code $? (continuing)"
done
