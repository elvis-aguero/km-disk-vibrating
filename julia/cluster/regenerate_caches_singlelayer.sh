#!/bin/bash
#SBATCH --job-name=km-dtn-sl-caches
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=elvis_vera@brown.edu
#SBATCH --output=output/tests/%x-%j.out
#SBATCH --error=output/tests/%x-%j.err
#
# Builds the single-layer DtN caches (the default method) for the registered domains and registers them next to the
# polar ones: 5x20 (nr 50), 25x200 and 50x100 (nr 2500). Each is a few seconds to ~25 s (24 s at nr = 2500, measured, 1 thread).
# Memory: G is (nr+3) x nr doubles = 50 MB at nr = 2500, a few copies. From the repo root:
#   sbatch julia/cluster/regenerate_caches_singlelayer.sh
set -euo pipefail
module load julia/1.12.4s-adme
cd "${SLURM_SUBMIT_DIR:?run with sbatch}"
for d in "5 20" "25 200" "50 100"; do
    julia -t "$SLURM_CPUS_PER_TASK" --project=julia/scripts julia/scripts/generate_dtn.jl $d single_layer
done
cat julia/data/dtn_cache/manifest.toml
