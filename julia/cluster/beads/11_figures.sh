#!/bin/bash
#SBATCH --job-name=bd-figs
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=output/beads/10_post/%x-%j.out
#SBATCH --error=output/beads/10_post/%x-%j.err
# Figures only, from a finished campaign (no re-processing): CAMP=output/beads/campaign_<...> sbatch ...
# A few hundred points: 1 core, 2G (ESTIMATE; the previous figure job peaked far below 4G).
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?run with sbatch}"
: "${CAMP:?set CAMP}"
module load py-matplotlib/3.10.3-2hcg
python3 figures/make_bead_figures.py "$CAMP"
