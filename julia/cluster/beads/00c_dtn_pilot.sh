#!/bin/bash
#SBATCH --job-name=bd-dtn-pilot
#SBATCH --account=mbessa-condo
#SBATCH --qos=mbessa-condo
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=output/beads/00_dtn/%x-%j.out
#SBATCH --error=output/beads/00_dtn/%x-%j.err
#
# Pilot of option D: DtN matrix as -Laplacian of a single-layer potential (no angular quadrature).
# Compares it with generate_dtn and with the exact Gaussian solution at nr = 100..400 (linear vs cubic interpolation, Richardson, error split).
# ESTIMATE, not measured: the reference generator is the slow part (row cost ~ nr^2; about 1 s per
# row at nr = 400, so ~7 min), the new method seconds. 1 thread on purpose: parallel threads
# in generate_dtn were measured to run at 1.00 busy core (see 00_dtn_D25Q200.sh).
# From the repo root:  sbatch julia/cluster/beads/00c_dtn_pilot.sh

set -euo pipefail
module load julia/1.12.4s-adme
cd "${SLURM_SUBMIT_DIR:?run this script with sbatch}"
echo "Host: $(hostname)   commit: $(git rev-parse --short HEAD)"
julia --startup-file=no julia/scripts/pilot_dtn_singlelayer.jl 100 200
