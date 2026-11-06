#!/bin/bash
#SBATCH --job-name=bd-dtn-spectrum
#SBATCH --account=mbessa-condo
#SBATCH --qos=mbessa-condo
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:40:00
#SBATCH --output=output/beads/00_dtn/%x-%j.out
#SBATCH --error=output/beads/00_dtn/%x-%j.err
#
# Error of the reference and of the single-layer matrix vs wavelength (wave packets, exact Hankel). ESTIMATE: reference build
# ~66 s at nr = 200 (measured earlier), exact integrals ~1-2 min. 1 thread.
# From the repo root:  sbatch julia/cluster/beads/00i_dtn_spectrum.sh

set -euo pipefail
module load julia/1.12.4s-adme
cd "${SLURM_SUBMIT_DIR:?run this script with sbatch}"
echo "Host: $(hostname)   commit: $(git rev-parse --short HEAD)"
julia --startup-file=no julia/scripts/pilot_dtn_spectrum.jl 200
