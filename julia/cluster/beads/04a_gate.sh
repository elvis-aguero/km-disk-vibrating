#!/bin/bash
#SBATCH --job-name=bd-gate
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:10:00
#SBATCH --output=output/beads/04a_gate/%x-%j.out
#SBATCH --error=output/beads/04a_gate/%x-%j.err
#
# Stage 04a: the gate. Decides, from stages 03 and 06, whether the PTFE gain / phase-lag
# curve depends on the acceleration. Pure Python (standard library): a few seconds.
# Exit code 0 = independent -> the main sweeps start. 3 = dependent, 4 = not enough good data
# -> the job FAILS on purpose, so the dependent jobs are cancelled and you get an email.
# Rules and thresholds: see the docstring of analyze_gate.py. The report is written to
# <campaign>/gate/gate_report.md.
#
#     CAMP=output/beads/campaign_<...> sbatch julia/cluster/beads/04a_gate.sh

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?run with sbatch}"
: "${CAMP:?set CAMP to the campaign directory}"
/usr/bin/python3 -I julia/cluster/beads/analyze_gate.py "$CAMP" --bead ptfe --ref-gamma "${MAIN_GAMMA:-0.5}" --main-gamma "${MAIN_GAMMA:-0.5}"
