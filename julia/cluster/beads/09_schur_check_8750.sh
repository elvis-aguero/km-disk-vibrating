#!/bin/bash
#SBATCH --job-name=schur-8750
#SBATCH --account=mbessa-condo
#SBATCH --qos=mbessa-condo
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=01:30:00
#SBATCH --output=output/beads/09_schur_check/%x-%j.out
#SBATCH --error=output/beads/09_schur_check/%x-%j.err
#
# Does the Schur step solver reproduce the dense LU solver at PRODUCTION size (25 x 700, nr = 8750)? Three cases that were
# already run with the LU solver (array 8482751): al16 rc 40 Hz, al32 rc 60 Hz, ptfe R 20 Hz (all quality good). Same DtN file,
# same parameters; only KMDISK_SOLVER=schur differs. Compare gain_phase.csv with the LU rows afterwards (compare script in the
# report). Memory ESTIMATE: 31 factorizations of 8751^2 = 19 GB + S0, G1 (1.2 GB) + transient Y (1.2 GB); budget = 0.7 x 64G = 45 GB >= estimate.
# Time ESTIMATE: LU phase ~1/8 of the 9-19 min per case of the dense path, so a few minutes per case; 1.5 h cap.
source "${SLURM_SUBMIT_DIR:?run with sbatch}/julia/cluster/beads/_common.sh"
make_run_dir 09_schur_check
DTN_FILE="$(ls output/beads/00_dtn_big/bd-dtn-big-al32rc_8459008_*/D25_Quant200_singlelayer_ord6_D700.jld2 | head -1)"
[[ -f "$DTN_FILE" ]] || { echo "DtN file not found" >&2; exit 1; }
for c in "al16 rc 40" "al32 rc 60" "ptfe R 20"; do
    read -r BEAD BASIS FREQ <<< "$c"
    echo "=== $BEAD $BASIS $FREQ Hz"
    KMDISK_SOLVER=schur KMDISK_DTN_FILE="$DTN_FILE" julia -t "$THREADS" "$BEADS_DIR/run_bead_sweep.jl" stage="09_schur_${BEAD}_${BASIS}" \
        outdir="$RUN_DIR/${BEAD}_${BASIS}" beads="$BEAD" basis="$BASIS" scheme=bdf2 gammas=0.5 freqs="$FREQ" domain=25,700 tres=60 \
        || echo "case $c exited with code $? (continuing)"
done
