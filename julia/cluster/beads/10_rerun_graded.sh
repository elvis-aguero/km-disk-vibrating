#!/bin/bash
#SBATCH --job-name=rerun-graded
#SBATCH --account=mbessa-condo
#SBATCH --qos=mbessa-condo
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=output/beads/10_rerun/%x-%j.out
#SBATCH --error=output/beads/10_rerun/%x-%j.err
#
# Re-run, with the current stack (graded grid: ppw 32, q 1.04, order 2, wall zone 2 + Schur solver via :auto), (1) the 51 rows of the
# deliverable CSV whose wall echo was not negligible a priori (25 production rows at 25 x 700, 26 gamma-scan/numerics rows at 25 x 200),
# each in a domain where the a-priori echo is < 1e-3 (R = 3.5 L_e, 350 a <= R <= 12000 a), and (2) 30 redundant rows where no change is
# expected (12 production 25 x 700, 12 production 25 x 200, 6 gamma-scan; fixed seed; al32 rc 50-70 Hz left out), same domain as the CSV row.
# Case list: output/beads/10_rerun/cases.tsv (written by the selection script). One case per line, run one after the other, no array.
# ESTIMATES (not measured): graded cases take 0.3-100 s each at nr 200-3500; the largest domains (R = 12000 a at 5 Hz) have nr of a few
# thousand and take minutes; Schur cache < 8 GB. Total about 1-2 h; 4 h cap.
source "${SLURM_SUBMIT_DIR:?run with sbatch}/julia/cluster/beads/_common.sh"
make_run_dir 10_rerun
CASES=output/beads/10_rerun/cases.tsv
[[ -f "$CASES" ]] || { echo "case list not found: $CASES" >&2; exit 1; }
tail -n +2 "$CASES" | while IFS=$'\t' read -r ID KIND BEAD BASIS GAMMA FREQ TRES D ROLE REFDOM; do
    echo "=== $ID $KIND $BEAD $BASIS gamma=$GAMMA f=$FREQ tres=$TRES domain=25x$D"
    julia -t "$THREADS" "$BEADS_DIR/run_bead_sweep.jl" stage="10_rerun_$ID" outdir="$RUN_DIR/$ID" beads="$BEAD" basis="$BASIS" \
        scheme=bdf2 gammas="$GAMMA" freqs="$FREQ" domain=25,"$D" tres="$TRES" grid=32,1.04,2 \
        || echo "case $ID exited with code $? (continuing)"
done
