#!/bin/bash
#SBATCH --job-name=bd-post-figs
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=output/beads/10_post/%x-%j.out
#SBATCH --error=output/beads/10_post/%x-%j.err
#
# Re-run the post-processing of a finished campaign with the corrected rules (rounding slack in
# the convergence check; new `subharmonic` label), re-run the gate on the corrected data, and
# make the figures. No simulation is re-run. One core and 4G: the Julia part loads the package
# (1.7 GB measured in the smoke test) and reads a few hundred small CSV files; the figures are a
# few hundred points. 30 min is a wide cap (package loading dominates).
#
#     CAMP=output/beads/campaign_<...> sbatch julia/cluster/beads/10_repostprocess_figures.sh
# The original files are kept as *.v1.csv and gate_v1_original/.

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?run with sbatch}"
: "${CAMP:?set CAMP to the campaign directory}"

module load julia/1.12.4s-adme

# Pre-flight: fail in a second on a syntax error in any Julia file of the campaign.
julia --startup-file=no -e '
    for f in ARGS
        ex = Meta.parseall(read(f, String))
        any(a -> a isa Expr && a.head in (:error, :incomplete), ex.args) && error("syntax error in $f")
    end
    println("pre-flight parse OK: ", length(ARGS), " files")
' julia/cluster/beads/*.jl

julia julia/cluster/beads/repostprocess.jl \
    "$CAMP/03_accel" "$CAMP/03b_accel_gt1" "$CAMP/04_sweep_R" "$CAMP/05_sweep_rc" \
    "$CAMP/06_numerics/euler_tres60" "$CAMP/06_numerics/bdf2_tres30"

[[ -d "$CAMP/gate_v1_original" ]] || cp -r "$CAMP/gate" "$CAMP/gate_v1_original"
/usr/bin/python3 -I julia/cluster/beads/analyze_gate.py "$CAMP" --out-dir "$CAMP/gate" || echo "gate on corrected data exited with code $? (3 = dependent, 4 = not enough good data)"

module load py-matplotlib/3.10.3-2hcg
python3 figures/make_bead_figures.py "$CAMP"
