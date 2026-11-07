#!/bin/bash
#SBATCH --job-name=bd-resim-s
#SBATCH --account=mbessa-condo
#SBATCH --qos=mbessa-condo
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=136G
#SBATCH --time=02:00:00
#SBATCH --array=0-95%12
#SBATCH --output=output/beads/08_resim/%x-%A_%a.out
#SBATCH --error=output/beads/08_resim/%x-%A_%a.err
#
# The 96 echo-contaminated rows (audit: wall-echo prediction >= 1% at 25 x 200), one case per array task, so each task is small enough to
# start in a gap. Domain 25 x 700, single-layer DtN (nr = 8750), Gamma 0.5, bdf2, 60 steps per period.
# Memory: one case needs ~76 GB measured (code estimate ~90 GB) for its 31 LU factorizations; the sweep budget is
# 0.7 x --mem = 95 GB, which is >= the estimate, so the LU path is chosen (GMRES fallback would log a warning).
# With one case per task the sweep runs concurrency 1 and gives all 8 BLAS threads to the LU phase.
# MEASURED (array 8482751): 9 to 19 min per case, peak 83 GB RSS. 2 h cap.
# %12 caps concurrent tasks: 12 x 136 G = 1.6 TB, under the QOS limit mem=2T.
source "${SLURM_SUBMIT_DIR:?run with sbatch}/julia/cluster/beads/_common.sh"
make_run_dir 08_resim

CASES=()
add() {   # bead basis freq...
    local b="$1" k="$2"; shift 2
    local f; for f in "$@"; do CASES+=("$b $k $f"); done
}
# group A
add al32  rc 5 10 15 20 25 30 35 40 45 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 75 80 85 90
add ptfe  R  5 10 15 20 25
add al16  R  5 10 15 20 25
add nylon R  5 10 15 20 25
# group B
add al32  R  5 10 15 20 25 30 35 40 45 50 51 52 53
add ptfe  rc 5 10 15 20 25 30 35 40 45
add al16  rc 5 10 15 20 25 30 35 40 45
add nylon rc 5 10 15 20 25 30 35 40 45 50 51 52 53 54 55 56
[[ ${#CASES[@]} -eq 96 ]] || { echo "expected 96 cases, got ${#CASES[@]}" >&2; exit 1; }

read -r BEAD BASIS FREQ <<< "${CASES[$SLURM_ARRAY_TASK_ID]}"
echo "Task ${SLURM_ARRAY_TASK_ID}: $BEAD $BASIS ${FREQ} Hz"

DTN_FILE="${DTN_FILE:-$(ls output/beads/00_dtn_big/bd-dtn-big-al32rc_8459008_*/D25_Quant200_singlelayer_ord6_D700.jld2 | head -1)}"
[[ -f "$DTN_FILE" ]] || { echo "DtN file not found: $DTN_FILE" >&2; exit 1; }
echo "DtN file: $DTN_FILE"

KMDISK_DTN_FILE="$DTN_FILE" julia -t "$THREADS" "$BEADS_DIR/run_bead_sweep.jl" stage="08_resim_${BEAD}_${BASIS}" \
    outdir="$RUN_DIR/${BEAD}_${BASIS}" beads="$BEAD" basis="$BASIS" scheme=bdf2 gammas=0.5 freqs="$FREQ" \
    domain=25,700 tres=60
