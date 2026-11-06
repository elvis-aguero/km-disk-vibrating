#!/bin/bash
#SBATCH --job-name=bd-resim
#SBATCH --account=mbessa-condo
#SBATCH --qos=mbessa-condo
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=780G
#SBATCH --time=06:00:00
#SBATCH --constraint=192core
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=elvis_vera@brown.edu
#SBATCH --output=output/beads/08_resim/%x-%j.out
#SBATCH --error=output/beads/08_resim/%x-%j.err
#
# Resimulate, in the LARGE domain 25 x 700, every production row (Gamma 0.5) whose wall-echo prediction was >= 1%
# (96 of 272 rows; audit: scratchpad audit.py). The DtN matrix is the single-layer construction (cubic, 6th-order FD)
# built in job 8459008 (nr = 8750); it matched production to <= 0.03% gain / 0.02 deg at 25 x 200 (jobs 8458894, 8459007).
# Two groups, run as two jobs (sbatch --export=ALL,GROUP=A ..., then GROUP=B), 128 threads / 780G each:
#   A: al32 rc (34 freqs), ptfe R, al16 R, nylon R (5 freqs each)        -> 49 cases
#   B: al32 R (13), ptfe rc (9), al16 rc (9), nylon rc (16)               -> 47 cases
# Resources (MEASURED in 8459008: 31 LU factorizations per case, 17,502^2 doubles each = ~76 GB per case, ~123 s per LU on
# 1 BLAS thread, stepping ~0.11 s/step; the sweep sets concurrency = min(threads, memory budget / LU cache) and
# BLAS threads = threads / concurrency). ESTIMATES: ~7-8 concurrent cases with ~16 BLAS threads each -> LU phase ~4-8 min
# per wave if LU scales (not measured), plus ~2-3 min stepping; 8 waves per group -> ~1.5-2 h per group. 6 h cap.
# A combination that fails does not stop the next one. Finished cases are skipped on resume (set RUN_DIR).
source "${SLURM_SUBMIT_DIR:?run with sbatch}/julia/cluster/beads/_common.sh"
make_run_dir 08_resim
GROUP="${GROUP:?set GROUP=A or GROUP=B}"
DTN_FILE="${DTN_FILE:-$(ls output/beads/00_dtn_big/bd-dtn-big-al32rc_8459008_*/D25_Quant200_singlelayer_ord6_D700.jld2 | head -1)}"
[[ -f "$DTN_FILE" ]] || { echo "DtN file not found: $DTN_FILE" >&2; exit 1; }
echo "DtN file: $DTN_FILE ($(stat -c %s "$DTN_FILE") bytes)   group: $GROUP"

run() {   # bead basis freqs
    KMDISK_DTN_FILE="$DTN_FILE" julia -t "$THREADS" "$BEADS_DIR/run_bead_sweep.jl" stage="08_resim_$1_$2" \
        outdir="$RUN_DIR/$1_$2" beads="$1" basis="$2" scheme=bdf2 gammas=0.5 freqs="$3" domain=25,700 tres=60 \
        || echo "combination $1/$2 exited with code $? (continuing)"
}
case "$GROUP" in
  A)
    run al32  rc 5,10,15,20,25,30,35,40,45,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,75,80,85,90
    run ptfe  R  5,10,15,20,25
    run al16  R  5,10,15,20,25
    run nylon R  5,10,15,20,25 ;;
  B)
    run al32  R  5,10,15,20,25,30,35,40,45,50,51,52,53
    run ptfe  rc 5,10,15,20,25,30,35,40,45
    run al16  rc 5,10,15,20,25,30,35,40,45
    run nylon rc 5,10,15,20,25,30,35,40,45,50,51,52,53,54,55,56 ;;
  *) echo "GROUP must be A or B" >&2; exit 1 ;;
esac
