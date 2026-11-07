#!/bin/bash
#SBATCH --job-name=bench-schur
#SBATCH --account=mbessa-condo
#SBATCH --qos=mbessa-condo
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=01:30:00
#SBATCH --output=output/tests/%x-%j.out
#SBATCH --error=output/tests/%x-%j.err
#
# Dense LU path vs Schur path on the same bath-driven case at nr = 2500 (5 periods at 30 Hz, Gamma 0.5, BDF2,
# 60 steps/period, 31 distinct matrices). Two separate processes so the peak RSS of each is its own. Prints time,
# peak memory and the max relative difference of the CoM traces. ESTIMATE of runtime: LU ~ 31 x ~1 min on 8 cores.
set -euo pipefail
MAIN="${SLURM_SUBMIT_DIR:?run with sbatch}"
module load julia/1.12.4s-adme
mkdir -p "$MAIN/output/tests"; cd "$MAIN"
OUT="$MAIN/output/tests/bench_schur_${SLURM_JOB_ID}"; mkdir -p "$OUT"
for s in lu schur; do
    julia -t "$SLURM_CPUS_PER_TASK" --project=julia/scripts julia/scripts/bench_schur.jl $s "$OUT/com_$s.csv"
done
julia --startup-file=no -e '
    rd(f) = [parse.(Float64, split(l, ",")) for l in readlines(f)]
    a = rd(ARGS[1]); b = rd(ARGS[2]); ca = [x[2] for x in a]; cb = [x[2] for x in b]
    println("max |CoM_schur - CoM_lu| / max|CoM_lu| = ", maximum(abs, ca .- cb) / maximum(abs, ca))
' "$OUT/com_lu.csv" "$OUT/com_schur.csv"
