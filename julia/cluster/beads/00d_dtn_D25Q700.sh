#!/bin/bash
#SBATCH --job-name=bd-dtn-D25Q700
#SBATCH --account=mbessa-condo
#SBATCH --qos=mbessa-condo
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=780G
#SBATCH --time=20:00:00
#SBATCH --constraint=192core
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=elvis_vera@brown.edu
#SBATCH --output=output/beads/00_dtn/%x-%j.out
#SBATCH --error=output/beads/00_dtn/%x-%j.err
#
# DtN operator for a LARGER domain: spatialResolution = 25 (same dr = 0.04 disk radii as the
# production cache), bathDiameter = 700 disk radii -> nr = 8750. Purpose: a physical domain
# radius of 350 radii (6.0 cm for Al 1/32" at rc, 1.71 cm today), so the wall echo is below
# ~1% down to ~35-40 Hz. The matrix is dimensionless: it serves every bead.
# From the repo root:  sbatch julia/cluster/beads/00d_dtn_D25Q700.sh
#
# Resources (ESTIMATES scaled from MEASUREMENTS of job 8264099 / 00_dtn_D25Q200.sh at nr = 2500:
# 42-58 s per row, 1.4 GiB per worker, 35.6 core-hours in all):
#   - Row cost ~ nr^2 (far-field loop ~ nr*refp radial nodes x pi*nr*refp angles): x12.25 ->
#     ~510-710 s per row x 8750 rows = ~1,460 core-hours. / 128 workers = ~11.5 h if scaling is
#     ideal (the 192-core Genoa nodes may be faster per core than the nodes of the benchmark).
#     Limit 20 h (2x margin). NO checkpoint: a time-limit kill loses everything. The progress log
#     (every 5% of rows, .err file) gives the ETA early.
#   - Memory ~ nr (the per-row temporaries scale with the angular count): 1.4 GiB x 3.5 =
#     ~4.9 GiB per worker x 128 = ~630 GiB, plus the master (matrix 0.6 GB, rows in flight).
#     780G is a ~24% margin. Check sacct MaxRSS afterwards.
#   - The script writes <name>.partial, reloads it, rechecks 5 rows (incl. the heaviest) and only
#     then renames and registers the cache (see regenerate_dtn.jl).

set -euo pipefail
module load julia/1.12.4s-adme
cd "${SLURM_SUBMIT_DIR:?run this script with sbatch}"
[[ -f julia/cluster/beads/regenerate_dtn.jl ]] || { echo "Error: run sbatch from the repo root." >&2; exit 1; }

echo "Host: $(hostname)   commit: $(git rev-parse --short HEAD)   workers: $SLURM_CPUS_PER_TASK"
julia --version
julia --startup-file=no -e '
    ex = Meta.parseall(read(ARGS[1], String))
    any(a -> a isa Expr && a.head in (:error, :incomplete), ex.args) && error("syntax error in $(ARGS[1])")
    println("pre-flight parse OK: ", ARGS[1])
' julia/cluster/beads/regenerate_dtn.jl

DTN_WORKERS="$SLURM_CPUS_PER_TASK" julia julia/cluster/beads/regenerate_dtn.jl 25 700

echo "Done. Cache: julia/data/dtn_cache/D25_Quant700_nr8750.jld2 (registered in manifest.toml)"
