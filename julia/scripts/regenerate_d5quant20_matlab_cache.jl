#!/usr/bin/env julia
#= One-off fix for a real bug: matlab/1_code/D5Quant20/DTNnew345nr50D5refp10.mat was
generated with the literal argument D=5, not this domain's nominal bathDiameter=20 (see
migrate_dtn_caches.jl's history for how this was discovered). solve_motion.m's own
"Machine-specific patch for D5Quant20" loaded that mismatched matrix regardless, meaning the
DTN operator and the bulk-domain operators (Laplacian, pressure integral) were built on
DIFFERENT physical domain sizes for the same simulation whenever this domain was invoked --
domain radius = bathDiameter/2, so D=5 vs D=20 is a 4x difference in physical domain size,
not a rounding-level discrepancy.

Confirmed via a real usage audit that D5Quant20 is never used by any production/experimental
sweep in either MATLAB or Julia (both sweeper.m and run_sweep.jl use bathDiameter=100,
spatialResolution=50) -- it exists solely as a cheap test/parity-check fixture, so
regenerating it changes no historical experimental result.

This script generates the DTN natively at the CORRECT D=20 and writes a MATLAB-loadable
.mat file with the standard naming convention (DTNnew345nr50D20refp10.mat), so
solve_motion.m's normal sprintf-based lookup finds it without any special case. Run once
from the repo root: `julia julia/scripts/regenerate_d5quant20_matlab_cache.jl`

Block comment, not a docstring: see _bootstrap.jl's header for why a bare top-level string
directly followed by `include(...)` here would fail to parse. =#

include("_bootstrap.jl")
using MAT

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const CACHE_DIR = joinpath(REPO_ROOT, "matlab", "1_code", "D5Quant20")
const OLD_MISLABELED_PATH = joinpath(CACHE_DIR, "DTNnew345nr50D5refp10.mat")
const NEW_PATH = joinpath(CACHE_DIR, "DTNnew345nr50D20refp10.mat")

nr = 50
D = 20.0

@info "generating corrected DTN" nr = nr D = D
DTN = generate_dtn(nr, D)
size(DTN) == (nr, nr) || error("expected $(nr)x$(nr), got $(size(DTN))")

MAT.matwrite(NEW_PATH, Dict("DTNnew345" => DTN))
@info "wrote corrected MATLAB-loadable cache" path = NEW_PATH

if isfile(OLD_MISLABELED_PATH)
    rm(OLD_MISLABELED_PATH)
    @info "removed old mislabeled cache" path = OLD_MISLABELED_PATH
end
