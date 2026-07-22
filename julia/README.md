# FaradayDisk.jl

A Julia port of the MATLAB disk-on-vibrating-bath simulator in `../matlab/`. `matlab/` is
left untouched as reference/legacy; this is an independent, from-scratch reimplementation
— not a transpiler output — with a cleaner module structure, a corrected early-stop
convergence check, structured logging, and real test coverage.

## Important: how this port was authored, and its current status

This port was originally written in an environment with **no way to install or run
Julia** (the official Julia binary host and package registry server were both blocked by
network policy, and no system package was available) — every file was written by careful
line-by-line translation of the MATLAB source and by reasoning through the math by hand,
with nothing actually executed before the first push.

It has since been run for real against GitHub Actions CI (`.github/workflows/julia-ci.yml`)
and iterated on fix-by-fix. **As of the latest CI run, the fast test tier is effectively
green** (104/105 tests passing) except for one known, understood, non-blocking gap
documented below. Bugs actually found and fixed this way included: a struct-constructor
signature collision, a docstring `$VAR` string-interpolation error, a wrong sign
convention documented (not implemented) for `fit_oscillation`, a test-tuning issue in the
convergence-check synthetic signal, an off-by-one in two tests' repo-root path
computation, three scripts missing `using Dates`/`using Printf` (a dependency's internal
`using` doesn't propagate to a separate script), a too-small test domain that let boundary
reflection contaminate results, `Krylov.GmresWorkspace`'s real constructor signature
(twice — first the positional args, then the storage-type argument), and a Julia
"world-age" bug from calling `include()` inside a running closure. None of that is
speculative anymore — it's what CI actually reported, fixed one push at a time.

One gap remains open:

- **`src/dtn/dtn_generator.jl`'s comparison against the legacy MATLAB cache
  (`test/integration/test_dtn_golden_small.jl`) is not fully resolved.** The test
  initially showed a 75% relative discrepancy at the domain size its folder name implies
  (`D5Quant20` → bathDiameter=20), but every formula/index/divisor in the generator was
  re-verified character-by-character against `DTNVectorized.m` with no error found. The
  actual cause turned out to be the legacy cache's *own* filename ambiguity — it's named
  `DTNnew345nr50D5refp10.mat` ("D5", not "D20"), and comparing against
  `generate_dtn(50, 5.0)` instead dropped the discrepancy to ~2.15%, strongly suggesting
  the legacy cache itself was generated with the "wrong" domain size relative to its
  folder name (a pre-existing MATLAB-side inconsistency, not something to fix in
  `matlab/`, which is left untouched). That remaining ~2.15% is larger than pure
  batched-vs-unbatched floating-point summation-order differences should produce (normally
  ~1e-10, not ~1e-2) and is **not** fully explained. See the module docstring in
  `dtn_generator.jl` and the test itself for what's been ruled out and the leading
  suspect (an `l_vals` angular-quadrature range edge effect). This is the one item left
  for whoever next has a working Julia install to dig into further — everything else
  driving actual simulation behavior (the solver, convergence check, sweep orchestration,
  I/O) is verified passing.
- There is deliberately no frozen/golden numeric regression fixture tier beyond the DTN
  comparison above — producing a frozen "run once, assert it never changes" reference
  would have required trusting a from-scratch, never-executed implementation as its own
  ground truth. Correctness for the solver, convergence check, and sweep pipeline instead
  comes from *physically-grounded* properties (does the center-of-mass oscillate at the
  driving frequency? does damping actually damp? do two independent linear solvers —
  LU-cached and GMRES — agree?) in `test/integration/`, all of which are now passing
  against real execution. See "Test coverage" below.

None of this should be read as "this code doesn't work" — the underlying math was traced
by hand against the MATLAB source and cross-checked internally (e.g. `materialize!` is
tested against a from-scratch matrix rebuild for consistency) wherever that was possible
without execution. It should be read as "verify this before trusting it," which is
exactly what CI (and the test suite) is set up to do on the very first push.

## Setup

```sh
cd julia
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Running things

**A single simulation:**
```julia
using FaradayDisk
# First, a DTN cache must exist for your (spatialResolution, bathDiameter) — either
# migrate an existing MATLAB cache (see below) or generate one natively:
#   julia --project=scripts scripts/generate_dtn.jl 5 20
params = SimulationParams(bathDiameter = 20.0, spatialResolution = 5.0, forceFrequency = 20.0)
result = run_simulation(params)  # resolves the DTN cache via the registry
```

**The (gamma, frequency) parameter sweep** (port of `sweeper.m`):
```sh
julia -t auto --project=julia julia/scripts/run_sweep.jl
```
Output goes to `julia/data/results/sweeps/sweep_<timestamp>/` — one CSV per case under
`cases/`, plus `summary.csv`.

**Generating a DTN cache natively** (port of `generate_dtn.m` / `DTNVectorized.m`):
```sh
julia -t auto --project=scripts julia/scripts/generate_dtn.jl <spatialResolution> <bathDiameter>
```
⚠️ The largest domains in this repo (`nr=2500`) took MATLAB's `parfor`-based generator on
the order of *hours*. Expect the same order of magnitude here regardless of language.

**Migrating the existing MATLAB `.mat` DTN caches to Julia-native JLD2** (one-off, run
once):
```sh
julia --project=scripts julia/scripts/migrate_dtn_caches.jl
```
This also cross-validates the small `D5Quant20` cache against the native generator (see
the migration script's docstring for why that particular cache's filename is ambiguous
enough to warrant a check).

**Validation overlay vs. experimental data** (port of `overlay_validation.py`, numeric
part only — see its docstring for why plotting was intentionally left out for now):
```sh
julia --project=scripts julia/scripts/run_validation_overlay.jl <sweep_dir>
```

**Tests:**
```sh
julia --project=julia julia/test/runtests.jl               # fast tier (seconds)
KMDISK_RUN_SLOW_TESTS=1 julia --project=julia julia/test/runtests.jl   # + slow tier (hours)
```

## Test coverage

- `test/unit/` — pure numerics: nondimensionalization, the (tridiagonal) Laplacian, the
  system-matrix assembly/patch consistency, the shared least-squares oscillation fit, and
  synthetic-signal tests for the corrected convergence check that directly encode the fix
  for the incident this port was commissioned after (see "Why the early-stop check was
  rewritten" below).
- `test/integration/` — the DTN-generator golden comparison (nr=50), an LU-cache
  hit-rate invariant, LU-vs-GMRES agreement, a small end-to-end sweep (schema + overlay
  sanity), and 9 physically-grounded checks (frequency-locking via a direct DFT,
  CoM-vs-wave amplitude sanity, static-equilibrium fixed-point, damping direction,
  linear-response amplitude scaling, no-blow-up, the divergence guard, and a
  boundary-reflection meta-check on the fixtures themselves) — all run on a small,
  fast domain.
- `test/slow/` — the nr=2500 DTN comparisons and the full 30-case MATLAB-baseline
  validation sweep. Gated behind `KMDISK_RUN_SLOW_TESTS=1`; never run in default CI.

## Why the early-stop check was rewritten

The MATLAB source's rolling-average convergence check stopped a simulation once
`std(rolling 1-period mean)/|mean|` dropped below a tolerance. Audited and found broken:
because the static-equilibrium initial condition removes most of the transient, that
ratio's denominator is dominated by a large, non-oscillating equilibrium offset, making
the ratio trivially small almost immediately — verified against a real 30-case sweep,
where every single case stopped at exactly the earliest allowed check point (3 forcing
periods), degrading the amplitude match to experimental data by roughly 10x and the phase
match by roughly 50x. `src/convergence.jl` replaces it with a check that fits oscillation
amplitude and phase (the same regression `overlay_validation.py` already used) over two
consecutive multi-period windows, and only declares convergence when both are stable
across them — see that file's docstring for the full story and
`test/unit/test_convergence_check.jl` for the synthetic-signal tests that pin the fix.

## Data layout

```
julia/data/
  dtn_cache/manifest.toml, *.jld2   — DTN operator caches, resolved by explicit registry
                                      lookup (src/dtn/dtn_registry.jl), not by
                                      reconstructing a filename (the MATLAB approach,
                                      which is already inconsistent for two of the three
                                      legacy caches — see migrate_dtn_caches.jl)
  results/runs/<run_id>/            — single-run output (gitignored; regenerate on demand)
  results/sweeps/<sweep_id>/        — sweep output (gitignored; regenerate on demand)
```

Large or pre-existing data (the `matlab/1_code/D*Quant*/*.mat` DTN caches, the
experimental measurement CSVs, the old sweep results used as MATLAB-baseline fixtures) is
read directly from `../matlab/0_data/...` / `../matlab/1_code/...` rather than duplicated
into `julia/` — see the path resolution in `scripts/migrate_dtn_caches.jl` and
`scripts/run_validation_overlay.jl`.

## Module-to-MATLAB cross-reference

| Julia | MATLAB |
|---|---|
| `src/domain.jl` | `simulation_code/domainMaker.m` |
| `src/dtn/dtn_generator.jl` | `simulation_code/DTNVectorized.m` |
| `src/dtn/dtn_registry.jl` | (new — replaces the filename-formatting lookup in `solve_motion.m`) |
| `src/parameters.jl` | `solve_motion.m`'s `arguments` block + nondimensionalization |
| `src/system_matrix.jl` | `simulation_code/buildSystemMatrix.m` |
| `src/solver/lu_cache.jl`, `solver_dispatch.jl` | the LU-caching logic in `advance_one_step.m` / `computeLU.m` |
| `src/solver/gmres_solver.jl` | the GMRES branch in `advance_one_step.m` |
| `src/solver/static_equilibrium.jl` | the `startStatic` block in `solve_motion.m` |
| `src/convergence.jl` | the (rewritten) early-stop check in `solve_motion.m`, plus the fit in `overlay_validation.py` |
| `src/state.jl` | `advance_one_step.m` |
| `src/simulate.jl` | `solve_motion.m` |
| `src/sweep.jl` | `sweeper.m` |
| `src/io.jl` | `results_saver` in `solve_motion.m` / the CSV writers in `sweeper.m` |
