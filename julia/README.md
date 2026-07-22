# FaradayDisk.jl

A Julia port of the MATLAB disk-on-vibrating-bath simulator in `../matlab/`. `matlab/` is
left untouched as reference/legacy; this is an independent, from-scratch reimplementation
— not a transpiler output — with a cleaner module structure, a corrected early-stop
convergence check, structured logging, and real test coverage.

## Important: how this port was authored

This port was written in an environment with **no way to install or run Julia** (the
official Julia binary host and package registry server were both blocked by network
policy, and no system package was available). Every file here was written by careful,
line-by-line translation of the MATLAB source and by reasoning through the math by hand
— **none of it has actually been executed**. `Pkg.instantiate()`, the test suite, and
every script need to be run for the first time by whoever picks this up next.

Two consequences worth knowing about before you start:

1. **The GMRES solver (`src/solver/gmres_solver.jl`) and the logging setup
   (`src/logging_setup.jl`) call `Krylov.jl`/`LoggingExtras.jl` APIs from memory.** If
   `Pkg.instantiate()` or the test suite reports a `MethodError`/`UndefKeywordError`
   originating in either file, that's almost certainly why — the fix is local to that one
   file. See the module docstring in each for specifics.
2. **`src/dtn/dtn_generator.jl` is the highest-risk file in the whole port** — a
   line-for-line transcription of `DTNVectorized.m`'s singular-integral quadrature (lots
   of hand-written polynomial coefficients, no way for me to check the arithmetic by
   running it). Its correctness is established entirely by
   `test/integration/test_dtn_golden_small.jl`, which compares its output against the
   legacy MATLAB `.mat` cache for the small `nr=50` domain. **If that test fails, this is
   the first and most likely place to look**, re-deriving from `DTNVectorized.m` by hand.
3. **There is no frozen/golden numeric regression fixture tier.** Normally a port like
   this would include "run once, freeze the answer, assert it never changes" tests — but
   producing that frozen answer requires actually running the corrected Julia code, which
   wasn't possible here. Instead, correctness for anything beyond the DTN generator is
   established through *physically-grounded* properties (does the center-of-mass
   oscillate at the driving frequency? does damping actually damp? do two independent
   linear solvers — LU-cached and GMRES — agree?) in `test/integration/`. See
   "Test coverage" below.

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
