# FaradayDisk.jl

A Julia port of the MATLAB disk-on-vibrating-bath simulator in `../matlab/`. `matlab/` is
left untouched as reference/legacy; this is an independent, from-scratch reimplementation
— not a transpiler output — with a cleaner module structure, a corrected early-stop
convergence check, structured logging, and real test coverage.

## Verification status

The full test suite (`julia/test/`) passes against both GitHub Actions CI
(`.github/workflows/julia-ci.yml`) and a local Julia install, including the DTN generator's
comparison against a real MATLAB-generated reference cache, which matches to machine
precision (~2e-15 relative difference — see `src/dtn/dtn_generator.jl` and
`scripts/migrate_dtn_caches.jl` for what that comparison actually checks and why).
There is deliberately no frozen/golden numeric regression tier beyond that DTN comparison
— correctness for the solver, convergence check, and sweep pipeline instead comes from
*physically-grounded* properties in `test/integration/` (does the center-of-mass oscillate
at the driving frequency? does damping actually damp? do two independent linear solvers —
LU-cached and GMRES — agree?). See "Test coverage" below.

One exception: **`scripts/live_plot.jl` (opt-in GLMakie live plotting) has never been run,
by CI or locally** — GLMakie needs a real display/GPU context that neither has. Treat it as
unverified until someone with local Julia and a display confirms it.

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

**Validation overlay vs. experimental data** (port of `overlay_validation.py`, including
plotting):
```sh
julia --project=scripts julia/scripts/run_validation_overlay.jl <sweep_dir>
```
Produces `val_amp_<sweep_name>.{png,pdf}` / `val_phase_<sweep_name>.{png,pdf}` in
`<sweep_dir>` (amplitude ratio and phase-difference overlays vs. the digitized
experimental data, one panel per, colored by gamma) via CairoMakie, a headless-safe
rasterizing/vector backend — this path is exercised by
`test/integration/test_plotting.jl` in the default (fast, CI) test tier.

**Opt-in live plotting during a run** (`scripts/live_plot.jl`, port of `solve_motion.m`'s
`drawnow` debug-plot block):
```julia
include("julia/scripts/live_plot.jl")   # bootstraps its own env (live/Project.toml) — see below

cb = make_live_plot_callback(problem)          # problem from build_problem(params, dtn)
run_simulation(params; dtn = dtn, on_step = cb)
```
Opens one GLMakie window and updates it in place every step (or every `refresh_every`
steps): the bath surface profile near the disk, the disk's own height, and the bath-drive
amplitude as a guide line. GLMakie needs a real display/GPU context, so this is **not
exercised by CI** (see "Verification status" above). It also lives in its own environment
(`scripts/live/Project.toml`), separate from the other scripts' shared one — it's the only
dependency in this repo that needs a real OpenGL/GLFW-capable display to even *install*,
so keeping it out of the environment every other script (migration, sweep, overlay) also
instantiates means a machine without one (a headless cluster node, say) never has to build
it just to run something that never touches it.
Never pass this callback into `run_sweep`/`run_sweep_case` — sweep cases run concurrently
via `Threads.@threads`, and touching a GLMakie window from more than one thread isn't safe.

**Tests:**
```sh
julia --project=julia julia/test/runtests.jl               # fast tier (seconds)
KMDISK_RUN_SLOW_TESTS=1 julia --project=julia julia/test/runtests.jl   # + slow tier (hours)
```

## Running the full validation sweep on a SLURM cluster

The full 30-case (gamma x frequency) sweep at the real `D50Quant100` domain (`nr=2500`) is
the genuinely slow, real-accuracy run — historically ~200-1500s/case in MATLAB — and is
what actually produces a comparison against the digitized experimental data. It's meant to
be submitted as a batch job, not run at a terminal:

```sh
sbatch julia/cluster/run_baseline_sweep.sh
```

Read the comments at the top of that script before your first submission — it has
placeholder `#SBATCH --partition`/`--account` lines you need to fill in for your cluster,
and a `module load julia` line that assumes an environment-modules setup (e.g. Brown's
Oscar); adjust or delete it if Julia is provided differently on yours. It:

1. Migrates the legacy MATLAB DTN `.mat` caches to JLD2 once (skipped on every later
   submission once `julia/data/dtn_cache/manifest.toml` exists — this never regenerates
   the operator, only loads the cached one, per the registry design in
   `src/dtn/dtn_registry.jl`).
2. Runs `scripts/run_matlab_baseline_check.jl`, multithreaded across the job's allocated
   cores (`Threads.@threads` over the sweep's Cartesian grid — see the design notes on why
   this port uses `Threads.jl` rather than `Distributed.jl`/MPI, so there's no benefit to
   requesting more than one node or task).
3. Writes everything into one folder per run, `output/<job_id>_<timestamp>/`, gitignored:
   ```
   output/<job_id>_<timestamp>/
     sweep/
       cases/gamma..._f...Hz.csv    — per-case (time_s, CoM_cm, eta_boundary_cm) traces
       summary.csv                    — one row/case: gamma, freq, amplitude, status, ...
       val_amp_<name>.{png,pdf}       — amplitude-ratio overlay vs. experiment
       val_phase_<name>.{png,pdf}     — phase-difference overlay vs. experiment
     logs/                          — structured .log + .jsonl logs for this run
   ```
   SLURM's own stdout/stderr (`<job-name>-<job_id>.out`/`.err`) land next to the script
   itself, per its `#SBATCH --output`/`--error` directives, so they're visible while the
   job is still queued or starting.

The `pass_primary`/`pass_incident_guard` fields the script prints (and logs) are exactly
the RMSE-vs-experiment thresholds described in "Why the early-stop check was rewritten"
below — a passing run means the corrected convergence check and solver reproduce the
known-good pre-regression MATLAB baseline (~0.10 amplitude RMSE, ~5.7° phase RMSE), not
the buggy sweep's (~1.09, ~285°).

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
| `scripts/run_validation_overlay.jl` | `overlay_validation.py` |
| `scripts/live_plot.jl` | the `drawnow` debug-plot block in `solve_motion.m` |
