# km-disk-vibrating

Simulates a disk floating on a vertically vibrating fluid bath (a Faraday-wave / walking-droplet
style setup): the bath surface is solved via a Dirichlet-to-Neumann (DTN) boundary operator, the
disk's center-of-mass motion is coupled to it, and the whole system is stepped forward in time
under sinusoidal forcing.

Two implementations live side by side:

- **`matlab/`** — the original MATLAB implementation. Left as-is; still the reference for what
  the experimental study's actual sweeps produced.
- **`julia/`** — an independent, from-scratch Julia port (not a transpiler output), with a
  cleaner module layout, a corrected convergence check, structured logging, and a real test
  suite. See `julia/README.md` for full setup/usage/testing instructions.

Both expose the same physical model through name-value-style APIs. This document is the API
reference for each, plus how they map onto each other.

## Repository layout

```
matlab/1_code/simulation_code/   solve_motion.m and its helpers (the core solver)
matlab/1_code/sweeper.m          (gamma, frequency) parameter sweep driver
matlab/1_code/tests/             MATLAB unit/integration/parity tests (matlab.unittest)
matlab/0_data/                   experimental data, sweep output (gitignored where large)
julia/src/                       FaradayDisk.jl package source
julia/scripts/                   CLI entry points (sweep runner, DTN generation, validation overlay, ...)
julia/test/                      unit/integration/slow test tiers
.github/workflows/                julia-ci.yml (Julia) and matlab-ci.yml (MATLAB + a live cross-language parity check)
```

## The physical model, in one place

Both implementations solve the same nondimensionalized problem:

- A disk of radius `diskRadius`, mass `diskMass`, floats on a circular bath of diameter
  `bathDiameter` (in disk radii), discretized with `spatialResolution` radial mesh points per
  disk radius (`nr = ceil(spatialResolution * bathDiameter / 2)` total radial nodes).
- The bath surface response is governed by a dense DTN operator matrix, precomputed once per
  domain (see "DTN operator" below) — this is the expensive part; solving it is a fast fixed-size
  linear system every time step.
- **Two independent forcing channels exist, and either or both can be active:**
  - **Disk forcing**: a sinusoidal force of amplitude `forceAmplitude` (dyne) and frequency
    `forceFrequency` (Hz) applied directly to the disk.
  - **Bath forcing**: the bath itself physically oscillates with amplitude `bathAmplitude` (cm)
    and frequency `bathFrequency` (Hz), `phaseDifference` (degrees) out of phase with the disk
    forcing.
  - **This experimental study's real convention is bath-driven**: `forceAmplitude = 0`, and
    `forceFrequency` is left at its default — it still sets the characteristic time unit
    (`T_unit = 1/forceFrequency`) even when it isn't a real force. `sweeper.m` / `run_sweep.jl`
    both use this convention. See "A note on forcing conventions" below — getting this backwards
    silently produces a plausible-looking but physically wrong result, and it happened once.
- `startStatic` solves for hydrostatic equilibrium once at t=0 as the initial condition (instead
  of starting from rest), matching how the real apparatus is set up.
- `earlyStop` (+ its `minPeriods`/`windowPeriods`/`checkEveryPeriods`/`amplitudeRelTol`/
  `phaseAbsTolDeg` knobs) stops the simulation once the center-of-mass oscillation's fitted
  amplitude and phase are stable across two consecutive multi-period windows, rather than running
  the full `simulationTime` — see "Why the convergence check looks the way it does" below.

## MATLAB API

### `solve_motion(...)` — `matlab/1_code/simulation_code/solve_motion.m`

```matlab
[t_s, CoM_cm, eta_history_cm] = solve_motion('name', value, ...)
```

All arguments are optional name-value pairs (a MATLAB `arguments` block); unless noted, units are
cgs.

| Name | Default | Units | Meaning |
|---|---|---|---|
| `diskRadius` | `0.2` | cm | disk radius |
| `diskMass` | `0.0283` | g | disk mass |
| `forceAmplitude` | `0` | dyne | disk-forcing amplitude |
| `forceFrequency` | `90` | Hz | disk-forcing frequency; also sets the characteristic time unit regardless of whether disk forcing is active |
| `bathAmplitude` | `0.01` | cm | bath-oscillation amplitude |
| `bathFrequency` | `90` | Hz | bath-oscillation frequency |
| `phaseDifference` | `-90` | degrees | phase of bath forcing relative to disk forcing |
| `bathDensity` | `1` | g/cm³ | fluid density |
| `bathSurfaceTension` | `72.20` | dyne/cm | fluid surface tension |
| `bathViscosity` | `0.978e-2` | Stokes | fluid kinematic viscosity |
| `g` | `981` | cm/s² | gravitational acceleration |
| `bathDiameter` | `100` | disk radii | bath domain diameter |
| `spatialResolution` | `50` | — | radial mesh points per disk radius (sets `nr`) |
| `temporalResolution` | `30` | — | time steps per dimensionless time unit |
| `simulationTime` | `10/90` | s | total time to simulate (upper bound if `earlyStop` fires first) |
| `debug_flag` | `true` | — | verbose logging + diary file |
| `solverType` | `"auto"` | — | `"auto"` (LU-cached if RAM allows, else GMRES), `"lu"`, or `"gmres"` |
| `gmresTolerance` | `1e-10` | — | GMRES convergence tolerance (when GMRES is active) |
| `startStatic` | `true` | — | start from hydrostatic equilibrium instead of rest |
| `earlyStop` | `true` | — | enable the convergence-based early stop |
| `minPeriods` | `10.0` | forcing periods | never declare convergence before this many periods |
| `windowPeriods` | `5.0` | forcing periods | length of each of the two compared convergence windows |
| `checkEveryPeriods` | `1.0` | forcing periods | convergence-check cadence |
| `amplitudeRelTol` | `0.02` | — | relative amplitude change tolerated between windows |
| `phaseAbsTolDeg` | `2.0` | degrees | absolute phase change tolerated between windows |

**Returns:**

| Value | Shape | Units | Meaning |
|---|---|---|---|
| `t_s` | `nsteps × 1` | s | simulation time at each recorded step |
| `CoM_cm` | `nsteps × 1` | cm | disk center-of-mass height at each recorded step |
| `eta_history_cm` | `nr × nsteps` | cm | bath surface elevation at each radial node, at each recorded step |

**Example** (bath-driven, matching `sweeper.m`'s convention):
```matlab
gamma = 0.2; freqHz = 30; omega = 2*pi*freqHz; g = 981;
[t_s, CoM_cm, eta] = solve_motion( ...
    'forceAmplitude', 0.0, 'forceFrequency', 90, ...
    'bathAmplitude', gamma*g/omega^2, 'bathFrequency', freqHz, ...
    'bathDiameter', 100, 'spatialResolution', 50, 'solverType', 'lu');
```

### `sweeper()` — `matlab/1_code/sweeper.m`

No arguments; edit the `sweep`/`fixed` structs at the top of the file directly, then run
`sweeper` from `matlab/1_code/`. Sweeps a Cartesian grid of `sweep.gamma × sweep.bathFrequency`
(dimensionless bath acceleration `Γ = A·ω²/g` at each swept frequency; converted to
`bathAmplitude = Γ·g/ω²` per case, bath-driven), running each case through `solve_motion` via
`parfor`.

Output, under `matlab/0_data/sweep_results/sweep_<timestamp>/`:
- One `CoM_gamma<g>_f<f>Hz.csv` per case: columns `time_s,CoM_cm,eta_boundary_cm` (the last is
  the max bath-surface elevation over the outer boundary region — a reflection-contamination
  indicator, not part of the physical result).
- One `summary_<timestamp>.csv`: columns `gamma,bathFrequency_Hz,bathAmplitude_cm,CoM_max_cm,
  CoM_rms_cm,elapsed_s,status`.

### Running the MATLAB tests

Requires MATLAB (`matlab.unittest`); there is no local-run substitute if you don't have it —
`.github/workflows/matlab-ci.yml` runs the full suite on every push touching `matlab/**` or
`julia/**` (auto-licensed for this public repo, no token needed). Locally, from a MATLAB session:
```matlab
addpath(genpath('matlab/1_code/simulation_code'));
results = runtests('matlab/1_code/tests');
```

## Julia API

### `SimulationParams` + `run_simulation` — `julia/src/parameters.jl`, `julia/src/simulate.jl`

```julia
using FaradayDisk
params = SimulationParams(; kwargs...)   # Base.@kwdef struct, all fields keyword/optional
result = run_simulation(params; dtn_registry=nothing, dtn=nothing, ram_budget_bytes=nothing, on_step=nothing)
```

`SimulationParams` fields mirror `solve_motion.m`'s arguments field-for-field, with three
differences: there is no `debug_flag` (use `julia/src/logging_setup.jl`'s `TeeLogger` instead —
see `julia/README.md`); `solverType` is a `Symbol` (`:auto`/`:lu`/`:gmres`) instead of a string;
and the five convergence knobs are grouped into a nested `convergence::ConvergenceOptions`
struct instead of five flat fields:

```julia
Base.@kwdef struct ConvergenceOptions
    enabled::Bool = true
    minPeriods::Float64 = 10.0
    windowPeriods::Float64 = 5.0
    checkEveryPeriods::Float64 = 1.0
    amplitudeRelTol::Float64 = 0.02
    phaseAbsTolDeg::Float64 = 2.0
end
```

| `SimulationParams` field | Default | Units |
|---|---|---|
| `diskRadius` | `0.2` | cm |
| `diskMass` | `0.0283` | g |
| `forceAmplitude` | `0.0` | dyne |
| `forceFrequency` | `90.0` | Hz |
| `bathAmplitude` | `0.01` | cm |
| `bathFrequency` | `90.0` | Hz |
| `phaseDifference` | `-90.0` | degrees |
| `bathDensity` | `1.0` | g/cm³ |
| `bathSurfaceTension` | `72.20` | dyne/cm |
| `bathViscosity` | `0.978e-2` | Stokes |
| `g` | `981.0` | cm/s² |
| `bathDiameter` | `100.0` | disk radii |
| `spatialResolution` | `50.0` | — |
| `temporalResolution` | `30` | steps per dimensionless time unit |
| `simulationTime` | `10/90` | s |
| `solverType` | `:auto` | `:auto`\|`:lu`\|`:gmres` |
| `gmresTolerance` | `1e-10` | — |
| `startStatic` | `true` | — |
| `earlyStop` | `true` | — |
| `convergence` | `ConvergenceOptions()` | — |

`with_overrides(params; kwargs...)` returns a copy of `params` with just the given fields
replaced — the idiomatic way to derive one case from a shared template (used throughout
`sweep.jl`) instead of re-specifying every field.

**`run_simulation` returns a `SimulationResult`:**

| Field | Type | Meaning |
|---|---|---|
| `params` | `SimulationParams` | the params actually used |
| `resolved_solver` | `Symbol` | `:lu` or `:gmres` — what `:auto` actually picked, and why (logged) |
| `status` | `Symbol` | `:ok` or `:diverged` |
| `t_s` | `Vector{Float64}` | simulation time at each recorded step (s) |
| `CoM_cm` | `Vector{Float64}` | disk center-of-mass height at each recorded step (cm) |
| `eta_history_cm` | `Matrix{Float64}` (`nr × nsteps`) | bath surface elevation (cm) |
| `convergence` | `ConvergenceResult` | the final convergence check outcome, recorded even if `earlyStop=false` |
| `elapsed_s` | `Float64` | wall-clock time |
| `lu_cache_hit_rate` | `Float64` | fraction of steps that reused a cached LU factorization |

A DTN operator matrix must be resolvable for `(spatialResolution, bathDiameter)` — either passed
directly via `dtn=`, or resolved from the on-disk cache (`dtn_registry.jl`; see
`julia/README.md`'s "Running things" for how to populate it). `resolve_dtn`/`load_dtn` are the
two functions involved if you want to manage this yourself.

**Example** (a cheap domain — `nr=50`; the real sweep's `bathDiameter=100`/`spatialResolution=50`
domain is `nr=2500` and native generation takes on the order of hours, see `julia/README.md`):
```julia
using FaradayDisk
gamma = 0.2; freqHz = 30.0; omega = 2*pi*freqHz
dtn = generate_dtn(50, 20.0)  # or resolve_dtn/load_dtn from an existing cache
params = SimulationParams(forceAmplitude=0.0, forceFrequency=90.0,
                           bathAmplitude=gamma*981.0/omega^2, bathFrequency=freqHz,
                           bathDiameter=20.0, spatialResolution=5.0, solverType=:lu)
result = run_simulation(params; dtn=dtn)
```

### `SweepSpec` + `run_sweep` — `julia/src/sweep.jl`

```julia
spec = SweepSpec(gammas=[...], bathFrequenciesHz=[...], fixed=SimulationParams(...))
results = run_sweep(spec, outdir; dtn_registry=nothing, dtn=nothing)  # Vector{SweepCaseResult}
```

Runs the Cartesian product of `gammas × bathFrequenciesHz` (bath-driven, same convention as
`sweeper.m`) across threads (`Threads.@threads`), writing the same CSV schema `sweeper.m` does
under `outdir/cases/` plus `outdir/summary.csv`, `outdir/sweep_metadata.toml`
(`is_bath_driven=true`, resolved solver, RAM budgeting), and structured `.log`/`.jsonl` logs.
`default_sweep_spec()` (`julia/scripts/run_sweep.jl`) reproduces `sweeper.m`'s exact fixed
parameters and grid. `SweepCaseResult` fields: `gamma`, `freqHz`, `bathAmplitude`, `CoM_max_cm`,
`CoM_rms_cm`, `elapsed_s`, `status` (`:ok`\|`:error`\|`:skipped`).

See `julia/README.md` for setup, running a sweep/single simulation/DTN generation from the
command line, the SLURM cluster driver, and the full test suite (`julia/test/`).

## MATLAB ↔ Julia parameter cross-reference

| Concept | MATLAB (`solve_motion.m`) | Julia (`SimulationParams`) |
|---|---|---|
| Physical parameters | flat name-value args | flat struct fields (identical names/units) |
| Solver choice | `solverType` (string) | `solverType` (`Symbol`) |
| Convergence knobs | 5 flat name-value args | nested `convergence::ConvergenceOptions` |
| Verbose logging | `debug_flag` (bool) | `Logging.jl` `TeeLogger` (see `julia/README.md`) |
| Sweep driver | `sweeper.m` | `julia/scripts/run_sweep.jl` / `SweepSpec` + `run_sweep` |
| DTN cache lookup | filename-formatting convention (`sprintf('DTNnew345nr%dD%drefp10.mat', ...)`) | explicit registry (`dtn_registry.jl`), no filename reconstruction |
| Output CSV schema | `time_s,CoM_cm,eta_boundary_cm` (case) / `gamma,bathFrequency_Hz,bathAmplitude_cm,CoM_max_cm,CoM_rms_cm,elapsed_s,status` (summary) | identical, plus a `sweep_metadata.toml` sidecar |

See `julia/README.md`'s "Module-to-MATLAB cross-reference" for the file-level (not just
API-level) mapping.


