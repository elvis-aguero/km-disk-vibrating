# Cluster run output

Each `julia/cluster/run_baseline_sweep.sh` submission writes one folder here,
`<slurm_job_id>_<timestamp>/`, with `sweep/` (per-case CSVs, `summary.csv`, overlay
plots vs. experiment) and `logs/`. New runs are gitignored.

Two runs are committed as reference evidence (full 30-case sweep, nr=2500):

| Run | Result vs. experiment | Status |
|---|---|---|
| `4193275_20260722_214351/` | amplitude RMSE 0.085, phase RMSE 6.0° — `pass_primary = true` | **reference: current good baseline** |
| `4185395_20260722_163519/` | amplitude RMSE 0.22, phase RMSE 276° — `pass_primary = false` | superseded failing run, kept for comparison |

Quick look: open `4193275_20260722_214351/sweep/val_amp_sweep.png` and
`val_phase_sweep.png` (simulation vs. experiment, one colour per Γ).
