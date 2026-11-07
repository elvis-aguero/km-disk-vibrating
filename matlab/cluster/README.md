# matlab/cluster

- `run_matlab_tests.sh`: runs all MATLAB tests (`matlab/1_code/tests`) in one Slurm job (8 cores, 32G, 1 h).
  Submit from the repo root: `sbatch matlab/cluster/run_matlab_tests.sh`. Output: `output/matlab/*.out`.
- Needs `module load matlab/R2025a-bcih julia/1.12.4s-adme` (the script does it). JuliaParityTest calls `julia`.
- Julia project for the parity reference: `julia/scripts` (override: `KM_JULIA_PROJECT=<dir>`; CI uses `julia/`).
- Never run MATLAB or Julia on a login node.
