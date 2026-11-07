#!/bin/bash
#SBATCH --job-name=km-matlab-tests
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT
#SBATCH --mail-user=elvis_vera@brown.edu
#SBATCH --output=output/matlab/%x-%j.out
#SBATCH --error=output/matlab/%x-%j.err
#
# MATLAB test suite (matlab/1_code/tests) on a compute node. JuliaParityTest shells out to julia, so both modules load.
#   CODE_DIR=<checkout> sbatch matlab/cluster/run_matlab_tests.sh      (default: the submit directory)
# Optional: TEST_FILTER=JuliaParityTest runs one test class only. Exit code is non-zero if any test fails.
# parity_reference.jl uses julia/scripts as project here (no GITHUB_ACTIONS); override with KM_JULIA_PROJECT.
set -euo pipefail
MAIN="${SLURM_SUBMIT_DIR:?run with sbatch}"
CODE_DIR="${CODE_DIR:-$MAIN}"
module load matlab/R2025a-bcih julia/1.12.4s-adme
mkdir -p "$MAIN/output/matlab"
cd "$CODE_DIR"
echo "code: $(git rev-parse --short HEAD)$(git status --porcelain | grep -q . && echo ' (dirty tree)') in $CODE_DIR"
export JULIA_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export TEST_FILTER="${TEST_FILTER:-}"
matlab -nodisplay -nosplash -nodesktop -batch "
addpath(genpath(fullfile(pwd,'matlab','1_code','simulation_code')));
import matlab.unittest.TestSuite; import matlab.unittest.TestRunner;
import matlab.unittest.plugins.TestReportPlugin;
suite = TestSuite.fromFolder(fullfile(pwd,'matlab','1_code','tests'));
f = getenv('TEST_FILTER'); if ~isempty(f), suite = suite(contains({suite.Name}, f)); end
runner = TestRunner.withTextOutput('OutputDetail', 3);
results = runner.run(suite);
disp(table(results));
if any([results.Failed] | [results.Incomplete]), exit(1); end
"
