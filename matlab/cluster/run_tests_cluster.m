function run_tests_cluster()
% Runs matlab/1_code/tests from the repo root (the current folder). TEST_FILTER (env) keeps tests whose name contains it.
% Exit status 1 if any test fails or is incomplete. Called by matlab/cluster/run_matlab_tests.sh with a one-line -batch command,
% because the cluster's matlab wrapper drops a multi-line quoted -batch argument.
root = pwd;
addpath(genpath(fullfile(root, 'matlab', '1_code', 'simulation_code')));
import matlab.unittest.TestSuite; import matlab.unittest.TestRunner;
suite = TestSuite.fromFolder(fullfile(root, 'matlab', '1_code', 'tests'));
f = getenv('TEST_FILTER');
if ~isempty(f), suite = suite(contains({suite.Name}, f)); end
runner = TestRunner.withTextOutput('OutputDetail', 3);
results = runner.run(suite);
disp(table(results));
if any([results.Failed] | [results.Incomplete]), exit(1); end
end
