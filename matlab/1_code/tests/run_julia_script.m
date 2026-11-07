function cmdout = run_julia_script(testCase, scriptName, varargin)
% RUN_JULIA_SCRIPT  Run julia/scripts/<scriptName> with args, assert exit 0, return stdout.
%   MATLAB's system() inherits MATLAB's LD_LIBRARY_PATH (bundled libs, crashes Julia): strip it with env -u.
repoRoot = fullfile(fileparts(mfilename('fullpath')), '..', '..', '..');
juliaScript = fullfile(repoRoot, 'julia', 'scripts', scriptName);
testCase.assertTrue(isfile(juliaScript), sprintf('not found: %s', juliaScript));
stderrFile = [tempname(), '.stderr.txt'];
cleanupStderr = onCleanup(@() delete_if_file(stderrFile)); %#ok<NASGU>
args = '';
for k = 1:numel(varargin)
    args = [args, ' "', varargin{k}, '"']; %#ok<AGROW>
end
cmd = sprintf('env -u LD_LIBRARY_PATH -u LD_PRELOAD julia "%s"%s 2>"%s"', juliaScript, args, stderrFile);
[status, cmdout] = system(cmd);
stderrText = '';
if isfile(stderrFile), stderrText = fileread(stderrFile); end
testCase.assertEqual(status, 0, sprintf('%s failed (exit %d)\nstdout: %s\nstderr: %s', scriptName, status, cmdout, stderrText));
end

function delete_if_file(f)
if isfile(f), delete(f); end
end
