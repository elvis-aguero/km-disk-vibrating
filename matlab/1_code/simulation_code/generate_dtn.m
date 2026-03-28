% Script to generate DtN matrix for spatialResolution=20, bathDiameter=20
nr = 200; % 20 points per radius * 20 radii / 2 = 200
D = 20;
fprintf('Starting DTN computation for nr=%d, D=%d...\n', nr, D);
tic;
DTNnew345 = DTNVectorized(nr, D);
time_elapsed = toc;
fprintf('DTN computation finished in %.2f minutes.\n', time_elapsed/60);

% Create folder and move the file
fold = fullfile('..', sprintf('D%dQuant%d', 20, 20));
if ~exist(fold, 'dir'); mkdir(fold); end
movefile(sprintf('DTNnew345nr%dD%drefp10.mat', nr, D), fold);
fprintf('Matrix saved to %s\n', fold);
exit;
