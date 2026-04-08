% sweeper.m
% Sweeps bath frequency at constant dimensionless bath acceleration
%   Gamma = A * omega^2 / g   (A = bath amplitude, omega = 2*pi*f)
%
% For each (Gamma, frequency) pair the bath amplitude is computed as
%   A = Gamma * g / omega^2
% so that Gamma is held fixed while frequency varies.
%
% Outputs: one CSV per run in sweep_results/, plus a summary CSV.
%
% Usage: sweeper  (run from the 1_code directory)

function sweeper()

%% ---- SWEEP PARAMETERS --------------------------------------------------
sweep.gamma         = [0.05, 0.2, 0.5];        % dimensionless bath acceleration Gamma = A*w^2/g
sweep.bathFrequency = 10:10:100;   % Hz

%% ---- FIXED PARAMETERS --------------------------------------------------
g_cgs                    = 981;       % cm/s^2

fixed.diskRadius         = 0.2;       % cm
fixed.diskMass           = 0.0283;    % g
fixed.forceAmplitude     = 0.0;       % dynes
fixed.forceFrequency     = 90;        % Hz  (sets time unit; bath drives motion)
fixed.phaseDifference    = -90;       % degrees
fixed.bathDensity        = 1.175;     % g/cm^3
fixed.bathSurfaceTension = 66.5;     % dynes/cm
fixed.bathViscosity      = 0.18/1.175;  % Stokes
fixed.bathDiameter       = 100;        % disk radii
fixed.spatialResolution  = 50;         % intervals per disk radius
fixed.temporalResolution = 60;        % steps per adimensional unit
fixed.simulationTime     = 30/90;      % s
fixed.debug_flag         = false;
fixed.solverType         = "auto";
%% -----------------------------------------------------------------------

% Locate simulation_code on the path
codeFolder = fullfile(fileparts(mfilename('fullpath')), 'simulation_code');
addpath(codeFolder);

% Output directory for CSVs (timestamped subfolder)
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
outDir = fullfile(fileparts(mfilename('fullpath')), '..', '0_data', 'sweep_results', ['sweep_', timestamp]);
if ~exist(outDir, 'dir'); mkdir(outDir); end

% Build Cartesian product of sweep axes
[GG, FF] = ndgrid(sweep.gamma, sweep.bathFrequency);
nCases   = numel(GG);

% Logging
logDir = fullfile(fileparts(mfilename('fullpath')), '..', '0_data', 'logs');
if ~exist(logDir, 'dir'); mkdir(logDir); end
logFile = fullfile(logDir, sprintf('sweeper_%s.txt', datestr(datetime(), 30)));
diary(logFile);
fprintf('Sweep started: %s\n', datestr(datetime()));
fprintf('Axes: gamma=[%s]  frequency=[%s] Hz\n', ...
    num2str(sweep.gamma), num2str(sweep.bathFrequency));
fprintf('Fixed params:\n'); disp(fixed);
fprintf('Total cases: %d\n\n', nCases);

% Summary table accumulators
summaryRows = cell(nCases, 1);

for ii = 1:nCases
    gamma_i = GG(ii);
    freq_i  = FF(ii);
    omega_i = 2 * pi * freq_i;
    
    % F / (m*g) = gamma  =>  F = gamma * m * g
    F_i = gamma_i * fixed.diskMass * g_cgs; 
    A_forcing_i = (gamma_i * g_cgs) / omega_i^2;

    fprintf('[%d/%d]  gamma_forcing=%.3f  f=%g Hz  -> F=%.4f dyn, A_f=%.4e cm\n', ...
        ii, nCases, gamma_i, freq_i, F_i, A_forcing_i);

    % CSV file for this run
    csvFile = fullfile(outDir, ...
        sprintf('CoM_gamma%.3f_f%gHz.csv', gamma_i, freq_i));

    if exist(csvFile, 'file')
        fprintf('  Skipping (CSV exists)\n');
        summaryRows{ii} = {gamma_i, freq_i, 0, NaN, NaN, NaN, 'skipped'};
        continue
    end

    t0 = tic;
    try
        [t_s, CoM_cm, eta_history_cm] = solve_motion( ...
            'diskRadius',         fixed.diskRadius, ...
            'diskMass',           fixed.diskMass, ...
            'forceAmplitude',     F_i, ...
            'forceFrequency',     freq_i, ...
            'bathAmplitude',      0.0, ...
            'bathFrequency',      freq_i, ...
            'phaseDifference',    fixed.phaseDifference, ...
            'bathDensity',        fixed.bathDensity, ...
            'bathSurfaceTension', fixed.bathSurfaceTension, ...
            'bathViscosity',      fixed.bathViscosity, ...
            'bathDiameter',       fixed.bathDiameter, ...
            'spatialResolution',  fixed.spatialResolution, ...
            'temporalResolution', fixed.temporalResolution, ...
            'simulationTime',     30/freq_i, ...
            'debug_flag',         fixed.debug_flag, ...
            'solverType',         fixed.solverType);

        % Post-processing: Calculate max elevation in outer boundary region
        % (Outer 1 disk diameter = 2 disk radii)
        boundary_nodes = ceil(2 * fixed.spatialResolution);
        eta_boundary_cm = max(abs(eta_history_cm(end-boundary_nodes+1:end, :)), [], 1);

        elapsed = toc(t0);

        % Write CoM CSV: header + data
        fid = fopen(csvFile, 'w');
        fprintf(fid, 'time_s,CoM_cm,eta_boundary_cm\n');
        fprintf(fid, '%.6e,%.6e,%.6e\n', [t_s(:), CoM_cm(:), eta_boundary_cm(:)]');
        fclose(fid);

        CoM_max = max(abs(CoM_cm));
        CoM_rms = sqrt(mean(CoM_cm.^2));
        fprintf('  Done in %.1f s  |CoM|_max=%.3e cm  |CoM|_rms=%.3e cm\n', ...
            elapsed, CoM_max, CoM_rms);
        summaryRows{ii} = {gamma_i, freq_i, A_forcing_i, CoM_max, CoM_rms, elapsed, 'ok'};

    catch ME
        elapsed = toc(t0);
        fprintf('  ERROR after %.1f s: %s\n', elapsed, ME.message);
        summaryRows{ii} = {gamma_i, freq_i, A_forcing_i, NaN, NaN, elapsed, 'error'};
    end
end

%% Write summary CSV
summaryFile = fullfile(outDir, sprintf('summary_%s.csv', datestr(datetime(), 30)));
fid = fopen(summaryFile, 'w');
fprintf(fid, 'gamma,bathFrequency_Hz,bathAmplitude_cm,CoM_max_cm,CoM_rms_cm,elapsed_s,status\n');
for ii = 1:nCases
    r = summaryRows{ii};
    fprintf(fid, '%.4f,%g,%.6e,%.4e,%.4e,%.1f,%s\n', r{:});
end
fclose(fid);

fprintf('\nSweep finished: %s\n', datestr(datetime()));
fprintf('Summary: %s\n', summaryFile);
diary off

end
