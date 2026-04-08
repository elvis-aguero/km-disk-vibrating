function plot_sweep_response(sweepDir)
% PLOT_SWEEP_RESPONSE  Normalised steady-state CoM amplitude vs bath frequency.
%
%   Reads all CoM_gamma*_f*Hz.csv files from a sweep_results directory.
%   Detects wave reflection by checking when eta_boundary_cm exceeds 10% 
%   of the bath amplitude. Truncates the data at that point.
%   Calculates the lab-frame CoM amplitude over the last few integer
%   periods before truncation.
%   Normalises by the bath amplitude A_bath = Gamma*g/omega^2, and overlays 
%   one curve per Gamma value.
%
%   Arguments:
%     sweepDir   - path to sweep_results folder
%                  (default: sweep_results/ next to this file)

if nargin < 1 || isempty(sweepDir)
    sweepDir = fullfile(fileparts(mfilename('fullpath')), 'sweep_results');
end

g_cgs = 981; % cm/s^2
phaseDifference = -90; % degrees
phase_rad = phaseDifference * pi / 180;

%% Discover files and parse (gamma, frequency) from filenames
files = dir(fullfile(sweepDir, 'CoM_gamma*_f*Hz.csv'));
if isempty(files)
    error('plot_sweep_response:noFiles', ...
          'No CoM CSV files found in:\n  %s', sweepDir);
end

nFiles = numel(files);
gammas = nan(nFiles, 1);
freqs  = nan(nFiles, 1);

for k = 1:nFiles
    tok = regexp(files(k).name, ...
                 'CoM_gamma([\d.]+)_f([\d.]+)Hz\.csv', 'tokens');
    if isempty(tok)
        warning('plot_sweep_response:badName', ...
                'Could not parse filename: %s — skipping.', files(k).name);
        continue
    end
    gammas(k) = str2double(tok{1}{1});
    freqs(k)  = str2double(tok{1}{2});
end

valid  = ~isnan(gammas);
files  = files(valid);
gammas = gammas(valid);
freqs  = freqs(valid);

%% Compute steady-state amplitude for each run
ampNorm = nan(numel(files), 1);

for k = 1:numel(files)
    % Try to read table to handle headers easily
    opts = detectImportOptions(fullfile(sweepDir, files(k).name));
    data = readtable(fullfile(sweepDir, files(k).name), opts);
    
    if height(data) < 10
        continue; % too short
    end
    
    t = data.time_s;
    CoM = data.CoM_cm;
    
    % Check if eta_boundary_cm is available
    if ismember('eta_boundary_cm', data.Properties.VariableNames)
        eta_boundary = data.eta_boundary_cm;
    else
        % Fallback for old files
        eta_boundary = zeros(size(t));
    end

    % Physical parameters
    omega_k = 2 * pi * freqs(k);
    A_bath = gammas(k) * g_cgs / omega_k^2;
    T_period = 1 / freqs(k);
    
    % Truncate data if wave reflection is detected
    wave_threshold = 0.10 * A_bath;
    wave_idx = find(abs(eta_boundary) > wave_threshold, 1);
    
    if ~isempty(wave_idx)
        if wave_idx < length(t)
            fprintf('  [f=%g Hz, gamma=%.3f] Wave reflection detected at t=%.3f s. Truncating data.\n', ...
                freqs(k), gammas(k), t(wave_idx));
        end
        % Truncate to just before the wave hits
        t = t(1:max(1, wave_idx-1));
        CoM = CoM(1:max(1, wave_idx-1));
    end
    
    % If there's barely any data left, skip
    total_time_simulated = t(end) - t(1);
    n_periods_available = floor(total_time_simulated / T_period);
    
    if n_periods_available < 2
        warning('plot_sweep_response:notEnoughData', ...
            'gamma=%.3f, f=%g Hz: Less than 2 periods of clean data available before reflection.', ...
            gammas(k), freqs(k));
        continue;
    end
    
    % Number of periods to extract for steady-state (e.g., last 3, or whatever is available if < 3)
    n_periods_eval = min(3, n_periods_available);
    
    % The time window for evaluating the amplitude
    t_start_eval = t(end) - n_periods_eval * T_period;
    
    eval_idx = t >= t_start_eval;
    t_eval = t(eval_idx);
    CoM_eval = CoM(eval_idx);
    
    % Convert to Lab Frame
    % The solver simulates the CoM relative to the bath.
    % z_bath(t) = A_bath * cos(omega * t + phase)
    z_bath = A_bath * cos(omega_k * t_eval + phase_rad);
    z_lab = CoM_eval + z_bath;
    
    % Compute Peak-to-Peak Amplitude / 2
    % We'll also check for drift by comparing the first and last half of the evaluation window
    half_idx = floor(length(z_lab)/2);
    ampFirst = (max(z_lab(1:half_idx)) - min(z_lab(1:half_idx))) / 2;
    ampLast  = (max(z_lab(half_idx+1:end)) - min(z_lab(half_idx+1:end))) / 2;
    
    ampSS = ampLast;
    
    % Warn if amplitude is still changing significantly (drift check)
    relDrift = abs(ampFirst - ampLast) / max(ampLast, eps);
    if relDrift > 0.05
        warning('plot_sweep_response:notSteady', ...
            'gamma=%.3g, f=%g Hz: amplitude drifted by %.1f%% in the final clean periods. Steady state might not be fully reached.', ...
            gammas(k), freqs(k), relDrift * 100);
    end
    
    % Normalise by bath amplitude
    ampNorm(k) = ampSS / A_bath;
end

%% Plot — one series per unique Gamma
uniqueGammas = unique(gammas);
nG = numel(uniqueGammas);
colors = lines(nG);

fig = figure;
ax  = axes(fig);
hold(ax, 'on');

for gi = 1:nG
    gval = uniqueGammas(gi);
    idx  = gammas == gval & ~isnan(ampNorm);

    f_g   = freqs(idx);
    amp_g = ampNorm(idx);

    % Sort by frequency so the line connects left-to-right
    [f_sorted, si] = sort(f_g);
    amp_sorted = amp_g(si);

    plot(ax, f_sorted, amp_sorted, 'o-', ...
         'Color',           colors(gi, :), ...
         'LineWidth',       1.5, ...
         'MarkerSize',      6, ...
         'MarkerFaceColor', colors(gi, :), ...
         'DisplayName',     sprintf('\\Gamma = %.3g', gval));
end

xlabel(ax, 'Bath frequency (Hz)');
ylabel(ax, 'A_{disk}^{lab} / A_{bath}  (dimensionless)');
title(ax, 'Steady-state Lab CoM amplitude vs Frequency');
legend(ax, 'Location', 'best');
grid(ax, 'on');
box(ax, 'on');
hold(ax, 'off');

% Set y-axis to start from 0 if appropriate, though let MATLAB decide limits mostly
% ylim(ax, [0, max(ampNorm)*1.1]);

end
