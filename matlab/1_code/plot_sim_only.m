function plot_sim_only(sweepDir)
    % PLOT_SIM_ONLY Plots simulation-only results (Amplitude and Phase) from a sweep.
    
    if nargin < 1 || isempty(sweepDir)
        baseDir = fullfile(fileparts(mfilename('fullpath')), '..', '0_data', 'sweep_results');
        d = dir(fullfile(baseDir, 'sweep_*'));
        if isempty(d); error('No sweep folders found.'); end
        [~, idx] = max([d.datenum]);
        sweepDir = fullfile(baseDir, d(idx).name);
    end
    [~, sweepName] = fileparts(sweepDir);
    fprintf('Plotting simulation results from: %s\n', sweepDir);

    g_cgs = 981;
    files = dir(fullfile(sweepDir, 'CoM_gamma*_f*Hz.csv'));
    res = struct('gamma', [], 'freq', [], 'ampNorm', [], 'phaseDiff', []);
    count = 0;

    % --- Mode Detection ---
    summaryFile = dir(fullfile(sweepDir, 'summary_*.csv'));
    isBathDriven = true; % Default
    if ~isempty(summaryFile)
        summaryData = readtable(fullfile(sweepDir, summaryFile(1).name));
        if all(summaryData.bathAmplitude_cm == 0)
            isBathDriven = false;
            fprintf('  Detected FORCING-driven sweep.\n');
        else
            fprintf('  Detected BATH-driven sweep.\n');
        end
    end
    phase_rad_bath = -90 * pi / 180;

    for k = 1:numel(files)
        % ... (loading logic unchanged) ...
        tok = regexp(files(k).name, 'CoM_gamma([\d.]+)_f([\d.]+)Hz\.csv', 'tokens');
        if isempty(tok); continue; end
        g_val = str2double(tok{1}{1});
        f_val = str2double(tok{1}{2});

        opts = detectImportOptions(fullfile(sweepDir, files(k).name));
        data = readtable(fullfile(sweepDir, files(k).name), opts);
        t = data.time_s; CoM = data.CoM_cm; eta_boundary = data.eta_boundary_cm;

        omega = 2 * pi * f_val;
        A_ref = g_val * g_cgs / omega^2;
        T_period = 1 / f_val;
        
        % Truncate at reflection (10% threshold)
        wave_idx = find(abs(eta_boundary) > 0.10 * A_ref, 1);
        if ~isempty(wave_idx); t = t(1:max(1, wave_idx-1)); CoM = CoM(1:max(1, wave_idx-1)); end

        n_periods_available = floor((t(end) - t(1)) / T_period);
        if n_periods_available >= 2
            count = count + 1;
            n_eval = min(3, n_periods_available); 
            t_eval = t(t >= (t(end) - n_eval * T_period));
            CoM_eval = CoM(t >= (t(end) - n_eval * T_period));
            
            % --- Physics Mode Adjustment ---
            if isBathDriven
                z_lab = CoM_eval + A_ref * cos(omega * t_eval + phase_rad_bath);
                phase_ref_rad = phase_rad_bath;
            else
                z_lab = CoM_eval;
                phase_ref_rad = 0; % Forcing is cos(omega*t)
            end
            
            % Robust LS Fit
            M = [cos(omega * t_eval), sin(omega * t_eval), ones(size(t_eval))];
            X = M \ z_lab; 
            ampSS = sqrt(X(1)^2 + X(2)^2);
            phi_disk = atan2(-X(2), X(1));
            dphi = mod(rad2deg(phase_ref_rad - phi_disk), 360); % Phase lag
            
            res(count).gamma = g_val; res(count).freq = f_val;
            res(count).ampNorm = ampSS / A_ref; res(count).phaseDiff = dphi;
        end
    end

    gammas = [res.gamma]; uG = unique(gammas); colors = lines(length(uG));

    % 1. Amplitude Plot
    figure('Name', ['Amplitude: ', sweepName], 'Color', 'w');
    hold on; grid on; box on;
    for i = 1:length(uG)
        idx = (gammas == uG(i));
        [f_s, si] = sort([res(idx).freq]); a_s = [res(idx).ampNorm]; a_s = a_s(si);
        plot(f_s, a_s, 'o-', 'Color', colors(i,:), 'LineWidth', 2, 'MarkerSize', 8, ...
            'MarkerFaceColor', colors(i,:), 'DisplayName', sprintf('\\Gamma = %.2f', uG(i)));
    end
    xlabel('Frequency (Hz)'); ylabel('A_{disk} / A_{forcing}');
    title(['Simulation Response: ', sweepName]); legend('Location', 'best');
    saveas(gcf, fullfile(fileparts(mfilename('fullpath')), 'Figures', sprintf('sim_amp_%s.png', sweepName)));

    % 2. Phase Plot
    figure('Name', ['Phase: ', sweepName], 'Color', 'w');
    hold on; grid on; box on;
    for i = 1:length(uG)
        idx = (gammas == uG(i));
        [f_s, si] = sort([res(idx).freq]); p_s = [res(idx).phaseDiff]; p_s = p_s(si);
        plot(f_s, p_s, 'o-', 'Color', colors(i,:), 'LineWidth', 2, 'MarkerSize', 8, ...
            'MarkerFaceColor', colors(i,:), 'DisplayName', sprintf('\\Gamma = %.2f', uG(i)));
    end
    xlabel('Frequency (Hz)'); ylabel('Phase Lag (deg)');
    title(['Simulation Phase: ', sweepName]); legend('Location', 'best');
    saveas(gcf, fullfile(fileparts(mfilename('fullpath')), 'Figures', sprintf('sim_phase_%s.png', sweepName)));
    
    fprintf('Simulation plots saved in 1_code/Figures/\n');
end
