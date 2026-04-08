function overlay_validation(sweepDir, ampFig, phaseFig)
    % OVERLAY_VALIDATION Overlays both amplitude and phase results on experimental figures.
    
    if nargin < 1 || isempty(sweepDir)
        baseDir = fullfile(fileparts(mfilename('fullpath')), '..', '0_data', 'sweep_results');
        % Automatically find the most recent 'sweep_*' subfolder
        d = dir(fullfile(baseDir, 'sweep_*'));
        if isempty(d)
            % Fallback to the base directory if no subfolders exist
            sweepDir = baseDir;
        else
            [~, idx] = max([d.datenum]);
            sweepDir = fullfile(baseDir, d(idx).name);
        end
    end
    fprintf('Using sweep results from: %s\n', sweepDir);

    if nargin < 2 || isempty(ampFig)
        ampFig = fullfile(fileparts(mfilename('fullpath')), 'Figures', 'ampiezza_solo_misure_3.fig');
    end
    if nargin < 3 || isempty(phaseFig)
        phaseFig = fullfile(fileparts(mfilename('fullpath')), 'Figures', 'fase_solo_misure_3.fig');
    end

    g_cgs = 981;
    phase_rad_bath = -90 * pi / 180;

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

    files = dir(fullfile(sweepDir, 'CoM_gamma*_f*Hz.csv'));
    
    res = struct('gamma', [], 'freq', [], 'ampNorm', [], 'phaseDiff', []);
    count = 0;

    for k = 1:numel(files)
        tok = regexp(files(k).name, 'CoM_gamma([\d.]+)_f([\d.]+)Hz\.csv', 'tokens');
        if isempty(tok); continue; end
        g_val = str2double(tok{1}{1});
        f_val = str2double(tok{1}{2});

        opts = detectImportOptions(fullfile(sweepDir, files(k).name));
        data = readtable(fullfile(sweepDir, files(k).name), opts);
        
        t = data.time_s;
        CoM = data.CoM_cm;
        eta_boundary = data.eta_boundary_cm;

        omega = 2 * pi * f_val;
        A_ref = g_val * g_cgs / omega^2; 
        T_period = 1 / f_val;
        
        % Truncate at reflection (10% threshold)
        wave_idx = find(abs(eta_boundary) > 0.10 * A_ref, 1);
        if ~isempty(wave_idx)
            t = t(1:max(1, wave_idx-1));
            CoM = CoM(1:max(1, wave_idx-1));
        end

        n_periods_available = floor((t(end) - t(1)) / T_period);
        if n_periods_available >= 2
            count = count + 1;
            
            % --- Use only the last 3 periods for the fit ---
            n_eval = min(3, n_periods_available); 
            eval_idx = t >= (t(end) - n_eval * T_period);
            t_eval = t(eval_idx);
            CoM_eval = CoM(eval_idx);
            
            % --- Physics Mode Adjustment ---
            if isBathDriven
                z_lab = CoM_eval + A_ref * cos(omega * t_eval + phase_rad_bath);
                phase_ref_rad = phase_rad_bath;
            else
                z_lab = CoM_eval;
                phase_ref_rad = 0; % Forcing is cos(omega*t)
            end
            
            % --- Robust Least Squares Fitting ---
            M = [cos(omega * t_eval), sin(omega * t_eval), ones(size(t_eval))];
            X = M \ z_lab; 
            
            % Amplitude A = sqrt(C1^2 + C2^2)
            ampSS = sqrt(X(1)^2 + X(2)^2);
            
            % Phase of the disk: z_lab = A*cos(omega*t + phi_disk)
            phi_disk = atan2(-X(2), X(1));
            
            % Phase difference (Convention: phi_ref - phi_disk for positive lag)
            dphi = mod(rad2deg(phase_ref_rad - phi_disk), 360); % Wrap to [0, 360]
            
            res(count).gamma = g_val;
            res(count).freq = f_val;
            res(count).ampNorm = ampSS / A_ref;
            res(count).phaseDiff = dphi;
        end
    end

    gammas = [res.gamma];
    uG = unique(gammas);
    colors = lines(length(uG));
    [~, sweepName] = fileparts(sweepDir);

    %% 2. Overlay Amplitude
    if exist(ampFig, 'file')
        fprintf('Updating Amplitude Figure: %s\n', ampFig);
        figAmp = openfig(ampFig, 'new', 'visible');
        ax = gca; hold on;
        h_sim = gobjects(length(uG), 1);
        for i = 1:length(uG)
            idx = (gammas == uG(i));
            f_s = [res(idx).freq];
            a_s = [res(idx).ampNorm];
            [f_s, si] = sort(f_s); a_s = a_s(si);
            h_sim(i) = plot(ax, f_s, a_s, 'o-', 'Color', colors(i,:), 'LineWidth', 2, ...
                'MarkerSize', 8, 'MarkerFaceColor', colors(i,:), 'DisplayName', sprintf('\\Gamma = %.2f', uG(i)));
        end
        legend(ax, h_sim, 'Location', 'best', 'FontSize', 10);
        set(ax, 'YTick', 0:0.2:2.0);
        outPath = fullfile(fileparts(mfilename('fullpath')), 'Figures', sprintf('val_amp_%s', sweepName));
        saveas(figAmp, [outPath, '.fig']);
        saveas(figAmp, [outPath, '.png']);
    end

    %% 3. Overlay Phase
    if exist(phaseFig, 'file')
        fprintf('Updating Phase Figure: %s\n', phaseFig);
        figPhase = openfig(phaseFig, 'new', 'visible');
        ax = gca; hold on;
        h_sim = gobjects(length(uG), 1);
        for i = 1:length(uG)
            idx = (gammas == uG(i));
            f_s = [res(idx).freq];
            p_s = [res(idx).phaseDiff];
            [f_s, si] = sort(f_s); p_s = p_s(si);
            h_sim(i) = plot(ax, f_s, p_s, 'o-', 'Color', colors(i,:), 'LineWidth', 2, ...
                'MarkerSize', 8, 'MarkerFaceColor', colors(i,:), 'DisplayName', sprintf('\\Gamma = %.2f', uG(i)));
        end
        legend(ax, h_sim, 'Location', 'best', 'FontSize', 10);
        ylabel(ax, 'Phase Difference (deg)');
        outPath = fullfile(fileparts(mfilename('fullpath')), 'Figures', sprintf('val_phase_%s', sweepName));
        saveas(figPhase, [outPath, '.fig']);
        saveas(figPhase, [outPath, '.png']);
    end
    
    fprintf('Validation complete. Figures saved in 1_code/Figures/\n');
end
