function check_early()
    sweepDir = 'sweep_results';
    g_cgs = 981;
    phase_rad = -90 * pi / 180;
    
    files = dir(fullfile(sweepDir, 'CoM_gamma*_f*Hz.csv'));
    for k = 1:numel(files)
        tok = regexp(files(k).name, 'CoM_gamma([\d.]+)_f([\d.]+)Hz\.csv', 'tokens');
        gamma = str2double(tok{1}{1});
        freq = str2double(tok{1}{2});
        
        opts = detectImportOptions(fullfile(sweepDir, files(k).name));
        data = readtable(fullfile(sweepDir, files(k).name), opts);
        
        t = data.time_s;
        CoM = data.CoM_cm;
        eta_boundary = data.eta_boundary_cm;
        
        omega_k = 2 * pi * freq;
        A_bath = gamma * g_cgs / omega_k^2;
        T_period = 1 / freq;
        
        wave_threshold = 0.10 * A_bath;
        wave_idx = find(abs(eta_boundary) > wave_threshold, 1);
        
        if ~isempty(wave_idx)
            periods_reached = t(wave_idx) / T_period;
            fprintf('f=%g Hz, gamma=%.3f: Wave reflection at t=%.3f s (%.1f periods)\n', freq, gamma, t(wave_idx), periods_reached);
            t = t(1:max(1, wave_idx-1));
            CoM = CoM(1:max(1, wave_idx-1));
        else
            periods_reached = t(end) / T_period;
            fprintf('f=%g Hz, gamma=%.3f: No reflection detected within %.1f periods\n', freq, gamma, periods_reached);
        end
        
        n_periods_available = floor((t(end) - t(1)) / T_period);
        n_periods_eval = min(3, n_periods_available);
        
        if n_periods_eval >= 2
            t_start_eval = t(end) - n_periods_eval * T_period;
            eval_idx = t >= t_start_eval;
            t_eval = t(eval_idx);
            z_lab = CoM(eval_idx) + A_bath * cos(omega_k * t_eval + phase_rad);
            
            ampSS = (max(z_lab) - min(z_lab)) / 2;
            ampNorm = ampSS / A_bath;
            fprintf('   -> A_disk/A_bath = %.4f\n', ampNorm);
        else
            fprintf('   -> Not enough periods for steady state.\n');
        end
    end
end
