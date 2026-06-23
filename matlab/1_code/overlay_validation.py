import os
import sys
import glob
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def run_overlay_validation(sweep_dir=None):
    # Determine base and sweep directories
    script_dir = os.path.dirname(os.path.abspath(__file__))
    base_dir = os.path.join(script_dir, '..', '0_data', 'sweep_results')
    
    if sweep_dir is None:
        # Find the most recent 'sweep_*' folder
        sweep_folders = glob.glob(os.path.join(base_dir, 'sweep_*'))
        if not sweep_folders:
            sweep_dir = base_dir
        else:
            sweep_dir = max(sweep_folders, key=os.path.getmtime)
            
    print(f"Using sweep results from: {sweep_dir}")
    
    # Modes & Physics constants
    g_cgs = 981.0
    phase_rad_bath = -90.0 * np.pi / 180.0
    
    # Detect forcing-driven vs bath-driven
    is_bath_driven = True
    summary_files = glob.glob(os.path.join(sweep_dir, 'summary_*.csv'))
    if summary_files:
        summary_data = pd.read_csv(summary_files[0])
        if 'bathAmplitude_cm' in summary_data.columns:
            if (summary_data['bathAmplitude_cm'] == 0).all():
                is_bath_driven = False
                print("  Detected FORCING-driven sweep.")
            else:
                print("  Detected BATH-driven sweep.")
                
    # Find CoM simulation result CSVs
    sim_files = glob.glob(os.path.join(sweep_dir, 'CoM_gamma*_f*Hz.csv'))
    if not sim_files:
        print(f"No simulation CSV files found in {sweep_dir}. Aborting validation overlay.")
        return
        
    results = []
    
    for fpath in sim_files:
        fname = os.path.basename(fpath)
        match = re.match(r'CoM_gamma([\d.]+)_f([\d.]+)Hz\.csv', fname)
        if not match:
            continue
            
        g_val = float(match.group(1))
        f_val = float(match.group(2))
        
        try:
            data = pd.read_csv(fpath)
        except Exception as e:
            print(f"Error reading {fname}: {e}")
            continue
            
        if len(data) < 10:
            continue
            
        t = data['time_s'].values
        CoM = data['CoM_cm'].values
        eta_boundary = data['eta_boundary_cm'].values if 'eta_boundary_cm' in data.columns else np.zeros_like(t)
        
        omega = 2.0 * np.pi * f_val
        A_ref = g_val * g_cgs / (omega**2)
        T_period = 1.0 / f_val
        
        # Truncate at reflection (10% threshold of A_ref)
        wave_idx = np.where(np.abs(eta_boundary) > 0.10 * A_ref)[0]
        if wave_idx.size > 0:
            idx_cutoff = max(1, wave_idx[0])
            t = t[:idx_cutoff]
            CoM = CoM[:idx_cutoff]
            
        total_time = t[-1] - t[0]
        n_periods_available = int(np.floor(total_time / T_period))
        
        if n_periods_available >= 2:
            n_eval = min(3, n_periods_available)
            t_start_eval = t[-1] - n_eval * T_period
            eval_idx = t >= t_start_eval
            
            t_eval = t[eval_idx]
            CoM_eval = CoM[eval_idx]
            
            if is_bath_driven:
                z_lab = CoM_eval + A_ref * np.cos(omega * t_eval + phase_rad_bath)
                phase_ref_rad = phase_rad_bath
            else:
                z_lab = CoM_eval
                phase_ref_rad = 0.0
                
            # Least squares fitting
            M = np.column_stack([np.cos(omega * t_eval), np.sin(omega * t_eval), np.ones_like(t_eval)])
            # Solve M * X = z_lab
            try:
                X, _, _, _ = np.linalg.lstsq(M, z_lab, rcond=None)
                ampSS = np.sqrt(X[0]**2 + X[1]**2)
                phi_disk = np.arctan2(-X[1], X[0])
                dphi = np.degrees(phase_ref_rad - phi_disk) % 360.0
                
                results.append({
                    'gamma': g_val,
                    'freq': f_val,
                    'ampNorm': ampSS / A_ref,
                    'phaseDiff': dphi
                })
            except Exception as e:
                print(f"Fitting failed for {fname}: {e}")
                
    if not results:
        print("No valid steady state periods found in simulation files.")
        return
        
    df_sim = pd.DataFrame(results)
    
    # Unique gammas in simulation results
    u_gammas = sorted(df_sim['gamma'].unique())
    # Color mapping for simulation curves (matplotlib equivalent of Matlab lines)
    cmap = plt.colormaps['tab10']
    sim_colors = {g: cmap(i) for i, g in enumerate(u_gammas)}
    
    # Load experimental measurement CSVs
    figures_dir = os.path.join(script_dir, 'Figures')
    exp_data_dir = os.path.abspath(os.path.join(script_dir, '..', '0_data', 'external', 'experimental_measurements'))
    amp_csv = os.path.join(exp_data_dir, 'ampiezza_solo_misure_3.csv')
    phase_csv = os.path.join(exp_data_dir, 'fase_solo_misure_3.csv')
    
    has_exp = os.path.exists(amp_csv) and os.path.exists(phase_csv)
    if has_exp:
        df_exp_amp = pd.read_csv(amp_csv)
        df_exp_phase = pd.read_csv(phase_csv)
        exp_colors = {
            0.05: (0.6, 1.0, 0.6),
            0.20: (0.35, 0.75, 0.35),
            0.50: (0.1, 0.5, 0.1)
        }
    else:
        print("Warning: Experimental measurement CSV files not found. Simulation curves will be plotted without experimental overlay.")
        
    sweep_name = os.path.basename(sweep_dir.rstrip('/'))
    
    # --- Overlay Plot 1: Amplitude ---
    plt.figure(figsize=(8, 6))
    ax = plt.gca()
    
    # Plot experimental data in the background
    if has_exp:
        for g in sorted(df_exp_amp['gamma'].unique()):
            df_g = df_exp_amp[df_exp_amp['gamma'] == g]
            plt.errorbar(
                df_g['frequency_Hz'], 
                df_g['value'], 
                yerr=[df_g['error_lower'], df_g['error_upper']], 
                fmt='o', 
                color='k', 
                mfc=exp_colors[g], 
                mec='k', 
                ms=7, 
                capsize=3, 
                elinewidth=1.0, 
                label=f'Exp $\\Gamma = {g:.2f}$',
                zorder=2
            )
            
    # Plot simulation overlay
    for g in u_gammas:
        df_g = df_sim[df_sim['gamma'] == g].sort_values(by='freq')
        plt.plot(
            df_g['freq'], 
            df_g['ampNorm'], 
            'o-', 
            color=sim_colors[g], 
            lw=2, 
            ms=6, 
            label=f'Sim $\\Gamma = {g:.2f}$',
            zorder=3
        )
        
    plt.xlabel(r'$f$ (Hz)', fontsize=12)
    plt.ylabel(r'$A_{\mathrm{disk}} / A_{\mathrm{base}}$', fontsize=12)
    plt.title(f'Amplitude Validation Overlay ({sweep_name})', fontsize=14, fontweight='bold', pad=15)
    plt.xlim(0, 105)
    plt.ylim(0, 1.2)
    plt.xticks([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
    plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2])
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend(loc='best', frameon=True, facecolor='white', framealpha=0.9, edgecolor='gray')
    ax.set_axisbelow(True)
    
    # Save Amplitude overlay
    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, f'val_amp_{sweep_name}.png'), dpi=300)
    plt.savefig(os.path.join(figures_dir, f'val_amp_{sweep_name}.pdf'))
    plt.close()
    
    # --- Overlay Plot 2: Phase ---
    plt.figure(figsize=(8, 6))
    ax = plt.gca()
    
    # Plot experimental data in the background
    if has_exp:
        for g in sorted(df_exp_phase['gamma'].unique()):
            df_g = df_exp_phase[df_exp_phase['gamma'] == g]
            plt.errorbar(
                df_g['frequency_Hz'], 
                df_g['value'], 
                yerr=[df_g['error_lower'], df_g['error_upper']], 
                fmt='o', 
                color='k', 
                mfc=exp_colors[g], 
                mec='k', 
                ms=7, 
                capsize=3, 
                elinewidth=1.0, 
                label=f'Exp $\\Gamma = {g:.2f}$',
                zorder=2
            )
            
    # Plot simulation overlay
    for g in u_gammas:
        df_g = df_sim[df_sim['gamma'] == g].sort_values(by='freq')
        plt.plot(
            df_g['freq'], 
            df_g['phaseDiff'], 
            'o-', 
            color=sim_colors[g], 
            lw=2, 
            ms=6, 
            label=f'Sim $\\Gamma = {g:.2f}$',
            zorder=3
        )
        
    plt.xlabel(r'$f$ (Hz)', fontsize=12)
    plt.ylabel(r'Phase Difference (deg)', fontsize=12)
    plt.title(f'Phase Validation Overlay ({sweep_name})', fontsize=14, fontweight='bold', pad=15)
    plt.xlim(0, 105)
    plt.ylim(0, 360)  # Wrapped to [0, 360]
    plt.xticks([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
    plt.yticks([0, 90, 180, 270, 360])
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend(loc='best', frameon=True, facecolor='white', framealpha=0.9, edgecolor='gray')
    ax.set_axisbelow(True)
    
    # Save Phase overlay
    plt.tight_layout()
    plt.savefig(os.path.join(figures_dir, f'val_phase_{sweep_name}.png'), dpi=300)
    plt.savefig(os.path.join(figures_dir, f'val_phase_{sweep_name}.pdf'))
    plt.close()
    
    print(f"Validation complete. Figures saved in Figures/ as val_amp_{sweep_name}.png/pdf and val_phase_{sweep_name}.png/pdf")

if __name__ == '__main__':
    sweep_dir_arg = sys.argv[1] if len(sys.argv) > 1 else None
    run_overlay_validation(sweep_dir_arg)
