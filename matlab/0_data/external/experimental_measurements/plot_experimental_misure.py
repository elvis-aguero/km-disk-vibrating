import os
import pandas as pd
import matplotlib.pyplot as plt

# Get the directory where this script is located
script_dir = os.path.dirname(os.path.abspath(__file__))

# Load the data
amp_csv = os.path.join(script_dir, 'ampiezza_solo_misure_3.csv')
phase_csv = os.path.join(script_dir, 'fase_solo_misure_3.csv')

if not os.path.exists(amp_csv) or not os.path.exists(phase_csv):
    raise FileNotFoundError("CSV files not found. Please run the data extraction first.")

df_amp = pd.read_csv(amp_csv)
df_phase = pd.read_csv(phase_csv)

# Define the colors matching the Matlab figure (RGB scaled to [0,1])
# Gamma 1 (0.05) = Light Green [0.6, 1.0, 0.6]
# Gamma 2 (0.20) = Medium Green [0.35, 0.75, 0.35]
# Gamma 3 (0.50) = Dark Green [0.1, 0.5, 0.1]
gamma_colors = {
    0.05: (0.6, 1.0, 0.6),
    0.20: (0.35, 0.75, 0.35),
    0.50: (0.1, 0.5, 0.1)
}

gamma_labels = {
    0.05: r'$\Gamma = 0.05$',
    0.20: r'$\Gamma = 0.20$',
    0.50: r'$\Gamma = 0.50$'
}

# --- Plot 1: Amplitude ---
plt.figure(figsize=(7, 5))
ax = plt.gca()

# Plot each gamma group
for g in sorted(df_amp['gamma'].unique()):
    df_g = df_amp[df_amp['gamma'] == g]
    plt.errorbar(
        df_g['frequency_Hz'], 
        df_g['value'], 
        yerr=[df_g['error_lower'], df_g['error_upper']], 
        fmt='o', 
        color='k', 
        mfc=gamma_colors[g], 
        mec='k', 
        ms=8, 
        capsize=3, 
        elinewidth=1, 
        label=gamma_labels[g]
    )

plt.xlabel(r'$f$ (Hz)', fontsize=12)
plt.ylabel(r'$A_{\mathrm{disk}} / A_{\mathrm{base}}$', fontsize=12)
plt.title('Experimental Amplitude Response', fontsize=14, fontweight='bold', pad=15)
plt.xlim(0, 105)
plt.ylim(0, 1.2)
plt.xticks([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
plt.yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2])
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend(loc='lower left', frameon=True, facecolor='white', framealpha=0.9, edgecolor='gray')
ax.set_axisbelow(True)

# Save figure
plt.tight_layout()
plt.savefig(os.path.join(script_dir, 'ampiezza_solo_misure_3.png'), dpi=300)
plt.savefig(os.path.join(script_dir, 'ampiezza_solo_misure_3.pdf'))
plt.close()

# --- Plot 2: Phase ---
plt.figure(figsize=(7, 5))
ax = plt.gca()

# Plot each gamma group
for g in sorted(df_phase['gamma'].unique()):
    df_g = df_phase[df_phase['gamma'] == g]
    plt.errorbar(
        df_g['frequency_Hz'], 
        df_g['value'], 
        yerr=[df_g['error_lower'], df_g['error_upper']], 
        fmt='o', 
        color='k', 
        mfc=gamma_colors[g], 
        mec='k', 
        ms=8, 
        capsize=3, 
        elinewidth=1, 
        label=gamma_labels[g]
    )

plt.xlabel(r'$f$ (Hz)', fontsize=12)
plt.ylabel(r'$\phi$ (deg)', fontsize=12)
plt.title('Experimental Phase Response', fontsize=14, fontweight='bold', pad=15)
plt.xlim(0, 105)
plt.ylim(0, 30)
plt.xticks([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
plt.yticks([0, 10, 20, 30])
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend(loc='lower right', frameon=True, facecolor='white', framealpha=0.9, edgecolor='gray')
ax.set_axisbelow(True)

# Save figure
plt.tight_layout()
plt.savefig(os.path.join(script_dir, 'fase_solo_misure_3.png'), dpi=300)
plt.savefig(os.path.join(script_dir, 'fase_solo_misure_3.pdf'))
plt.close()

print("Plotting complete. PNG and PDF figures saved in the Figures/ directory.")
