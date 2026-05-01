import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

os.makedirs('figures', exist_ok=True)

# Clinical data (Table A.2)
t_covid = np.array([10, 11, 12, 13])
C_covid = np.array([1.585, 0.100, 0.016, 0.001])
slope_c, _, r_c, _, _ = linregress(t_covid, np.log(C_covid))
lam_c = -slope_c
t12_c = np.log(2) / lam_c * 24
R2_c = r_c**2

t_ebola = np.array([9, 11, 13, 15])
C_ebola = np.array([0.165, 0.040, 0.0034, 0.00017])
slope_e, _, r_e, _, _ = linregress(t_ebola, np.log(C_ebola))
lam_e = -slope_e
t12_e = np.log(2) / lam_e * 24
R2_e = r_e**2

fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))

# Panel A: decay rates
ax = axes[0]
bars = ax.bar(['SARS-CoV-2', 'Ebola'], [lam_c, lam_e],
              color=['#e74c3c', '#3498db'], edgecolor='black', linewidth=1.2)
ax.set_ylabel('Decay rate $\\lambda$ (day$^{-1}$)', fontsize=12)
ax.set_title('Post-peak decay rates', fontsize=12)
ax.grid(axis='y', alpha=0.3)
for bar in bars:
    h = bar.get_height()
    ax.text(bar.get_x() + bar.get_width() / 2.0, h, f'{h:.2f}',
            ha='center', va='bottom', fontsize=11)

# Panel B: half-lives
ax = axes[1]
bars = ax.bar(['SARS-CoV-2', 'Ebola'], [t12_c, t12_e],
              color=['#c0392b', '#2980b9'], edgecolor='black', linewidth=1.2)
ax.set_ylabel('Half-life (hours)', fontsize=12)
ax.set_title('Viral clearance half-life', fontsize=12)
ax.grid(axis='y', alpha=0.3)
for bar in bars:
    h = bar.get_height()
    ax.text(bar.get_x() + bar.get_width() / 2.0, h, f'{h:.1f}',
            ha='center', va='bottom', fontsize=11)

plt.tight_layout()
plt.savefig('figures/Comparative_Decay_Rates.eps', format='eps', dpi=300)
plt.savefig('figures/Comparative_Decay_Rates.png', format='png', dpi=300)
plt.close()
print("[run_clinical] Figure 14 generated.")
