import os
import numpy as np
import matplotlib.pyplot as plt
from src.model import *
from src.experiments import run_experiment

os.makedirs('figures', exist_ok=True)

# ±30% parameter ranges
taus   = [1.4, 2.0, 2.6]
gammas = [0.18, 0.25, 0.32]
rhos   = [0.20, 0.35, 0.50]

immune_set = [
    ('Linear',      lambda C: f_linear(C)),
    ('MM',          lambda C: f_mm(C, 0.5)),
    ('Hill',        lambda C: f_hill(C, 0.5, 4)),
    ('Switch',      lambda C: f_switch(C, 10.0, 0.5))
]

exposure = lambda t, C: exposure_constant(t, C)
data = {name: [] for name, _ in immune_set}

for tau in taus:
    for g in gammas:
        for r_val in rhos:
            for name, f_imm in immune_set:
                p = {'r': 2.8, 'gamma': g, 'rho': r_val}
                _, _, _, m = run_experiment(f_imm, exposure, p, tau=tau, T=100.0, dt=0.01)
                lam = m['lambda']
                if not np.isnan(lam):
                    data[name].append(lam)

fig, ax = plt.subplots(figsize=(7, 4.5))
bp = ax.boxplot([data[name] for name, _ in immune_set],
                positions=[1, 2, 3, 4], widths=0.6, patch_artist=True,
                medianprops=dict(color='black', linewidth=1.5))
colors = ['#0072BD', '#D95319', '#EDB120', '#7E2F8E']
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)

ax.set_xticklabels([name for name, _ in immune_set])
ax.set_ylabel('$\\lambda_{\\mathrm{decay}}$ (day$^{-1}$)', fontsize=12)
ax.set_title('Parameter sensitivity: stability of decay hierarchy', fontsize=12)
ax.grid(axis='y', alpha=0.3)
plt.tight_layout()
plt.savefig('figures/Sensitivity_Lambda_Distribution.eps', format='eps', dpi=300)
plt.savefig('figures/Sensitivity_Lambda_Distribution.png', format='png', dpi=300)
plt.close()
print("[run_sensitivity] Figure 13 generated.")
