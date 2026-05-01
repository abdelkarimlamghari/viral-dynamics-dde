import os
import numpy as np
import matplotlib.pyplot as plt
from src.model import *
from src.experiments import run_experiment

os.makedirs('figures', exist_ok=True)

# -------------------------------------------------
# Baseline parameters (Table 3 of manuscript)
# -------------------------------------------------
params = {'r': 2.8, 'gamma': 0.25, 'rho': 0.35}
tau = 2.0
T = 100.0
dt = 0.01

immune_funcs = [
    ('Linear',        lambda C: f_linear(C)),
    ('Michaelis--Menten', lambda C: f_mm(C, K=0.5)),
    ('Hill-type',     lambda C: f_hill(C, K=0.5, n=4)),
    ('Switch-like',   lambda C: f_switch(C, kappa=10.0, C_star=0.5))
]

exposure_funcs = [
    ('Constant',    lambda t, C: exposure_constant(t, C)),
    ('Periodic',    lambda t, C: exposure_periodic(t, C)),
    ('Event-driven',lambda t, C: exposure_event(t, C)),
    ('Adaptive',    lambda t, C: exposure_adaptive(t, C))
]

# -------------------------------------------------
# Run factorial sweep (16 scenarios)
# -------------------------------------------------
results = {k: np.zeros((4, 4)) for k in ['C_max','t_max','Lag','AUC','Amp',
                                          'C_final','N_peaks','T_rel','lambda','R2']}
results['regime'] = [['' for _ in range(4)] for _ in range(4)]
trajectories = {}

for i, (ilab, f_imm) in enumerate(immune_funcs):
    for j, (elab, f_exp) in enumerate(exposure_funcs):
        t, C, Ca, m = run_experiment(f_imm, f_exp, params, tau=tau, T=T, dt=dt)
        trajectories[(i, j)] = (t, C, Ca)
        for k in results:
            if k != 'regime':
                results[k][i, j] = m[k]
            else:
                results[k][i][j] = m[k]

# -------------------------------------------------
# Figure 1: immune & exposure panels
# -------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))

C_vals = np.linspace(0, 1.5, 200)
for lab, func in immune_funcs:
    ax1.plot(C_vals, func(C_vals), lw=2.5, label=lab)
ax1.set_xlabel('Viral load $C$', fontsize=12)
ax1.set_ylabel('Immune response $f(C)$', fontsize=12)
ax1.set_title('(A) Immune response functions', fontsize=12)
ax1.legend(frameon=True, fontsize=10)
ax1.set_xlim([0, 1.5])
ax1.set_ylim([0, 1.2])
ax1.grid(True, alpha=0.3)

t_vals = np.linspace(0, 50, 1000)
# Representative adaptive profile (computed along a trajectory)
t_adap, C_adap, _, _ = run_experiment(lambda C: f_linear(C),
                                      lambda t, C: exposure_adaptive(t, C),
                                      params, tau=tau, T=50, dt=dt)
Ca_adap = np.array([exposure_adaptive(ti, ci) for ti, ci in zip(t_adap, C_adap)])

ax2.plot(t_vals, [exposure_constant(ti, 0) for ti in t_vals], lw=2.5, label='Constant')
ax2.plot(t_vals, [exposure_periodic(ti, 0) for ti in t_vals], lw=2.5, label='Periodic')
ax2.plot(t_vals, [exposure_event(ti, 0) for ti in t_vals], lw=2.5, label='Event-driven')
ax2.plot(t_adap, Ca_adap, lw=2.5, ls='--', label='Adaptive (representative)')
ax2.set_xlabel('Time $t$ (days)', fontsize=12)
ax2.set_ylabel('Exposure $C_a(t)$', fontsize=12)
ax2.set_title('(B) Exposure profiles', fontsize=12)
ax2.legend(frameon=True, fontsize=10)
ax2.set_xlim([0, 50])
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('figures/immune_exposure_panels.eps', format='eps', dpi=300)
plt.savefig('figures/immune_exposure_panels.png', format='png', dpi=300)
plt.close()

# -------------------------------------------------
# Figure 2: representative time series
# -------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(11, 8))
panel_cfg = [((0, 0), '(A)'), ((1, 1), '(B)'), ((2, 2), '(C)'), ((3, 3), '(D)')]
positions = [(0, 0), (0, 1), (1, 0), (1, 1)]

for (i, j, panel), (row, col) in zip(panel_cfg, positions):
    t, C, _, _ = trajectories[(i, j)]
    ax = axes[row, col]
    ax.plot(t, C, 'k-', lw=1.5)
    ax.set_xlabel('Time (days)', fontsize=11)
    ax.set_ylabel('Viral load $C(t)$', fontsize=11)
    ax.set_title(f"{panel} {immune_funcs[i][0]} + {exposure_funcs[j][0]}", fontsize=11)
    ax.set_xlim([0, T])
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('figures/viral_time_series_panels.eps', format='eps', dpi=300)
plt.savefig('figures/viral_time_series_panels.png', format='png', dpi=300)
plt.close()

# -------------------------------------------------
# Helper: grouped bar plots
# -------------------------------------------------
def grouped_bar(matrix, ylabel, fname, fmt='.3f'):
    fig, ax = plt.subplots(figsize=(9, 5.5))
    x = np.arange(4)
    width = 0.18
    imm_labels = ['Linear', 'MM', 'Hill', 'Switch']
    exp_labels = ['Constant', 'Periodic', 'Event-driven', 'Adaptive']
    colors = ['#0072BD', '#D95319', '#EDB120', '#7E2F8E']

    for i in range(4):
        offset = width * (i - 1.5)
        vals = matrix[i]
        bars = ax.bar(x + offset, vals, width, label=imm_labels[i],
                      color=colors[i], edgecolor='black', linewidth=0.5)
        for bar, val in zip(bars, vals):
            if not (np.isnan(val) or np.isinf(val)):
                ax.text(bar.get_x() + bar.get_width() / 2.0, val,
                        f'{val:{fmt}}', ha='center', va='bottom',
                        fontsize=7, rotation=90)

    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_xticks(x)
    ax.set_xticklabels(exp_labels, fontsize=11)
    ax.legend(fontsize=10)
    ax.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'figures/{fname}.eps', format='eps', dpi=300)
    plt.savefig(f'figures/{fname}.png', format='png', dpi=300)
    plt.close()

# -------------------------------------------------
# Figures 3–10: continuous metrics
# -------------------------------------------------
metric_specs = [
    ('C_max',   'Peak viral load $C_{\\max}$',                'Cmax_tau2',        '.3f'),
    ('t_max',   'Time to peak $t_{\\max}$ (days)',            'tmax_tau2',        '.2f'),
    ('Lag',     'Response lag $\\mathcal{L}$ (days)',         'Lag_tau2',         '.3f'),
    ('AUC',     'Cumulative viral burden (AUC)',              'AUC_tau2',         '.1f'),
    ('Amp',     'Asymptotic amplitude',                       'Amp_tau2',         '.4f'),
    ('C_final', 'Terminal viral load $C_{\\mathrm{final}}$',  'Cfinal_tau2',      '.3f'),
    ('N_peaks', 'Number of local maxima $N_{\\mathrm{peaks}}$','Npeaks_tau2',     '.0f'),
    ('lambda',  'Post-peak decay rate $\\lambda_{\\mathrm{decay}}$','lambda_decay_tau2','.4f')
]

for key, ylab, fname, fmt in metric_specs:
    grouped_bar(results[key], ylab, fname, fmt)

# -------------------------------------------------
# Figure 11: Dynamical regime heatmap
# -------------------------------------------------
regime_map = {'stable': 0, 'equilibrium': 1, 'limit_cycle': 2,
              'non_convergent': 3, 'unstable': 4}
regime_txt = {0: 'STABLE', 1: 'EQ', 2: 'LIMIT', 3: 'NC', 4: 'UNSTABLE'}

fig, ax = plt.subplots(figsize=(9, 6))
data_num = np.array([[regime_map.get(results['regime'][i][j], 3)
                      for j in range(4)] for i in range(4)])
im = ax.imshow(data_num, cmap='RdYlGn_r', aspect='auto', vmin=0, vmax=4)

ax.set_xticks(range(4))
ax.set_yticks(range(4))
ax.set_xticklabels(['Constant', 'Periodic', 'Event-driven', 'Adaptive'],
                   rotation=45, ha='right')
ax.set_yticklabels(['Linear', 'Michaelis--Menten', 'Hill-type', 'Switch-like'])
ax.set_xlabel('Exposure function $C_a(t)$', fontsize=12)
ax.set_ylabel('Immune response function $f(C)$', fontsize=12)

for i in range(4):
    for j in range(4):
        val = data_num[i, j]
        txt = regime_txt[val]
        col = 'white' if val >= 3 else 'black'
        ax.text(j, i, txt, ha='center', va='center', color=col,
                fontsize=11, weight='bold')

cbar = plt.colorbar(im, ax=ax, ticks=range(5))
cbar.ax.set_yticklabels(['Stable', 'Equilibrium', 'Limit Cycle',
                         'Non-convergent', 'Unstable'])
cbar.set_label('Dynamical Regime', fontsize=11)
plt.tight_layout()
plt.savefig('figures/Dynamical_Regime_tau2.eps', format='eps', dpi=300)
plt.savefig('figures/Dynamical_Regime_tau2.png', format='png', dpi=300)
plt.close()

print("[run_main] Figures 1–11 generated in figures/")
