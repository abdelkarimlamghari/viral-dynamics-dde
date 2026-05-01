import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from src.model import *
from src.solver import simulate_dde

os.makedirs('figures', exist_ok=True)

# Hopf validation parameters (Section 5.1)
r = 1.5
gamma = 3.5
rho = 0.005
K = 0.5
n = 4
kappa = 10.0
Cstar = 0.5

models = ['Linear', 'MM', 'Hill', 'Switch']
colors = ['#0072BD', '#D95319', '#EDB120', '#7E2F8E']

def equilibrium_and_tauc(model):
    def F_eq(C):
        if model == 'Linear':
            fC = C
        elif model == 'MM':
            fC = C / (K + C)
        elif model == 'Hill':
            fC = C**n / (K**n + C**n)
        elif model == 'Switch':
            fC = 1.0 / (1.0 + np.exp(-kappa * (C - Cstar)))
        return r * C * (1.0 - C) - gamma * C * fC + rho

    C_eq = fsolve(F_eq, 0.3)[0]

    if model == 'Linear':
        f_star, df_star = C_eq, 1.0
    elif model == 'MM':
        f_star = C_eq / (K + C_eq)
        df_star = K / (K + C_eq)**2
    elif model == 'Hill':
        f_star = C_eq**n / (K**n + C_eq**n)
        df_star = (n * K**n * C_eq**(n - 1)) / (K**n + C_eq**n)**2
    elif model == 'Switch':
        f_star = 1.0 / (1.0 + np.exp(-kappa * (C_eq - Cstar)))
        df_star = (kappa * np.exp(-kappa * (C_eq - Cstar))) / (1.0 + np.exp(-kappa * (C_eq - Cstar)))**2

    a = r * (1.0 - 2.0 * C_eq) - gamma * f_star
    b = -gamma * C_eq * df_star

    tau_c = np.nan
    if b**2 > a**2:
        omega = np.sqrt(b**2 - a**2)
        tau_c = np.arccos(-a / b) / omega
    return C_eq, tau_c

tau_list = np.linspace(0.1, 3.5, 50)
T = 400.0
dt = 0.01
transient = 300.0

fig, ax = plt.subplots(figsize=(8.5, 5))
handles = []

for m, model in enumerate(models):
    C_eq, tau_c = equilibrium_and_tauc(model)

    if model == 'Linear':
        f_imm = lambda C: f_linear(C)
    elif model == 'MM':
        f_imm = lambda C: f_mm(C, K)
    elif model == 'Hill':
        f_imm = lambda C: f_hill(C, K, n)
    elif model == 'Switch':
        f_imm = lambda C: f_switch(C, kappa, Cstar)

    amps = []
    for tau in tau_list:
        t = np.arange(0.0, T + dt, dt)
        C = simulate_dde(t, dt, tau, f_imm,
                         lambda ti, Ci: 1.0,  # constant exposure Ca=1
                         {'r': r, 'gamma': gamma, 'rho': rho},
                         C0=C_eq * 1.01)
        idx = t > transient
        if np.any(idx):
            a = np.max(C[idx]) - np.min(C[idx])
            amps.append(a if a > 1e-4 else 0.0)
        else:
            amps.append(0.0)

    line, = ax.plot(tau_list, amps, lw=2.5, color=colors[m], label=model)
    handles.append(line)
    if not np.isnan(tau_c):
        ax.axvline(tau_c, color=colors[m], ls='--', lw=1.5, alpha=0.7)

ax.set_xlabel('Immune delay $\\tau$ (days)', fontsize=12)
ax.set_ylabel('Limit cycle amplitude', fontsize=12)
ax.set_xlim([0.1, 3.5])
ax.legend(handles, models, loc='upper left', fontsize=10)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('figures/Bif_Hopf.eps', format='eps', dpi=300)
plt.savefig('figures/Bif_Hopf.png', format='png', dpi=300)
plt.close()
print("[run_hopf] Figure 12 generated.")
