import numpy as np
from .solver import simulate_dde
from .metrics import compute_all_metrics

def run_experiment(f, exposure, params, tau=2.0, T=100.0, dt=0.01, C0=0.01):
    """
    Run a single factorial scenario and return trajectory + metrics.
    """
    t = np.arange(0.0, T + dt, dt)
    C = simulate_dde(t, dt, tau, f, exposure, params, C0)
    Ca = np.array([exposure(ti, ci) for ti, ci in zip(t, C)])
    metrics = compute_all_metrics(t, C, Ca, dt)
    return t, C, Ca, metrics
