import numpy as np

def simulate_dde(t, dt, tau, f, exposure, params, C0=0.01):
    """
    Explicit Euler integration of the delayed viral dynamics model.

    dC/dt = r*C*(1-C) - gamma*C*f(C(t-tau)) + rho*Ca(t)

    Parameters
    ----------
    t        : 1D array, uniform time grid
    dt       : scalar, time step (must equal t[1]-t[0])
    tau      : scalar, immune delay
    f        : callable, immune response function f(C)
    exposure : callable, exposure function Ca(t, C)
    params   : dict with keys 'r', 'gamma', 'rho'
    C0       : float, initial viral load / history constant

    Returns
    -------
    C        : 1D array, viral load trajectory
    """
    n_steps = len(t)
    C = np.zeros(n_steps)

    # history on [-tau, 0]
    delay_steps = int(np.ceil(tau / dt))
    C[:delay_steps] = C0

    for k in range(delay_steps, n_steps):
        # linear interpolation for delayed state
        t_delayed = t[k] - tau
        idx = int(t_delayed / dt)
        frac = (t_delayed - idx * dt) / dt
        idx = min(idx, n_steps - 2)
        C_tau = C[idx] + frac * (C[idx + 1] - C[idx])

        # exposure at previous step (explicit)
        Ca = exposure(t[k], C[k - 1])

        dC = (params["r"] * C[k - 1] * (1.0 - C[k - 1])
              - params["gamma"] * C[k - 1] * f(C_tau)
              + params["rho"] * Ca)

        C[k] = max(C[k - 1] + dt * dC, 0.0)

    return C
