import numpy as np
from scipy.signal import find_peaks
from scipy.stats import linregress
from scipy.interpolate import interp1d

def find_peak_interpolated(t_arr, C_arr):
    """Cubic-spline interpolated peak."""
    idx_approx = np.argmax(C_arr)
    window = min(10, idx_approx, len(C_arr) - idx_approx - 1)
    if window < 2:
        return t_arr[idx_approx], C_arr[idx_approx]
    t_win = t_arr[idx_approx - window:idx_approx + window + 1]
    C_win = C_arr[idx_approx - window:idx_approx + window + 1]
    f_interp = interp1d(t_win, C_win, kind='cubic', fill_value='extrapolate')
    t_fine = np.linspace(t_win[0], t_win[-1], 1000)
    C_fine = f_interp(t_fine)
    i_max = np.argmax(C_fine)
    return t_fine[i_max], C_fine[i_max]

def find_lag_interpolated(t_arr, C_arr, Ca_arr, dt):
    """Cross-correlation lag with cubic refinement."""
    corr = np.correlate(C_arr - np.mean(C_arr), Ca_arr - np.mean(Ca_arr), mode='full')
    lag_idx = np.arange(-len(t_arr) + 1, len(t_arr))
    idx_approx = np.argmax(np.abs(corr))
    window = min(5, idx_approx, len(corr) - idx_approx - 1)
    if window < 2:
        return lag_idx[idx_approx] * dt
    lag_win = lag_idx[idx_approx - window:idx_approx + window + 1] * dt
    corr_win = corr[idx_approx - window:idx_approx + window + 1]
    f_interp = interp1d(lag_win, corr_win, kind='cubic', fill_value='extrapolate')
    lag_fine = np.linspace(lag_win[0], lag_win[-1], 1000)
    return lag_fine[np.argmax(np.abs(f_interp(lag_fine)))]

def compute_lambda_decay(t, C, dt=0.01, tail_fraction=0.2, min_points=8):
    """
    Publication-ready post-peak decay rate.
    Returns (lambda_decay, R2, regime).
    """
    n_tail = int(len(C) * tail_fraction)
    if n_tail < 5:
        return np.nan, np.nan, 'non_convergent'

    tail = C[-n_tail:]
    peaks, _ = find_peaks(tail, prominence=0.01 * np.max(C))
    troughs, _ = find_peaks(-tail, prominence=0.01 * np.max(C))

    if len(peaks) >= 2 and len(troughs) >= 2:
        env_max = np.mean(tail[peaks])
        env_min = np.mean(tail[troughs])
        C_star = (env_max + env_min) / 2.0
    else:
        C_star = np.mean(tail)

    A = np.abs(C - C_star)
    max_idx = np.argmax(C)
    t_post = t[max_idx:]
    A_post = A[max_idx:]

    valid = A_post > 1e-10
    if np.sum(valid) < min_points:
        return np.nan, np.nan, 'non_convergent'

    t_v = t_post[valid]
    A_v = A_post[valid]

    distance = max(1, int(0.5 / dt))
    peaks_post, _ = find_peaks(A_v, distance=distance)

    if len(peaks_post) >= 3:
        t_fit = t_v[peaks_post]
        A_fit = A_v[peaks_post]
    else:
        t_fit = t_v
        A_fit = A_v

    if len(t_fit) < 5:
        return np.nan, np.nan, 'non_convergent'

    logA = np.log(A_fit + 1e-10)
    slope, intercept, r_value, p_value, std_err = linregress(t_fit, logA)
    lambda_decay = -slope
    r2 = r_value ** 2

    std_tail = np.std(tail)
    C_max = np.max(C)

    if r2 < 0.6:
        regime = 'non_convergent'
    elif lambda_decay < -1e-3:
        regime = 'unstable'
    elif np.abs(lambda_decay) < 1e-3:
        regime = 'equilibrium' if std_tail < 0.01 * C_max else 'limit_cycle'
    elif lambda_decay > 0:
        regime = 'stable'
    else:
        regime = 'non_convergent'

    return lambda_decay, r2, regime

def compute_all_metrics(t, C, Ca, dt):
    """Compute the full metric suite for a single trajectory."""
    C_max = np.max(C)
    t_peak, _ = find_peak_interpolated(t, C)
    L = max(0.0, abs(find_lag_interpolated(t, C, Ca, dt)))
    AUC = np.trapz(C, t)

    tail_idx = int(0.8 * len(t))
    Amp = np.max(C[tail_idx:]) - np.min(C[tail_idx:])
    C_final = C[-1]

    peaks, _ = find_peaks(C, prominence=0.01)
    N_peaks = len(peaks)

    # Relative clearance time (10% of peak)
    peak_idx = np.argmax(C)
    post_C = C[peak_idx:]
    post_t = t[peak_idx:]
    thresh = 0.1 * C_max
    idx_below = np.where(post_C <= thresh)[0]
    T_rel = post_t[idx_below[0]] if len(idx_below) > 0 else post_t[np.argmin(post_C)]

    lam, r2, regime = compute_lambda_decay(t, C, dt)

    return {
        'C_max': C_max,
        't_max': t_peak,
        'Lag': L,
        'AUC': AUC,
        'Amp': Amp,
        'C_final': C_final,
        'N_peaks': N_peaks,
        'T_rel': T_rel,
        'lambda': lam,
        'R2': r2,
        'regime': regime
    }
