import numpy as np

# ------------------------------
# Immune response architectures
# ------------------------------
def f_linear(C):
    """Proportional immune clearance."""
    return C

def f_mm(C, K=0.5):
    """Michaelis--Menten saturating response."""
    return C / (K + C)

def f_hill(C, K=0.5, n=4):
    """Hill-type ultrasensitive response."""
    return C**n / (K**n + C**n)

def f_switch(C, kappa=10.0, C_star=0.5):
    """Switch-like threshold response."""
    return 1.0 / (1.0 + np.exp(-kappa * (C - C_star)))

# ------------------------------
# Exposure profiles
# ------------------------------
def exposure_constant(t, C, C0=1.0):
    """Persistent background exposure."""
    return C0

def exposure_periodic(t, C, C0=1.0, omega=2.0 * np.pi / 10.0):
    """Seasonal / rhythmic contact patterns."""
    return C0 * (1.5 + 0.5 * np.sin(omega * t)) / 2.0

def exposure_event(t, C, centers=(10.0, 30.0), sigma=1.0):
    """Discrete exposure events (Gaussian pulses)."""
    return sum(np.exp(-((t - c) ** 2) / (2.0 * sigma ** 2)) for c in centers)

def exposure_adaptive(t, C, C0=1.0, alpha=0.5):
    """Behavioural feedback reducing exposure with viral load."""
    return C0 / (1.0 + alpha * C)
