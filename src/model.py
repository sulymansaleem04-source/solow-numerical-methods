import numpy as np


def f(k, s, alpha, n, delta, g):
    """Solow ODE: dk/dt = s*k^alpha - (n + delta + g)*k"""
    return s * k**alpha - (n + delta + g) * k


def steady_state(s, alpha, n, delta, g):
    """Steady-state capital per effective worker, k*."""
    return (s / (n + delta + g)) ** (1 / (1 - alpha))


def convergence_speed(alpha, n, delta, g):
    """lambda = (1 - alpha)(n + delta + g)."""
    return (1 - alpha) * (n + delta + g)


def exact(t, k0, s, alpha, n, delta, g):
    """Exact solution k(t) of the Solow ODE (Cobb-Douglas case)."""
    kstar = steady_state(s, alpha, n, delta, g)
    lam = convergence_speed(alpha, n, delta, g)
    a = 1 - alpha
    u = kstar**a + (k0**a - kstar**a) * np.exp(-lam * t)
    return u ** (1 / a)
