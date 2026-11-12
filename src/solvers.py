import numpy as np


def euler(f, k0, t_end, h, args=()):
    """Solve dk/dt = f(k) with the Euler method.

    f      : the slope function, called as f(k, *args)
    k0     : starting value of k at t = 0
    t_end  : final time (years)
    h      : step size (years)
    args   : extra parameters passed to f (s, alpha, n, delta, g)
    """
    n_steps = int(round(t_end / h))
    t = np.linspace(0, n_steps * h, n_steps + 1)
    k = np.empty(n_steps + 1)
    k[0] = k0
    for i in range(n_steps):
        k[i + 1] = k[i] + h * f(k[i], *args)
    return t, k
