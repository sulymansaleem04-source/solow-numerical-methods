import numpy as np
from src.model import f, steady_state, convergence_speed, exact
from src.solvers import euler, heun, rk4

s, alpha, n, delta, g = 0.21365, 0.68792, 0.02588, 0.03015, 0.02
par = (s, alpha, n, delta, g)
kstar = steady_state(*par)
lam = convergence_speed(alpha, n, delta, g)
k0 = 0.2 * kstar
T = 100
prev = {}
for h in (1.0, 0.5):
    print("h =", h)
    for name, solver in (("Euler", euler), ("Heun", heun), ("RK4", rk4)):
        t, k = solver(f, k0, T, h, args=par)
        err = np.max(np.abs(k - exact(t, k0, *par)))
        ratio = prev[name] / err if name in prev else float("nan")
        prev[name] = err
        print(name, "max error =", err, "ratio =", ratio)
