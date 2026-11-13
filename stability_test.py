import numpy as np, warnings
warnings.filterwarnings("ignore")
from src.model import f, steady_state, convergence_speed
from src.solvers import euler, heun, rk4


def converges(solver, par, h, steps=3000):
    kstar = steady_state(*par)
    k0 = 0.9 * kstar
    t, k = solver(f, k0, steps * h, h, args=par)
    return np.isfinite(k[-1]) and abs(k[-1] - kstar) < abs(k0 - kstar) * 0.01


def threshold(solver, par, lo, hi):
    for _ in range(40):
        mid = (lo + hi) / 2
        if converges(solver, par, mid):
            lo = mid
        else:
            hi = mid
    return lo


for label, alpha in (("Nigeria", 0.68792), ("Textbook", 0.33)):
    par = (0.21365, alpha, 0.02588, 0.03015, 0.02)
    lam = convergence_speed(alpha, *par[2:])
    print(label, "lambda =", lam)
    for name, sol in (("Euler", euler), ("Heun", heun), ("RK4", rk4)):
        print(" ", name, "limit h =", threshold(sol, par, 0.01, 500))
