"""Reproduce every table and figure in the paper.

Run from the repository root:  python run_all.py
"""
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from src.model import f, steady_state, convergence_speed, exact
from src.solvers import euler, heun, rk4

warnings.filterwarnings("ignore")
plt.rcParams.update({"font.family": "serif", "font.size": 11,
                     "axes.grid": True, "grid.alpha": 0.3})

SOLVERS = (("Euler", euler), ("Heun", heun), ("RK4", rk4))

# ---------- 1. Calibration from the data ----------
d = pd.read_csv("data/nigeria_solow_data.csv").set_index("year").loc[1990:2023]
G = 0.02  # technology growth, fixed at the common literature value
calib = {"s": d["s"].mean(), "alpha": d["alpha"].mean(),
         "n": d["n"].mean(), "delta": d["delta"].mean(), "g": G}
pd.DataFrame({"parameter": list(calib), "value": [round(v, 5) for v in calib.values()]}
             ).to_csv("results/calibration.csv", index=False)
par = (calib["s"], calib["alpha"], calib["n"], calib["delta"], calib["g"])
kstar = steady_state(*par)
lam = convergence_speed(calib["alpha"], calib["n"], calib["delta"], calib["g"])
k0, T = 0.2 * kstar, 100
print("k* = %.4f, lambda = %.5f" % (kstar, lam))

# ---------- 2. Error table ----------
hs = [0.05, 0.1, 0.25, 0.5, 1, 2, 5, 10]
rows = []
for h in hs:
    row = {"h_years": h}
    for name, sol in SOLVERS:
        t, k = sol(f, k0, T, h, args=par)
        row[name] = np.max(np.abs(k - exact(t, k0, *par)))
    rows.append(row)
err = pd.DataFrame(rows)
err.to_csv("results/error_table.csv", index=False, float_format="%.3e")
orders = {name: np.polyfit(np.log(err.h_years[:6]), np.log(err[name][:6]), 1)[0]
          for name, _ in SOLVERS}
print("fitted orders:", {k: round(v, 2) for k, v in orders.items()})

# ---------- 3. Stability limits ----------
def converges(solver, p, h, steps=3000):
    ks = steady_state(*p)
    t, k = solver(f, 0.9 * ks, steps * h, h, args=p)
    return np.isfinite(k[-1]) and abs(k[-1] - ks) < abs(0.1 * ks) * 0.01

def threshold(solver, p, lo=0.01, hi=500):
    for _ in range(40):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if converges(solver, p, mid) else (lo, mid)
    return lo

rows = []
for label, a in (("Nigeria", calib["alpha"]), ("Textbook", 0.33)):
    p = (calib["s"], a, calib["n"], calib["delta"], calib["g"])
    l = convergence_speed(a, calib["n"], calib["delta"], calib["g"])
    rows.append({"case": label, "alpha": round(a, 4), "lambda": round(l, 5),
                 "Euler_theory": round(2 / l, 1),
                 "Euler": round(threshold(euler, p), 1),
                 "Heun": round(threshold(heun, p), 1),
                 "RK4": round(threshold(rk4, p), 1)})
pd.DataFrame(rows).to_csv("results/stability_limits.csv", index=False)

# ---------- 4. Figures ----------
h = 10.0
tt = np.linspace(0, T, 400)
fig, ax = plt.subplots(figsize=(7, 4.5))
ax.plot(tt, exact(tt, k0, *par), "k-", lw=2.2, label="Exact solution")
for (name, sol), st in zip(SOLVERS, ("o--", "s--", "^--")):
    t, k = sol(f, k0, T, h, args=par)
    ax.plot(t, k, st, ms=4.5, lw=1.2, label=name)
ax.axhline(kstar, color="grey", ls=":", lw=1)
ax.text(2, kstar - 1.3, "steady state $k^*$ = %.1f" % kstar, fontsize=9, color="grey")
ax.set_xlabel("Time (years)"); ax.set_ylabel("Capital per effective worker, $k$")
ax.set_title("Solow paths for Nigeria, step size h = 10 years")
ax.legend(loc="lower right"); fig.tight_layout()
fig.savefig("figures/fig1_solution_paths.png", dpi=200)

fig, ax = plt.subplots(figsize=(7, 4.5))
for (name, sol), st in zip(SOLVERS, ("o-", "s-", "^-")):
    ax.loglog(err.h_years, err[name], st, label=name)
ax.set_xlabel("Step size h (years)"); ax.set_ylabel("Maximum absolute error")
ax.set_title("Error against step size (log-log scale)")
ax.legend(); fig.tight_layout()
fig.savefig("figures/fig2_error_vs_step.png", dpi=200)

fig, ax = plt.subplots(figsize=(7, 4.5))
for z, st, lab in ((0.5, "o-", "h$\\lambda$ = 0.5 (h = 21 yrs): smooth decay"),
                   (1.5, "s-", "h$\\lambda$ = 1.5 (h = 63 yrs): oscillating, decaying"),
                   (2.5, "^-", "h$\\lambda$ = 2.5 (h = 105 yrs): oscillating, growing")):
    hh = z / lam
    t, k = euler(f, 0.9 * kstar, 6 * hh, hh, args=par)
    ax.plot(np.arange(len(k)), (k - kstar) / kstar, st, label=lab)
ax.axhline(0, color="k", lw=1)
ax.set_xlabel("Step number"); ax.set_ylabel("Gap from steady state, $(k-k^*)/k^*$")
ax.set_title("Euler stability regimes (Nigeria, $\\lambda$ = %.4f)" % lam)
ax.legend(fontsize=9); fig.tight_layout()
fig.savefig("figures/fig3_euler_stability.png", dpi=200)
print("Done: tables in results/, figures in figures/")
