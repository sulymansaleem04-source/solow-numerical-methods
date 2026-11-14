# Numerical Solutions and Stability Analysis of the Solow Growth Model

A comparative study of the Euler, Heun (RK2) and classical Runge-Kutta (RK4) methods, applied to the Solow growth model calibrated to Nigeria (1990 to 2023).

**Author:** Saleem Sulyman

## Summary

The Solow growth model reduces to one nonlinear differential equation for capital per effective worker. With Cobb-Douglas production it has a closed-form solution, so the error of every numerical method can be measured exactly. This repository holds the paper, code, data and results.

## Main results

- Measured convergence orders: Euler 1.00, Heun 2.00, RK4 3.97.
- At a one-year step, the maximum error is 5.0e-02 (Euler), 6.1e-04 (Heun) and 1.9e-08 (RK4).
- Nigerian calibration: capital share 0.688, convergence speed 0.0237 per year, steady-state capital per effective worker 27.40.
- Largest stable step: about 84 years (Euler and Heun) and 116 years (RK4). With a textbook capital share of 0.33 these fall to about 39 and 54 years.

## Repository structure

- `paper/` : the paper (PDF and Word)
- `src/` : the model and the three solvers
- `data/` : Nigerian calibration data (`nigeria_solow_data.csv`)
- `results/` : calibration, error and stability tables
- `figures/` : figures used in the paper
- `run_all.py` : reproduces every table and figure
- `stability_test.py`, `first_test.py` : standalone checks

## How to reproduce
