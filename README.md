# Numerical Solutions and Stability Analysis of the Solow Growth Model

[

![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24213594.svg)

](https://doi.org/10.5281/zenodo.24213594)

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

```
pip install -r requirements.txt
python run_all.py
```

This rebuilds the tables in `results/` and the figures in `figures/` from the data file.

## Data sources

- Population growth: World Bank, World Development Indicators.
- Investment share, labour share and depreciation rate: Penn World Table 11.0 (Feenstra et al., 2015).
- Technology growth is assumed to be 0.02.

## How to cite

Sulyman, S. (2026). *Numerical solutions and stability analysis of the Solow growth model: A comparative study of Euler and Runge-Kutta methods* [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.24213594

## Licence

MIT. See the `LICENSE` file.
