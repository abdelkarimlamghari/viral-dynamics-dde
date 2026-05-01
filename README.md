# Viral Dynamics with Delayed Immune Response

Numerical implementation for the manuscript:

> "Delayed Immune Responses and Heterogeneous Exposure Shape Within-Host Viral Dynamics"

## Reproducibility

All figures can be reproduced by running:

```bash
pip install -r requirements.txt
python scripts/generate_all.py
```

Individual analyses are also available as separate scripts in `scripts/`:
- `run_main.py` — Figures 1–11 (immune/exposure panels, time series, metrics)
- `run_hopf.py` — Figure 12 (Hopf bifurcation validation)
- `run_sensitivity.py` — Figure 13 (parameter sensitivity)
- `run_clinical.py` — Figure 14 (clinical decay-rate validation)

## Repository Structure

- `src/` — model definitions, DDE solver, metrics, and experiment runner
- `scripts/` — figure generation scripts
- `figures/` — output directory (created automatically)

## Numerical Scheme

The DDE is integrated with an explicit forward-Euler method of steps
(dt = 0.01 days). Delayed terms are evaluated by linear interpolation.
This matches the numerical methodology described in the manuscript.

## Requirements

- Python >= 3.8
- NumPy, Matplotlib, SciPy
