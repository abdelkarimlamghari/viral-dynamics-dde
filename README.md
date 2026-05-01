# Viral Dynamics with Delayed Immune Response

Reproducible computational implementation accompanying the manuscript:

> *Delayed Immune Responses and Heterogeneous Exposure Shape Within-Host Viral Dynamics*

This repository provides a complete numerical framework for simulating delayed within-host viral dynamics under heterogeneous exposure. It includes all simulations, bifurcation analyses, sensitivity experiments, and figure generation used in the manuscript.

---

## Reproducibility

All figures reported in the manuscript can be reproduced by running:

```bash
pip install -r requirements.txt
python scripts/generate_all.py

All generated figures are saved in the figures/ directory.

Individual analyses are also available as separate scripts in scripts/:

run_main.py — Figures 1–11 (immune/exposure panels, time series, metrics)
run_hopf.py — Figure 12 (Hopf bifurcation validation)
run_sensitivity.py — Figure 13 (parameter sensitivity)
run_clinical.py — Figure 14 (clinical decay-rate validation)
Repository Structure
src/ — model definitions, DDE solver, metrics, and experiment runner
scripts/ — figure generation scripts
figures/ — output directory (created automatically)
Numerical Scheme

The DDE is integrated using a method-of-steps approach with an explicit forward Euler scheme (time step $\Delta t = 0.01$ days).
Delayed terms are evaluated by linear interpolation between grid points.

This implementation is consistent with the numerical methodology described in the manuscript.

Requirements
Python ≥ 3.8
NumPy
SciPy
Matplotlib

Install dependencies with:

pip install -r requirements.txt
Citation

If you use this code, please cite the associated manuscript.
