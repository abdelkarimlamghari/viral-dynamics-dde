# Viral Dynamics with Delayed Immune Response

Reproducible computational implementation accompanying the manuscript:

&gt; *Delayed Immune Responses and Heterogeneous Exposure Shape Within-Host Viral Dynamics*

This repository provides a complete numerical framework for simulating delayed within-host viral dynamics under heterogeneous exposure. It includes all simulations, bifurcation analyses, sensitivity experiments, and figure generation used in the manuscript.

---

## Reproducibility

All figures reported in the manuscript can be reproduced by running:

```bash
pip install -r requirements.txt
python scripts/generate_all.py

---
## Individual analyses are also available as separate scripts in scripts/:
| Script               | Output                                                                               |
| -------------------- | ------------------------------------------------------------------------------------ |
| `run_main.py`        | Figures 1–11 (immune/exposure panels, time series, metric bar plots, regime heatmap) |
| `run_hopf.py`        | Figure 12 (Hopf bifurcation validation)                                              |
| `run_sensitivity.py` | Figure 13 (parameter sensitivity analysis)                                           |
| `run_clinical.py`    | Figure 14 (clinical decay-rate validation)                                           |

---
## Repository Structure
viral-dynamics-dde/
├── requirements.txt          # Python dependencies
├── README.md                 # This file
├── src/                      # Core library
│   ├── __init__.py
│   ├── model.py              # Immune response & exposure functions
│   ├── solver.py             # DDE solver (explicit Euler + linear interpolation)
│   ├── metrics.py            # Post-peak decay rate, regime classification, etc.
│   └── experiments.py        # Experiment runner
├── scripts/                  # Figure generation
│   ├── generate_all.py       # Master script: runs all analyses
│   ├── run_main.py           # Main factorial simulations (16 scenarios)
│   ├── run_hopf.py           # Bifurcation diagram
│   ├── run_sensitivity.py    # Parameter robustness
│   └── run_clinical.py       # Clinical validation (SARS-CoV-2 & Ebola)
└── figures/                  # Generated outputs (auto-created)

---
## Numerical Scheme
The delayed differential equation is integrated with an explicit forward-Euler method of steps (dt = 0.01 days). Delayed terms C(t − τ) are evaluated by linear interpolation of previously computed solution values. This scheme matches the numerical methodology described in Section 4 of the manuscript.
Baseline parameters (Table 3): r = 2.8, γ = 0.25, ρ = 0.35, τ = 2.0 days.

---
## Requirements
Python ≥ 3.8
NumPy ≥ 1.21.0
SciPy ≥ 1.7.0
Matplotlib ≥ 3.5.0
Install dependencies with:

```bash
pip install -r requirements.txt

---
## Citation
If you use this code, please cite the associated manuscript.
