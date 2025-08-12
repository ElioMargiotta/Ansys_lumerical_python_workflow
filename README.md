# FDTD Optical Stack — Angle & Hemispherical Transmittance

> Angle-resolved and hemispherical transmittance of a layered optical stack, automated with **Ansys Lumerical FDTD** and **Python**.

## Overview
This repository builds a layered optical model, runs an **incidence-angle sweep**, and analyzes the results to produce:
- angle-resolved transmittance \(T(\theta)\),
- optional spectral transmittance \(T(\lambda)\),
- a cosine-weighted hemispherical metric \(\langle T\rangle_{\text{hem}}\).

The workflow is fully scripted around Ansys Lumerical’s Python API.

## Repository Structure
The current repo contains the following entry points (at repo root):
```
app.py                # Orchestrates the pipeline (build → analyze)
build_structure.py    # Creates geometry, source/monitors, sweep; saves .fsp
analyse_results.py    # Runs sweep, computes metrics, saves plots/CSV
```
You may later introduce folders (`scripts/`, `configs/`, `results/`) if preferred; the scripts work from the root as-is.

## Requirements
- **OS:** Windows 10/11 (recommended for Lumerical; API path below)
- **Ansys Lumerical FDTD** v242 (or close) installed
- **Python:** ≥ 3.9
- **Lumerical Python API:** available at something like  
  `C:\Program Files\Lumerical\v242\api\python`

Python packages:
```
numpy
scipy
pandas
matplotlib
PyYAML
```

Install:
```bash
pip install -r requirements.txt
# or install individually:
pip install numpy scipy pandas matplotlib PyYAML
```

## Configure the Lumerical API Path
Ensure the scripts can import `lumapi`. Inside the Python files, adjust:
```python
import sys, os
sys.path.append(r"C:\Program Files\Lumerical\v242\api\python")  # Update if needed
import lumapi
```

## Configuration (`params.yaml` expected)
Create a `params.yaml` at the repo root (or adjust the script path) with, e.g.:
```yaml
wavelength:
  lambda_min: 0.45e-6
  lambda_max: 0.70e-6
  lambda_analysis: 0.55e-6

FDTD:
  mesh_accuracy: 2
  x_span: 20e-6
  y_span: 10e-6
  boundary_x: PML
  boundary_y: PML

layers:
  - name: substrate
    material: SiO2 (Glass)
    thickness: 1.0e-3
    type: flat

  - name: coating
    material: Polymer (n=1.50)
    thickness: 2.0e-6
    type: flat
    # optional roughness (if used by your build script):
    sigma_rms: 10e-9
    corr_length_x: 200e-9
    corr_length_y: 200e-9

sweep:
  angle_min_deg: 0
  angle_max_deg: 80
  points: 7
  result_monitor: T_monitor_top
  result_quantity: T
```

## How to Run
From the repo root:
```bash
python app.py
```
This will:
1) run `build_structure.py` to create the FDTD model and sweep (saving a `.fsp`), then  
2) run `analyse_results.py` to execute the sweep, produce plots, and log metrics.

You can also run scripts individually:
```bash
python build_structure.py
python analyse_results.py
```

## Outputs
- **Plots**: angle-resolved \(T(\theta)\) at \(\lambda_{\text{analysis}}\) (e.g., `results/plots/T_theta.png` if your script writes there)
- **Optional**: spectral \(T(\lambda)\) at normal incidence
- **CSV log**: cosine-weighted hemispherical averages (simulation vs. analytic baseline), e.g., `results/logs/res.csv`

> If your current scripts write to different folders (e.g., `C:\ansys\results\`), either keep that convention or edit the output paths in the scripts.

## Tips & Good Practices
- **Normalization checks:** At normal incidence, compare against Fresnel/TMM to validate monitors and sign/orientation.
- **Large-angle behavior:** Strong deviations vs. Fresnel may indicate absorption, roughness-induced scattering, or insufficient PML padding.
- **Mesh & boundaries:** Refine mesh near high-contrast interfaces; keep scatterers well-separated from PML.

## Troubleshooting
- `ModuleNotFoundError: lumapi` → Fix the API path (`sys.path.append(...)`) to match your Lumerical install.
- No sweep data → Confirm sweep name (`angle_sweep`), parameter key (e.g., `::model::source::angle theta`), and result mapping (e.g., `::model::T_monitor_top::T`).
- Lumerical not launching → Verify license and that the GUI/API runs once outside Python.

## License
Add your preferred license (e.g., MIT).

## Citation
If this work contributes to a publication, please cite **Ansys Lumerical FDTD** and this repository.
