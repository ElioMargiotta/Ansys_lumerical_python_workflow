# FDTD — Angle-Resolved & Hemispherical Extraction in a Planar Optical Stack

> **Ansys Lumerical FDTD** + **Python** pipeline to study the escape-cone transmission and the impact of **surface roughness** in an extraction layer. Full-wave results are compared against a **TMM** (Transfer Matrix Method) baseline (and an optional **Fresnel** shortcut for quick checks).

---

## 1) Project Goals

This repository builds a planar multilayer stack (substrate → extraction/“anti-glare” film → air), runs an **incidence-angle sweep**, and reports:

- angle-resolved transmittance $T(\theta)$,
- optional spectral transmittance $T(\lambda)$,
- a cosine-weighted hemispherical average $\langle T\rangle_{\mathrm{hem}}$,
- an optional **beyond–escape-cone** fraction (how much light exits above the critical angle due to roughness/grating effects).

All geometry, source, boundaries, and roughness statistics come from a single YAML config.

---

## 2) Repository Layout

```
build_structure.py    # constructs geometry (flat or rough), boundaries, source, monitors + angle sweep; saves structure.fsp
analyse_results.py    # runs the sweep, fetches T(θ), computes TMM (and optional Fresnel), hemispherical metrics, and plots
app.py                # orchestrates the pipeline (build → analyze)
configs/params.yaml   # single source of truth (indices, thicknesses, roughness, spans, boundaries, source)
results/              # plots (.png) and logs (.csv)
base.fsp              # base Lumerical project (materials/units) loaded at startup
structure.fsp         # generated FDTD project
```

---

## 3) Requirements & Setup

- **Ansys Lumerical FDTD** v242 (or close)
- **Python** ≥ 3.9
- **Lumerical Python API** on your `PYTHONPATH` (e.g., `C:\Program Files\Lumerical\v242\api\python`)

Python packages:
```
numpy
scipy
pandas
matplotlib
PyYAML
python-dotenv
```

Install:
```bash
pip install -r requirements.txt
# or
pip install numpy scipy pandas matplotlib PyYAML python-dotenv
```

Configure the Lumerical API path via `.env`:
```bash
cp .env.example .env
# edit .env to point to your Lumerical API folder
```

---

## 4) Configuration — `configs/params.yaml` (example)

```yaml
# Layers (bottom → top)
layers:
  - name: Layer1
    material: "1.5"       # must exist in base.fsp
    thickness: 0.1e-6

  - name: layer2
    group: true           # true → two phase surface witn np
    material: "1.6"       # matrix material
    nanoparticle_material: "1.5"
    thickness: 0.1e-6
    nanoparticle:
      radius: 5e-9
      density: 0.1        # filling factor (50%)

  - name: Layer3
    surface: true         # true → rough surface generated & imported; false → flat rectangle
    material: "1.7"
    thickness: 5e-6
    sigma_rms: 0.07e-6    # RMS roughness height (m)
    corr_length_x: 0.1e-6 # correlation length in x (m)
    corr_length_y: 0.1e-6 # correlation length in y (m)
    seed_process: 6       # random seed for reproducibility
    delta: 0.005e-6       # sampling step of height map Z(x,y) (m)

# Lateral spans (Bloch period in x)
x_thickness: 15e-7        # x period for Bloch; periodizes roughness realization
z_thickness: 10e-9        # unused in 2D

# FDTD region & boundaries (y span is derived from layer sum + buffer)
FDTD:
  dimension: "2D"
  x: 0
  z: 0.0
  x_span: 1e-6
  z_span: 0
  Boundary:
    x_min: "Bloch"
    x_max: "Bloch"
    y_min: "PML"
    y_max: "PML"
    z_min: "PML"
    z_max: "PML"

# Narrowband source (single wavelength)
Source:
  lambda_min: 5e-07
  lambda_max: 5e-07
```

> **Tip.** With Bloch in $x$, any “random” rough surface becomes **periodic** over `x_thickness`. This opens **grating orders** beyond the critical angle, which is expected for a periodic model. For diffuse studies, switch to **PML in $x$**, increase lateral span, and average several seeds.

---

## 5) Running the Pipeline

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

**Outputs**
- `results/plots/*.png` — curves $T(\theta)$ (FDTD vs **TMM**; Fresnel optional), $\theta_c$ markers, etc.
- `results/res.csv` — hemispherical averages and beyond-cone fractions.
- `structure.fsp` — generated FDTD project.

---

## 6) Simulation Setup (used in this study)

![Simulation setup](images/setup.png)

**Geometry & materials**
- Layers from YAML (flat or rough top interface with $\sigma_{\mathrm{rms}}, L_x=L_y, \delta$).
- Superstrate: air $n=1$.

**Boundaries & mesh**
- 2D FDTD; **Bloch** in $x$, **PML** in $y$.

**Sources & monitors**
- **Angle sweep (baseline):** plane wave along $+y$; top **Y-normal power** monitor in air returns $T(\theta_1)$.

---

## 7) What the Angle Sweep Computes (FDTD)

Time-domain Maxwell (solved iteratively):
```math
\nabla\times \mathbf{E} = -\,\mu\,\frac{\partial \mathbf{H}}{\partial t},
\qquad
\nabla\times \mathbf{H} = \varepsilon\,\frac{\partial \mathbf{E}}{\partial t}.
```

Per angle, the top monitor yields the **transmittance per period**:
```math
T_{\mathrm{FDTD}}(\theta_1)=\frac{P_{\mathrm{top}}(\theta_1)}{P_{\mathrm{inc}}(\theta_1)},
```
with energy balance (when converged):
```math
R(\theta_1)+T(\theta_1)+A(\theta_1)=1.
```

**Escape-cone edge**:
```math
\theta_c \;=\; \arcsin\!\left(\frac{n_{\mathrm{exit}}}{n_{\mathrm{inc}}}\right)\approx \arcsin\!\left(\frac{1}{\Re\,n_1}\right).
```

---

## 8) Theory Baseline Used by `analyse_results.py` (STACK / TMM)

This repository’s analysis script uses Ansys’ **STACK** optical solver (i.e., **TMM**) via `fdtd.stackrt(...)` to compute the multilayer baseline and overlay it on the FDTD sweep.

### 8.1 How inputs are assembled
- **Frequency grid** $f$: if `fmin==fmax`, a single-frequency array; else `linspace(fmin, fmax, 100)`.
- **Angles** $\theta_{\mathrm{th}}$: `np.arange(0, 81, 1)` (degrees in the incident medium).
- **Thickness list** `tickness_list`: append `0.0` for the **exit half-space** (semi-infinite).
- **Refractive indices** `nd_list`: append **air** $n=1+0j$ as the last medium.
- Build the dispersive index matrix $N \in \mathbb{C}^{(n_\text{layers})\times (n_f)}$ by broadcasting scalars or copying arrays.

### 8.2 Running `stackrt` and extracting $T(\theta)$
```python
RT = fdtd.stackrt(N, tickness_list, f, theta_th)  # returns 'theta','Ts','Tp','lambda',...
theta_th_plot = np.ravel(RT["theta"])
Ts = np.ravel(RT["Ts"])
Tp = np.ravel(RT["Tp"])
Tavg = 0.5 * (Ts + Tp)   # unpolarized average
```

**What TMM includes**
- **Specular** transmission/reflection and **film interference** (Fabry–Pérot) for planar, laterally uniform stacks.

**What TMM does not include**
- Scattering from **rough** or **textured** interfaces. Thus, **post-critical** transmission in FDTD (e.g., grating orders under Bloch) will **not** appear in TMM.

> **Optional Fresnel shortcut.** A two-interface Fresnel product $T_{12}(\theta_1)\,T_{23}(\theta_2)$ can be plotted for quick checks, but TMM is the main baseline here.

---

## 9) Why Roughness Creates Post-Critical Transmission (with Bloch)

Under Bloch BCs the rough top surface $y=h(x)$ is **periodized** with period $\Lambda$. Its lateral spectrum has discrete components at $k_x=mG$, with
```math
G=\frac{2\pi}{\Lambda}, \qquad m\in\mathbb{Z}\ \ (\text{Floquet orders}).
```

Radiation into air follows grating kinematics:
```math
k_{x,\mathrm{out}} \;=\; n_1 k_0 \sin\theta_1 + m\,G .
```
Order $m$ is **propagating** in air iff
```math
\bigl|k_{x,\mathrm{out}}\bigr| \le k_0
\ \Longleftrightarrow\
\left|\,n_1\sin\theta_1 + m\,\frac{\lambda}{\Lambda}\,\right| \le 1.
```

Therefore, several negative orders $m<0$ can open just above $\theta_c$, producing **discrete lobes** in $T(\theta_1)$ even when the specular order ($m=0$) is closed.  
If $\sigma=0$ (perfectly flat), all $|m|\ge 1$ channels vanish → **no** post-critical transmission (FDTD ≈ TMM).

---

## 10) Hemispherical Average & Beyond-Cone Metric

We report a (3D-equivalent) flux-weighted average (simulation itself is 2D):
```math
\langle T\rangle_{\mathrm{hem}}=
\frac{\displaystyle\int_0^{\pi/2} T(\theta)\,\cos\theta\,\sin\theta\,\mathrm{d}\theta}
     {\displaystyle\int_0^{\pi/2} \cos\theta\,\sin\theta\,\mathrm{d}\theta}.
```

Restricting the numerator to $\theta>\theta_c$ gives the **beyond-cone** fraction:
```math
\langle T\rangle_{\mathrm{out\_cone}}=
\frac{\displaystyle\int_{\theta>\theta_c} T(\theta)\,\cos\theta\,\sin\theta\,\mathrm{d}\theta}
     {\displaystyle\int_0^{\pi/2} \cos\theta\,\sin\theta\,\mathrm{d}\theta}.
```

> **Normalization note (/4 factor).** Some workflows divide the metric by 4 to reflect that the model represents a **quarter** of the full hemispherical phase space (e.g., 2D symmetry, per-period normalization). Use **1** for a full-hemisphere metric; keep **/4** only if it matches your physical normalization.

---

## 11) Roughness Parameters (from YAML) and How They’re Used

```yaml
- name: Layer1
  surface: true/false
  material: "1.7"
  thickness: 5e-6
  nd: 1.7
  sigma_rms: 0.07e-6
  corr_length_x: 0.1e-6
  corr_length_y: 0.1e-6
  seed_process: 6
  delta: 0.005e-6
```

| Key | Meaning | Used when |
|---|---|---|
| `name` | Layer label | Always (object naming) |
| `surface` | `false` → flat rectangle; `true` → rough surface via `sroughness` + Surface-import | Selects build branch |
| `material` | Material name/id (FDTD uses this database entry) | Geometry |
| `thickness` | Film thickness (m). For flat: `y span`; for imported: extrusion | Geometry |
| `nd` | Scalar index for theory (TMM); does **not** override `material` in FDTD | Analysis |
| `sigma_rms` | RMS height (m) of roughness | Only if `surface: true` |
| `corr_length_x/y` | Correlation lengths (m) | Only if `surface: true` |
| `seed_process` | Random seed for reproducibility | Only if `surface: true` |
| `delta` | Sampling step (m) of height map $Z(x,y)$ | Only if `surface: true` |

**Modeling notes.** The generated surface is periodic over the chosen span; with Bloch in $x$ the model behaves like a **grating** (discrete diffraction orders). For diffuse studies use **PML in $x$**, enlarge span $\gg L$, and average several seeds.

**Practical tips.** Choose $\Delta\lesssim \min\!\big(L/6,\ \sigma_{\mathrm{rms}}/3,\ \lambda/(10\,n_{\max})\big)$; refine mesh near the rough interface; verify energy balance $R+T(+A)\approx 1$ per angle; increase PML thickness for grazing angles.

---

## 12) Numerical Hygiene

- **Units:** convert monitor outputs to **fractions** for math and to **percent** for display.  
- **Post-critical filtering:** optionally ignore samples with $\theta>\theta_c$ and $T>1.2$ to avoid interpolation/normalization outliers.  
- **Interpolation:** prefer **PCHIP** or clip to $[0,1]$ to avoid cubic-spline overshoot above 1.  
- **Reproducibility:** fix `seed_process`; for statistics, average multiple seeds.

---

## 13) References (Ansys Documentation & Related)

- Bloch boundary conditions in FDTD and MODE — Ansys Optics Help  
  https://optics.ansys.com/hc/en-us/articles/360034382714-Bloch-boundary-conditions-in-FDTD-and-MODE
- Periodic boundary conditions in FDTD and MODE — Ansys Optics Help  
  https://optics.ansys.com/hc/en-us/articles/360034382734-Periodic-boundary-conditions-in-FDTD-and-MODE
- Understanding injection angles in broadband simulations  
  https://optics.ansys.com/hc/en-us/articles/360034382894-Understanding-injection-angles-in-broadband-simulations
- PML boundary conditions in FDTD and MODE  
  https://optics.ansys.com/hc/en-us/articles/360034382674-PML-boundary-conditions-in-FDTD-and-MODE
- Far-field projections in FDTD — overview  
  https://optics.ansys.com/hc/en-us/articles/360034914713-Far-field-projections-in-FDTD-overview
- Far-field from a box of monitors  
  https://optics.ansys.com/hc/en-us/articles/360034915613-Far-field-projections-from-a-box-of-monitors
- TFSF source — tips & best practices  
  https://optics.ansys.com/hc/en-us/articles/360034382934-Tips-and-best-practices-when-using-the-FDTD-TFSF-source
- TFSF source — simulation object  
  https://optics.ansys.com/hc/en-us/articles/360034902093-Total-Field-Scattered-Field-TFSF-source-Simulation-object
- `sroughness` — script command (Gaussian-correlated rough surface)  
  https://optics.ansys.com/hc/en-us/articles/360034926193-sroughness-Script-command
- Grating projections in FDTD — overview  
  https://optics.ansys.com/hc/en-us/articles/360034394354-Grating-projections-in-FDTD-overview
- `gratingangle` / `gratingpolar` — script commands  
  https://optics.ansys.com/hc/en-us/articles/360034927273-gratingangle-Script-command  
  https://optics.ansys.com/hc/en-us/articles/360034407034-gratingpolar-Script-command
- STACK optical solver (TMM) — overview  
  https://optics.ansys.com/hc/en-us/articles/360034914653-STACK-Optical-Solver-Overview
- `stackrt` — script command (TMM R/T)  
  https://optics.ansys.com/hc/en-us/articles/360034406254-stackrt-Script-command

---

## 14) License / Citation

- **License**: CSEM (adapt as needed).  
- **Citation**: Please cite **Ansys Lumerical FDTD** and this repository if you use these results.
