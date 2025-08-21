# FDTD Optical Stack — Angle & Hemispherical Transmittance

> Angle-resolved and hemispherical transmittance of a layered optical stack, automated with **Ansys Lumerical FDTD** and **Python**.

## Overview
This repository builds a layered optical model, runs an **incidence-angle sweep**, and analyzes the results to produce:
- angle-resolved transmittance $(T(\theta))$,
- optional spectral transmittance $(T(\lambda))$,
- a cosine-weighted hemispherical metric $(\langle T\rangle_{\text{hem}})$.

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
Ensure the scripts can import `lumapi`. Insert path inside `.env.example`.
```bash
cp .env.example .env
```

Optionnal :
resore `.env.example`:

```bash
git checkout -- env.example
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
python scripts\build_structure.py
python scripts\analyse_results.py
```

## Outputs
- **Plots**: angle-resolved \(T(\theta)\) at \(\lambda_{\text{analysis}}\) (e.g., `results/plots/T_theta.png` if your script writes there)
- **Optional**: spectral \(T(\lambda)\) at normal incidence
- **CSV log**: cosine-weighted hemispherical averages (simulation vs. analytic baseline), e.g., `results/logs/res.csv`

> If your current scripts write to different folders (e.g., `C:\ansys\results\`), either keep that convention or edit the output paths in the scripts.

## Troubleshooting
- `ModuleNotFoundError: lumapi` → Fix the API path (`sys.path.append(...)`) to match your Lumerical install.
- Lumerical not launching → Verify license and that the GUI/API runs once outside Python.

---

# Theory Behind the Angle‑Sweep Analysis

This document summarizes the **physics** and **math** used by the analysis script, tailored to the current configuration:

## 0 Simulation setup (used in this study)

![Simulation setup](images/setup.png)

**Geometry & materials**  
- Downshifter: \(n_1=1.5\), thickness \(d_1=100\,\mathrm{nm}\).  
- Anti-glare / extraction film: \(n_2=1.7\), thickness \(d_2=5\,\mu\mathrm{m}\). Optional rough top surface with \(\sigma_{\mathrm{rms}}=70\,\mathrm{nm}\) and correlation length \(L_x=L_y=100\,\mathrm{nm}\).  
- Superstrate: air \(n_3=1\).

**Boundaries & mesh**  
- 2D FDTD. Bloch in \(x\) with period \(\Lambda=1\,\mu\mathrm{m}\). PML in \(y\). Powers are per period.  
- Mesh override around interfaces/roughness; PML spacing chosen to reduce grazing-angle reflection.

**Two source/monitor configurations used in this repository**  
1. **Angle-sweep (specular baseline & grating physics):** plane wave injected along \(+y\) with Bloch phase corresponding to the internal angle \(\theta_1\). Top **power** monitor in air returns \(T(\theta_1)\). This is the configuration compared against Fresnel/TMM in Sections 2–3.  
2. **Re-emission (point source):** a compact **point electric dipole** (or equivalent compact broadband source) **directed toward \(+y\)** is placed in the downshifter near the interface to model spontaneous re-emission/photoluminescence out-coupling. The **detector** is a Y-normal **power** monitor at (or slightly above) the top interface in air measuring the extracted power. Optional far-field projection can be used to obtain angle-resolved escape distributions.

**Modulable layers**  
- All indices, thicknesses, and roughness statistics are read from `configs/params.yaml`, so the stack can be changed without editing the scripts.
"""
## 1 What the FDTD sweep computes

Maxwell (time‑domain, linear isotropic media):
$$
\nabla\times \mathbf{E} = -\,\mu\,\frac{\partial \mathbf{H}}{\partial t},\qquad
\nabla\times \mathbf{H} = \varepsilon\,\frac{\partial \mathbf{E}}{\partial t}.
$$

Bloch periodicity in \(x\) (per‑period fields and powers):
$$
\mathbf{F}(x+\Lambda,y)=\mathbf{F}(x,y)\,e^{i k_x \Lambda},\qquad
k_x = n_1 k_0 \sin\theta_1,\quad k_0=\frac{2\pi}{\lambda}.
$$

Per angle, the script uses power monitors to compute the **transmittance per period**
$$
T_{\mathrm{FDTD}}(\theta_1) \;=\; \frac{P_{\mathrm{top}}(\theta_1)}{P_{\mathrm{inc}}(\theta_1)},
$$
and with consistent normalization/convergence the per‑period energy balance holds
$$
R(\theta_1)+T(\theta_1)+A(\theta_1)=1.
$$

**Critical angle (escape‑cone edge):**
$$
\theta_c=\arcsin\!\left(\frac{n_3}{n_1}\right)
$$
---

## 2 Smooth references used in the plots

### 2.1 Fresnel (two interfaces, no thickness/phase)

At an interface \(n_i\!\to n_t\) with Snell \(n_i\sin\theta_i=n_t\sin\theta_t\), the **power** transmittance for TE/TM is
$$
T_s=\frac{n_t\cos\theta_t}{n_i\cos\theta_i}\,\lvert t_s\rvert^2,\qquad
T_p=\frac{n_t\cos\theta_t}{n_i\cos\theta_i}\,\lvert t_p\rvert^2.
$$
In the three‑medium stack, the simple (no‑multiple‑bounce) model is
$$
T_{\mathrm{Fresnel}}(\theta_1)\;\approx\; T_{12}(\theta_1)\,T_{23}(\theta_2),
$$
and we average TE/TM for unpolarized light. This explains the **sub‑unity baseline** at small \(\theta_1\) and \(T\!\to\!0\) beyond \(\theta_c\).

### 2.2 TMM (Transfer Matrix Method: film thickness + phase)

A film \(n_2,d_2\) is represented by the characteristic matrix
$$
M_j=\begin{pmatrix}
\cos\delta_j & \dfrac{i}{Y_j}\sin\delta_j\\[6pt]
i\,Y_j\sin\delta_j & \cos\delta_j
\end{pmatrix},\qquad
\delta_j=k_0 n_j d_j \cos\theta_j,
$$
with admittance \(Y_j=n_j\cos\theta_j\) (TE) or \(Y_j=n_j/\cos\theta_j\) (TM).  
For the stack \(M=\prod_j M_j\), the transmission into the exit medium is
$$
t=\frac{2Y_{\mathrm{in}}}{Y_{\mathrm{in}}M_{11}+Y_{\mathrm{in}}Y_{\mathrm{out}}M_{12}+M_{21}+Y_{\mathrm{out}}M_{22}},\qquad
T(\theta_1)=\frac{\Re Y_{\mathrm{out}}}{\Re Y_{\mathrm{in}}}\,\lvert t\rvert^2,
$$
again averaged over TE/TM. TMM adds Fabry–Pérot phase while still respecting the escape cone.

---

## 3) Why roughness creates post‑critical transmission (with Bloch)

The rough top surface \(y=h(x)\) is random (RMS \(\sigma\), correlation \(L\)) but **periodized** with period \(\Lambda\). Its lateral spectrum therefore contains discrete components at \(k_x=mG\) with the grating vector
$$
G=\frac{2\pi}{\Lambda},\qquad m\in\mathbb{Z}\quad (\text{Floquet harmonics}).
$$

Radiation into air follows **grating kinematics**:
$$
k_{x,\mathrm{out}} \;=\; n_1 k_0 \sin\theta_1 + m\,G.
$$
Order \(m\) is **propagating** in air iff
$$
\lvert k_{x,\mathrm{out}}\rvert \le k_0 \;\;\Longleftrightarrow\;\;
\bigl\lvert\,n_1\sin\theta_1 + m\,\tfrac{\lambda}{\Lambda}\,\bigr\rvert \le 1.
$$
With \(\lambda/\Lambda=0.5\), several negative orders \(m<0\) are already **open immediately above \(\theta_c\)**, yielding **discrete lobes** in \(T(\theta_1)\) even though the specular order \(m=0\) is closed. Different random seeds change the roughness Fourier amplitudes near \(mG\) (lobe **heights** vary), while the **angles** follow the grating condition.

A thick film (\(d_2=5\,\mu\mathrm{m}\), \(n_2=1.7\)) supports many slab modes; roughness + grating can phase‑match guided/leaky modes to radiate into particular orders (narrow, seed‑dependent peaks).

> If \(\sigma=0\) (plane surface), all \(|m|\ge 1\) channels vanish ⇒ no post‑critical transmission, matching Fresnel/TMM.

---

## 4) “Hemispherical” average and “beyond‑cone” metric

The script reports a **3D‑equivalent** flux‑weighted average (the simulation itself is 2D):
$$
\langle T\rangle \;=\;
\frac{\displaystyle\int T(\theta)\,\cos\theta\,\sin\theta\,\mathrm{d}\theta}
{\displaystyle\int \cos\theta\,\sin\theta\,\mathrm{d}\theta},
$$
and the **escape‑cone** contribution by restricting the numerator to \(\theta>\theta_c\). The weight \(\cos\theta\,\sin\theta\) is the standard hemispherical radiance factor.

---

## 5) Numerical hygiene in the script

- **Units:** helper functions convert monitor output to **fractions** for math and to **percent** for plots.  
- **Post‑critical filter:** samples with \(\theta>\theta_c\) and \(T>1.2\) may be ignored to avoid normalization/interpolation artifacts (physics enforces \(0\le T\le 1\) for passive media).  
- **Interpolation:** use a **monotone** interpolant (PCHIP) or clip to \([0,1]\) to avoid cubic‑spline overshoot above 1.

---

## References (Ansys documentation & related)

1. **Bloch boundary conditions in FDTD and MODE.** Ansys Optics Help.  
   https://optics.ansys.com/hc/en-us/articles/360034382714-Bloch-boundary-conditions-in-FDTD-and-MODE  
2. **Understanding injection angles in broadband simulations.**  
   https://optics.ansys.com/hc/en-us/articles/360034382894-Understanding-injection-angles-in-broadband-simulations  
3. **PML boundary conditions in FDTD and MODE.**  
   https://optics.ansys.com/hc/en-us/articles/360034382674-PML-boundary-conditions-in-FDTD-and-MODE  
4. **Far field projections in FDTD — overview.**  
   https://optics.ansys.com/hc/en-us/articles/360034914713-Far-field-projections-in-FDTD-overview  
5. **Far field projections from a box of monitors.**  
   https://optics.ansys.com/hc/en-us/articles/360034915613-Far-field-projections-from-a-box-of-monitors  
6. **TFSF source — tips and best practices.**  
   https://optics.ansys.com/hc/en-us/articles/360034382934-Tips-and-best-practices-when-using-the-FDTD-TFSF-source  
7. **TFSF source — simulation object.**  
   https://optics.ansys.com/hc/en-us/articles/360034902093-Total-Field-Scattered-Field-TFSF-source-Simulation-object  
8. **sroughness — script command.** (Generates a rough, periodized surface.)  
   https://optics.ansys.com/hc/en-us/articles/360034926193-sroughness-Script-command  
9. **Grating projections in FDTD — overview.**  
   https://optics.ansys.com/hc/en-us/articles/360034394354-Grating-projections-in-FDTD-overview  
10. **gratingangle / gratingpolar — script commands.**  
    https://optics.ansys.com/hc/en-us/articles/360034927273-gratingangle-Script-command  
    https://optics.ansys.com/hc/en-us/articles/360034407034-gratingpolar-Script-command  
11. **STACK optical solver — overview (TMM).**  
    https://optics.ansys.com/hc/en-us/articles/360034914653-STACK-Optical-Solver-Overview  
12. **stackrt — script command (TMM reflection/transmission).**  
    https://optics.ansys.com/hc/en-us/articles/360034406254-stackrt-Script-command  
13. **Periodic boundary conditions in FDTD and MODE.**  
    https://optics.ansys.com/hc/en-us/articles/360034382734-Periodic-boundary-conditions-in-FDTD-and-MODE

---

**TL;DR** — Plane (Fresnel/TMM) respects Snell and cuts off at \(\theta_c\); periodic roughness opens Floquet orders and creates post‑critical lobes. The script compares FDTD to Fresnel/TMM, computes a hemispherical average, and reports the beyond‑cone fraction to quantify the roughness effect.
---

# Roughness Parameters in the YAML (and how the build script uses them)

This note explains each key in your layer block and maps it to what happens in **Lumerical FDTD** via your `build_structure.py`. It also includes practical guidance and links to the **official Ansys Optics** help pages.

```yaml
- name: Anti-glare
  surface: false
  material: "1.7"
  thickness: 5e-6
  nd: 1.7
  sigma_rms: 0.07e-6
  corr_length_x: 0.1e-6
  corr_length_y: 0.1e-6
  seed_process: 6
  delta: 0.005e-6
```

---

## 1) Meaning of each parameter

| Key | What it means | When it is used |
|---|---|---|
| `name` | Human‑readable label of the layer. | Always (for naming the created object). |
| `surface` | **Switch** for how the layer is built: `false` → **flat** rectangle; `true` → **rough surface** generated with `sroughness` and imported as a “Surface import” object. | Controls the whole branch in the script. |
| `material` | Material name/id to assign (here a constant‑index material “1.7”). | Assigned to the rectangle or to the imported‑surface object. |
| `thickness` | Film thickness (m). For the rectangle it is the **y‑span**; for the imported surface it is the extrusion distance defining the solid. | Geometry creation. |
| `nd` | Refractive index used by your **post‑processing** (Fresnel/TMM theory curves). It does **not** override the FDTD material if `material` is set. | Analysis only. |
| `sigma_rms` | **RMS height** (standard deviation) of the rough profile in meters. Larger → stronger scattering. | Only when `surface: true`; passed to `sroughness`. |
| `corr_length_x`, `corr_length_y` | **Correlation lengths** in x/y (m). They set the lateral scale of roughness; larger → smoother undulations. | Only when `surface: true`; passed to `sroughness`. |
| `seed_process` | Random **seed** for the generator; same seed reproduces the same surface. | Only when `surface: true`; passed to `sroughness`. |
| `delta` | **Sampling step** (m) used to discretize the height map Z(x,y) before importing. Smaller `delta` resolves finer features at higher memory/mesh cost. | Only when `surface: true`; used to build the x/y sample vectors and grid size. |

> **Units:** All distances are meters in the YAML (`5e-6` = 5 μm, etc.).

---

## 2) How the script uses the keys (mapping to code)

There are two branches in `build_structure.py`:

### A) `surface: false` → flat film (rectangle)
```python
fdtd.addrect()
fdtd.set("name", name)
fdtd.set("material", layer["material"])
fdtd.set("y", y_center)
fdtd.set("y span", thickness)
fdtd.set("x", 0); fdtd.set("x span", params["x_thickness"])
fdtd.set("z", 0); fdtd.set("z span", params["z_thickness"])
```
This creates a **flat** layer. All roughness keys are ignored.

### B) `surface: true` → rough film (surface import)
1) **Build the sampling grid** from your `x_thickness`, `thickness`, and `delta`:
```python
x_span = float(params["x_thickness"])
y_span = float(thickness)
delta  = float(layer["delta"])

Nx = int(x_span/delta) + 1
Ny = int(y_span/delta) + 1

x_w = np.linspace(-x_span/2,  x_span/2,  Nx)   # center x at 0
z_w = np.linspace(-thickness/2, thickness/2, Ny)  # provisional axis before rotation
```
2) **Generate the rough surface** height map using Lumerical’s **`sroughness`** (Gaussian‑correlated model):
```python
H = fdtd.sroughness(x_w, z_w,
                    sigma_rms, corr_length_x, corr_length_y, seed_process)
```
3) **Import the surface** as geometry, assign material, place & rotate:
```python
fdtd.addimport()
fdtd.set("name", f"{name}_rough")
fdtd.set("material", layer["material"])
fdtd.importsurface(H_path, 1, "m", 0, 1, thickness, 0)   # import Z(x,y) in meters
# center & spans
fdtd.set("x", 0); fdtd.set("x span", x_span)
fdtd.set("y", y_center); fdtd.set("y span", 0)
fdtd.set("z", 0); fdtd.set("z span", thickness)
# orient the imported surface so that the height variation becomes the film top
fdtd.set("first axis", "x"); fdtd.set("rotation 1", 90)
fdtd.set("second axis","z"); fdtd.set("rotation 2",180)
```
This path builds a **solid** whose top interface is the imported rough profile Z(x,y).

---

## 3) Modeling notes specific to Lumerical

### 3.1 Roughness is periodic over the span (important)
`sroughness` operates in **Fourier space** over your chosen span and returns an **IFFT** realization. The generated surface is therefore **periodic** with that span. If the simulation also uses **Bloch** or **periodic** x‑boundaries, the model behaves like a **grating** (discrete diffraction orders beyond the critical angle). For truly **diffuse** scattering studies, switch x‑boundaries to **PML**, enlarge the lateral span (≫ correlation length), and average over multiple seeds.

- Rough surface generator (`sroughness`):  
  https://optics.ansys.com/hc/en-us/articles/360034926193-sroughness-Script-command
- Bloch/Periodic BCs and when to use them:  
  https://optics.ansys.com/hc/en-us/articles/360034382714-Bloch-boundary-conditions-in-FDTD-and-MODE  
  https://optics.ansys.com/hc/en-us/articles/360034382734-Periodic-boundary-conditions-in-FDTD-and-MODE

### 3.2 Surface import object
The **Surface import** geometry takes a sampled height map Z(x,y) with units, and extrudes/places it according to the object settings you set after `importsurface(...)`. Use sufficient sampling density (`delta`) to capture the smallest features you want to model (see tips below).  
https://optics.ansys.com/hc/en-us/articles/360034915553-Surface-Import-Object

---

## 4) Practical tips for choosing the parameters

- **Sampling step `delta`:** choose it to resolve the **correlation length** and the **RMS amplitude**:
  - good starting rule: `delta ≲ min(L/6, sigma_rms/3, λ/(10 n_max))`
  - finer `delta` → more mesh and memory; do mesh override near the rough interface.
- **Span vs. correlation length:** make sure the lateral **x_span** contains **many** correlation lengths (e.g., ≥ 10× L) if you want a representative realization (especially with PML in x).
- **Seed:** keep a fixed seed for reproducibility; for **statistics** (diffuse tails) average over several seeds.
- **Boundaries:** PML in y (open), and either Bloch (for controlled grating physics/angle sweep) or PML (for diffuse). Increase PML thickness/strength for grazing angles.
- **Normalization:** compute transmittance as `P_top / P_inc` per angle; check **R+T(+A) ≈ 1** as a convergence test.

---

## 5) Where `nd` fits (theory curves only)

Your analysis scripts compute Fresnel/TMM predictions using the **scalar index** `nd` from the YAML. In FDTD, the **actual material** is the one named in `material`. For multilayer theory, Ansys provides the **STACK** (TMM) solver and the `stackrt` scripting API:

- STACK (Optical Solver Overview):  
  https://optics.ansys.com/hc/en-us/articles/360034914653-STACK-Optical-Solver-Overview  
- `stackrt` (R/T calculation):  
  https://optics.ansys.com/hc/en-us/articles/360034406254-stackrt-Script-command

---

## 6) Common pitfalls and how to avoid them

- **Post‑critical peaks in a “rough” periodic run:** those are **grating orders** (not numerical noise). Confirm with grating‑projection tools (`gratingangle`, `gratingpolar`).  
  https://optics.ansys.com/hc/en-us/articles/360034394354-Grating-projections-in-FDTD-overview
- **T > 1 spikes:** often due to interpolation overshoot or slight normalization mismatch at grazing angles. Use a **monotone** interpolant (PCHIP) for plotting and clip to [0,1] for display; verify energy balance.
- **Too coarse `delta`:** under‑resolves the surface; scattering is underestimated and results depend on the sampling grid.
- **Too small span with PML in x:** one or two correlation lengths across → strong realization‑to‑realization variability; increase span and average seeds.

---

### TL;DR

- Set `surface: true` to generate and import a **rough** top interface with stats controlled by `sigma_rms`, `corr_length_x/y`, `seed_process`, and representation resolution `delta`.
- With **Bloch** x‑boundaries the surface is **periodic** and acts like a **grating**; with **PML** x‑boundaries and a large span it approximates **diffuse** roughness.
- Use `nd` only for theory curves; the FDTD material comes from `material`.
---

## License
CSEM

## Citation
If this work contributes to a publication, please cite **Ansys Lumerical FDTD** and this repository.
