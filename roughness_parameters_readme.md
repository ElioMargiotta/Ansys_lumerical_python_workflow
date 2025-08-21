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
