# Theory Behind the Angle‑Sweep Analysis

This document summarizes the **physics** and **math** used by the analysis script, tailored to the current configuration:

---
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
