import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import datetime
import yaml
import matplotlib.animation as animation
from scipy.ndimage import gaussian_filter1d
from scipy.integrate import simpson
from scipy.interpolate import make_interp_spline
from build_structure import nd_list # Import the refractive indices from the build_structure file
from dotenv import load_dotenv

# -----------------------------
# Helper: robust unit handling
# -----------------------------
def _to_fraction(T_array):
    """Ensure transmission is a fraction [0..1] even if input is in %."""
    T_array = np.asarray(T_array, dtype=float)
    if np.nanmax(T_array) > 1.000001:  # looks like percent
        return T_array / 100.0
    return T_array

def _for_plot_percent(T_array):
    """Return transmission scaled to percent for plotting."""
    T_array = np.asarray(T_array, dtype=float)
    if np.nanmax(T_array) <= 1.000001:  # looks like fraction
        return 1* T_array
    return T_array

# -----------------------------
# Fresnel (double interface) helpers
# -----------------------------
def fresnel_T_s(n_i, n_t, theta_i, theta_t):
    num = 4 * n_i * n_t * np.cos(theta_i) * np.cos(theta_t)
    den = (n_i * np.cos(theta_i) + n_t * np.cos(theta_t)) ** 2
    return num / den

def fresnel_T_p(n_i, n_t, theta_i, theta_t):
    num = 4 * n_i * n_t * np.cos(theta_i) * np.cos(theta_t)
    den = (n_i * np.cos(theta_t) + n_t * np.cos(theta_i)) ** 2
    return num / den

# -----------------------------
# TMM (Transfer Matrix Method)
# -----------------------------
def _cos_theta_in_layer(n0, theta0_deg, nj):
    # Snell: n0 sinθ0 = nj sinθj; allow complex (evanescent) via scimath.sqrt
    theta0 = np.radians(theta0_deg)
    sin0 = np.sin(theta0)
    arg = 1 - (n0 * sin0 / nj) ** 2
    return np.lib.scimath.sqrt(arg)

def _char_matrix_layer(k0, nj, dj, cos_tj, pol):
    # Admittance Y_j (TE=s, TM=p)
    if pol == "s":
        Yj = nj * cos_tj
    else:
        # Avoid divide-by-zero at grazing
        Yj = nj / np.where(np.abs(cos_tj) == 0, 1e-300, cos_tj)
    delta = k0 * nj * dj * cos_tj
    c, s = np.cos(delta), 1j * np.sin(delta)
    M11 = c
    M12 = s / Yj
    M21 = s * Yj
    M22 = c
    return M11, M12, M21, M22

def tmm_multilayer_T(theta_deg, n_in, n_list, d_list, n_out, wavelength_m, pol="unpolarized"):
    """
    Planar multilayer TMM with full interference.
    n_in, n_list, n_out can be complex (supports absorption).
    theta_deg: incidence angles in degrees in the INCIDENT medium n_in.
    pol: 's', 'p', or 'unpolarized'.
    Returns array of power transmittance into n_out.
    """
    theta_deg = np.atleast_1d(theta_deg).astype(float)
    k0 = 2 * np.pi / wavelength_m

    def _T_pol(polcode):
        # Global characteristic matrix initialised to identity
        M11 = np.ones_like(theta_deg, dtype=complex)
        M12 = np.zeros_like(theta_deg, dtype=complex)
        M21 = np.zeros_like(theta_deg, dtype=complex)
        M22 = np.ones_like(theta_deg, dtype=complex)

        cos_in = _cos_theta_in_layer(n_in, theta_deg, n_in)

        for nj, dj in zip(n_list, d_list):
            cosj = _cos_theta_in_layer(n_in, theta_deg, nj)
            m11, m12, m21, m22 = _char_matrix_layer(k0, nj, dj, cosj, polcode)
            # M <- M @ m  (element-wise for each angle)
            M11, M12, M21, M22 = (
                M11 * m11 + M12 * m21,
                M11 * m12 + M12 * m22,
                M21 * m11 + M22 * m21,
                M21 * m12 + M22 * m22,
            )

        cos_out = _cos_theta_in_layer(n_in, theta_deg, n_out)
        if polcode == "s":
            Yin = n_in * cos_in
            Yout = n_out * cos_out
        else:
            Yin = n_in / np.where(np.abs(cos_in) == 0, 1e-300, cos_in)
            Yout = n_out / np.where(np.abs(cos_out) == 0, 1e-300, cos_out)

        denom = (Yin * M11 + Yin * Yout * M12 + M21 + Yout * M22)
        t = 2 * Yin / denom
        T = (np.real(Yout) / np.real(Yin)) * np.abs(t) ** 2
        return np.real_if_close(T)

    if pol == "s":
        return _T_pol("s")
    elif pol == "p":
        return _T_pol("p")
    else:
        Ts = _T_pol("s")
        Tp = _T_pol("p")
        return 0.5 * (Ts + Tp)
    
# 1) Load .env from the current working directory (or pass a path to your .env)
load_dotenv()  # same as load_dotenv(dotenv_path=".env")

# 2) Grab the path and import
LUMAPI_PATH = os.getenv("LUMAPI_PATH")  # e.g. C:\Program Files\Lumerical\v242\api\python
if LUMAPI_PATH and LUMAPI_PATH not in sys.path:
    sys.path.append(LUMAPI_PATH)

import lumapi  # type: ignore
# Chargement des paramètres utilisateur
with open("./configs/params.yaml", "r", encoding="utf‑8") as f:
    params = yaml.safe_load(f)

# Convertit la liste `layers` en dictionnaire indexé par le champ "name"
layers = {
    layer.get("name", f"layer_{i}"): layer for i, layer in enumerate(params.get("layers", []))
}
print(nd_list)

# Ouverture du fichier simulé
fdtd = lumapi.FDTD(hide=True)
fdtd.load("./structure.fsp")

# Exécution du sweep d'angle
fdtd.runsweep("angle_sweep")
theta_raw = np.array(fdtd.getsweepdata("angle_sweep", "incidence angle"))
transmission_raw = np.array(fdtd.getsweepdata("angle_sweep", "Transmittance"))

theta = theta_raw.flatten()               # [deg]
transmission = transmission_raw.flatten()  # [%]

# monotonic interpolator
x_smooth = np.linspace(theta.min(), theta.max(), 300)
spline = make_interp_spline(theta, transmission, k=3) #spline cubique lissé
y_smooth = np.clip(spline(x_smooth), 0, None)

# Fraction transmise (simulation) — intégrale hémisphérique
cos_theta_all = np.cos(np.radians(theta))
sin_theta_all = np.sin(np.radians(theta))

num_sim = simpson(transmission * cos_theta_all * sin_theta_all, x=theta)
den_sim = simpson(cos_theta_all * sin_theta_all, x=theta)  # = 0.5 si maille fine

pourcentage_transmis = num_sim / den_sim * 100  # en %

# ╔══════════════════════════════════════════════════════════════════════╗
# ║                      ANALYSE THÉORIQUE – DOUBLE INTERFACE            ║
# ╚══════════════════════════════════════════════════════════════════════╝

n1 = float(layers["downshifters"]["nd"])  # indice down‑shifter (layer 1)
n2 = float(layers["Anti-glare"]["nd"])  # couche de couplage (layer 2, 3 µm)
d2 = float(layers["Anti-glare"]["thickness"])   # meters
n3 = 1.0   # air
lam_m = float(params["Source"]["lambda_min"])    # meters (single wavelength)
lam_nm = lam_m * 1e9
# → Angle en radians
theta_rad = np.radians(theta)
sin_theta1 = np.sin(theta_rad)

# → 1) refraction layer 1 → layer 2
sin_theta2 = n1 / n2 * sin_theta1
valid12 = sin_theta2 <= 1.0

# → 2) refraction layer 2 → air  (n1 sinθ1 = n3 sinθ3)
sin_theta3 = n1 / n3 * sin_theta1
valid23 = sin_theta3 <= 1.0

valid = valid12 & valid23                 # le rayon traverse les deux interfaces

# Pré‑allocation
ts_tot = np.zeros_like(theta_rad)
tp_tot = np.zeros_like(theta_rad)

# ----- Fonctions Fresnel -----

def fresnel_T_s(n_i, n_t, theta_i, theta_t):
    """Transmittance de puissance (pol s)."""
    num = 4 * n_i * n_t * np.cos(theta_i) * np.cos(theta_t)
    den = (n_i * np.cos(theta_i) + n_t * np.cos(theta_t)) ** 2
    return num / den


def fresnel_T_p(n_i, n_t, theta_i, theta_t):
    """Transmittance de puissance (pol p)."""
    num = 4 * n_i * n_t * np.cos(theta_i) * np.cos(theta_t)
    den = (n_i * np.cos(theta_t) + n_t * np.cos(theta_i)) ** 2
    return num / den

# ----- Calcul couche par couche (seulement pour les angles valides) -----
if np.any(valid):
    theta2 = np.empty_like(theta_rad)
    theta3 = np.empty_like(theta_rad)
    theta2[valid] = np.arcsin(sin_theta2[valid])
    theta3[valid] = np.arcsin(sin_theta3[valid])

    # Interface 1→2
    Ts12 = fresnel_T_s(n1, n2, theta_rad[valid], theta2[valid])
    Tp12 = fresnel_T_p(n1, n2, theta_rad[valid], theta2[valid])

    # Interface 2→3
    Ts23 = fresnel_T_s(n2, n3, theta2[valid], theta3[valid])
    Tp23 = fresnel_T_p(n2, n3, theta2[valid], theta3[valid])

    # Transmittance totale
    ts_tot[valid] = Ts12 * Ts23
    tp_tot[valid] = Tp12 * Tp23

# → Lumière non polarisée
T_analytique = 0.5 * (ts_tot + tp_tot)

# -----------------------------
# TMM (film with thickness d2)
# -----------------------------
T_tmm_frac = tmm_multilayer_T(theta, n_in=n1, n_list=[n2], d_list=[d2], n_out=n3,
                              wavelength_m=lam_m, pol="unpolarized")
T_tmm_plot = _for_plot_percent(T_tmm_frac)


# ----- Fraction pondérée (cosθ sinθ) -----
cos_theta = np.cos(theta_rad)
sin_theta = np.sin(theta_rad)

a_num = simpson(T_analytique * cos_theta * sin_theta, x=theta)
a_den = simpson(cos_theta * sin_theta, x=theta)
fraction_transmise_ponderee_th = a_num / a_den
frac_tmm = simpson(T_tmm_frac * cos_theta_all * sin_theta_all, x=theta) / a_den
theta_c = np.degrees(np.arcsin(np.real(1.0 / n1))) if np.real(n1) > 1 else 90.0

mask_out = theta > theta_c

def _escape_fraction(Tfrac):
    if not np.any(mask_out):
        return 0.0
    return simpson(Tfrac[mask_out] * cos_theta_all[mask_out] * sin_theta_all[mask_out],
                   x=theta[mask_out]) / a_den

escape_tmm = _escape_fraction(T_tmm_frac)


print(f"Fraction transmise pondérée (simulation): {pourcentage_transmis:.2f} % du maximum")
print(f"Fraction transmise pondérée (théorie): {fraction_transmise_ponderee_th:.2%}")

# ╔══════════════════════════════════════════════════════════════════════╗
# ║                              GRAPHIQUES                              ║
# ╚══════════════════════════════════════════════════════════════════════╝

timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
figure_label = f"run_{timestamp}"

plt.figure(figure_label)
plt.plot(theta, transmission, "+", label="Simulation FDTD")
plt.plot(theta, T_analytique, "--", label=f"Fresnel (n₁={n1}→{n2}→1)")
plt.plot(x_smooth, y_smooth, "-", label="Spline lissée cubique")
plt.plot(theta, T_tmm_plot, "-", label=f"TMM (film {d2:.1f} µm, n₂={np.real(n2)})")
plt.axvline(theta_c, linestyle=":", label=r"$\theta_c$")
plt.xlabel("Angle d'incidence θ₁ (°)")
plt.ylabel(f"Transmission (%) à {params["Source"]["lambda_min"]} nm")
plt.title("T(θ₁) simulation vs. modèle analytique")
plt.grid(True)
plt.legend()
plt.tight_layout()

# Logging résultats
folder_results = "./results"
folder_plots = os.path.join(folder_results, "plots")
os.makedirs(folder_plots, exist_ok=True)

plot_path = os.path.join(folder_plots, f"{figure_label}.png")
plt.savefig(plot_path)

# Log CSV 
log_path = os.path.join(folder_results, "res.csv")
col_sim = f"%T_sim_{figure_label}"
col_th = f"%T_th_{figure_label}"

if os.path.exists(log_path):
    log_df = pd.read_csv(log_path)
else:
    log_df = pd.DataFrame()

log_df[col_sim] = [pourcentage_transmis]
log_df[col_th] = [fraction_transmise_ponderee_th * 100]

log_df.to_csv(log_path, index=False)

# Clean up Lumerical instance
fdtd.close()

print(f"Figure enregistrée : {plot_path}")
print(f"Log CSV mis à jour : {log_path}")
