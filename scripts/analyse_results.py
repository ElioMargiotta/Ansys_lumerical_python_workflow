import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import datetime
import yaml
import json
import matplotlib.animation as animation
from scipy.ndimage import gaussian_filter1d
from scipy.integrate import simpson
from scipy.interpolate import make_interp_spline
from build_structure import nd_list, fmax, fmin, sweep_param, tickness_list  # Import the refractive indices and sweep parameters from the build_structure file
from dotenv import load_dotenv


# 1) Load .env from the current working directory (or pass a path to your .env)
load_dotenv()  # same as load_dotenv(dotenv_path=".env")

# 2) Grab the path and import
LUMAPI_PATH = os.getenv("LUMAPI_PATH")  # e.g. C:\Program Files\Lumerical\v242\api\python
if LUMAPI_PATH and LUMAPI_PATH not in sys.path:
    sys.path.append(LUMAPI_PATH)

import lumapi  # type: ignore
import atexit

# ----------------------------
# Helper for JSON serialization
# ----------------------------
def _to_serializable(obj):
    """Best-effort converter for numpy/scalars/complex to JSON-safe types."""
    if isinstance(obj, (np.integer, np.floating)):
        return obj.item()
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, complex):
        return {"real": obj.real, "imag": obj.imag}
    if isinstance(obj, (datetime.datetime, datetime.date)):
        return obj.isoformat()
    try:
        # Fallback: attempt float conversion
        return float(obj)
    except Exception:
        return str(obj)

# Chargement des paramètres utilisateur
with open("./configs/params.yaml", "r", encoding="utf-8") as f:
    params = yaml.safe_load(f)

# Convertit la liste `layers` en dictionnaire indexé par le champ "name"
layers = {
    layer.get("name", f"layer_{i}"): layer for i, layer in enumerate(params.get("layers", []))
}

# Ouverture du fichier simulé
fdtd = lumapi.FDTD(hide=True)
fdtd.load("./structure.fsp")

# Ensure FDTD process is closed even if an exception occurs

def _cleanup():
    try:
        fdtd.close()
    except Exception:
        pass
atexit.register(_cleanup)

# ╔══════════════════════════════════════════════════════════════════════╗
# ║                         FDTD / sweep Angle                            ║
# ╚══════════════════════════════════════════════════════════════════════╝

# Exécution du sweep d'angle
fdtd.runsweep("angle_sweep")
theta_raw = np.array(fdtd.getsweepdata("angle_sweep", "incidence angle"))
transmission_raw = np.array(fdtd.getsweepdata("angle_sweep", "Transmittance"))

theta = theta_raw.flatten()               # [deg]
transmission = transmission_raw.flatten()  # [%] (may be in % or fraction depending on monitor settings)

# monotonic interpolator
x_smooth = np.linspace(theta.min(), theta.max(), 300)
spline = make_interp_spline(theta, transmission, k=3) # spline cubique lissé
y_smooth = np.clip(spline(x_smooth), 0, None)

# Fraction transmise (simulation) — intégrale hémisphérique
cos_theta_all = np.cos(np.radians(theta))
sin_theta_all = np.sin(np.radians(theta))

num_sim = simpson(transmission * cos_theta_all * sin_theta_all, x=theta)
den_sim = simpson(cos_theta_all * sin_theta_all, x=theta)  # = 0.5 si maille fine

pourcentage_transmis = num_sim / den_sim * 100  # en %

# ╔══════════════════════════════════════════════════════════════════════╗
# ║                     Theorique analysis (TMM) / using Ansys           ║
# ╚══════════════════════════════════════════════════════════════════════╝
# compute theta critical angle

# --- robust scalar critical angle (theta_c) ---
# n1 may come as a 0-D/1-D ndarray; extract the scalar safely
n1 = complex(np.asarray(nd_list[0]).ravel()[0])
# domain clamp to avoid tiny numerical overflow of arcsin
arg = np.real(1.0 / n1)
arg = float(np.clip(arg, -1.0, 1.0))
theta_c = float(np.degrees(np.arcsin(arg)))

# construct matrix
if fmin == fmax:
    f = np.array([fmin])  # Fréquence de la source
else:
    f = np.linspace(fmin, fmax, 100)  # Fréquence de la source

start = 0
stop = 80
theta_th = np.arange(float(start), float(stop) + 1.0, 1.0)  # Angle d'incidence de 0 à 80° par pas de 1°

# (your original list)
tickness_list = list(tickness_list)
# add exit half-space thickness
tickness_list.append(0.0)
# ensure ndarray of floats
tickness_list = np.asarray(tickness_list, dtype=float)
nd_list = list(nd_list)

# append AIR as a (1,1) complex array (physically n=1+0j)
nd_list.append(np.array([[1+0j]], dtype=complex))

nf = len(f)   # number of frequency samples
nd = len(tickness_list)  # number of layers (including incident and exit half-spaces)

N = np.zeros((nd, nf), dtype=complex)

for i in range(nd):
    n_i = nd_list[i]
    if np.isscalar(n_i):
        N[i, :] = n_i
    else:
        n_i = np.asarray(n_i)
        if n_i.size != nf:
            raise ValueError(f"nd_list[{i}] must have length {nf}, got {n_i.size}.")
        N[i, :] = n_i

# Ansys STACK/TMM (angle sweep, dispersive indices)
RT = fdtd.stackrt(N, tickness_list, f, theta_th)  # returns dict-like with 'theta','Ts','Tp','lambda',...

# unpack and flatten to 1D (single wavelength case: lam constant)
theta_th_plot = np.ravel(RT["theta"])        # (M,)
Ts = np.ravel(RT["Ts"])                      # (M,)
Tp = np.ravel(RT["Tp"])                      # (M,)
lam = float(np.ravel(RT["lambda"])[0])

# unpolarized average
Tavg = 0.5 * (Ts + Tp)

# ╔══════════════════════════════════════════════════════════════════════╗
# ║                          Outcoupling fraction                        ║
# ╚══════════════════════════════════════════════════════════════════════╝

def _to_fraction(arr):
    arr = np.asarray(arr, dtype=float)
    if np.nanmax(arr) > 1.000001:  # appears to be in %
        return arr / 1
    return arr

# --- FDTD (uses the actual sweep angles 'theta') ---
T_sim_frac = _to_fraction(transmission)
# --- TMM/STACK (uses 'theta_th_plot') ---
T_tmm_frac = _to_fraction(Tavg)

# angles in degrees in your arrays; build radian versions for integration
rad = np.radians(theta)           # FDTD angles
rad_t = np.radians(theta_th_plot) # TMM angles

cos_sim = np.cos(rad);  sin_sim = np.sin(rad)
cos_t   = np.cos(rad_t); sin_t   = np.sin(rad_t)

# --- FDTD: escape fraction (one front surface) ---
outcoupling_fdtd = 0.5 * simpson(T_sim_frac * cos_sim * sin_sim, x=rad)

# beyond-cone (θ > θc). This captures scattering/diffraction-assisted escape.
mask_out_fdtd = theta > theta_c
outcone_fdtd = (0.5 * simpson(T_sim_frac[mask_out_fdtd] *
                              cos_sim[mask_out_fdtd] *
                              sin_sim[mask_out_fdtd],
                              x=rad[mask_out_fdtd])
               ) if np.any(mask_out_fdtd) else 0.0

# --- TMM/STACK: same formulas ---
outcoupling_tmm = 0.5 * simpson(T_tmm_frac * cos_t * sin_t, x=rad_t)

mask_out_tmm = theta_th_plot > theta_c
outcone_tmm = (0.5 * simpson(T_tmm_frac[mask_out_tmm] *
                             cos_t[mask_out_tmm] *
                             sin_t[mask_out_tmm],
                             x=rad_t[mask_out_tmm])
              ) if np.any(mask_out_tmm) else 0.0

print(f"[Outcoupling] Hemispherical (cos(theta)*sin(theta)) - FDTD: {100*outcoupling_fdtd:.2f}%   TMM: {100*outcoupling_tmm:.2f}%")
print(f"[Outcoupling] Beyond escape cone (theta > theta_c ~= {theta_c:.2f} deg) - FDTD: {100*outcone_fdtd:.2f}%   TMM: {100*outcone_tmm:.2f}%")

# ╔══════════════════════════════════════════════════════════════════════╗
# ║                              GRAPHIQUES                              ║
# ╚══════════════════════════════════════════════════════════════════════╝

timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
figure_label = f"run_{timestamp}"

plt.figure(figure_label)
plt.plot(theta, transmission, "+", label="Simulation FDTD")
plt.plot(theta_th_plot, Tavg, linestyle="-", label="T_TMM (unpolarized)")
plt.xlabel("θ (deg)")
plt.plot(x_smooth, y_smooth, "-", label="Spline lissée cubique")
plt.axvline(theta_c, linestyle=":", label="θ_c")
plt.xlabel("θ (deg)")
plt.ylabel("Transmittance")
plt.title(f"T vs θ @ λ = {lam:.0e} m")
plt.grid(True)
plt.legend()
plt.tight_layout()

# -----------------
# Results directories
# -----------------
folder_results = "./results"
folder_plots = os.path.join(folder_results, "plots")
os.makedirs(folder_plots, exist_ok=True)

# NEW: Dedicated logs directory at project root
folder_logs = "./logs"
os.makedirs(folder_logs, exist_ok=True)

# Save plot
plot_path = os.path.join(folder_plots, f"{figure_label}.png")
plt.savefig(plot_path)

# Log CSV 
log_path = os.path.join(folder_results, "res.csv")
col_fdtd = f"%T_hemis_FDTD_{figure_label}"
col_tmm = f"%T_hemis_TMM_{figure_label}"
col_fdtd_out = f"%T_outcone_FDTD_{figure_label}"
col_tmm_out = f"%T_outcone_TMM_{figure_label}"

if os.path.exists(log_path):
    log_df = pd.read_csv(log_path)
else:
    log_df = pd.DataFrame()

# Save hemispherical and out-of-cone metrics (in %)
log_df[col_fdtd] = [100*outcoupling_fdtd]
log_df[col_tmm] = [100*outcoupling_tmm]
log_df[col_fdtd_out] = [100*outcone_fdtd]
log_df[col_tmm_out] = [100*outcone_tmm]

log_df.to_csv(log_path, index=False)

# ╔══════════════════════════════════════════════════════════════════════╗
# ║                                LOGS JSON                             ║
# ╚══════════════════════════════════════════════════════════════════════╝
# Snapshot of configuration + key results, saved with the SAME timestamp
json_log_path = os.path.join(folder_logs, f"{figure_label}.json")

config_snapshot = {
    "timestamp": timestamp,
    "figure_label": figure_label,
    "env": {"LUMAPI_PATH": LUMAPI_PATH},
    "config_file": "./configs/params.yaml",
    "params": params,  # full YAML parameters used
    "sweep": {
        "fmin": fmin,
        "fmax": fmax,
        "theta_start_deg": float(0),
        "theta_stop_deg": float(80),
        "theta_step_deg": 1.0,
    },
    "structure": {
        "thickness_list": tickness_list.tolist(),
        "n_layers": int(nd),
        "critical_angle_deg": theta_c,
        "wavelength_m": lam,
    },
    "results": {
        "outcoupling_fdtd_percent": 100 * outcoupling_fdtd,
        "outcoupling_tmm_percent": 100 * outcoupling_tmm,
        "outcone_fdtd_percent": 100 * outcone_fdtd,
        "outcone_tmm_percent": 100 * outcone_tmm,
        "plot_path": plot_path,
        "csv_log_path": log_path,
    },
}

with open(json_log_path, "w", encoding="utf-8") as jf:
    json.dump(config_snapshot, jf, ensure_ascii=False, indent=2, default=_to_serializable)

# Clean up Lumerical instance
fdtd.close()

print(f"Figure enregistrée : {plot_path}")
print(f"Log CSV mis à jour : {log_path}")
print(f"Config snapshot JSON sauvegardé : {json_log_path}")
