import sys
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import datetime
import yaml
from scipy.integrate import simpson
from scipy.interpolate import make_interp_spline
from dotenv import load_dotenv

# Optional: try to set UTF-8 output; strings below are ASCII-safe anyway
try:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="backslashreplace")
except Exception:
    pass

# 1) Load .env and add Lumerical Python API to path if provided
load_dotenv()
LUMAPI_PATH = os.getenv("LUMAPI_PATH")  # e.g. C:\\Program Files\\Lumerical\\v242\\api\\python
if LUMAPI_PATH and LUMAPI_PATH not in sys.path:
    sys.path.append(LUMAPI_PATH)

import lumapi  # type: ignore

# -----------------------------
# Helpers
# -----------------------------
def _to_fraction(T_array):
    """Ensure transmission is a fraction [0..1] even if input is in %."""
    T_array = np.asarray(T_array, dtype=float)
    if np.nanmax(T_array) > 1.000001:  # looks like percent
        return T_array / 1.0
    return T_array

def _for_plot_percent(T_array):
    """Return transmission scaled to percent for plotting."""
    T_array = np.asarray(T_array, dtype=float)
    if np.nanmax(T_array) <= 1.000001:  # looks like fraction
        return 1 * T_array
    return T_array

def _apply_post_crit_filter(theta_deg, Tfrac, theta_c_deg, thresh=1.01):
    """
    Remove values above 'thresh' for theta > theta_c.
    Returns (theta_valid, Tfrac_valid, mask_valid, n_removed).
    """
    theta_deg = np.asarray(theta_deg, dtype=float)
    Tfrac = np.asarray(Tfrac, dtype=float)
    mask = np.ones_like(theta_deg, dtype=bool)
    post_crit = theta_deg > theta_c_deg
    too_high = Tfrac > thresh
    bad = post_crit & too_high
    mask[bad] = False
    return theta_deg[mask], Tfrac[mask], mask, int(np.count_nonzero(bad))

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
    # Snell: n0 sin(theta0) = nj sin(thetaj); allow complex (evanescent)
    theta0 = np.radians(theta0_deg)
    sin0 = np.sin(theta0)
    arg = 1 - (n0 * sin0 / nj) ** 2
    return np.lib.scimath.sqrt(arg)

def _char_matrix_layer(k0, nj, dj, cos_tj, pol):
    # Admittance Y_j (TE=s, TM=p)
    if pol == "s":
        Yj = nj * cos_tj
    else:  # "p"
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

# -----------------------------
# Load parameters and stack
# -----------------------------
with open("./configs/params.yaml", "r", encoding="utf-8") as f:
    params = yaml.safe_load(f)

layers_list = params.get("layers", [])
layers = {layer.get("name", f"layer_{i}"): layer for i, layer in enumerate(layers_list)}

# Refractive indices (allow complex) and thicknesses from YAML
n1 = complex(layers["downshifters"]["nd"])      # incident medium (semi-infinite)
n2 = complex(layers["Anti-glare"]["nd"])        # film
d2 = float(layers["Anti-glare"]["thickness"])   # meters
n3 = complex(1.0)                                # air

lam_m = float(params["Source"]["lambda_min"])    # meters (single wavelength)
lam_nm = lam_m * 1e9

# -----------------------------
# Open FDTD project and run sweep
# -----------------------------
fdtd = lumapi.FDTD(hide=True)
fdtd.load("./structure.fsp")
fdtd.runsweep("angle_sweep")

theta_raw = np.array(fdtd.getsweepdata("angle_sweep", "incidence angle"))
transmission_raw = np.array(fdtd.getsweepdata("angle_sweep", "Transmittance"))

theta = theta_raw.flatten()           # degrees
T_sim_frac = _to_fraction(transmission_raw.flatten())  # fraction for math

# -----------------------------
# Critical angle and post-critical filtering
# -----------------------------
theta_c = float(np.degrees(np.arcsin(float(np.real(1.0 / n1))))) if np.real(n1) > 1 else 90.0
theta_valid, T_sim_valid, mask_valid, n_removed = _apply_post_crit_filter(theta, T_sim_frac, theta_c, thresh=1.2)

if n_removed > 0:
    print(f"[INFO] Removed {n_removed} post-critical samples with T > 1.2 (fraction units)")

# For plotting: percent scale
T_sim_plot_all = _for_plot_percent(T_sim_frac)
T_sim_plot_valid = _for_plot_percent(T_sim_valid)

# Smooth only over valid points
if len(theta_valid) >= 4:
    x_smooth = np.linspace(theta_valid.min(), theta_valid.max(), 400)
    spline = make_interp_spline(theta_valid, T_sim_valid, k=3)
    y_smooth_frac = np.clip(spline(x_smooth), 0, None)
    y_smooth_plot = _for_plot_percent(y_smooth_frac)
else:
    x_smooth = theta_valid
    y_smooth_plot = T_sim_plot_valid

# -----------------------------
# Fresnel (two interfaces, no interference)
# -----------------------------
theta_rad = np.radians(theta)
sin_theta1 = np.sin(theta_rad)

sin_theta2 = (n1 / n2) * sin_theta1
sin_theta3 = (n1 / n3) * sin_theta1

valid12 = np.abs(sin_theta2) <= 1.0
valid23 = np.abs(sin_theta3) <= 1.0
valid = valid12 & valid23

ts_tot = np.zeros_like(theta_rad, dtype=float)
tp_tot = np.zeros_like(theta_rad, dtype=float)

if np.any(valid):
    theta2 = np.empty_like(theta_rad, dtype=complex)
    theta3 = np.empty_like(theta_rad, dtype=complex)
    theta2[valid] = np.arcsin(sin_theta2[valid])
    theta3[valid] = np.arcsin(sin_theta3[valid])

    Ts12 = fresnel_T_s(n1, n2, theta_rad[valid], theta2[valid])
    Tp12 = fresnel_T_p(n1, n2, theta_rad[valid], theta2[valid])
    Ts23 = fresnel_T_s(n2, n3, theta2[valid], theta3[valid])
    Tp23 = fresnel_T_p(n2, n3, theta2[valid], theta3[valid])

    ts_tot[valid] = np.real(Ts12 * Ts23)
    tp_tot[valid] = np.real(Tp12 * Tp23)

T_fresnel_frac = 0.5 * (ts_tot + tp_tot)
T_fresnel_plot = _for_plot_percent(T_fresnel_frac)

# -----------------------------
# TMM (film with thickness d2)
# -----------------------------
T_tmm_frac = tmm_multilayer_T(theta, n_in=n1, n_list=[n2], d_list=[d2], n_out=n3,
                              wavelength_m=lam_m, pol="unpolarized")
T_tmm_plot = _for_plot_percent(T_tmm_frac)

# -----------------------------
# Hemispherical-weighted fractions and escape beyond cone
# (use only valid simulation samples for integration)
# -----------------------------
cos_theta_all = np.cos(np.radians(theta))
sin_theta_all = np.sin(np.radians(theta))

# For integrals, drop the invalid sim samples
theta_int = theta[mask_valid]
cos_int = cos_theta_all[mask_valid]
sin_int = sin_theta_all[mask_valid]
T_sim_int = T_sim_frac[mask_valid]

den = simpson(cos_int * sin_int, x=theta_int)

def _weighted_fraction(theta_arr, Tfrac_arr, cos_arr, sin_arr):
    return simpson(Tfrac_arr * cos_arr * sin_arr, x=theta_arr) / den

frac_sim = _weighted_fraction(theta_int, T_sim_int, cos_int, sin_int)
frac_fres = _weighted_fraction(theta, T_fresnel_frac, cos_theta_all, sin_theta_all)
frac_tmm = _weighted_fraction(theta, T_tmm_frac, cos_theta_all, sin_theta_all)

mask_out_all = theta > theta_c
mask_out_int = mask_out_all & mask_valid

def _escape_fraction(theta_arr, Tfrac_arr, cos_arr, sin_arr, mask_out):
    if not np.any(mask_out):
        return 0.0
    return simpson(Tfrac_arr[mask_out] * cos_arr[mask_out] * sin_arr[mask_out],
                   x=theta_arr[mask_out]) / den

escape_sim = _escape_fraction(theta, T_sim_frac, cos_theta_all, sin_theta_all, mask_out_int)
escape_fres = _escape_fraction(theta, T_fresnel_frac, cos_theta_all, sin_theta_all, mask_out_all)
escape_tmm = _escape_fraction(theta, T_tmm_frac, cos_theta_all, sin_theta_all, mask_out_all)

print(f"Weighted transmitted fraction (FDTD, filtered): {100*frac_sim:.2f} %")
print(f"Weighted transmitted fraction (Fresnel): {100*frac_fres:.2f} %")
print(f"Weighted transmitted fraction (TMM): {100*frac_tmm:.2f} %")
print(f"Beyond the cone (theta > theta_c ~= {theta_c:.2f} deg):")
print(f"  FDTD (filtered): {100*escape_sim:.2f} %")
print(f"  Fresnel:         {100*escape_fres:.2f} %")
print(f"  TMM:             {100*escape_tmm:.2f} %")

# -----------------------------
# Plot
# -----------------------------
timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
figure_label = f"run_{timestamp}"

plt.figure(figure_label)
plt.plot(theta_valid, T_sim_plot_valid, "+", ms=3, label="FDTD (kept)")
plt.plot(x_smooth, y_smooth_plot, "-", label="FDTD spline (kept)")
plt.plot(theta, T_fresnel_plot, "--", label=f"Fresnel (n1={np.real(n1)} -> {np.real(n2)} -> 1)")
plt.plot(theta, T_tmm_plot, "-", label=f"TMM (film {d2*1e6:.1f} um, n2={np.real(n2)})")
plt.axvline(theta_c, linestyle=":", label="theta_c")
plt.xlabel("Incidence angle theta1 (deg)")
plt.ylabel(f"Transmission (%) at {lam_nm:.0f} nm")
plt.title("T(theta1): FDTD vs Fresnel vs TMM (post-critical outliers filtered)")
plt.grid(True)
plt.legend()
plt.tight_layout()

# -----------------------------
# Save artifacts
# -----------------------------
folder_results = "./results"
folder_plots = os.path.join(folder_results, "plots")
os.makedirs(folder_plots, exist_ok=True)

plot_path = os.path.join(folder_plots, f"{figure_label}.png")
plt.savefig(plot_path)

# Log CSV (append columns with this run's summaries)
log_path = os.path.join(folder_results, "res.csv")
col_sim = f"%T_sim_filtered_{figure_label}"
col_fres = f"%T_fres_{figure_label}"
col_tmm = f"%T_tmm_{figure_label}"
col_sim_escape = f"%T_sim_outcone_filtered_{figure_label}"
col_tmm_escape = f"%T_tmm_outcone_{figure_label}"

if os.path.exists(log_path):
    log_df = pd.read_csv(log_path)
else:
    log_df = pd.DataFrame()

log_df[col_sim] = [100*frac_sim]
log_df[col_fres] = [100*frac_fres]
log_df[col_tmm] = [100*frac_tmm]
log_df[col_sim_escape] = [100*escape_sim]
log_df[col_tmm_escape] = [100*escape_tmm]
log_df.to_csv(log_path, index=False)

# Clean up
fdtd.close()

print(f"Saved figure: {plot_path}")
print(f"Updated CSV log: {log_path}")
