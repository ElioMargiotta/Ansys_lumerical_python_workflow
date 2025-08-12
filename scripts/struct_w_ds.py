import sys
import os
import yaml
import numpy as np
from dotenv import load_dotenv

# 1) Load .env from the current working directory (or pass a path to your .env)
load_dotenv()  # same as load_dotenv(dotenv_path=".env")

# 2) Grab the path and import
LUMAPI_PATH = os.getenv("LUMAPI_PATH")  # e.g. C:\Program Files\Lumerical\v242\api\python
if LUMAPI_PATH and LUMAPI_PATH not in sys.path:
    sys.path.append(LUMAPI_PATH)

import lumapi  # type: ignore

# === Load parameters ===
with open("./configs/params.yaml", "r") as f:
    params = yaml.safe_load(f)
    
layers = params["layers"] #structure
fdtd_conf = params["FDTD"] #x dimension
bounds = fdtd_conf["Boundary"] #boundary conditions (PML, periodic, etc)
source= params["Source"]

# 1) Start a new FDTD session & project
fdtd = lumapi.FDTD(hide=True)
fdtd.load("C:/ansys/base.fsp")
fdtd.eval("selectall; delete;")



# -- Build the layer stack in 2D (XY plane) --
y_prev = None
t_prev = None
tickness_list = []

for i, layer in enumerate(layers):
    name = layer["name"]
    thickness = float(layer["thickness"])
    tickness_list.append(thickness)

    if i == 0:
        y_center = 0
    else:
        y_center = y_prev + t_prev / 2 + thickness / 2

    # ==============
    # Case 1: Simple Layer (rectangle)
    # ==============
    if not layer.get("group", False):
        mat = layer["material"]
        fdtd.addrect()
        fdtd.set("name", name)
        fdtd.set("material", mat)
        fdtd.set("y", y_center)
        fdtd.set("y span", thickness)
        fdtd.set("x", 0)
        fdtd.set("x span", float(params["x_thickness"]))
        fdtd.set("z", 0)
        fdtd.set("z span", float(params["z_thickness"]))

    # ==============
    # Case 2: Group Layer (downshifter)
    # ==============
    else:

        # Create the group
        fdtd.addstructuregroup()
        fdtd.set("name", f"{name}_group")
        
        # --- Select the group structure ---
        fdtd.select(f"{name}_group")

        # Add matrix rectangle inside the group
        fdtd.addrect()
        fdtd.set("name", f"{name}_matrix")
        fdtd.set("material", layer["matrix_material"])
        fdtd.set("y", y_center)
        fdtd.set("y span", thickness)
        fdtd.set("x", 0)
        fdtd.set("x span", float(params["x_thickness"]))
        fdtd.set("z", 0)
        fdtd.set("z span", float(params["z_thickness"]))
        fdtd.addtogroup(f"{name}_group")


        # Compute number of nanoparticles
        particle_radius = layer["nanoparticle"]["radius"]
        particle_density = layer["nanoparticle"]["density"]
        volume_layer = float(params["x_thickness"]) * float(params["z_thickness"]) * thickness
        volume_particle = (4/3) * np.pi * float(particle_radius)**3
        N_particles = int(particle_density * volume_layer / volume_particle)

        # Add particles
        for j in range(N_particles):
            fdtd.addsphere()
            fdtd.set("name", f"{name}_particle_{j}")
            fdtd.set("material", layer["nanoparticle_material"])
            fdtd.set("radius", float(particle_radius))
            fdtd.set("x", (np.random.rand() - 0.5) * float(params["x_thickness"]))
            fdtd.set("y", (np.random.rand() - 0.5) * thickness + y_center)
            fdtd.set("z", (np.random.rand() - 0.5) * float(params["z_thickness"]))
            fdtd.addtogroup(f"{name}_group")

    # Update stacking for next layer
    y_prev = y_center
    t_prev = thickness

# 2) Create & configure the FDTD region
fdtd.addfdtd(dimension= str(fdtd_conf["dimension"]))
fdtd.set("mesh accuracy", 2)  # Sets global mesh accuracy

# Compute total thickness + buffer
total_y_thickness = sum(tickness_list)
buffer = 0.75 * total_y_thickness  # e.g. 20% margin
y_span_FDTD = total_y_thickness + buffer
y0_FDTD =  - (tickness_list[0])/2 + (total_y_thickness) / 2  # center of the region
x_FDTD =  abs(float(params["x_thickness"]))
z_FDTD = abs(float(params["z_thickness"]))

# Position & size
fdtd.set("x", float(fdtd_conf["x"]))
fdtd.set("y", y0_FDTD)
fdtd.set("z", float(fdtd_conf["z"]))

fdtd.set("x span", x_FDTD)
fdtd.set("y span", y_span_FDTD)
fdtd.set("z span", z_FDTD)  # 0 for 2D

# Boundary conditions
fdtd.set("x min bc", bounds["x_min"])
fdtd.set("x max bc", bounds["x_max"])
fdtd.set("y min bc", bounds["y_min"])
fdtd.set("y max bc", bounds["y_max"])


# === 3) Add a plane-wave source near the BOTTOM, injecting downward (negative Y) ===
lambda_min = float(source["lambda_min"])
lambda_max = float(source["lambda_max"])

# Calculate top/bottom edges of FDTD
top_y = y0_FDTD + (y_span_FDTD / 2)
bot_y = y0_FDTD - (y_span_FDTD / 2)


# Source
fdtd.addplane()
fdtd.set("name", "source")
fdtd.set("injection axis", "y")
fdtd.set("direction", "Forward")  # negative Y direction → downward
fdtd.set("y", float(tickness_list[0]/4))  # just inside the lower boundary
fdtd.set("x", 0)
fdtd.set("x span", x_FDTD)
fdtd.set("z", 0)
fdtd.set("z span", 0)
fdtd.set("wavelength start", lambda_min)
fdtd.set("wavelength stop", lambda_max)


# === 4) Add monitors for transmittance at top & bottom ===
# Top monitor
fdtd.addpower()
fdtd.set("name", "T_monitor_top")
fdtd.set("monitor type", "2D Y-normal")  # line along X
fdtd.set("y", top_y - 0.2*buffer)
fdtd.set("x", 0)
fdtd.set("x span", x_FDTD)
fdtd.set("z", 0)
fdtd.set("z span", 0)
# ✅ Override settings to force spectral resolution
fdtd.set("override global monitor settings", True)
fdtd.set("use wavelength spacing", True)
fdtd.set("frequency points", 10)


# Bottom monitor
fdtd.addpower()
fdtd.set("name", "T_monitor_bottom")
fdtd.set("monitor type", "2D Y-normal")
fdtd.set("y", bot_y + 0.2*buffer)
fdtd.set("x", 0)
fdtd.set("x span", x_FDTD)
fdtd.set("z", 0)
fdtd.set("z span", 0)
# ✅ Override settings to force spectral resolution
fdtd.set("override global monitor settings", True)
fdtd.set("use wavelength spacing", True)
fdtd.set("frequency points", 10)



# === Add absorbed power monitor ===
fdtd.addpower()
fdtd.set("name","P_abs")
fdtd.set("monitor type", "2D Z-normal")
fdtd.set("x", 0)
fdtd.set("x span", x_FDTD)
fdtd.set("y", y0_FDTD)
fdtd.set("y span", y_span_FDTD)
fdtd.set("z", 0)
fdtd.set("override global monitor settings", True)
fdtd.set("use wavelength spacing", True)
fdtd.set("frequency points", 10)

# Directly extract thicknesses
t_silicon = tickness_list[0]
t_encap = tickness_list[1]
t_downshifter = tickness_list[2]

# Compute y_center
y_center_downshifter = 0.5*float(t_silicon) + float(t_encap) + 0.5 * float(t_downshifter)
n_dipoles = 10
lambda_emit = 600e-9


output_path = "./structure.fsp"
fdtd.save(output_path)

