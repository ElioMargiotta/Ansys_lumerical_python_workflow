import sys
import os
import yaml
import numpy as np
from scipy.constants import c
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
fdtd.load("./base.fsp")
fdtd.eval("selectall; delete;")

# -- Build the layer stack in 2D (XY plane) --
y_prev = None
t_prev = None
tickness_list = []
fmin=c* (float(params["Source"]["lambda_min"]))**(-1)  # Hz
f=fmax=fmin # source max and min frequency are the same (monochromatic source)
nd_list = []

for i, layer in enumerate(layers):
    name = layer["name"]
    thickness = float(layer["thickness"])
    tickness_list.append(thickness)

    if i == 0:
        y_center = 0
    else:
        y_center = y_prev + t_prev / 2 + thickness / 2

    # ==============
    # Case 1: Group Layer (with nanoparticles)
    # ==============
    if  layer.get("group", False):

    # ==============
    # Create the group
        fdtd.addstructuregroup()
        fdtd.set("name", f"{name}_group")
        
        # --- Select the group structure ---
        fdtd.select(f"{name}_group")

        # Add matrix rectangle inside the group
        fdtd.addrect()
        fdtd.set("name", f"{name}_matrix")
        fdtd.set("material", layer["material"])
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
            ni1=fdtd.getfdtdindex(f"{layer["material"]}", f, fmin, fmax)  # Ensure the object is created in the FDTD model
            ni2=fdtd.getfdtdindex(f"{layer["nanoparticle_material"]}", f, fmin, fmax)  # Ensure the object is created in the FDTD model
            nimean= ni1 + particle_density*(ni2-ni1) # Mean refractive index for the group 
            nd_list.append(ni)  # Store the refractive index for this layer
    
    # ==============
    # Case 2: Surface
    # ==============    
    elif layer.get("surface", False): 
        x_span = float(params["x_thickness"])
        y_span = float(thickness)
        sigma_rms = float(layer["sigma_rms"])
        corr_length_x = float(layer["corr_length_x"])
        corr_length_y = float(layer["corr_length_y"])
        seed_process = float(layer["seed_process"])
        delta = float(layer["delta"])  # delta:              sampling resolution of the surface.
#                     For example, delta=10nm with a 1000nm span means the surface will have 100 sample points 
        
        Nx = int(x_span/delta)+1
        Ny = int(y_span/delta)+1
        
        x_w = np.linspace(-x_span/2,  x_span/2,  Nx)    # garde x centré sur 0
        z_w = np.linspace(-thickness/2, thickness/2, Ny)  # provisoire, avant pivot
        H = fdtd.sroughness( x_w, z_w, sigma_rms, corr_length_x, corr_length_y, seed_process )
        folder_results = "C:/ansys/results"
        H_path = os.path.join(folder_results, "rough_surface_H.txt")
        np.savetxt(H_path, H,
                fmt="%.6e",                # scientific, 6 decimals
                delimiter="\t")
        
        fdtd.addimport()
        fdtd.set("name", f"{name}_rough")
        fdtd.set("material", layer["material"])
        fdtd.importsurface(f"{H_path}",1,"m",0,1,thickness,0)   # top face
        fdtd.set("x", 0)  # center on x=0
        fdtd.set("z", 0)  # center on z=0
        fdtd.set("y", y_center)
        fdtd.set("x span", x_span)
        fdtd.set("z span", thickness)
        fdtd.set("y span", 0)
        fdtd.set("first axis",  "x")             # Rotations-tab equivalent
        fdtd.set("rotation 1",  90)              # +90 ° around x ⇒ z→y
        fdtd.set("second axis",  "z")             # Rotations-tab equivalent
        fdtd.set("rotation 2",  180)
        ni=fdtd.getfdtdindex(f"{layer["material"]}", f, fmin, fmax)  # Ensure the object is created in the FDTD model
        nd_list.append(ni)  # Store the refractive index for this layer
       
       
    # ==============
    # Case 3: Simple Layer (rectangle)
    # ============== 
    else :
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
        ni=fdtd.getfdtdindex(f"{layer["material"]}",f, fmin, fmax)  # Ensure the object is created in the FDTD model
        nd_list.append(ni)  # Store the refractive index for this layer 
            
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


# Source (to represente the relfection of the downshifter layer)
fdtd.addplane()
fdtd.set("name", "source")
fdtd.set("injection axis", "y")
fdtd.set("direction", "Forward")  # positive Y direction → upward
fdtd.set("y", 0)  # just inside the lower boundary
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

sweep_param = {
    "Name": "incidence angle",
    "Parameter": "::model::source::angle theta",
    "Type": "Number",
    "Start": 0,
    "Stop": 80,
}

sweep_results = {
    "Name": "Transmittance",
    "Result": "::model::T_monitor_top::T",
}

fdtd.addsweep()
fdtd.setsweep("sweep", "name", "angle_sweep")
fdtd.setsweep("angle_sweep", "type", "Ranges")
fdtd.setsweep("angle_sweep", "number of points", 30)
fdtd.addsweepparameter("angle_sweep", sweep_param)
fdtd.addsweepresult("angle_sweep", sweep_results)

output_path = "./structure.fsp"
fdtd.save(output_path)