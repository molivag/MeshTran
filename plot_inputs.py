import numpy as np
import matplotlib.pyplot as plt

# ---------- CARGA DATOS ----------
# UTM
dem_utm   = np.loadtxt("./input_data/pyplots/dem_utm_km.dat")
sites_utm = np.loadtxt("./input_data/pyplots/sites_utm_km.dat", usecols=(1,2), ndmin=2)


# Maya / FEMTIC local
dem_maya   = np.loadtxt("./input_data/topography_for.dat")
sites_maya = np.loadtxt("./input_data/sites_coord_elev.dat", usecols=(1,2,3), ndmin=2)
domain     = np.loadtxt("./input_data/analysis_domain.dat")
coast      = np.loadtxt("./input_data/coast_line.dat", skiprows=1)


resample = 1            # <-- aquí cambias 3, 5, 10, etc.
ms = 8                  # marker size sites
# ---------------------------



# ---------- FIGURA ----------
fig, axs = plt.subplots(1, 2, figsize=(12, 5))

# ===============================
# SUBPLOT 1 — UTM
# ===============================
ax = axs[0]
ax.scatter(dem_utm[::resample,0],
           dem_utm[::resample,1],
           s=0.5)

ax.scatter(sites_utm[:,0],
           sites_utm[:,1],
           s=ms, marker="v", c="m")

ax.set_title("UTM — DEM + Sites")
ax.set_xlabel("X (km)")
ax.set_ylabel("Y (km)")
ax.set_aspect("equal", adjustable="box")

# ===============================
# SUBPLOT 2 — MAYA / FEMTIC
# ===============================
ax = axs[1]
ax.scatter(dem_maya[::resample,0],
           dem_maya[::resample,1],
           s=0.5)

ax.scatter(sites_maya[:,0],
           sites_maya[:,1],
           s=ms, marker="v", c="m")

# Analysis domain (rectángulo)
xmin, xmax = domain[0,0], domain[0,1]
ymin, ymax = domain[1,0], domain[1,1]
ax.plot([xmin, xmax, xmax, xmin, xmin],
        [ymin, ymin, ymax, ymax, ymin], 'k-')

# Coastline
ax.plot(coast[:,0], coast[:,1], 'r')

ax.set_title("Maya / FEMTIC — DEM + Sites + Domain + Coast")
ax.set_xlabel("X (km)")
ax.set_ylabel("Y (km)")
ax.set_aspect("equal", adjustable="box")

plt.tight_layout()
plt.show()
