import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import os

# ---------- CARGA DATOS ----------
# UTM
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
dem_path = os.path.join(BASE_DIR,"preprocessing","PlotWithPython","dem_utm_km.dat")
sites_path = os.path.join(BASE_DIR,"preprocessing","PlotWithPython","sites_utm_km.dat")
maya_path = os.path.join(BASE_DIR,"preprocessing","geometry","topography_for.dat")
sites_maya_path = os.path.join(BASE_DIR,"preprocessing","geometry","sites_coord_elev.dat")
domain_path = os.path.join(BASE_DIR,"preprocessing","geometry","analysis_domain.dat")
coast_path = os.path.join(BASE_DIR,"preprocessing","geometry","coast_line.dat")




dem_utm   = np.loadtxt(dem_path)
sites_utm = np.loadtxt(sites_path, usecols=(1,2), ndmin=2)

# Maya / FEMTIC local
dem_maya   = np.loadtxt(maya_path)
sites_maya = np.loadtxt(sites_maya_path, usecols=(1,2,3), ndmin=2)
domain     = np.loadtxt(domain_path)
coast      = np.loadtxt(coast_path, skiprows=1)


resample = 1            # <-- aquí cambias 3, 5, 10, etc.
ms = 8                  # marker size sites
# ---------------------------



# ---------- FIGURA ----------
fig, axs = plt.subplots(1, 2, figsize=(12, 5))

# ===============================
# SUBPLOT 1 — UTM
# ===============================
ax = axs[0]
ax.scatter(dem_utm[::resample,1], dem_utm[::resample,0],  s=0.5)

ax.scatter(sites_utm[:,1],
           sites_utm[:,0],
           s=ms, marker="v", c="m")

ax.set_title("UTM — DEM + Sites")
ax.set_xlabel("X (km)")
ax.set_ylabel("Y (km)")
ax.set_aspect("equal", adjustable="box")

# ===============================
# SUBPLOT 2 — MALLA / FEMTIC
# ===============================
ax = axs[1]
# ax.scatter(dem_maya[::resample,0],
#            dem_maya[::resample,1],
#            s=0.5)
#
sc = ax.scatter(dem_maya[::resample,1],
                dem_maya[::resample,0],
                c=dem_maya[::resample,2],   # <-- z
                cmap='terrain',
                vmin=-10,
                vmax=10,
                s=1)

fig.colorbar(sc, ax=ax, label="Elevation (m)")





ax.scatter(sites_maya[:,1], sites_maya[:,0], s=ms, marker="v", c="m")


ax.set_title("Analysis domain + Sites in mesh coordinate")
ax.set_xlabel("Y (km)")
ax.set_ylabel("X (km)")
ax.set_aspect("equal", adjustable="box")

# Analysis domain (rectángulo)
xmin, xmax = domain[1,1], domain[1,0]
ymin, ymax = domain[0,1], domain[0,0]
ax.plot([xmin, xmax, xmax, xmin, xmin],
        [ymin, ymin, ymax, ymax, ymin], 'k-')

# Coastline
ax.plot(coast[:,1], coast[:,0], 'r')

ax.set_title("Domain + Sites + Coast")
ax.set_xlabel("Y (km)")
ax.set_ylabel("X (km)")
ax.set_aspect("equal", adjustable="box")


legend_elements = [
    Line2D([0], [0], color='k', linestyle='-',
           label='Analysis Domain '),
    Line2D([0], [0], color='r', linestyle='-',
           label='Coast Line ')
]

fig.legend(
    handles=legend_elements,
    loc='upper center',
    ncol=2,
    frameon=False,
    bbox_to_anchor=(0.86, 0.07)
)

x0 = ax.get_xlim()[0] + 5
y0 = ax.get_ylim()[0] + 5

ax.annotate("N", xy=(x0, y0+15), xytext=(x0, y0),
            arrowprops=dict(arrowstyle="->", lw=1.5),
            ha="center")
ax.annotate("E", xy=(x0+15, y0), xytext=(x0+2, y0),
            arrowprops=dict(arrowstyle="->", lw=1.5),
            va="center")
plt.tight_layout()
plt.show()
