# MeshTran-Femtic

Fortran-based preprocessing and meshing pipeline for FEMTIC.

This repository provides a robust wrapper around the FEMTIC preprocessing workflow, including:

- DEM and bathymetry handling
- Coastline and analysis domain generation
- Automated execution of makeTetraMesh (steps 1–4)
- Region tagging and TetGen execution
- Validation and visualization helpers

Designed for large-scale MT inversions and HPC workflows.

## Requirements
- FEMTIC
- makeTetraMesh
- TetGen
- ifort/ifx (or compatible)

## Status
Active development.

### [1] Leer EDI → lat/lon/elev
### [2] latlon_to_utm → siteX_m, siteY_m
### [3] data_m_to_km
#####   ├─ genera siteXkm, siteYkm
#####   ├─ genera demXkm, demYkm
#####   └─ escribe archivos SOLO para plot (diagnóstico)

### [4] compute_center(siteXkm, siteYkm) → x0, y0

### [5] utm_to_mesh_coords
#####  ├─ siteXmesh = siteXkm - x0
#####  ├─ siteYmesh = siteYkm - y0
#####  ├─ demXmesh  = demXkm  - x0
#####  └─ demYmesh  = demYkm  - y0

[6] define_analysis_domain (YA EN COORDENADAS DE MALLA)

[7] check_domain (DEM vs dominio + padding)

[8] write FEMTIC input files
        ├─ topography.dat
        ├─ bathymetry.dat
        ├─ coastline.dat (si aplica)
        ├─ observing_site.dat
        └─ analysis_domain.dat




cosas que faltan, simplificar la linea de costa:
polygons_simplfied = []
for i in range(len(polygons)):
        cl = polygons[i]
        mask = rdp(cl, algo="iter", return_mask=True, epsilon=3) # chech rdp repo for details, epsilon controls the simplification
        plt.plot(cl[:,0],cl[:,1], 'k-',lw=1)
        clm = cl[mask]
        plt.plot(clm[:,0],clm[:,1], 'r-',lw=1)
        polygons_simplfied.append(clm)

plt.xlim(-200,200)
plt.ylim(-200,200)
