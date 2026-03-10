*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-
* 
* MeshTran-Femtic configuration file
*
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-

MESH_NATURE = native

DEM_FILE = TOTPO_KAL_BBMT_FLAT_475m.xyz
DEM_UNITS = kilometers 

TOPO_FILE = topography_for.dat 
BATHY_FILE = bathymetry_for.dat
COSLI_FILE = coast_line.dat

XMIN_DOM = -40.0
XMAX_DOM =  40.0
YMIN_DOM = -50.0
YMAX_DOM =  50.0
ZMIN_DOM = -50.0
ZMAX_DOM =  50.0

PAD_X = 10.0
PAD_Y = 10.0

HAS_SEA = NO
SEA_LEVEL = 0.0

## => Meshing for sites and ellipsoides at the surface
NUM_ELIPSES = 10
NUM_ESFERAS = 5
ROTATION = 0.0
minRADIO = 0.1              ! 0.3 1.0 3.0 5.0
maxRADIO = 5.0
minEDGES = 0.10             !0.05 0.10 0.30 0.50

SITEpadding = 30.0         ! Margen extra sobre el radio de los sites (km)
FARelemSIZE = 50.0        ! Tamaño máximo de triángulo en el borde (km)
SURF_RESOLUTION = 0.8     ! len en el centro (km)
GROWTH = 2.5              ! Cuánto aumenta 'len' en cada elipse hacia afuera



## => Refinement for 3d mesh and sites
REFI_ELLIPSES = 10
HIGH_RESOL_LAYER = 1.6    !thikness of high resolution from surface = HRL*TargDep


NEAR_FIELD_RESOL = 5      !min elem size at sites (core_resolution)
FAR_FIELD_RESOL =  50    !max elem size 

ITER_TET_REFI = 6
BKGRD_RHO = 100
F_MIN_HZ = 0.001

!esto controla radio a en makeMTR
H_PADDING = 15             !horizontal extension beyond last site, a1 = max_distance_site x H_PADDING 
TARGET_DEPTH = 20

!target depth controlled by following 3
SITE_ELLIPSES = 6
V_PADDING = 2.8             !vertical extension beyond the refinement


!esto seria para malla de parametros en makeMTR no influye
FRAC_SKIN_DEPTH = 5.55    !No quiero que el tamaño de elemento (o parámetro) sea mayor que el 25% de la skin depth
## => Regions Atributes in the model
ELLIPSES = 9 
REGIONS = 2
LOCATION = 1.0 2.0 -39.0 0.0 0.0 39.0 	!xyz xyz
ID_REGION = 10 20			!reg1 reg2
RHO_REGIONS = 1.0e9 1.0e2		!rho1 rho2
REP_PARTITION = -1 9
FIX_RESISTIVITY = 1 0

## => Refinement for mesh parameter based on spheres
PARAM_ESFER = 2
PARAM_RADIOS = 3.0 5.0 
PARAM_EDGES = 2.0 3.0 


