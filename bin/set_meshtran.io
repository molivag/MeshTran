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
minEDGES = 0.02             !0.05 0.10 0.30 0.50

SITEpadding = 5.0         ! Margen extra sobre el radio de los sites (km)
FARelemSIZE = 25.0        ! Tamaño máximo de triángulo en el borde (km)
SURF_RESOLUTION = 0.5     ! len en el centro (km)
GROWTH = 1.8              ! Cuánto aumenta 'len' en cada elipse hacia afuera





## => Refinement for sites at the volumetric mesh iterative
SITE_ELLIPSES = 6
LEN_ELLIPSE =  0.5  1.0  1.5  2.0  3.0  5.0
MAX_EDGE_LEN = 0.10 0.20 0.30 0.50 1.00 2.00
TET_REFINEMENT= 2

## => Global Mesh Refinement based 
##    on ellipses growing from the center of the model




## => Regions Atributes in the model
REGIONS = 2
LOCATION = 1.0 2.0 -39.0 0.0 0.0 39.0 	!xyz xyz
ID_REGION = 20 40			!reg1 reg2
RHO_REGIONS = 1.0e9 1.0e2		!rho1 rho2
REP_PARTITION = -1 9
FIX_RESISTIVITY = 1 0

## => Refinement for mesh parameter based on spheres
PARAM_ESFER = 2
PARAM_RADIOS = 3.0 5.0 
PARAM_EDGES = 2.0 3.0 


