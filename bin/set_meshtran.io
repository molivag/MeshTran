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

XMIN_DOM = -50.0
XMAX_DOM =  50.0
YMIN_DOM = -40.0
YMAX_DOM =  40.0
ZMIN_DOM = -40.0
ZMAX_DOM =  40.0

PAD_X = 10.0
PAD_Y = 10.0

HAS_SEA = NO
SEA_LEVEL = 0.0

## => Refinement for sites at the surface
ESFERAS = 5
RADIOS = 0.1 0.3 1.0 3.0 5.0
EDGES = 0.02 0.05 0.10 0.30 0.50

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

## => Refinement for sites at the volumetric mesh
SITE_ELLIPSES = 6
LEN_ELLIPSE =  0.5  1.0  1.5  2.0  3.0  5.0
MAX_EDGE_LEN = 0.10 0.20 0.30 0.50 1.00 2.00

## => Global Mesh Refinement based 
##    on ellipses growing from the center of the model
ROTATION=0.0
NUM_ELIPSES=10

## => Iterative tetgen refinement
TET_REFINEMENT= 2

