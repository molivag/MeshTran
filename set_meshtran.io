*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
*                                                                 *
* MeshTran configuration file                                     *
*                                                                 *
*                                        MAOG    Bcn, Mzo. 2026   *
*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

MESH_NATURE = native
THREADS = 18

DEM_FILE = MET_50m_mosaic_valles_Lcosta_cropped.xyz
__DEM_FILE = TOTPO_KAL_BBMT_FLAT_475m.xyz
DEM_LatLon = no 
DEM_UNITS = meters                    !meters or kilometers
assumming elevation in meters

TOPO_FILE = topography_for.dat 
BATHY_FILE = bathymetry_for.dat
COSLI_FILE = coast_line.dat

XMIN_DOM = -32.0
XMAX_DOM =  20.0
YMIN_DOM = -25.0
YMAX_DOM =  30.0
ZMIN_DOM = -50.0
ZMAX_DOM =  50.0

PAD_X = 5.0
PAD_Y = 5.0

HAS_SEA = yes
SEA_LEVEL = 0.0

## =>  Observation data
IN_diag_ERR = 0.10
OFF_diag_ERR = 0.05
SUBSAMPLING = 0

## => If subsampling is one, select the frequency range
FMIN = 00
FMAX = 00



## =>  Ellipsoids Surface Mesh
NUM_ELIPSES = 5
SITEpadding = 4.0        ! Margen extra sobre el radio de los sites (km)
FARelemSIZE = 15.0        ! Tamaño máximo de triángulo en el borde (km) OBJparamRefi%sizeBoundary
SURF_RESOLUTION = 5.0     ! len en el centro (km) OBJparamRefi%coreResol
GROWTH = 1.0              ! Cuánto aumenta 'len' en cada elipse hacia afuera

## => Meshing for sites at the surface
NUM_ESFERAS = 10
ROTATION = 0.0
minRADIO = 0.01        ! Si este num se reduce, mejor formart de esfera y menos elementos pero no se si menos resol, parece que no
maxRADIO = 10.0
minEDGES = 0.01        ! Esto controla fuertemente la densidad de elementos en los sites



## => Refinement for 3D mesh and sites
REFI_ELLIPSES = 5
HIGH_RESOL_LAYER = 1.6    !thikness of high resolution from surface = HRL*TargDep
V_PADDING = 2.3           !vertical extension beyond the refinement


NEAR_FIELD_RESOL = 4      !min elem size at sites (core_resolution)
FAR_FIELD_RESOL =  10     !max elem size 

ITER_TET_REFI = 7

!esto controla radio a en makeMTR
H_PADDING = 30             !horizontal extension (in deep) beyond last site, a1 = max_distance_site x H_PADDING 
TARGET_DEPTH = 15






!esto seria para malla de parametros en makeMTR no influye
## => Para el calculo del skin depth en el mallado de parametros
BKGRD_RHO = 100
F_MIN_HZ = 0.001

FRAC_SKIN_DEPTH = 5.55    !No quiero que el tamaño de elemento (o parámetro) sea mayor que el 25% de la skin depth


## => Regions Atributes in the model
ELLIPSES = 9 
REGIONS = 3

LOCATION = 0.0 0.0 -20.0    -15.0 15.0  2.0    0.0 0.0 45.0       	  !xyz xyz  xyz
ID_REGION =     10               20                 30	    		         !reg1 reg2
RHO_REGIONS =   1.0e9          0.01                1.0e2		             !rho1 rho2
REP_PARTITION = -1              -1                   9           
FIX_RESISTIVITY = 1              1                   0           

## => Refinement for mesh parameter based on spheres
PARAM_ESFER = 2
PARAM_RADIOS = 3.0 5.0 
PARAM_EDGES = 2.0 3.0 


