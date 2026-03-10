module mesh_entities
   use mesh_config, only: dp

   implicit none
   ! integer, parameter :: dp = kind(1.0d0)

   type :: SiteSet
      character(len=256),  allocatable :: id(:)
      real(dp),            allocatable :: x(:), y(:), z(:)
      real(dp),            allocatable :: z_fix(:)
      integer                          :: n = 0
   end type SiteSet

   type :: GlobalRefinement
      ! integer               :: ellipsForSite = 0 !
      integer     :: levels = 0                !REFI_ELLIPSES
      integer     :: n_ellipses_site = 0       !SITE_ELLIPSES
      integer     :: n_iterative_refi = 0
      real(dp)    :: target_depth
      real(dp)    :: depth_min_factor
      real(dp)    :: background_rho
      real(dp)    :: core_resolution
      real(dp)    :: farfield_resolution
      real(dp)    :: horizontal_padding
      real(dp)    :: vertical_padding_factor
      real(dp)    :: frac_skin_depth
      real(dp)    :: fmin_hz

      ! Internos (no TOML por ahora)
      real(dp)    :: a_ratio           = 1.3_dp       !refinamiento más alargado horizontalmente
      real(dp)    :: len_growth        = 1.45_dp       !factor de crecimiento del tamaño de elemento
      real(dp)    :: vol_len_factor    = 1.6_dp       !volumen máximo permitido del tetraedro
      real(dp)    :: air_c_over_a_core = 0.30_dp      !casquete elipsoidal encima del core
      real(dp)    :: fh_core_vol       = 0.30_dp      !factor de refinamiento volumétrico dentro del core.
   end type GlobalRefinement

   type :: ParamRefinement
      integer :: Nsph = 0           !site
      integer :: Nelipses = 0
      real(dp) :: rotation
      real(dp) :: minEDGE           !site
      real(dp) :: minRAD            !site
      real(dp) :: maxRAD            !site
      real(dp) :: paddingRefi
      real(dp) :: sizeBoundary
      real(dp) :: coreResol
      real(dp) :: growthFactor
   end type ParamRefinement

   type :: ModelRegion
      integer :: Nregions = 0
      integer :: NparamEsfer = 0
      integer :: NpEllipses = 0
      integer, allocatable  :: ID(:)
      integer, allocatable  :: repeatPartition(:)
      integer, allocatable  :: isRHOfix(:)
      real(dp), allocatable :: rho(:)
      real(dp), allocatable :: coord(:, :)
      real(dp), allocatable :: radiusForEsfer(:)
      real(dp), allocatable :: edgesForEsfer(:)
   end type ModelRegion

end module mesh_entities
