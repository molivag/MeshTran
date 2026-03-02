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
      integer               :: ellipsForSite = 0
      integer               :: n_ellipses_site = 0
      integer               :: n_iterative_refi = 0
      real(dp)              :: rotation = 0.0_dp
      real(dp), allocatable :: maxSiteEdge(:)
      real(dp), allocatable :: lenEllipseSite(:)
   end type GlobalRefinement

   type :: ParamRefinement
      integer :: Nsph = 0
      integer :: Nelipses = 0
      real(dp) :: rotation
      real(dp) :: minEDGE
      real(dp) :: minRAD
      real(dp) :: maxRAD
      real(dp) :: paddingRefi
      real(dp) :: sizeBoundary
      real(dp) :: coreResol
      real(dp) :: growthFactor
   end type ParamRefinement

   type :: ModelRegion
      integer :: Nregions = 0
      integer :: NparamEsfer = 0
      integer, allocatable  :: ID(:)
      integer, allocatable  :: repeatPartition(:)
      integer, allocatable  :: isRHOfix(:)
      real(dp), allocatable :: rho(:)
      real(dp), allocatable :: coord(:, :)
      real(dp), allocatable :: radiusForEsfer(:)
      real(dp), allocatable :: edgesForEsfer(:)
   end type ModelRegion

end module mesh_entities
