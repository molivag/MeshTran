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
      real(dp), allocatable :: maxSiteEdge(:)
      real(dp), allocatable :: lenEllipseSite(:)
      integer               :: n_ellipses_site = 0
      real(dp)              :: rotation = 0.0_dp
      integer               :: ellipsForSite = 0
   end type GlobalRefinement

   type :: ParamRefinement
      integer :: Nsph = 0
      real(dp), allocatable :: edges(:)
      real(dp), allocatable :: radius(:)
   end type ParamRefinement

   type :: ModelRegion
      integer :: Nregions = 0
      integer :: NparamEsfer = 0
      integer, allocatable  :: ID(:)
      real(dp), allocatable :: rho(:)
      integer, allocatable  :: repeatPartition(:)
      integer, allocatable  :: isRHOfix(:)
      real(dp), allocatable :: coord(:, :)
      
      real(dp), allocatable :: radiusForEsfer(:)
      real(dp), allocatable :: edgesForEsfer(:)
   end type ModelRegion

end module mesh_entities
