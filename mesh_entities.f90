module mesh_entities
   use mesh_config
   implicit none
   ! integer, parameter :: dp = kind(1.0d0)

   type :: MeshIO
      character(len=256) :: outdir
   end type MeshIO

   type :: SiteSet
      character(len=256),  allocatable :: id(:)
      real(dp),            allocatable :: x(:), y(:), z(:)
      real(dp),            allocatable :: z_fix(:)
      integer                          :: n = 0
   end type SiteSet

   type :: GlobalRefinement
      real(dp), allocatable :: maxSiteEdge(:)
      real(dp), allocatable :: lenEllipseSite(:)
      real(dp)              :: rotation = 0.0_dp
      integer               :: n_ellipses = 0
      integer               :: ellipsForSite = 0
   end type GlobalRefinement

   type :: ParamRefinement
      integer :: n = 0
      real(dp), allocatable :: edges(:)
      real(dp), allocatable :: radius(:)
   end type ParamRefinement

   type :: RegionModel
      integer :: n = 0
      integer, allocatable :: id(:)
      real(dp), allocatable :: rho(:)
      integer, allocatable :: repeatPartition(:)
      integer, allocatable :: isRHOfix(:)
      real(dp), allocatable :: coord(:, :)
   end type RegionModel

end module mesh_entities
