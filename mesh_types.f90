module mesh_types
   implicit none
   integer, parameter :: dp = kind(1.0d0)

   type :: MeshSettings

      ! ----- Archivos -----
      character(len=156) :: dem_file
      character(len=156) :: dem_units
      character(len=156) :: outdir
      character(len=156) :: topography_file
      character(len=156) :: bathymetry_file
      character(len=156) :: coastLine_file

      ! ----- Dominio -----
      logical :: has_sea
      real(dp) :: sea_level
      real(dp) :: xminDOM, xmaxDOM
      real(dp) :: yminDOM, ymaxDOM
      real(dp) :: zminDOM, zmaxDOM
      real(dp) :: pad_x, pad_y
      real(dp) :: rotation

      ! ----- Refinamiento global -----
      integer :: Nsph
      real(dp), allocatable :: radius(:)
      real(dp), allocatable :: edges(:)

      ! ----- Regiones -----
      integer :: Nregions
      integer, allocatable :: regionsID(:)
      real(dp), allocatable :: coord_regions(:, :)
      real(dp), allocatable :: rho_regions(:)

   end type MeshSettings

end module mesh_types
