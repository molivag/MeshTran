module mesh_config

   implicit none
   integer, parameter :: dp = kind(1.0d0)
   CHARACTER(len=50),PARAMETER :: outdir ='preprocessing/geometry/'

   type :: MeshSettings

      ! ----- Archivos -----
      character(len=150) :: dem_file
      character(len=50) :: dem_units
      character(len=20) :: mesh_nature
      ! character(len=156) :: outdir
      character(len=100) :: topography_file
      character(len=100) :: bathymetry_file
      character(len=100) :: coastLine_file

      ! ----- Dominio -----
      logical :: has_sea
      real(dp) :: sea_level
      real(dp) :: xminDOM, xmaxDOM
      real(dp) :: yminDOM, ymaxDOM
      real(dp) :: zminDOM, zmaxDOM
      real(dp) :: pad_x, pad_y
      real(dp) :: rotation


   end type MeshSettings

end module mesh_config
