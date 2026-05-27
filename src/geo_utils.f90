module geo_utils
  
  use mesh_config, only: dp
  implicit none
  
  
  contains

!=========================================================
!=======
!=========================================================trendline
   subroutine lat_long_to_UTM_km(lat, lon, n_points, east_km, north_km)

      implicit none
      integer, intent(in)     :: n_points
      real(dp), intent(in)    :: lat(n_points), lon(n_points)             ! grados decimales
      ! real(dp), intent(out)   :: east_km(n_points), north_km(n_points)    ! ya en kilometros
      real(dp), allocatable, intent(out)   :: east_km(:), north_km(:)    ! ya en kilometros

      integer:: i, zone
      real(dp)               :: x(n_points), y(n_points)                 ! metros
      real(dp):: lon_mean

      ! --- constantes---
      real(dp), parameter:: semieje = 6378137.0d0
      real(dp), parameter:: f = 1.0d0/298.257223563d0
      real(dp), parameter:: k0 = 0.9996d0
      real(dp), parameter:: pi = 3.141592653589793d0

      real(dp):: e2, ep2
      real(dp):: latr(n_points), lonr(n_points), lon0r
      real(dp):: N(n_points), T(n_points), C(n_points), A(n_points), M(n_points)
      real(dp):: lon0

      ! --- elipsoide---
      e2 = f*(2.0d0 - f)
      ep2 = e2/(1.0d0 - e2)

      ! --- zona UTM---
      lon_mean = sum(lon)/dble(n_points)
      zone = int((lon_mean + 180.0d0)/6.0d0) + 1
      lon0 = (zone - 1)*6.0d0 - 180.0d0 + 3.0d0

      do i = 1, n_points
         ! --- radianes---
         latr(i) = lat(i)*pi/180.0d0
         lonr(i) = lon(i)*pi/180.0d0
         lon0r = lon0*pi/180.0d0

         ! --- términos auxiliares---
         N(i) = semieje/sqrt(1.0d0 - e2*sin(latr(i))**2)
         T(i) = tan(latr(i))**2
         C(i) = ep2*cos(latr(i))**2
         A(i) = cos(latr(i))*(lonr(i) - lon0r)

         ! --- arco meridional---
         M(i) = semieje*((1.0d0 - e2/4.0d0 - 3.0d0*e2**2/64.0d0 - 5.0d0*e2**3/256.0d0)*latr(i) &
                         - (3.0d0*e2/8.0d0 + 3.0d0*e2**2/32.0d0 + 45.0d0*e2**3/1024.0d0)*sin(2.0d0*latr(i)) &
                         + (15.0d0*e2**2/256.0d0 + 45.0d0*e2**3/1024.0d0)*sin(4.0d0*latr(i)) &
                         - (35.0d0*e2**3/3072.0d0)*sin(6.0d0*latr(i)))

         ! --- coordenadas UTM---
         x(i) = k0*N(i)*(A(i) + (1.0d0 - T(i) + C(i))*A(i)**3/6.0d0 &
                         + (5.0d0 - 18.0d0*T(i) + T(i)**2 + 72.0d0*C(i) - 58.0d0*ep2)*A(i)**5/120.0d0) + 500000.0d0

         y(i) = k0*(M(i) + N(i)*tan(latr(i))*(A(i)**2/2.0d0 &
                                              + (5.0d0 - T(i) + 9.0d0*C(i) + 4.0d0*C(i)**2)*A(i)**4/24.0d0 &
                                              + (61.0d0 - 58.0d0*T(i) + T(i)**2 + 600.0d0*C(i) - 330.0d0*ep2)*A(i)**6/720.0d0))

         ! hemisferio sur
         !to avoind negative coordinates in northing
         if (lat(i) < 0.0d0) y(i) = y(i) + 10000000.0d0
      end do

      !Covierto UTM en metros a UTM en kilometros
      ! x = x/1000.0d0
      ! y = y/1000.0d0
      ALLOCATE(east_km(n_points), north_km(n_points))

      north_km = y/1000.0d0    ! Y = Norte
      east_km = x/1000.0d0    ! X = Este

      !Lo ideal es que esta conversion se haga fuera e la rutina, al igual en read_dem o read_EDIs
      ! x = north_km   ! X = Norte
      ! y = east_km    ! Y = Este

   end subroutine lat_long_to_UTM_km
!=========================================================
!=======
!=========================================================
   subroutine debug_coastline(east, north, n)
      !
      ! This routine receive a east and north variables in km UTM to plot the DEM and Coast Line together
      !

      implicit none

      integer, intent(in) :: n
      real(8), intent(in) :: east(n)
      real(8), intent(in) :: north(n)

      integer :: i, stat
      character(len=20) :: answer, filename
      character(len=256) :: cmd
      character(len=1) :: q
      q = "'"
      !----------------------------------------
      ! write temporary coastline
      !----------------------------------------

      filename = "coast_debug.xyz"
      open(unit=99, file='./preprocessing/geometry/'//trim(filename), status='replace')

      do i = 1, n
         write(99,*) east(i), north(i)
      end do

      close(99)

      !----------------------------------------
      ! launch gnuplot
      !----------------------------------------

      
      write(cmd,'(A)') &
      'gnuplot -persist -e " ' // &
      'set size ratio -1; ' // &
      'plot ' // &
      q // './preprocessing/DEM/UTM_gebco.xyz' // q // &
      ' u 1:2:3 w dots palette notitle, ' // &
      q // './preprocessing/geometry/coast_debug.xyz' // q // &
      ' u 1:2 w lp lw 1 lc rgb ' // q // 'black' // q // ' title ' // q // 'Closed Coast' // q // &
      '"'

call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
      ! write (cmd, '(A)') gnuplot -persist -e 'plot "coast_debug.dat" w lp' '
      ! write(cmd,'(A)') 'gnuplot -persist -e "plot ''coast_debug.xy'' w lp"'
      ! write(cmd,'(A)') 'gnuplot -persist -e "plot './preprocessing/geometry/'//trim(filename)" w lp"'

! write(cmd,'(A)') 'gnuplot -persist -e "plot ' // q //'./preprocessing/geometry/' // trim(filename) // q // ' w lp"'
         ! call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
         if (stat /= 0) error stop 'ERROR executing gnuplot'

      ! call execute_command_line( &
      ! "gnuplot -persist -e " &
      ! // 'set title ''MeshTran CoastLine Debug''; ' &
      ! // 'set size ratio -1; ' &
      ! // 'plot ''coast_debug.xy'' w lp pt 7 lw 2' &
      ! // '" '&
      ! )
      !
      !----------------------------------------
      ! ask user validation
      !----------------------------------------

      print *
      print *, 'Is coastline correct? (yes/no)'
      read(*,*) answer

      if (trim(adjustl(answer)) /= 'yes') then

         print *
         print *, 'ERROR: coastline validation failed'
         stop

      end if

   end subroutine
end module geo_utils

