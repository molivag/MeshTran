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
!=========================================================
!=======
!=========================================================
subroutine ray_casting(coast_x, coast_y, n_coast, dem_x, dem_y, dem_z, n_dem, is_land)

   implicit none

   integer,  intent(in) :: n_coast
   real(dp), intent(in) :: coast_x(n_coast), coast_y(n_coast)
   integer,          intent(in)  :: n_dem
   real(dp),         intent(in)  :: dem_x(n_dem)   ! Norte (km)
   real(dp),         intent(in)  :: dem_y(n_dem)   ! Este  (km)
   real(dp),         intent(in)  :: dem_z(n_dem)   ! Elevacion (km)
   logical,          intent(out) :: is_land(n_dem)

   integer  :: i, j, cnt
   real(dp) :: px, py, xi, yi, xj, yj
   real(dp), parameter :: eps = 1.0d-10

   print*, " "
   print*, " ---> Running Ray Casting algorithm to clasify DEM"
   print*, " "

   do i = 1, n_dem
      px  = dem_x(i)
      py  = dem_y(i)
      cnt = 0

      do j = 1, n_coast - 1
         xi = coast_x(j)
         yi = coast_y(j)
         xj = coast_x(j+1)
         yj = coast_y(j+1)

         if (ray_intersects_seg(px, py, xi, yi, xj, yj, eps)) &
            cnt = cnt + 1
      end do

      is_land(i) = (mod(cnt, 2) /= 0)
   end do

contains

   function ray_intersects_seg(px, py, ax, ay, bx, by, eps) result(intersect)
      implicit none
      real(dp), intent(in) :: px, py, ax, ay, bx, by, eps
      logical :: intersect

      real(dp) :: ppy, seg_ax, seg_ay, seg_bx, seg_by
      real(dp) :: m_red, m_blue

      ! A debe ser el punto con menor x (Norte)
      if (ax > bx) then
         seg_ax = bx;  seg_ay = by
         seg_bx = ax;  seg_by = ay
      else
         seg_ax = ax;  seg_ay = ay
         seg_bx = bx;  seg_by = by
      end if

      ! Evitar rayo sobre vertice
      ppy = px
      if (ppy == seg_ax .or. ppy == seg_bx) ppy = ppy + eps

      intersect = .false.

      if (ppy < seg_ax .or. ppy > seg_bx) return
      if (py  > max(seg_ay, seg_by))       return

      if (py < min(seg_ay, seg_by)) then
         intersect = .true.
         return
      end if

      if (abs(seg_bx - seg_ax) > tiny(1.0_dp)) then
         m_red = (seg_by - seg_ay) / (seg_bx - seg_ax)
      else
         m_red = huge(1.0_dp)
      end if

      if (abs(px - seg_ax) > tiny(1.0_dp)) then
         m_blue = (ppy - seg_ay) / (py - seg_ay)
      else
         m_blue = huge(1.0_dp)
      end if

      intersect = (m_blue >= m_red)

   end function ray_intersects_seg

end subroutine ray_casting
!=========================================================
!=======
!=========================================================
subroutine write_edi_site_tmp(tmpfile, nf, freq, zr, zi, zerr, mask)
   !
   ! Escribe archivo temporal con datos del sitio para gnuplot
   ! Columnas: freq Re_Zxy Im_Zxy SE_xy Re_Zyx Im_Zyx SE_yx
   !                Re_Zxx Im_Zxx SE_xx Re_Zyy Im_Zyy SE_yy
   !
   implicit none
   character(len=*), intent(in) :: tmpfile
   integer,          intent(in) :: nf
   real(dp),         intent(in) :: freq(nf)
   real(dp),         intent(in) :: zr(nf,4), zi(nf,4), zerr(nf,4)
   logical,          intent(in) :: mask(nf)

   integer            :: j, unit_out, ios

   open(newunit=unit_out, file=trim(tmpfile), &
        status='replace', action='write', iostat=ios)
   if (ios /= 0) then
      print*, 'ERROR: No se pudo abrir archivo temporal ', trim(tmpfile)
      error stop
   end if

   write(unit_out,'(A)') '# col1:freq  col2:Re_Zxy  col3:Im_Zxy  col4:SE_xy  &
                          &col5:Re_Zyx  col6:Im_Zyx  col7:SE_yx   &
                          &col8:Re_Zxx  col9:Im_Zxx  col10:SE_xx  &
                          &col11:Re_Zyy col12:Im_Zyy col13:SE_yy'

   do j = 1, nf
      if (.not. mask(j)) cycle
      write(unit_out, '(13(1x, ES14.6))') &
         freq(j),                       &  ! col 1
         zr(j,2), zi(j,2), zerr(j,2),  &  ! col 2,3,4   XY
         zr(j,3), zi(j,3), zerr(j,3),  &  ! col 5,6,7   YX
         zr(j,1), zi(j,1), zerr(j,1),  &  ! col 8,9,10  XX
         zr(j,4), zi(j,4), zerr(j,4)      ! col 11,12,13 YY
   end do

   close(unit_out)

end subroutine write_edi_site_tmp
!=========================================================
!=======
!=========================================================
! subroutine plot_edi_site(site_name)
!    implicit none
!    character(len=*), intent(in) :: site_name
!
!    integer             :: stat, iu
!    character(len=512)  :: tmpfile, pdffile, gpfile
!
!    tmpfile = '/tmp/'//trim(site_name)//'.dat'
!    pdffile = 'preprocessing/edi_sites/'//trim(site_name)//'.pdf'
!    gpfile  = '/tmp/'//trim(site_name)//'.gp'
!
!    call execute_command_line('mkdir -p preprocessing/edi_sites', &
!                               wait=.true., exitstat=stat)
!
!    !----------------------------------------------------------
!    ! Escribir script gnuplot a archivo temporal
!    !----------------------------------------------------------
!    open(newunit=iu, file=trim(gpfile), status='replace')
!
!    write(iu,'(A)') 'set terminal pdfcairo enhanced color size 12,8 font "Helvetica,10"'
!    write(iu,'(A)') 'set output "'//trim(pdffile)//'"'
!    write(iu,'(A)') 'mu0 = 1.2566370614e-6'
!    write(iu,'(A)') 'rho(f,Re,Im) = (1.0/(2.0*pi*f*mu0)) * (Re**2 + Im**2)'
!    write(iu,'(A)') 'phase(Re,Im) = atan2(Im,Re) * 180.0/pi'
!       !para las error bars
! write(iu,'(A)') 'drho(f,Re,Im,SE) = (1.0/(pi*f*mu0)) * sqrt(Re**2+Im**2) * SE'
! write(iu,'(A)') 'dphase(Re,Im,SE) = (SE / sqrt(Re**2+Im**2)) * 180.0/pi'
!
!    write(iu,'(A)') 'set multiplot layout 2,2 title "'//trim(site_name)//'"'
!    write(iu,'(A)') 'set logscale x'
!    write(iu,'(A)') 'set yrange [1e3:1e-3]'
!    write(iu,'(A)') 'set xlabel "Frequency (Hz)"'
!    write(iu,'(A)') 'set format x "10^{%T}"'
!
!    ! Panel (1,1): rho_a off-diagonal
!    write(iu,'(A)') 'set logscale y'
!    write(iu,'(A)') 'set yrange [0.5:3500]'
!    write(iu,'(A)') 'set ylabel "$rho_a$ (Ohm·m)"'
!    write(iu,'(A)') 'set title "Off-diagonal: rho_a XY & YX"'
!
!    ! write(iu,'(A)') 'plot "'//trim(tmpfile)//'" u 1:(rho($1,$2,$3)) w lp pt 4 lc rgb "purple" title "rho XY", \'
!    ! write(iu,'(A)') '     "'//trim(tmpfile)//'" u 1:(rho($1,$5,$6)) w lp pt 6 lc rgb "red" title "rho YX"'
! write(iu,'(A)') 'plot "'//trim(tmpfile)//'" u 1:(rho($1,$2,$3)):(drho($1,$2,$3,$4)) w yerrorbars pt 4 lc rgb "purple" title "rho XY", \'
! write(iu,'(A)') '     "'//trim(tmpfile)//'" u 1:(rho($1,$5,$6)):(drho($1,$5,$6,$7)) w yerrorbars pt 6 lc rgb "red" title "rho YX"'
!
!    ! Panel (1,2): rho_a diagonal
!    write(iu,'(A)') 'set title "Diagonal: rho_a XX & YY"'
!    ! write(iu,'(A)') 'plot "'//trim(tmpfile)//'" u 1:(rho($1,$8,$9)) w lp pt 4 lc rgb "purple" title "rho XX", \'
!    ! write(iu,'(A)') '     "'//trim(tmpfile)//'" u 1:(rho($1,$11,$12)) w lp pt 6 lc rgb "red" title "rho YY"'
! write(iu,'(A)') 'plot "'//trim(tmpfile)//'" u 1:(rho($1,$8,$9)):(drho($1,$8,$9,$10)) w yerrorbars pt 4 lc rgb "purple" title "rho XX", \'
! write(iu,'(A)') '     "'//trim(tmpfile)//'" u 1:(rho($1,$11,$12)):(drho($1,$11,$12,$13)) w yerrorbars pt 6 lc rgb "red" title "rho YY"'
!
!    ! Panel (2,1): fase off-diagonal
!    write(iu,'(A)') 'unset logscale y'
!    write(iu,'(A)') 'set yrange [-180:180]'
!    write(iu,'(A)') 'set ytics 45'
!    write(iu,'(A)') 'set ylabel "Phase (deg)"'
!    write(iu,'(A)') 'set title "Off-diagonal: Phase XY & YX"'
!    ! write(iu,'(A)') 'plot "'//trim(tmpfile)//'" u 1:(phase($2,$3)) w lp pt 4 lc rgb "purple" title "Phase XY", \'
!    ! write(iu,'(A)') '     "'//trim(tmpfile)//'" u 1:(phase($5,$6)) w lp pt 6 lc rgb "red" title "Phase YX"'
! write(iu,'(A)') 'plot "'//trim(tmpfile)//'" u 1:(phase($2,$3)):(dphase($2,$3,$4)) w yerrorbars pt 4 lc rgb "purple" title "Phase XY", \'
! write(iu,'(A)') '     "'//trim(tmpfile)//'" u 1:(phase($5,$6)):(dphase($5,$6,$7)) w yerrorbars pt 6 lc rgb "red" title "Phase YX"'
!
!    ! Panel (2,2): fase diagonal
!    write(iu,'(A)') 'set title "Diagonal: Phase XX & YY"'
!    ! write(iu,'(A)') 'plot "'//trim(tmpfile)//'" u 1:(phase($8,$9)) w lp pt 4 lc rgb "purple" title "Phase XX", \'
!    ! write(iu,'(A)') '     "'//trim(tmpfile)//'" u 1:(phase($11,$12)) w lp pt 6 lc rgb "red" title "Phase YY"'
! write(iu,'(A)') 'plot "'//trim(tmpfile)//'" u 1:(phase($8,$9)):(dphase($8,$9,$10)) w yerrorbars pt 4 lc rgb "purple" title "Phase XX", \'
! write(iu,'(A)') '     "'//trim(tmpfile)//'" u 1:(phase($11,$12)):(dphase($11,$12,$13)) w yerrorbars pt 6 lc rgb "red" title "Phase YY"'
!
!    write(iu,'(A)') 'unset multiplot'
!    close(iu)
!
!    !----------------------------------------------------------
!    ! Ejecutar gnuplot
!    !----------------------------------------------------------
!    call execute_command_line('gnuplot '//trim(gpfile), &
!                               wait=.true., exitstat=stat)
!    if (stat /= 0) then
!       print*, 'ERROR: gnuplot fallo para sitio ', trim(site_name)
!       error stop
!    end if
!
!    !----------------------------------------------------------
!    ! Borrar temporales
!    !----------------------------------------------------------
!    call execute_command_line('rm -f '//trim(tmpfile)//' '//trim(gpfile), &
!                               wait=.true., exitstat=stat)
!
!    ! print*, '✅ PDF generado: ', trim(pdffile)
!
! end subroutine plot_edi_site
subroutine plot_edi_site(site_name)
   implicit none
   character(len=*), intent(in) :: site_name

   integer             :: stat, iu
   character(len=512)  :: tmpfile, pdffile, gpfile

   tmpfile = '/tmp/'//trim(site_name)//'.dat'
   pdffile = 'preprocessing/edi_sites/'//trim(site_name)//'.pdf'
   gpfile  = '/tmp/'//trim(site_name)//'.gp'

   call execute_command_line('mkdir -p preprocessing/edi_sites', &
                              wait=.true., exitstat=stat)

   open(newunit=iu, file=trim(gpfile), status='replace')

   ! Terminal
   write(iu,'(A)') 'set terminal pdfcairo enhanced color size 12,8 font "Helvetica,10"'
   write(iu,'(A)') 'set output "'//trim(pdffile)//'"'

   ! Funciones
   write(iu,'(A)') 'mu0 = 1.2566370614e-6'
   write(iu,'(A)') 'rho(f,Re,Im) = (1.0/(2.0*pi*f*mu0)) * (Re**2 + Im**2)'
   write(iu,'(A)') 'phase(Re,Im) = atan2(-Im,Re) * 180.0/pi'
   write(iu,'(A)') 'drho(f,Re,Im,SE) = (1.0/(pi*f*mu0)) * sqrt(Re**2+Im**2) * SE'
   write(iu,'(A)') 'dphase(Re,Im,SE) = (SE / sqrt(Re**2+Im**2)) * 180.0/pi'

   ! Multiplot
   write(iu,'(A)') 'set multiplot layout 2,2 title "'//trim(site_name)//'"'

   ! Eje X comun
   write(iu,'(A)') 'set logscale x'
   ! write(iu,'(A)') 'set xrange [1e-3:1e3]'
   write(iu,'(A)') 'set xrange [1e4:1e-4]'      !de altas a bajas freq
   write(iu,'(A)') 'set mxtics 10'
   write(iu,'(A)') 'set xlabel "Frequency (Hz)"'
   write(iu,'(A)') 'set format x "10^{%T}"'

   ! Grid punteado
   write(iu,'(A)') 'set grid xtics ytics lt 0 lw 0.5 lc rgb "gray70"'

   ! ============================================================
   ! Panel (1,1): rho_a off-diagonal
   ! ============================================================
   write(iu,'(A)') 'set logscale y'
   write(iu,'(A)') 'set yrange [0.5:10000]'
   write(iu,'(A)') 'set ylabel "rho_a (Ohm·m)"'
   write(iu,'(A)') 'set title "Off-diagonal: rho_a XY & YX"'
   write(iu,'(A)') 'plot "'//trim(tmpfile)// &
      '" u 1:(rho($1,$2,$3)):(drho($1,$2,$3,$4)) w yerrorbars pt 5 ps 0.7 lc rgb "dark-turquoise" title "rho XY", \'
   write(iu,'(A)') '     "'//trim(tmpfile)// &
      '" u 1:(rho($1,$5,$6)):(drho($1,$5,$6,$7)) w yerrorbars pt 7 ps 0.7 lc rgb "dark-orange" title "rho YX"'

   ! ============================================================
   ! Panel (1,2): rho_a diagonal
   ! ============================================================
   write(iu,'(A)') 'set title "Diagonal: rho_a XX & YY"'
   write(iu,'(A)') 'plot "'//trim(tmpfile)// &
      '" u 1:(rho($1,$8,$9)):(drho($1,$8,$9,$10)) w yerrorbars pt 5 ps 0.7 lc rgb "dark-turquoise" title "rho XX", \'
   write(iu,'(A)') '     "'//trim(tmpfile)// &
      '" u 1:(rho($1,$11,$12)):(drho($1,$11,$12,$13)) w yerrorbars pt 7 ps 0.7 lc rgb "dark-orange" title "rho YY"'

   ! ============================================================
   ! Panel (2,1): fase off-diagonal
   ! ============================================================
   write(iu,'(A)') 'unset logscale y'
   write(iu,'(A)') 'set yrange [-200:200]'
   write(iu,'(A)') 'set ytics 45'
      write(iu,'(A)') 'set grid xtics noytics lt 0 lw 0.5 lc rgb "gray70"'

   write(iu,'(A)') 'set ylabel "Phase (deg)"'
   write(iu,'(A)') 'set title "Off-diagonal: Phase XY & YX"'
   ! Linea de referencia en fase = 0
   write(iu,'(A)') 'set arrow 1 from 1e-3,0 to 1e3,0 nohead lt 0 lw 0.8 lc rgb "gray50"'
   write(iu,'(A)') 'plot "'//trim(tmpfile)// &
      '" u 1:(phase($2,$3)):(dphase($2,$3,$4)) w yerrorbars pt 5 ps 0.7 lc rgb "dark-turquoise" title "Phase XY", \'
   write(iu,'(A)') '     "'//trim(tmpfile)// &
      '" u 1:(phase($5,$6)):(dphase($5,$6,$7)) w yerrorbars pt 7 ps 0.7 lc rgb "dark-orange" title "Phase YX"'

   ! ============================================================
   ! Panel (2,2): fase diagonal
   ! ============================================================
   write(iu,'(A)') 'set title "Diagonal: Phase XX & YY"'
   write(iu,'(A)') 'plot "'//trim(tmpfile)// &
      '" u 1:(phase($8,$9)):(dphase($8,$9,$10)) w yerrorbars pt 5 ps 0.7 lc rgb "dark-turquoise" title "Phase XX", \'
   write(iu,'(A)') '     "'//trim(tmpfile)// &
      '" u 1:(phase($11,$12)):(dphase($11,$12,$13)) w yerrorbars pt 7 ps 0.7 lc rgb "dark-orange" title "Phase YY"'

   ! Limpiar arrow
   write(iu,'(A)') 'unset arrow 1'
   write(iu,'(A)') 'unset multiplot'
   close(iu)

   ! Ejecutar gnuplot
   call execute_command_line('gnuplot '//trim(gpfile), &
                              wait=.true., exitstat=stat)
   if (stat /= 0) then
      print*, 'ERROR: gnuplot fallo para sitio ', trim(site_name)
      error stop
   end if

   ! Borrar temporales
   call execute_command_line('rm -f '//trim(tmpfile)//' '//trim(gpfile), &
                              wait=.true., exitstat=stat)

   ! print*, '✅ PDF generado: ', trim(pdffile)

end subroutine plot_edi_site
end module geo_utils

