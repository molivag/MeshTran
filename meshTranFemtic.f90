program femtic_mesh_driver

   use mesh_config
   use mesh_entities
   implicit none
   ! integer, parameter:: dp = kind(1.0d0)

   TYPE(MeshSettings)::settings
   type(MeshIO) :: io
   type(SiteSet) :: sites
   type(GlobalRefinement) :: globRefi
   type(ParamRefinement) :: paramRefi
   type(RegionModel) :: regions

   ! Control parameters
   character(len=256):: dem_file, outdir, dem_units, topography_file, bathymetry_file, coast_line_file
   character(len=256), allocatable:: edi_files(:), edi_id(:)
   logical:: has_sea
   real(dp):: sea_level
   real(dp):: xminDOM, xmaxDOM, yminDOM, ymaxDOM, zminDOM, zmaxDOM, x0, y0
   real(dp), allocatable:: edi_lat(:), edi_lon(:), edi_elev(:)
   real(dp):: pad_x, pad_y, rotation

   ! DEM data
   integer:: n_dem, zone, n_edi_files, Nsph, n_sites, n_elipses, ellipsForSite, NparamEsfer, Nregions, refi_tetgen
   real(dp), allocatable:: site_x(:), site_y(:), site_z(:), dem_x(:), dem_y(:), dem_z(:), fix_elev_site(:)
   real(dp), allocatable:: site_x_km(:), site_y_km(:), dem_x_km(:), dem_y_km(:), radius(:), edges(:)
   real(dp), allocatable:: maxSiteEdge(:), lenEllipseSite(:), edgesForEsfer(:), radiusForEsfer(:)
   real(dp), allocatable:: coordinate_regions(:, :), rhoOFregions(:)
   integer, allocatable :: repeatPartition(:), isRHOfix(:), IDregions(:)


   

   !---------------------------------------------------
   !     Read input configutation file
   !---------------------------------------------------
   call read_set_femtic("set_meshtran.io", settings)
   !---------------------------------------------------
   !     Read Digital Elevation Model
   !---------------------------------------------------
   call read_dem(dem_file, dem_units, dem_x, dem_y, dem_z, n_dem)
   allocate (dem_x_km(n_dem), dem_y_km(n_dem))

   !---------------------------------------------------
   !     Read EDI files
   !---------------------------------------------------
   call get_edi_file_list(edi_files, n_edi_files)
   n_sites = n_edi_files
   allocate (edi_id(n_edi_files), edi_lat(n_edi_files), edi_lon(n_edi_files), edi_elev(n_edi_files))
   allocate (site_x(n_edi_files), site_y(n_edi_files), site_z(n_edi_files), fix_elev_site(n_sites))
   ALLOCATE (site_x_km(n_edi_files), site_y_km(n_edi_files))
   call read_edi_files(edi_files, n_edi_files, edi_id, edi_lat, edi_lon, edi_elev)

   !---------------------------------------------------
   !     Convert Lat-Long to UTm
   !---------------------------------------------------
   !                                                esto deberia llamarse edi_x y edi_y
   call edi_to_utm(edi_lat, edi_lon, n_edi_files, site_x, site_y, zone)

   call write_dem_sites_utm(edi_id, site_x, site_y, n_edi_files, dem_x, &
                            dem_y, dem_z, n_dem, site_x_km, site_y_km, dem_x_km, dem_y_km)

   dem_x_km = dem_x
   dem_y_km = dem_y
   site_y_km = site_y
   site_x_km = site_x

   !---------------------------------------------------
   !     Compute anchor point (central point) of analysis domain based on sites locations
   !---------------------------------------------------
   call compute_center(site_x_km, site_y_km, n_edi_files, x0, y0)

   !---------------------------------------------------
   !     Recenter
   !---------------------------------------------------
   call recenter_all(n_edi_files, n_dem, x0, y0, site_x_km, site_y_km, dem_x_km, dem_y_km)
   !a parti de aqui coordenadas de malla

   call snap_sites_to_dem(edi_id, site_x_km, site_y_km, n_sites, dem_x, dem_y, dem_z, n_dem, fix_elev_site)

   ! !Check if DEM cover whole analysis domain+padding area
   call check_domain(dem_x_km, dem_y_km, xminDOM, xmaxDOM, yminDOM, ymaxDOM, pad_x, pad_y)

   ! !Write files in mesh coordinates centered at anchor point
   call define_analysis_domain(xminDOM, xmaxDOM, yminDOM, ymaxDOM, zminDOM, zmaxDOM)
   call write_topography(topography_file, dem_x_km, dem_y_km, dem_z, n_dem, sea_level)
   call write_bathymetry(bathymetry_file, dem_x_km, dem_y_km, dem_z, n_dem, sea_level, outdir)
   call write_coast_line(coast_line_file, xminDOM, xmaxDOM, yminDOM, ymaxDOM, outdir)
   call write_observing_sites(site_x_km, site_y_km, n_edi_files, outdir, Nsph, edges, radius)

   call setGlobalMeshRefinement(outdir, 0.0d0, 0.0d0, n_elipses, rotation, site_x_km, site_y_km, &
                                n_edi_files, ellipsForSite, maxSiteEdge, lenEllipseSite)
   ! call setGlobalMeshRefinement(io, sites, globRefi)

   call resistivitty_attribute(outdir, n_sites, site_x_km, site_y_km, fix_elev_site, NparamEsfer, edgesForEsfer, &
                               radiusForEsfer, Nregions, IDregions, rhoOFregions, repeatPartition, isRHOfix)

   ! call run_makeTetraMesh_and_assign_regions(Nregions,coordinate_regions,IDregions)
   ! call run_TETGEN_and_refine_mesh(refi_tetgen)
   call run_TetGen2Femtic(refi_tetgen)

contains
   !=========================================================
   !=======
   !=========================================================
   subroutine read_set_femtic(fname, settings)

      use mesh_config
      implicit none

      character(len=*), intent(in) :: fname
      type(MeshSettings), intent(inout) :: settings

      character(len=256) :: line, key, val
      integer :: iu
      real(dp), allocatable :: temp_vec(:)

      open (newunit=iu, file=fname, status='old')

      do
         read (iu, '(A)', end=100) line
         if (index(line, '=') == 0) cycle

         key = adjustl(trim(line(:index(line, '=') - 1)))
         val = adjustl(trim(line(index(line, '=') + 1:)))

         select case (trim(key))

         case ('DEM_FILE')
            settings%dem_file = trim(val)

         case ('DEM_UNITS')
            settings%dem_units = trim(val)

         case ('OUTDIR')
            settings%outdir = trim(val)

         case ('SEA_LEVEL')
            read (val, *) settings%sea_level

         case ('XMIN_DOM')
            read (val, *) settings%xminDOM

         case ('XMAX_DOM')
            read (val, *) settings%xmaxDOM

         case ('YMIN_DOM')
            read (val, *) settings%yminDOM

         case ('YMAX_DOM')
            read (val, *) settings%ymaxDOM

         case ('ZMIN_DOM')
            read (val, *) settings%zminDOM

         case ('ZMAX_DOM')
            read (val, *) settings%zmaxDOM

         case ('PAD_X')
            read (val, *) settings%pad_x

         case ('PAD_Y')
            read (val, *) settings%pad_y

         case ('ROTATION')
            read (val, *) settings%rotation

         case ('ESFERAS')
            read (val, *) settings%Nsph

            if (allocated(settings%radius)) deallocate (settings%radius)
            if (allocated(settings%edges)) deallocate (settings%edges)

            allocate (settings%radius(settings%Nsph))
            allocate (settings%edges(settings%Nsph))

         case ('RADIOS')
            read (val, *) settings%radius

         case ('EDGES')
            read (val, *) settings%edges

         case ('REGIONS')
            read (val, *) settings%Nregions

            if (allocated(settings%regionsID)) deallocate (settings%regionsID)
            if (allocated(settings%rho_regions)) deallocate (settings%rho_regions)
            if (allocated(settings%coord_regions)) deallocate (settings%coord_regions)

            allocate (settings%regionsID(settings%Nregions))
            allocate (settings%rho_regions(settings%Nregions))
            allocate (settings%coord_regions(3, settings%Nregions))
            allocate (temp_vec(3*settings%Nregions))

         case ('LOCATION')
            read (val, *) temp_vec
            settings%coord_regions = reshape(temp_vec, (/3, settings%Nregions/))
            deallocate (temp_vec)

         case ('ID_REGION')
            read (val, *) settings%regionsID

         case ('RHO_REGIONS')
            read (val, *) settings%rho_regions

         end select

      end do

100   close (iu)

   end subroutine read_set_femtic
!=========================================================
!=======
!=========================================================
   subroutine read_dem(fname, units, x, y, z, n)
      implicit none
      character(len=*), intent(in):: fname, units
      real(dp), allocatable, intent(out):: x(:), y(:), z(:)
      integer, intent(out):: n

      integer:: iu, i
      real(dp):: xx, yy, zz
      real(dp):: scale

      scale = 1.0_dp
      if (units == 'meters') scale = 1.0d-3

      open (newunit=iu, file=fname, status='old')
      n = 0
      do
         read (iu, *, end=10)
         n = n + 1
      end do
10    rewind (iu)

      allocate (x(n), y(n), z(n))
      do i = 1, n
         read (iu, *) xx, yy, zz
         x(i) = xx*scale
         y(i) = yy*scale
         z(i) = zz*scale
      end do
      close (iu)
   end subroutine
!=========================================================
!=======
!=========================================================
   subroutine check_domain(x, y, xmin, xmax, ymin, ymax, pad_x, pad_y)
      implicit none
      real(dp), intent(in):: x(:), y(:)
      real(dp), intent(in):: xmin, xmax, ymin, ymax, pad_x, pad_y
      real(dp)            :: xminDEM, yminDEM, xmaxDEM, ymaxDEM, dem_size_x, dem_size_y

      xminDEM = minval(x)
      xmaxDEM = maxval(x)
      yminDEM = minval(y)
      ymaxDEM = maxval(y)
      dem_size_x = xmaxDEM - xminDEM
      dem_size_y = ymaxDEM - yminDEM

      print *, 'DEM X range      :', xminDEM, xmaxDEM
      print *, 'DEM Y range      :', yminDEM, ymaxDEM
      print *, 'Domain X (+pad)  :', xminDOM - pad_x, xmaxDOM + pad_x
      print *, 'Domain Y (+pad)  :', yminDOM - pad_y, ymaxDOM + pad_y
      print *, 'DEM size (km)    :', dem_size_x, dem_size_y
      ! print*, 'Domain+pad (km)  :', dom_pad_size_x, dom_pad_size_y

      if (minval(x) > xmin - pad_x .or. maxval(x) < xmax + pad_x) stop 'DEM X does not cover domain+padding'
      if (minval(y) > ymin - pad_y .or. maxval(y) < ymax + pad_y) stop 'DEM Y does not cover domain+padding'
   end subroutine
!=========================================================
!=======
!=========================================================
   real function dms_to_decimal(deg, min, sec)
      real:: deg, min, sec
      real:: sign

      sign = 1.0
      if (deg < 0.0) sign = -1.0

      dms_to_decimal = sign*(abs(deg) + min/60.0 + sec/3600.0)
   end function
!=========================================================
!=======
!=========================================================
! subroutine parse_dms_line(line, value)
!   character(len=*), intent(in):: line
!   real(8), intent(out):: value

!   integer:: p_eq, p1, p2
!   real:: deg, min, sec
!   character(len = 64):: tmp

!   p_eq = index(line, "=")
!   tmp = adjustl(line(p_eq+1:))

!   read(tmp, '(f6.0, ":", f6.0, ":")') deg, min

!   value = dms_to_decimal(deg, min, sec)
! end subroutine
   subroutine parse_ref_value(line, value)
      implicit none
      character(len=*), intent(in)  :: line
      real(8), intent(out)          :: value

      character(len=:), allocatable:: s
      integer:: ipos, p1, p2
      real(8):: deg, min, sec
      real(8):: sign

      !----------------------------------------
      ! 1) Extraer texto después del último '='
      !----------------------------------------
      ipos = scan(line, "=", back=.true.)
      if (ipos == 0) then
         write (*, *) 'ERROR: No "=" found in line: ', trim(line)
         stop
      end if

      s = adjustl(trim(line(ipos + 1:)))

      !----------------------------------------
      ! 2) Detectar si es D:M:S o escalar
      !----------------------------------------
      if (index(s, ":") == 0) then
         ! Caso simple: número escalar
         read (s, *) value
         return
      end if

      !----------------------------------------
      ! 3) Conversión D:M:S → decimal
      !----------------------------------------
      sign = 1.0d0
      if (s(1:1) == '-') then
         sign = -1.0d0
         s = s(2:)
      end if

      p1 = index(s, ":")
      p2 = index(s(p1 + 1:), ":") + p1

      read (s(1:p1 - 1), *) deg
      read (s(p1 + 1:p2 - 1), *) min
      read (s(p2 + 1:), *) sec

      value = sign*(abs(deg) + min/60.0d0 + sec/3600.0d0)

   end subroutine parse_ref_value
!=========================================================
!=======
!=========================================================
   subroutine get_edi_file_list(edi_files, n_files)

      use iso_fortran_env, only: error_unit
      implicit none

      character(len=256), allocatable, intent(out):: edi_files(:)
      integer, intent(out):: n_files

      integer:: unit, ios
      character(len=256):: line
      character(len=*), parameter:: tmpfile = "edi_list.tmp"

      !------------------------------------------------------------
      ! 1) Listar archivos EDI usando el sistema
      !------------------------------------------------------------
      call execute_command_line( &
         "ls ./input_data/edi_files/*.edi > "//tmpfile, &
         exitstat=ios)

      if (ios /= 0) then
         write (error_unit, *) "ERROR: failed to list ./input_data/edi_files/*.edi"
         stop
      end if

      !------------------------------------------------------------
      ! 2) Contar archivos
      !------------------------------------------------------------
      n_files = 0
      open (newunit=unit, file=tmpfile, status="old", action="read")
      do
         read (unit, '(A)', iostat=ios) line
         if (ios /= 0) exit
         n_files = n_files + 1
      end do
      close (unit)

      if (n_files == 0) then
         write (error_unit, *) "ERROR: No EDI files found in ./input_data/edi_files/"
         stop
      end if

      !------------------------------------------------------------
      ! 3) Alocar vector de archivos
      !------------------------------------------------------------
      allocate (edi_files(n_files))

      !------------------------------------------------------------
      ! 4) Leer nombres de archivo
      !------------------------------------------------------------
      open (newunit=unit, file=tmpfile, status="old", action="read")
      do ios = 1, n_files
         read (unit, '(A)') edi_files(ios)
      end do
      close (unit)

      !------------------------------------------------------------
      ! 5) Limpiar archivo temporal
      !------------------------------------------------------------
      call execute_command_line("rm  -f "//tmpfile)

   end subroutine get_edi_file_list
!=========================================================
!=======
!=========================================================
   subroutine read_edi_files(edi_files, n, edi_id, edi_lat, edi_lon, edi_elev)
      implicit none

      integer, intent(in):: n
      character(len=*), intent(in):: edi_files(n)
      real(8), intent(out)         :: edi_lat(n), edi_lon(n), edi_elev(n)
      character(len=*), intent(out)         :: edi_id(n)

      integer:: i, unit, ios, p1, p2
      character(len=512):: line

      logical:: found_id, found_lat, found_lon, found_elev

      do i = 1, n

         found_id = .false.
         found_lat = .false.
         found_lon = .false.
         found_elev = .false.

         open (newunit=unit, file=edi_files(i), status="old", action="read")

         do
            read (unit, '(A)', iostat=ios) line
            if (ios /= 0) exit

            if (index(line, 'DATAID=') > 0) then
               p1 = index(line, '"')
               p2 = index(line(p1 + 1:), '"') + p1
               edi_id(i) = line(p1 + 1:p2 - 1)
               found_id = .true.
            end if

            if (index(line, 'REFLAT=') > 0) then
               call parse_ref_value(line, edi_lat(i))
               found_lat = .true.
            end if

            if (index(line, 'REFLONG=') > 0) then
               call parse_ref_value(line, edi_lon(i))
               found_lon = .true.
            end if

            if (index(line, 'REFELEV=') > 0) then
               read (line(index(line, '=') + 1:), *) edi_elev(i)
               found_elev = .true.
            end if

            if (found_id .and. found_lat .and. found_lon .and. found_elev) exit
         end do

         close (unit)

         if (.not. (found_lat .and. found_lon .and. found_elev)) then
            print *, "ERROR: Missing coordinates in ", trim(edi_files(i))
            stop
         end if
      end do
   end subroutine read_edi_files

!=========================================================
!=======
!=========================================================
   subroutine edi_to_utm(lat, lon, n_files, x, y, zone)
      implicit none
      integer, intent(in)  :: n_files
      real(8), intent(in)  :: lat(n_files), lon(n_files)      ! grados decimales
      real(8), intent(out):: x(n_files), y(n_files)           ! metros
      integer, intent(out):: zone

      integer:: i
      real(8):: lon_mean

      ! --- constantes---
      real(8), parameter:: semieje = 6378137.0d0
      real(8), parameter:: f = 1.0d0/298.257223563d0
      real(8), parameter:: k0 = 0.9996d0
      real(8), parameter:: pi = 3.141592653589793d0

      real(8):: e2, ep2
      real(8):: latr(n_files), lonr(n_files), lon0r
      real(8):: N(n_files), T(n_files), C(n_files), A(n_files), M(n_files)
      real(8):: lon0

      ! --- elipsoide---
      e2 = f*(2.0d0 - f)
      ep2 = e2/(1.0d0 - e2)

      ! --- zona UTM---
      lon_mean = sum(lon)/dble(n_files)
      zone = int((lon_mean + 180.0d0)/6.0d0) + 1
      lon0 = (zone - 1)*6.0d0 - 180.0d0 + 3.0d0

      do i = 1, n_files
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
         if (lat(i) < 0.0d0) y(i) = y(i) + 10000000.0d0
      end do

      !Covierto UTM en metros a UTM en kilometros
      x = x/1000.0d0
      y = y/1000.0d0

   end subroutine edi_to_utm
!=========================================================
!=======
!=========================================================
   subroutine compute_center(x, y, n, x0, y0)
      real(8), intent(in)  :: x(n), y(n)
      real(8), intent(out):: x0, y0
      integer, intent(in)  :: n

      x0 = sum(x)/dble(n)
      y0 = sum(y)/dble(n)
   end subroutine
!=========================================================
!=======
!=========================================================
   subroutine recenter_all(n_sites, n_dem, x0, y0, site_x, site_y, dem_x, dem_y)
      real(8), intent(inout):: site_x(n_sites), site_y(n_sites), dem_x(n_dem), dem_y(n_dem)
      real(8), intent(in)    :: x0, y0
      integer, intent(in)    :: n_sites, n_dem

      site_x(:) = site_x(:) - x0
      site_y(:) = site_y(:) - y0

      dem_x(:) = dem_x(:) - x0
      dem_y(:) = dem_y(:) - y0
   end subroutine
!=========================================================
!=======
!=========================================================
   subroutine write_dem_sites_utm(ediID, siteXm, siteYm, nn_sites, demXmts, demYmts, demZmts, &
                                  nn_dem, siteXkm, siteYkm, demXkm, demYkm)
      implicit none
      real(8), intent(in)  :: siteXm(nn_sites), siteYm(nn_sites), demXmts(nn_dem), demYmts(nn_dem), demZmts(nn_dem)
      character(len=*), intent(in)  :: ediID(nn_sites)
      integer, intent(in)  :: nn_sites, nn_dem
      real(8), intent(out):: siteXkm(nn_sites), siteYkm(nn_sites), demXkm(nn_dem), demYkm(nn_dem)
      character(len=11) :: dir
      integer:: i, iu

      siteXkm = siteXm
      siteYkm = siteYm

      demXkm = demXmts
      demYkm = demYmts

      dir = 'input_data/'
      ! ojo aqui cuando meshtran se pase a la estructura de folders de mtif entonces la ruta
      !debera ser:
      !     dir = '../preprocessing'
      !y esto es asi porque meshtran se ejecutara en bin/ y debe subir un nivel in entrar
      !en preprocessing
   
      ! ==========================
      ! Write SITE coordinates
      ! ==========================
      open (newunit=iu, file=trim(dir)//'PlotWithPython/sites_utm_km.dat', status='replace', action='write')

      write (iu, '(A)') '# x_km   y_km'
      do i = 1, nn_sites
         write (iu, '(A12, F15.5, 1X, F15.5)') ediID(i), siteXkm(i), siteYkm(i)
      end do
      close (iu)

      ! ==========================
      ! Write DEM coordinates
      ! ==========================
      open (newunit=iu, file=trim(dir)//'PlotWithPython/dem_utm_km.dat', status='replace', action='write')

      write (iu, '(A)') '# x_km   y_km   z_m'
      do i = 1, nn_dem
         write (iu, '(F15.5, 1X, F15.5, 1X, F15.5)') demXkm(i), demYkm(i), demZmts(i)
      end do
      close (iu)

   end subroutine
!=========================================================
!=======
!=========================================================
   subroutine snap_sites_to_dem(ediID, site_x, site_y, n_sites, dem_x, dem_y, dem_z, n_dem, fix_elev_site)

      implicit none
      integer, intent(in):: n_sites, n_dem
      real(8), intent(in):: site_x(n_sites), site_y(n_sites)
      real(8), intent(in):: dem_x(n_dem), dem_y(n_dem), dem_z(n_dem)
      real(8), intent(out):: fix_elev_site(n_sites)
      character(len=*) :: ediID(n_sites)

      character(len=512):: fname
      integer:: i, j, jmin, iu
      real(8):: dx, dy, d2, d2min
      real(8):: zmin

      ! Protección mínima: elevación mínima del DEM
      zmin = minval(dem_z)
      if (zmin < 1.0d0) zmin = 1.0d0   ! nunca 0 ni negativa

      do i = 1, n_sites
         d2min = huge(1.0d0)
         jmin = 1

         do j = 1, n_dem
            dx = dem_x(j) - site_x(i)
            dy = dem_y(j) - site_y(i)
            d2 = dx*dx + dy*dy

            if (d2 < d2min) then
               d2min = d2
               jmin = j
            end if
         end do

         fix_elev_site(i) = dem_z(jmin)

         ! seguridad extra
         if (fix_elev_site(i) < zmin) fix_elev_site(i) = zmin
      end do

      ! Writing final coordinate site files
      fname = trim(outdir)//'/sites_coord_elev.dat'
      open (newunit=iu, file=fname, status='replace', action='write')
      do i = 1, n_sites
         write (iu, '(A12, 3f15.5)') ediID(i), site_x(i), site_y(i), fix_elev_site(i)/1000.0d0
      end do

      close (iu)

   end subroutine
!=========================================================
!=======
!=========================================================
   subroutine write_observing_sites(site_x, site_y, n_sites, outdir, Nsph, edges, radius)
      implicit none

      ! ======================
      ! Inputs
      ! ======================
      integer, intent(in):: n_sites, Nsph
      real(8), intent(in):: site_x(n_sites), site_y(n_sites), edges(Nsph), radius(Nsph)
      character(len=*), intent(in):: outdir

      ! ======================
      ! Local variables
      ! ======================
      integer:: iu, i, k
      character(len=512):: fname

      ! ======================
      ! Define spheres (km)
      ! ======================
      ! radius = (/ 0.1d0, 0.3d0, 1.0d0, 3.0d0, 5.0d0 /)
      ! edge   = (/ 0.02d0, 0.05d0, 0.10d0, 0.30d0, 0.50d0 /)

      ! ======================
      ! Output file
      ! ======================
      fname = trim(outdir)//'/observing_site.dat'
      open (newunit=iu, file=fname, status='replace', action='write', form='formatted')

      ! ======================
      ! Write number of sites
      ! ======================
      write (iu, '(I9)') n_sites

      ! ======================
      ! Write each observation site
      ! ======================
      do i = 1, n_sites

         ! X  Y  Nsph
         write (iu, '(2F15.6, 1X, I3, 2x)', advance='no') site_x(i), site_y(i), Nsph

         ! (Ri, Li) pairs
         do k = 1, Nsph
            write (iu, '(F5.2, F5.2, 1x)', advance='no') radius(k), edges(k)
         end do

         ! New line
         write (iu, *)

      end do

      close (iu)

   end subroutine write_observing_sites

!=========================================================
!=======
!=========================================================
   subroutine define_analysis_domain(xmin, xmax, ymin, ymax, zmin, zmax)
      use mesh_config
      implicit none
      real(8), intent(in):: xmin, xmax, ymin, ymax, zmin, zmax
      integer:: iu

      ! xmin =  - pad_x
      ! xmax = maxval(x) + pad_x
      ! ymin = minval(y) - pad_y
      ! ymax = maxval(y) + pad_y

      open (newunit=iu, file="./"//trim(config%outdir)//"analysis_domain.dat", status='replace', action='write')
      write (iu, '(2F10.2)') xmin, xmax
      write (iu, '(2F10.2)') ymin, ymax
      write (iu, '(2F10.2)') zmin, zmax
      close (iu)

   end subroutine
!=========================================================
!=======
!=========================================================
   subroutine write_topography(topo_file, x, y, z, n, sea_level)!, outdir)
      use mesh_config
      implicit none
      character(len=*), intent(in):: topo_file
      real(dp), intent(in):: x(:), y(:)
      real(dp)   :: z(:)
      integer, intent(in):: n
      real(dp), intent(in):: sea_level
      ! character(len=*), intent(in):: outdir
      integer:: i, iu

      z = z/1000.0d0

      open (newunit=iu, file=trim(config%outdir)//topo_file, status='replace')
      do i = 1, n
         if (.not. has_sea) then
            ! Caso SIN mar: todo es tierra
            write (iu, '(3F15.5)') x(i), y(i), z(i)
         else
            ! Caso CON mar
            if (z(i) > sea_level) then
               write (iu, '(3F15.5)') x(i), y(i), z(i)
            else
               write (iu, '(3F15.5)') x(i), y(i), -1.0_dp
            end if
         end if
      end do
      close (iu)
   end subroutine
!=========================================================
!=======
!=========================================================
   subroutine write_bathymetry(bathy_file, x, y, z, n, sea_level, outdir)
      implicit none
      character(len=*), intent(in):: bathy_file
      real(dp), intent(in):: x(:), y(:)
      real(dp)    :: z(:)
      integer, intent(in):: n
      real(dp), intent(in):: sea_level
      character(len=*), intent(in):: outdir
      integer:: i, iu

      z = z/1000.0d0

      open (newunit=iu, file=trim(outdir)//bathy_file, status='replace')
      do i = 1, n
         if (.not. has_sea) then
            ! Caso SIN mar: no hay batimetría
            write (iu, '(3F15.5)') x(i), y(i), -1.0_dp/1000.0d0
         else
            ! Caso CON mar
            if (z(i) < sea_level) then
               write (iu, '(3F15.5)') x(i), y(i), -z(i)
            else
               write (iu, '(3F15.5)') x(i), y(i), -1.0_dp/1000.0d0
            end if
         end if
      end do
      close (iu)
   end subroutine
!=========================================================
!=======
!=========================================================
   subroutine write_coast_line(coast_line_file, xmin, xmax, ymin, ymax, outdir)
      implicit none
      character(len=*), intent(in):: coast_line_file
      real(dp), intent(in):: xmin, xmax, ymin, ymax
      character(len=*), intent(in):: outdir
      integer:: iu

      if (.not. has_sea) then
         open (newunit=iu, file=trim(outdir)//coast_line_file, status='replace')
         write (iu, *) 1
         write (iu, '(2(F21.15, 1X), 2I2)') xmin - 5.0d0, ymin - 5.0d0, 0, 0
         write (iu, '(2(F21.15, 1X), 2I2)') xmax + 5.0d0, ymin - 5.0d0, 0, 0
         write (iu, '(2(F21.15, 1X), 2I2)') xmax + 5.0d0, ymax + 5.0d0, 0, 0
         write (iu, '(2(F21.15, 1X), 2I2)') xmin - 5.0d0, ymax + 5.0d0, 1, 0
         close (iu)
      else
         stop 'Sea case not implemented yet'
      end if
   end subroutine
!=========================================================
!=======
!=========================================================
   subroutine run_makeTetraMesh_and_assign_regions(N_regions, coord_regions, ID_regions)
      implicit none

      integer, intent(in) :: N_regions, ID_regions(N_regions)
      real(8), intent(in) :: coord_regions(3, N_regions)

      integer :: stat, iu_in, iu_out, k
      logical :: ex, found_part4
      character(len=512) :: line

      ! -----------------------------
      ! 1. Ir a buildMesh
      ! -----------------------------
      ! call execute_command_line('cd buildMesh', wait=.true., exitstat=stat)
      ! if (stat /= 0) stop 'ERROR: cannot cd to buildMesh'

      ! -----------------------------
      ! 2. Copiar geometría
      ! -----------------------------
      call execute_command_line('echo " "')
      call execute_command_line('echo "--Running makeTetraMesh..."')
      call execute_command_line('cd buildMesh && cp ../input_data/geometry/* .', &
                                wait=.true., exitstat=stat)
      if (stat /= 0) stop 'ERROR: copying geometry failed'
      call execute_command_line('echo " "')

      ! -----------------------------
      ! 3. Ejecutar steps 1–4
      ! -----------------------------
      call execute_command_line('echo "step 1 --> Defining Land/Sea Boundary"')
      call execute_command_line('cd buildMesh && makeTetraMesh -stp 1', wait=.true., exitstat=stat)
      if (stat /= 0) then
         error stop 'ERROR: makeTetraMesh step 1 failed'
      end if
      call execute_command_line('echo "done..." && sleep 1')

      call execute_command_line('echo "step 2 --> Building 2D mesh..."')
      call execute_command_line('cd buildMesh && makeTetraMesh -stp 2', wait=.true., exitstat=stat)
      if (stat /= 0) then
         error stop 'ERROR: makeTetraMesh step 2 failed'
      end if
      call execute_command_line('echo "done..." && sleep 1')

      call execute_command_line('echo "step 3 --> Interpolating altitudes"')
      call execute_command_line('cd buildMesh && makeTetraMesh -stp 3', wait=.true., exitstat=stat)
      if (stat /= 0) then
         error stop 'ERROR: makeTetraMesh step 3 failed'
      end if
      call execute_command_line('echo "done..." && sleep 1')

      call execute_command_line('echo "step 4 --> Making surface mesh"')
      call execute_command_line('cd buildMesh && makeTetraMesh -stp 4', wait=.true., exitstat=stat)
      if (stat /= 0) then
         error stop 'ERROR: makeTetraMesh step 4 failed'
      end if
      call execute_command_line('echo "done..." && sleep 1')

      ! -----------------------------
      ! 4. Verificar output.poly
      ! -----------------------------
      inquire (file='buildMesh/output.poly', exist=ex)
      if (.not. ex) stop 'ERROR: output.poly was not generated'

      ! -----------------------------
      ! 5. Parchear regiones
      ! -----------------------------

      ! after_part4 = .false.
      found_part4 = .false.

      call execute_command_line('echo " " ')
      call execute_command_line('echo "Assigning regions in output.poly file"')
      open (newunit=iu_in, file='buildMesh/output.poly', status='old')
      open (newunit=iu_out, file='buildMesh/output.poly.tmp', status='replace')

      ! do
      !   read(iu_in, '(A)', end=100) line

      !   write(iu_out, '(A)') trim(line)

      !   if (index(line, '# Part 4') > 0) then
      !     read(iu_in, '(A)') line   ! esta es la línea "0"
      !     write(iu_out, '(I0)') 2   ! número de regiones

      !     ! ---- REGIONES ----
      !     do k =1,Nregions
      !       write(iu_out,'(I3,3F10.3,I4,1PE12.4)') k , coord_regions(:,k) , ID_regions(k), 1.0e9
      !     enddo
      !   else
      !     write(*,*) ' There is no #Part 4 content on output.poly where it assign regions'
      !     write(*,*) 'Aborting tetgen execution'
      !     stop
      !   end if
      ! end do

      do
         read (iu_in, '(A)', end=100) line

         write (iu_out, '(A)') trim(line)

         if (index(line, '# Part 4') > 0) then
            found_part4 = .true.

            read (iu_in, '(A)') line   ! leer el "0"
            write (iu_out, '(I0)') Nregions

            do k = 1, Nregions
               write (iu_out, '(I3,3F10.3,I4,1PE12.4)') k, coord_regions(:, k), ID_regions(k), 1.0e9
            end do
         end if

      end do

100   continue

      if (.not. found_part4) then
         write (*, *) 'There is no # Part 4 section in output.poly'
         error stop
      end if

      close (iu_in)
      close (iu_out)

      call execute_command_line('mv buildMesh/output.poly.tmp buildMesh/output.poly')

      write (*, *) 'Ok: ✅ makeTetraMesh steps 1–4 done and regions added to output.poly'

   end subroutine run_makeTetraMesh_and_assign_regions
!=========================================================
!=======
!=========================================================
   subroutine setGlobalMeshRefinement(outdir_refi, xcenter, ycenter, Nglob_ellipses, rot_glob, corXsite, corYsite, &
                                      nSites, n_site_ellipses, max_edge_len_within_ellipse, len_alon_x_axis)

      implicit none
      character(len=*), intent(in):: outdir_refi
      integer, intent(in)   :: Nglob_ellipses
      real(dp), intent(in)  :: ycenter, xcenter, rot_glob
      character(len=512)    :: fname
      real(dp)              :: z0, oblatness_on_ZXplane
      integer               :: iu, i, nSites, n_site_ellipses, k
      ! Parámetros hardcodeados (pueden ser arreglos en el futuro)
      real(dp) :: a(Nglob_ellipses), lengths(Nglob_ellipses), fh(Nglob_ellipses), fvp(Nglob_ellipses), fvm(Nglob_ellipses)
      real(dp), intent(in) :: len_alon_x_axis(n_site_ellipses), max_edge_len_within_ellipse(n_site_ellipses)
      real(dp), intent(in) :: corXsite(nSites), corYsite(nSites)

      z0 = 0.0d0

      ! Datos del ejemplo del autor
      a = (/40.0, 45.0, 50.0, 60.0, 80.0, 100.0, 200.0, 300.0, 400.0, 500.0/)
      lengths = (/1.0, 1.5, 3.0, 5.0, 8.0, 10.0, 15.0, 20.0, 30.0, 45.0/)
      fh = (/0.5, 0.5, 0.5, 0.3, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0/)
      fvp = (/0.7, 0.5, 0.4, 0.3, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0/)
      fvm = (/0.9, 0.7, 0.7, 0.5, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0/)

      a = a*1.0d0

      fname = trim(outdir_refi)//'makeMtr.param'
      open (newunit=iu, file=fname, status='replace', action='write')

      ! 1. Coordenadas del centro (Y, X, Z)
      write (iu, '(3F12.5)') xcenter, ycenter, z0

      ! 2. Ángulo de rotación
      write (iu, '(F12.2)') rot_glob

      ! 3. Número de elipsoides
      write (iu, '(I5)') Nglob_ellipses

      ! 4. Bloque de elipsoides
      do i = 1, min(Nglob_ellipses, 10)
         write (iu, '(F6.1, F6.1, 1x, 3F6.2)') a(i), lengths(i), fh(i), fvp(i), fvm(i)
      end do

      close (iu)
      print *, 'Archivo makeMtr.param creado exitosamente en buildMesh/'

      oblatness_on_ZXplane = 0.3

      fname = trim(outdir_refi)//'/obs_site.dat'
      open (newunit=iu, file=fname, status='replace', action='write', form='formatted')

      ! ======================
      ! Write number of sites
      ! ======================
      write (iu, '(I4)') nSites

      ! ======================
      ! Write each observation site
      ! ======================
      do i = 1, nSites
         ! X  Y  z0
         write (iu, '(3F15.6 )') corXsite(i), corYsite(i), z0
         write (iu, '(I1)') n_site_ellipses

         ! (Ri, Li) pairs
         do k = 1, n_site_ellipses
            write (iu, '(3F5.2 )') len_alon_x_axis(k), max_edge_len_within_ellipse(k), oblatness_on_zxplane
         end do
      end do

      close (iu)
      print *, 'Archivo obs_site.dat creado exitosamente en buildMesh/'

   end subroutine setGlobalMeshRefinement
!=========================================================
!=======
!=========================================================
   subroutine resistivitty_attribute(out_dir, Nsites, coorXsite, coorYsite, coorZsite, &
                                     N_paramEsfer, esferEdges, radiusEsfer, N_regions, ID_regions, rhoRegions, repRegion, fixed)

      implicit none

      ! ======================
      ! Inputs
      ! ======================
      character(len=*), intent(in):: out_dir
      integer, intent(in) :: Nsites, N_regions, N_paramEsfer
      real(8), intent(in) :: coorXsite(Nsites), coorYsite(Nsites), coorZsite(Nsites)
      real(8), intent(in) :: esferEdges(N_paramEsfer), radiusEsfer(N_paramEsfer)
      real(8), intent(in) :: rhoRegions(N_regions)
      integer, intent(in) :: repRegion(N_regions), fixed(N_regions), ID_regions(N_regions)

      ! ======================
      ! Local variables
      ! ======================
      integer :: iu, i, j
      character(len=512) :: fname

      ! Parámetros para la Parte 2 (Elipsoides Regionales)
      integer, parameter :: Ne = 9
      real(8) :: a_reg(Ne), len_reg(Ne), fh_reg(Ne), fv_reg(Ne)

      ! Datos del ejemplo de Usui para la Parte 2
      a_reg = (/40.0, 45.0, 50.0, 60.0, 100.0, 200.0, 300.0, 500.0, 1000.0/)
      len_reg = (/2.0, 3.0, 5.0, 10.0, 100.0, 200.0, 300.0, 500.0, 1000.0/)
      fh_reg = (/0.5, 0.5, 0.4, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0/)
      fv_reg = (/0.7, 0.7, 0.7, 0.6, 0.5, 0.3, 0.2, 0.1, 0.0/)

      fname = trim(out_dir)//'resistivity_attr.dat'
      open (newunit=iu, file=fname, status='replace', action='write')

      ! -----------------------------------------------------------
      ! PARTE 1: Atributos de Región (Nreg = 3: Aire, Mar, Subsuelo)
      ! Formato: ID  Resistividad  RepeatNumber  FixFlag
      ! -----------------------------------------------------------
      write (iu, '(I2)') N_regions
      do j = 1, N_regions
         write (iu, '(I0, 2X, ES10.3, 1X, I3, 1X, I3)') ID_regions(j), rhoRegions(j), repRegion(j), fixed(j)
      end do

      ! -----------------------------------------------------------
      ! PARTE 2: Elipsoides de Control Regional
      ! -----------------------------------------------------------
      write (iu, '(3F12.4)') 0.0, 0.0, 0.0  ! Centro (Y, X, Z) en km
      write (iu, '(F12.4)') 0.0            ! Rotación (grados)
      write (iu, '(I0)') Ne                ! Número de elipsoides
      do i = 1, Ne
         write (iu, '(F10.1, F10.1, 2F8.2)') a_reg(i), len_reg(i), fh_reg(i), fv_reg(i)
      end do

      ! -----------------------------------------------------------
      ! PARTE 3: Refinamiento Local por Puntos (Estaciones)
      ! -----------------------------------------------------------
      write (iu, '(I0)') Nsites           ! Número de estaciones
      do i = 1, Nsites
         ! Coordenadas (X, Y, Z) en km. Nota: Z viene en metros en tu código.
         write (iu, '(3(F0.6,2x))') coorXsite(i), coorYsite(i), coorZsite(i)

         ! Número de esferas por sitio (Ejemplo: 2 esferas)
         write (iu, '(I0)') N_paramEsfer
         ! Formato: Radio_km  len_km
         do j = 1, N_paramEsfer
            write (iu, '(F5.1, 1X, F5.1)') radiusEsfer(j), esferEdges(j)
         end do

      end do

      close (iu)
      print *, 'OK: resistivity_attr.dat escrito exitosamente'

   end subroutine
!=========================================================
!=======
!=========================================================
   subroutine run_TETGEN_and_refine_mesh(Nref)
      implicit none

      integer, intent(in) :: Nref
      integer :: r, stat
      character(len=512) :: cmd

      ! -----------------------------
      ! Runing tetgen on output.poly file
      ! -----------------------------
      call execute_command_line('echo " "')
      call execute_command_line('echo "running tetgen --> Building 2D Mesh including topography"')
      call execute_command_line('cd buildMesh && tetgen -nVpYAakq3.0/0 output.poly', wait=.true., exitstat=stat)
      if (stat /= 0) then
         error stop 'ERROR: tetgen execution failed'
      end if
      call execute_command_line(' cd buildMesh && head output.1.ele')
      call execute_command_line('echo "done..." && sleep 3')
      call execute_command_line('echo " " ')
      call execute_command_line('echo " Preparing mesh refinement " && cd buildMesh && mkdir -p refinement', &
                                wait=.true., exitstat=stat)
      if (stat /= 0) error stop 'ERROR: could not create refinement directory'

      ! inquire(file='buildMesh/refinement', exist=dir_exists)
      ! if (dir_exists) then
      !     write(*,*) 'Refinement directory successfully created.'
      ! else
      !     error stop 'ERROR: refinement directory does not exist.'
      ! endif

      call execute_command_line('cd buildMesh && cp output.1* refinement')
      pause
      call execute_command_line('cd buildMesh && cp makeMtr.param obs_site.dat refinement', wait=.true., exitstat=stat)
      if (stat /= 0) stop 'ERROR: there is no files content refinement parameters'

      ! --------------------------------------------------
      ! Iterative refinement
      ! for i in 1 2 ... Nref
      ! --------------------------------------------------
      do r = 1, Nref - 1

         write (cmd, '(A,I0,A)') 'cd buildMesh/refinement && head output.', r, '.ele'
         call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
         write (*, '(A,I0)') 'Refinement iteration: ', r
         pause
         ! makeMtr output.$r
         write (cmd, '(A,I0)') 'cd buildMesh/refinement && makeMtr output.', r
         call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
         if (stat /= 0) error stop 'ERROR in makeMtr'
         ! tetgen ... output.$r
         write (cmd, '(A,I0)') 'cd buildMesh/refinement && tetgen -nmpYVrAakq3.0/0 output.', r
         call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
         if (stat /= 0) error stop 'ERROR in tetgen refinement'

      end do

      call execute_command_line('echo "Refinement done..." && sleep 1')

   end subroutine run_TETGEN_and_refine_mesh
!=========================================================
!=======
!=========================================================
   subroutine run_TetGen2Femtic(Nref)
      implicit none

      integer, intent(in) :: Nref
      character(len=256) :: cmd, fname
      integer :: stat
      logical :: ex1, ex2, ex3

      call execute_command_line('echo " "')
      call execute_command_line('echo "Final step: tetgen2femtic execution"')
      call execute_command_line('cd buildMesh && mkdir -p tetgenTOfemtic', wait=.true., exitstat=stat)
      if (stat /= 0) then
         error stop 'ERROR: could not be created tetgenTOfemtic directory'
      end if
      ! Copiar el último refinement
      write (cmd, '(A,I0,A)') 'cd buildMesh && cp refinement/output.', Nref, '* tetgenTOfemtic'
      call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
      if (stat /= 0) error stop 'ERROR copying final refinement to tetgenTOfemtic folder'

      call execute_command_line('cd buildMesh && cp resistivity_attr.dat tetgenTOfemtic')

      ! Ejecutar TetGen2Femtic con el último número
      write (cmd, '(A,I0)') 'cd buildMesh/tetgenTOfemtic && TetGen2Femtic output.', Nref
      call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
      if (stat /= 0) error stop 'ERROR running TetGen2Femtic'

      call execute_command_line('echo "done..." && sleep 1')

      write (cmd, '(A,I0,A)') 'cp buildMesh/tetgenTOfemtic/&
    &      {mesh.dat,resistivity_block_iter0.dat,output.', Nref, '.femtic.vtk} input_data/inv'
      call execute_command_line(trim(cmd), wait=.true., exitstat=stat)

      ! Verificación
      inquire (file='input_data/inv/mesh.dat', exist=ex1)
      inquire (file='input_data/inv/resistivity_block_iter0.dat', exist=ex2)

      write (fname, '(A,I0,A)') 'input_data/inv/output.', Nref, '.femtic.vtk'
      inquire (file=fname, exist=ex3)

      if (ex1 .and. ex2 .and. ex3) then
         write (*, *) 'Input files for running Femtic are successfully copied to input_data/inv'
      else
         write (*, *) 'ERROR: Missing files after copy.'
         write (*, *) 'mesh.dat exists? ', ex1
         write (*, *) 'resistivity_block_iter0.dat exists? ', ex2
         write (*, *) 'vtk exists? ', ex3
         error stop
      end if

      call execute_command_line('echo " " ')
      call execute_command_line('echo " 🏁 Finishing meshTran execution..." && sleep 3')

   end subroutine run_TetGen2Femtic

end program femtic_mesh_driver
