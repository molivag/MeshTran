program femtic_mesh_driver

   use mesh_config
   use mesh_entities
   implicit none

   type(GlobalRefinement) :: globRefi
   type(ParamRefinement) :: paramRefi
   TYPE(MeshSettings)::settings
   type(ModelRegion) :: regions
   ! type(SiteSet) :: sites

   ! Control parameters
   character(len=256), allocatable:: edi_files(:), edi_id(:)
   ! logical:: has_sea
   ! real(dp):: xminDOM, xmaxDOM, yminDOM, ymaxDOM!, zminDOM, zmaxDOM, pad_x, pad_y,
   real(dp), allocatable:: edi_lat(:), edi_lon(:), edi_elev(:), site_x_km(:), site_y_km(:), dem_x_km(:), dem_y_km(:)
   real(dp)::  sea_level, x0, y0

   ! DEM data
   integer:: i, n_dem, zone, n_edi_files, n_sites!, Nsph, n_elipses, ellipsForSite, NparamEsfer, Nregions, refi_tetgen
   real(dp), allocatable:: site_x(:), site_y(:), site_z(:), dem_x(:), dem_y(:), dem_z(:), fix_elev_site(:)

   !---------------------------------------------------
   !     Read input configutation file
   !---------------------------------------------------
   call read_set_femtic('bin/set_meshtran.io', settings, paramRefi, globRefi, regions)
   !---------------------------------------------------
   !     Check if Read input configutation file
   !---------------------------------------------------
   call mesh_convertion(settings)
   !---------------------------------------------------
   !     Read Digital Elevation Model
   !---------------------------------------------------
   call read_dem(settings, dem_x, dem_y, dem_z, n_dem)
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
   ! call check_domain(dem_x_km, dem_y_km, xminDOM, xmaxDOM, yminDOM, ymaxDOM, pad_x, pad_y)
   call check_domain(dem_x_km, dem_y_km, settings)

   ! !Write files in mesh coordinates centered at anchor point

   call generate_observe_dat(edi_files, n_edi_files, site_x_km, site_y_km)
   call define_analysis_domain(settings)
   call write_topography(settings, dem_x_km, dem_y_km, dem_z, n_dem)
   call write_bathymetry(settings, dem_x_km, dem_y_km, dem_z, n_dem)
   call write_coast_line(settings)
   call write_observing_sites(site_x_km, site_y_km, n_edi_files, paramRefi)
   call control_mesh(x0, y0, site_x_km, site_y_km, n_edi_files, settings, paramRefi)

   call write_makeMtr_param(0.0d0, 0.0d0, site_x_km, site_y_km, n_edi_files, paramRefi, globRefi)
   call write_obs_sites(site_x_km, site_y_km, fix_elev_site, n_edi_files, paramRefi)

   call resistivitty_attribute(n_sites, site_x_km, site_y_km, fix_elev_site, settings, globRefi, regions)

   call run_makeTetraMesh_and_assign_regions(regions)
   call run_TETGEN_and_refine_mesh(globRefi)
   call run_TetGen2Femtic(settings, globRefi)

contains
   !=========================================================
   !=======
   !=========================================================
   subroutine read_set_femtic(fname, OBJsettings, OBJparamRefi, OBJglobRefi, OBJmodReg)

      use mesh_entities
      implicit none

      character(len=*), intent(in) :: fname
      type(MeshSettings), intent(inout) :: OBJsettings
      type(ParamRefinement), INTENT(INout) :: OBJparamRefi
      type(GlobalRefinement), INTENT(INOUT) :: OBJglobRefi
      type(ModelRegion), INTENT(INOUT) :: OBJmodReg

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

         case ('MESH_NATURE')
            OBJsettings%mesh_nature = trim(val)

         case ('DEM_FILE')
            OBJsettings%dem_file = trim(val)

         case ('TOPO_FILE')
            OBJsettings%topography_file = trim(val)

         case ('BATHY_FILE')
            OBJsettings%bathymetry_file = trim(val)

         case ('COSLI_FILE')
            OBJsettings%coastLine_file = trim(val)

         case ('DEM_UNITS')
            OBJsettings%dem_units = trim(val)

         case ('SEA_LEVEL')
            read (val, *) OBJsettings%sea_level

         case ('XMIN_DOM')
            read (val, *) OBJsettings%xminDOM

         case ('XMAX_DOM')
            read (val, *) OBJsettings%xmaxDOM

         case ('YMIN_DOM')
            read (val, *) OBJsettings%yminDOM

         case ('YMAX_DOM')
            read (val, *) OBJsettings%ymaxDOM

         case ('ZMIN_DOM')
            read (val, *) OBJsettings%zminDOM

         case ('ZMAX_DOM')
            read (val, *) OBJsettings%zmaxDOM

         case ('PAD_X')
            read (val, *) OBJsettings%pad_x

         case ('PAD_Y')
            read (val, *) OBJsettings%pad_y

            ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
            ! = = = = refinement for sites at the surface and control = = = =
            ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
         case ('NUM_ELIPSES')
            read (val, *) OBJparamRefi%Nelipses

         case ('NUM_ESFERAS')
            read (val, *) OBJparamRefi%Nsph

         case ('ROTATION')
            read (val, *) OBJparamRefi%rotation

         case ('minRADIO')
            read (val, *) OBJparamRefi%minRAD

         case ('maxRADIO')
            read (val, *) OBJparamRefi%maxRAD

         case ('minEDGES')
            read (val, *) OBJparamRefi%minEDGE

         ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
         ! = = = = To build the surface mesh (cascaroin) = = = = = = = = = = 
         ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
         case ('SITEpadding')
            read (val, *) OBJparamRefi%paddingRefi

         case ('FARelemSIZE')
            read (val, *) OBJparamRefi%sizeBoundary

         case ('SURF_RESOLUTION')
            read (val, *) OBJparamRefi%coreResol

         case ('GROWTH')
            read (val, *) OBJparamRefi%growthFactor

            ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
            ! = = = = Refinement for sites at the volumetric mesh = = = = = =
            ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
         case ('REFI_ELLIPSES')
            read (val, *) OBJglobRefi%levels

         case ('HIGH_RESOL_LAYER')
            read (val, *) OBJglobRefi%depth_min_factor

         case ('NEAR_FIELD_RESOL')
            read (val, *) OBJglobRefi%core_resolution

         case ('FAR_FIELD_RESOL')
            read (val, *) OBJglobRefi%farfield_resolution

         case ('FRAC_SKIN_DEPTH')
            read (val, *) OBJglobRefi%frac_skin_depth

         case ('ITER_TET_REFI')
            read (val, *) OBJglobRefi%n_iterative_refi

         case ('BKGRD_RHO')
            read (val, *) OBJglobRefi%background_rho

         case ('H_PADDING')
            read (val, *) OBJglobRefi%horizontal_padding

         case ('F_MIN_HZ')
            read (val, *) OBJglobRefi%fmin_hz

         case ('SITE_ELLIPSES')
            read (val, *) OBJglobRefi%n_ellipses_site

         case ('TARGET_DEPTH')
            read (val, *) OBJglobRefi%target_depth

         case ('V_PADDING')
            read (val, *) OBJglobRefi%vertical_padding_factor

            ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
            ! = = = = Refinement for mesh parameter based on spheres  = = = =
            ! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
         case ('ELLIPSES')
            read (val, *) OBJmodReg%NpEllipses

         case ('REGIONS')
            read (val, *) OBJmodReg%Nregions

            if (allocated(OBJmodReg%ID)) deallocate (OBJmodReg%ID)
            if (allocated(OBJmodReg%rho)) deallocate (OBJmodReg%rho)
            if (allocated(OBJmodReg%coord)) deallocate (OBJmodReg%coord)

            allocate (OBJmodReg%ID(OBJmodReg%Nregions))
            allocate (OBJmodReg%rho(OBJmodReg%Nregions))
            allocate (OBJmodReg%coord(3, OBJmodReg%Nregions))
            allocate (OBJmodReg%repeatPartition(OBJmodReg%Nregions))
            allocate (OBJmodReg%isRHOfix(OBJmodReg%Nregions))
            allocate (temp_vec(3*OBJmodReg%Nregions))

         case ('LOCATION')
            read (val, *) temp_vec
            OBJmodReg%coord = reshape(temp_vec, (/3, OBJmodReg%Nregions/))
            deallocate (temp_vec)

         case ('ID_REGION')
            read (val, *) OBJmodReg%ID

         case ('RHO_REGIONS')
            if (.not. allocated(OBJmodReg%rho)) then
               write (*, *) 'ERROR: RHO_REGIONS defined before REGIONS'
               stop
            end if

            read (val, *) OBJmodReg%rho

         case ('REP_PARTITION')
            if (.not. allocated(OBJmodReg%repeatPartition)) then
               write (*, *) 'ERROR: REP_PARTITION defined before REGIONS'
               stop
            end if
            read (val, *) OBJmodReg%repeatPartition

         case ('FIX_RESISTIVITY')
            if (.not. allocated(OBJmodReg%isRHOfix)) then
               write (*, *) 'ERROR: FIX_RESISTIVITY defined before REGIONS'
               stop
            end if
            read (val, *) OBJmodReg%isRHOfix

         case ('PARAM_ESFER')
            read (val, *) OBJmodReg%NparamEsfer

            if (allocated(OBJmodReg%radiusForEsfer)) deallocate (OBJmodReg%radiusForEsfer)
            if (allocated(OBJmodReg%edgesForEsfer)) deallocate (OBJmodReg%edgesForEsfer)
            allocate (OBJmodReg%radiusForEsfer(OBJmodReg%NparamEsfer))
            allocate (OBJmodReg%edgesForEsfer(OBJmodReg%NparamEsfer))

         case ('PARAM_RADIOS')
            if (.not. allocated(OBJmodReg%radiusForEsfer)) then
               write (*, *) 'ERROR: EDGES defined before paramESFERAS'
               stop
            end if
            read (val, *) OBJmodReg%radiusForEsfer
         case ('PARAM_EDGES')
            if (.not. allocated(OBJmodReg%edgesForEsfer)) then
               write (*, *) 'ERROR: EDGES defined before paramESFERAS'
               stop
            end if
            read (val, *) OBJmodReg%edgesForEsfer

            !## => Regions Atributes in the model

         case ('ELIPSES')
            read (val, *) OBJmodReg%NpEllipses
         end select

      end do

100   close (iu)

   end subroutine read_set_femtic
!=========================================================
!=======
!=========================================================
   subroutine generate_observe_dat(edi_files, n_files, site_x, site_y)
      implicit none
      integer, intent(in) :: n_files
      character(len=*), intent(in) :: edi_files(n_files)
      real(dp), intent(in) :: site_x(n_files), site_y(n_files)

      integer :: i, unit_in, unit_out, ios, nf, nf_final, j
      character(len=512) :: line

      ! Arrays dinámicos para los datos de cada estación
      real(dp), allocatable :: freq(:), zr(:, :), zi(:, :), zvar(:,:), zerr(:, :), sdSI(:,:)
      ! El tensor Z en EDI se lee por componentes: XX, XY, YX, YY
      ! Usaremos un array zr(nf, 4) donde 1=XX, 2=XY, 3=YX, 4=YY
      logical :: ex1
      logical, allocatable :: mask(:)  ! Para marcar qué frecuencias se quedan

      ! Variables de control (Estas vendrían de tu archivo de conf)
      integer :: subsampling_mode
      real(dp) :: fmin,noise_floor
      real(dp) :: fmax

      noise_floor = 1.0d-15 ! Epsilon de seguridad para evitar división por cero
      subsampling_mode = 0   ! 0=No, 1=Sí (Rango fmin/fmax)
      fmin = 0.01d0          ! Frecuencia mínima deseada
      fmax = 10.0d0          ! Frecuencia máxima deseada (evita las muy altas)

      open (newunit=unit_out, file='computing/observe.dat', status='replace')

      ! 1. Cabecera Global de FEMTIC
      write (unit_out, '(A, I6)') 'MT', n_files

      do i = 1, n_files
         open (newunit=unit_in, file=edi_files(i), status='old', action='read')

         ! 2. Encontrar NFREQ general de la estación
         nf = 0
         do
            read (unit_in, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line, 'NFREQ=') > 0) then
               read (line(index(line, '=') + 1:), *) nf
               exit
            end if
         end do

         if (nf <= 0) then
            print *, "Error: NFREQ no encontrado en ", trim(edi_files(i))
            stop
         end if

         ! Dimensionar vectores para esta estación específica
         allocate (freq(nf), zr(nf, 4), zi(nf, 4), zerr(nf, 4), zvar(nf, 4), sdSI(nf,4))
         allocate (mask(nf))

         ! 3. Leer Frecuencias (>FREQ)
         call read_edi_block(unit_in, '>FREQ', nf, freq)

         ! 4. Leer Componentes del Tensor (Real, Imag, Varianza)
         ! XX
         call read_edi_block(unit_in, '>ZXXR', nf, zr(:, 1))
         call read_edi_block(unit_in, '>ZXXI', nf, zi(:, 1))
         call read_edi_block(unit_in, '>ZXX.VAR', nf, zvar(:, 1))
         ! XY
         call read_edi_block(unit_in, '>ZXYR', nf, zr(:, 2))
         call read_edi_block(unit_in, '>ZXYI', nf, zi(:, 2))
         call read_edi_block(unit_in, '>ZXY.VAR', nf, zvar(:, 2))
         ! YX
         call read_edi_block(unit_in, '>ZYXR', nf, zr(:, 3))
         call read_edi_block(unit_in, '>ZYXI', nf, zi(:, 3))
         call read_edi_block(unit_in, '>ZYX.VAR', nf, zvar(:, 3))
         ! YY
         call read_edi_block(unit_in, '>ZYYR', nf, zr(:, 4))
         call read_edi_block(unit_in, '>ZYYI', nf, zi(:, 4))
         call read_edi_block(unit_in, '>ZYY.VAR', nf, zvar(:, 4))


         ! --- TRANSFORMACIONES CRÍTICAS ---
         ! 5. Convertir a unidades SI (V/m / A/m)
         call convert_EDI_file_units(nf, zvar, zr, zi, sdSI)
         
         ! 6. Apply conjugate to accomplish FEMTIC time dependence exp(-iωt) 
         call apply_complex_conjugate(nf,zi)
         ! 6. Calcular Error Floor (basado en el 5% y 20%)
         call apply_error_floor(nf, zr, zi, sdSI, zerr)
         ! ---------------------------------

         ! 7. Lógica de Subsampling
         mask = .true.    ! Por defecto todas se quedan
         nf_final = nf ! Por defecto el número final es el del EDI

         if (subsampling_mode == 1) then
            nf_final = 0
            do j = 1, nf
               ! Si la frecuencia está fuera del rango, la descartamos
               if (freq(j) < fmin .or. freq(j) > fmax) then
                  mask(j) = .false.
               else
                  nf_final = nf_final + 1
               end if
            end do
         end if

         ! 8. Escritura condicionada
         ! Si nf_final es 0 porque el rango fue muy estricto, avisar
         if (nf_final == 0) then
            print *, "Aviso: Estación ", trim(edi_files(i)), " no tiene frecuencias en el rango."
         end if

         ! 9. Escribir bloque de estación en observe.dat
         ! Formato: ID_est, ID_H_field, X, Y (en km o m según tu malla)
         write (unit_out, '(I6, I6, 1x, 2F12.5)') i, i, site_x(i), site_y(i)
         write (unit_out, '(I6)') nf_final

         do j = 1, nf
            if (mask(j)) then
               ! FEMTIC: Freq, ZxxR, ZxxI, ZxyR, ZxyI, ZyxR, ZyxI, ZyyR, ZyyI, + 8 errores
               ! Error = sqrt(Varianza)
               write (unit_out, '(F12.5, 1x, 16(1x, ES12.5))') &
                  freq(j), &
                  zr(j, 1), zi(j, 1), zr(j, 2), zi(j, 2), &
                  zr(j, 3), zi(j, 3), zr(j, 4), zi(j, 4), &
                  zerr(j, 1), zerr(j, 1), & ! Error XX (R e I)
                  zerr(j, 2), zerr(j, 2), & ! Error XY
                  zerr(j, 3), zerr(j, 3), & ! Error YX
                  zerr(j, 4), zerr(j, 4)   ! Error YY
            end if
         end do

         ! deallocate (freq, zr, zi, zvar, mask)
         deallocate (freq, zr, zi, zvar, zerr, sdSI, mask)
         close (unit_in)

      end do
      write(unit_out, '(A)') 'END'

      close (unit_out)
      inquire (file='computing/observe.dat', exist=ex1)

      if (ex1) then
         write (*, *) 'Observe.dat file successfully generated into computing/'
      else
         write (*, *) 'ERROR: Missing observe.dat file in computing/'
         write (*, *) 'observe.dat exists? ', ex1
         error stop
      end if
   end subroutine generate_observe_dat
!=========================================================
!=======
!=========================================================
subroutine apply_complex_conjugate(nf, zi)
    implicit none
    integer, intent(in) :: nf
    real(dp), intent(inout) :: zi(nf, 4)
    
    ! Multiplicación por -1 de la parte imaginaria
    ! Esto rota la fase de +45 (EDI) a -45 (FEMTIC)
    zi = -1.0d0 * zi
end subroutine apply_complex_conjugate
!=========================================================
!=======
!=========================================================
   subroutine convert_EDI_file_units(nf, zvar, zr, zi, sdSI)
      real(dp), parameter :: PI = 4 * atan(1.0_16)
      integer, intent(in) :: nf
      real(dp), intent(inout) :: zr(nf, 4), zi(nf, 4)
      real(dp), intent(inout) :: zvar(nf,4)
      real(dp), INTENT(OUT) :: sdSI(nf,4)
      real(dp) :: mu0, scale_factor, sd_sigma(nf,4)

      sd_sigma = sqrt(zvar)

      mu0 = 4.0d0 * PI * 1.0E-7 
      scale_factor = mu0 * 1.0E3  ! De mV/km/nT a V/m / A/m

      zr = zr * scale_factor
      zi = zi * scale_factor
      sdSI = sd_sigma * scale_factor

   end subroutine convert_EDI_file_units
!=========================================================
!=======
!=========================================================
   subroutine apply_error_floor(nf, zr, zi, SD, zerr)
      integer,   intent(in) :: nf
      real(dp),  intent(in) :: zr(nf, 4), zi(nf, 4), SD(nf,4)
      real(dp)              :: amp_xy, amp_yx, amp_ref
      real(dp),  intent(out) :: zerr(nf, 4)
      integer  :: j

      do j = 1, nf
         ! 1. Calculamos la amplitud individual de las componentes principales
         amp_xy = sqrt(zr(j, 2)**2 + zi(j, 2)**2)
         amp_yx = sqrt(zr(j, 3)**2 + zi(j, 3)**2)

         ! 2. Media geométrica de las amplitudes principales
         amp_ref = sqrt(amp_xy * amp_yx)

         ! Aplicamos los porcentajes de acuerdo al manual/femticPy
         zerr(j, 1) = max(SD(j,1), amp_ref * 0.20d0)  ! INdiag  ZXX (20%)
         zerr(j, 2) = max(SD(j,2), amp_ref * 0.05d0)  ! OFFdiag ZXY (5%)
         zerr(j, 3) = max(SD(j,3), amp_ref * 0.20d0)  ! OFFdiag ZYX (5%)
         zerr(j, 4) = max(SD(j,4), amp_ref * 0.05d0)  ! INdiag  ZYY (20%)
      end do
   end subroutine apply_error_floor
!=========================================================
!=======
!=========================================================
   subroutine read_edi_block(unit, tag, n, values)
      integer, intent(in) :: unit, n
      character(len=*), intent(in) :: tag
      real(dp), intent(out) :: values(n)

      character(len=512) :: line
      integer :: ios

      rewind (unit) ! Volver al inicio para buscar el tag
      do
         read (unit, '(A)', iostat=ios) line
         if (ios /= 0) return
         if (index(line, trim(tag)) > 0) exit
      end do

      ! Leer n valores. Fortran automáticamente saltará líneas
      ! hasta encontrar la cantidad de números solicitada por el array 'values'
      read (unit, *, iostat=ios) values
   end subroutine read_edi_block
!=========================================================
!=======
!=========================================================
   subroutine control_mesh(x0, y0, siteX, siteY, n_sites, OBJsettings, OBJparamRefi)

      implicit none

      integer :: iu, ios
      integer, intent(in) :: n_sites
      real(dp), intent(in) :: x0, y0
      real(dp), intent(in) :: siteX(n_sites), siteY(n_sites)
      type(MeshSettings), intent(in) :: OBJsettings
      type(ParamRefinement), INTENT(IN) :: OBJparamRefi

      !-----------------------------------------
      ! Abrir archivo
      !-----------------------------------------
      open (newunit=iu, file=trim(outdir)//'control.dat', status='replace', action='write', iostat=ios)
      if (ios /= 0) stop "Error opening control.dat"

      !-----------------------------------------
      ! Escribir encabezado
      !-----------------------------------------
      write (iu, 101) 'CENTER'
      write (iu, '(3(F3.1, 1x))') 0.0d0, 0.0d0, 0.0d0

      write (iu, 101) 'ROTATION'
      write (iu, '(F3.1)') OBJparamRefi%rotation

      write (iu, 101) 'NUM_THREADS'
      write (iu, '(I0)') 8

      write (iu, 101) 'SURF_MESH'
      write (iu, 101) 'ELLIPSOIDS'

      !-----------------------------------------
      ! Escribir elipsoides
      !-----------------------------------------
      call surface_ellipsoids(x0, y0, siteX, siteY, n_sites, OBJsettings, OBJparamRefi, iu)

      !-----------------------------------------
      ! Bloque INTERPOLATE y otros
      !-----------------------------------------
      write (iu, 101) 'INTERPOLATE'
      write (iu, '(F0.2)') 10.0d0
      write (iu, '(I0)') 4
      write (iu, '(ES0.4)') 1.0d-6

      write (iu, 101) 'ALTITUDE'
      write (iu, 101) OBJsettings%topography_file
      write (iu, '(F5.3)') 0.0d0
      write (iu, '(ES0.4)') 1.0d10

      write (iu, 101) 'SEA_DEPTH'
      write (iu, 101) OBJsettings%bathymetry_file
      write (iu, '(F5.3)') 0.0d0
      write (iu, '(ES0.4)') 1.0d10

      write (iu, 101) 'END'

      !-----------------------------------------
      ! Cerrar archivo
      !-----------------------------------------
      close (iu)
101   format(A)

   end subroutine control_mesh
!=========================================================
!=======
!=========================================================
   subroutine surface_ellipsoids(x_0, y_0, sitex, sitey, Nsites, OBJsettings, OBJrefiParam, unit)

      implicit none

      type(ParamRefinement), intent(in) :: OBJrefiParam
      type(MeshSettings), INTENT(IN) :: OBJsettings
      ! Inputs
      integer, intent(in) :: Nsites, unit
      real(dp), intent(in) :: x_0, y_0
      real(dp), intent(in) :: sitex(Nsites), sitey(Nsites)

      ! Local variables
      integer :: ki
      real(dp) :: max_r, dist, current_len, domain_r
      real(dp) :: a1, a_last, a_ratio, current_a, fh, fvp, fvm, t
      real(dp) :: fh_core, fvp_core, fvm_core, fh_final, fvp_final, fvm_final

      ! ----------------------------------------------------------
      ! 1) Calcular radio máximo desde el centro a las estaciones
      ! ----------------------------------------------------------

      print *, 'esto es siteX y siteY antes del loop'
      do ki = 1, 4
         print'(2f15.5)', sitex(ki), sitey(ki)
      end do
      max_r = 0.0d0
      do ki = 1, n_sites
         dist = sqrt((sitex(ki))**2 + (sitey(ki))**2)
         if (dist > max_r) max_r = dist
      end do
      print'(A,f15.5)', 'esto es DIST', dist

      ! ----------------------------------------------------------
      ! 2) Definir primer y último radio de elipsoide
      ! ----------------------------------------------------------

      ! Radio del core (zona survey + padding)
      print'(A,f15.5)', 'esto es max_r', max_r
      a1 = max_r + OBJrefiParam%paddingRefi

      ! Radio del dominio (distancia máxima desde centro al borde)
      domain_r = max(abs(OBJsettings%xminDOM - x_0), abs(OBJsettings%xmaxDOM - x_0), &
                     abs(OBJsettings%yminDOM - y_0), abs(OBJsettings%ymaxDOM - y_0))

      print'(A,f15.5)', 'esto es a1', a1
      print *, ' '
      a_last = domain_r
      print'(A,f15.5)', ' a_last', domain_r

      ! Evitar degeneración si levels = 1
      if (OBJrefiParam%Nelipses > 1) then
         a_ratio = (a_last/a1)**(1.0d0/(OBJrefiParam%Nelipses - 1))
      else
         a_ratio = 1.0d0
      end if

      ! ----------------------------------------------------------
      ! 3) Inicializar tamaño de elemento (len)
      ! ----------------------------------------------------------
      current_len = OBJrefiParam%coreResol

      ! ----------------------------------------------------------
      ! 4) Generar elipsoides jerárquicos
      ! ----------------------------------------------------------
      current_a = a1
      fh_core = 0.5d0
      fvp_core = 0.7d0
      fvm_core = 0.9d0

      fh_final = 0.0d0
      fvp_final = 0.0d0
      fvm_final = 0.0d0

      write (unit, '(I0)') OBJrefiParam%Nelipses
      a_ratio = 1.3
      do ki = 1, OBJrefiParam%Nelipses
         t = real(ki - 1, dp)/real(OBJrefiParam%Nelipses - 1, dp)

         fh = (1.0d0 - t)*fh_core + t*fh_final
         fvp = (1.0d0 - t)*fvp_core + t*fvp_final
         fvm = (1.0d0 - t)*fvm_core + t*fvm_final

         write (unit, 10) current_a, current_len, fh, fvp, fvm

         ! Actualizar radio geométrico
         current_a = current_a*a_ratio

         ! Actualizar tamaño de elemento con cap
         current_len = min(current_len*OBJrefiParam%growthFactor, &
                           OBJrefiParam%sizeBoundary)

      end do

10    format(F7.2, 1X, F6.2, 1X, F6.3, 1X, F6.3, 1X, F6.3)
   end subroutine surface_ellipsoids

   ! subroutine surface_ellipsoids(x0, y0, site_x, site_y, n_sites, config)
   !     ! ... declaraciones de variables ...
   !     real :: max_r, current_a, current_len
   !     integer :: i
   !
   !     ! 1. Calcular el radio máximo desde el centro (x0, y0) a la estación más lejana
   !     max_r = 0.0
   !     do i = 1, n_sites
   !         dist = sqrt((site_x(i)-x0)**2 + (site_y(i)-y0)**2)
   !         if (dist > max_r) max_r = dist
   !     end do
   !
   !     ! 2. Radio inicial de la elipse core
   !     current_a = max_r + config%core_radius_padding
   !     current_len = config%core_resolution
   !
   !     ! Escribir encabezado del control.dat (puntos de rotación, etc)
   !     write(unit, *) 0.0  ! Rotation
   !     write(unit, *) config%levels ! Número de elipses
   !
   !     ! 3. Bucle para generar las elipses (Matrioshka de superficie)
   !     do i = 1, config%levels
   !         ! Escribimos: a, len, fh, fvp, fvm
   !         ! Usamos fvm=0.99 para que sea un disco plano en superficie
   !         write(unit, '(F10.2, F10.2, F10.2, F10.2, F10.2)') &
   !             current_a, current_len, 0.0, 0.0, 0.99
   !
   !         ! Actualizamos para la siguiente capa hacia afuera
   !         current_a = current_a * 1.5   ! Expandimos el radio
   !         current_len = min(current_len * config%growth_factor, config%boundary_resolution)
   !     end do
   ! end subroutine
!=========================================================
!=======
!=========================================================
   subroutine read_dem(OBJsettings, x, y, z, n)

      use mesh_config
      implicit none

      type(MeshSettings), intent(in) :: OBJsettings
      real(dp), allocatable, intent(out):: x(:), y(:), z(:)
      integer, intent(out):: n

      integer:: iu, i
      real(dp):: xx, yy, zz
      real(dp):: scale

      scale = 1.0_dp
      if (OBJsettings%dem_units == 'meters') scale = 1.0d-3

      open (newunit=iu, file='preprocessing/DEM/'//OBJsettings%dem_file, status='old')
      n = 0
      do
         read (iu, *, end=10)
         n = n + 1
      end do
10    rewind (iu)

      allocate (x(n), y(n), z(n))
      do i = 1, n
         read (iu, *) xx, yy, zz
         ! x(i) = xx*scale
         ! y(i) = yy*scale
         x(i) = yy*scale   ! North → X ; applying MT convention
         y(i) = xx*scale   ! East  → Y ; applying MT convention
         z(i) = zz*scale
      end do
      close (iu)
   end subroutine
!=========================================================
!=======
!=========================================================
   subroutine check_domain(x, y, OBJsettings)!xmin, xmax, ymin, ymax, pad_x, pad_y)
      implicit none
      type(MeshSettings), INTENT(IN) :: OBJsettings
      real(dp), intent(in):: x(:), y(:)
      ! real(dp), intent(in):: xmin, xmax, ymin, ymax, pad_x, pad_y
      real(dp)            :: xminDEM, yminDEM, xmaxDEM, ymaxDEM, dem_size_x, dem_size_y
      real(dp)            :: condition_1, condition_2, condition_3, condition_4, x1, x2, y1, y2

      xminDEM = minval(x)
      xmaxDEM = maxval(x)
      yminDEM = minval(y)
      ymaxDEM = maxval(y)
      dem_size_x = xmaxDEM - xminDEM
      dem_size_y = ymaxDEM - yminDEM

      x1 = OBJsettings%xminDOM - OBJsettings%pad_x
      x2 = OBJsettings%xmaxDOM + OBJsettings%pad_x
      y1 = OBJsettings%yminDOM - OBJsettings%pad_y
      y2 = OBJsettings%ymaxDOM + OBJsettings%pad_y

      print *, 'DEM X range      :', xminDEM, xmaxDEM
      print *, 'DEM Y range      :', yminDEM, ymaxDEM
      print *, 'Domain X (+pad)  :', x1, x2
      print *, 'Domain Y (+pad)  :', y1, y2
      print *, 'DEM size (km)    :', dem_size_x, dem_size_y

      condition_1 = OBJsettings%xminDOM - OBJsettings%pad_x
      condition_2 = OBJsettings%xmaxDOM + OBJsettings%pad_x
      condition_3 = OBJsettings%yminDOM - OBJsettings%pad_y
      condition_4 = OBJsettings%ymaxDOM + OBJsettings%pad_y

      if (minval(x) > condition_1 .or. maxval(x) < condition_2) stop 'DEM X does not cover domain+padding'
      if (minval(y) > condition_3 .or. maxval(y) < condition_4) stop 'DEM Y does not cover domain+padding'
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
      real(dp), intent(out)          :: value

      character(len=:), allocatable:: s
      integer:: ipos, p1, p2
      real(dp):: deg, min, sec
      real(dp):: sign

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
         "ls preprocessing/edi_files/*.edi > "//tmpfile, &
         exitstat=ios)

      if (ios /= 0) then
         write (error_unit, *) "ERROR: failed to list preprocessing/edi_files/*.edi"
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
      real(dp), intent(out)         :: edi_lat(n), edi_lon(n), edi_elev(n)
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
      real(dp), intent(in)  :: lat(n_files), lon(n_files)      ! grados decimales
      real(dp), intent(out) :: x(n_files), y(n_files)           ! metros
      integer, intent(out) :: zone

      integer:: i
      real(dp):: lon_mean, east_km(n_files), north_km(n_files)

      ! --- constantes---
      real(dp), parameter:: semieje = 6378137.0d0
      real(dp), parameter:: f = 1.0d0/298.257223563d0
      real(dp), parameter:: k0 = 0.9996d0
      real(dp), parameter:: pi = 3.141592653589793d0

      real(dp):: e2, ep2
      real(dp):: latr(n_files), lonr(n_files), lon0r
      real(dp):: N(n_files), T(n_files), C(n_files), A(n_files), M(n_files)
      real(dp):: lon0

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
      ! x = x/1000.0d0
      ! y = y/1000.0d0
      !
      east_km = x/1000.0d0
      north_km = y/1000.0d0

      x = north_km   ! X = Norte
      y = east_km    ! Y = Este

   end subroutine edi_to_utm
!=========================================================
!=======
!=========================================================
   subroutine compute_center(x, y, n, x0, y0)
      real(dp), intent(in)  :: x(n), y(n)
      real(dp), intent(out):: x0, y0
      integer, intent(in)  :: n

      x0 = sum(x)/dble(n)
      y0 = sum(y)/dble(n)
   end subroutine
!=========================================================
!=======
!=========================================================
   subroutine recenter_all(n_sites, n_dem, x0, y0, site_x, site_y, dem_x, dem_y)
      real(dp), intent(inout):: site_x(n_sites), site_y(n_sites), dem_x(n_dem), dem_y(n_dem)
      real(dp), intent(in)    :: x0, y0
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
      real(dp), intent(in)  :: siteXm(nn_sites), siteYm(nn_sites), demXmts(nn_dem), demYmts(nn_dem), demZmts(nn_dem)
      character(len=*), intent(in)  :: ediID(nn_sites)
      integer, intent(in)  :: nn_sites, nn_dem
      real(dp), intent(out):: siteXkm(nn_sites), siteYkm(nn_sites), demXkm(nn_dem), demYkm(nn_dem)
      character(len=20) :: dir
      integer:: i, iu

      siteXkm = siteXm
      siteYkm = siteYm

      demXkm = demXmts
      demYkm = demYmts

      dir = 'preprocessing/'
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
   subroutine snap_sites_to_dem(ediID, siteX, siteY, Nsites, DEMcorX, DEMcorY, DEMcorZ, NDEM, siteZ)

      use mesh_config
      implicit none
      integer, intent(in):: Nsites, NDEM
      real(dp), intent(in):: siteX(Nsites), siteY(Nsites)
      real(dp), intent(in):: DEMcorX(NDEM), DEMcorY(NDEM), DEMcorZ(NDEM)
      real(dp), intent(out):: siteZ(Nsites)
      character(len=*) :: ediID(Nsites)

      character(len=512):: fname
      integer:: ii, j, jmin, iu
      real(dp):: dx, dy, d2, d2min
      real(dp):: zmin

      ! Protección mínima: elevación mínima del DEM
      zmin = minval(DEMcorZ)
      if (zmin < 1.0d0) zmin = 1.0d0   ! nunca 0 ni negativa

      do ii = 1, Nsites
         d2min = huge(1.0d0)
         jmin = 1

         do j = 1, NDEM
            dx = DEMcorX(j) - siteX(ii)
            dy = DEMcorY(j) - siteY(ii)
            d2 = dx*dx + dy*dy

            if (d2 < d2min) then
               d2min = d2
               jmin = j
            end if
         end do

         siteZ(ii) = DEMcorZ(jmin)

         ! seguridad extra
         if (siteZ(ii) < zmin) siteZ(ii) = zmin
      end do

      ! if (OBJsetting%dem_units == 'kilometer') then
      siteZ = siteZ/1000.0d0
      ! end if
      ! Writing final coordinate site files
      fname = trim(outdir)//'sites_coord_elev.dat'
      open (newunit=iu, file=fname, status='replace', action='write')
      do ii = 1, Nsites
         write (iu, '(A12, 3f15.5)') ediID(ii), siteX(ii), siteY(ii), siteZ(ii)
      end do

      close (iu)

   end subroutine
!=========================================================
!=======
!=========================================================
   subroutine write_observing_sites(site_x, site_y, n_sites, OBJparamRefi)

      implicit none

      type(ParamRefinement), intent(in) :: OBJparamRefi
      integer, intent(in) :: n_sites
      real(dp), intent(in) :: site_x(n_sites), site_y(n_sites)

      integer :: iu, i, k
      real(dp) :: r_ratio, e_ratio, current_r, current_e
      character(len=512) :: fname

      ! 1. Calculamos los factores de crecimiento (Geometric Progression)
      ! Si n=5, calculamos el ratio para llegar de min a max en n-1 pasos
      r_ratio = (OBJparamRefi%maxRAD/OBJparamRefi%minRAD)**(1.0d0/(OBJparamRefi%Nsph - 1))

      ! El e_ratio asegura que la última esfera tenga el tamaño de la malla general
      e_ratio = (OBJparamRefi%coreResol/OBJparamRefi%minEDGE)**(1.0d0/(OBJparamRefi%Nsph - 1))

      fname = trim(outdir)//'/observing_site.dat'
      open (newunit=iu, file=fname, status='replace', action='write')

      write (iu, '(I9)') n_sites

      do i = 1, n_sites
         ! Escribimos coordenadas y número de esferas
         write (iu, '(2(F12.6, 1X))', advance='yes') site_x(i), site_y(i)
         write (iu, '(I0)') OBJparamRefi%Nsph

         current_r = OBJparamRefi%minRAD
         current_e = OBJparamRefi%minEDGE

         do k = 1, OBJparamRefi%Nsph
            write (iu, '(F7.3, 2x, F7.4)', advance='yes') current_r, current_e
            ! Actualizamos para la siguiente capa
            current_r = current_r*r_ratio
            current_e = current_e*e_ratio
         end do
         write (iu, *)
      end do

      close (iu)
   end subroutine write_observing_sites

!=========================================================
!=======
!=========================================================
   subroutine define_analysis_domain(settings)
      implicit none
      type(MeshSettings), INTENT(IN) :: settings

      integer:: iu

      ! xmin =  - pad_x
      ! xmax = maxval(x) + pad_x
      ! ymin = minval(y) - pad_y
      ! ymax = maxval(y) + pad_y

      open (newunit=iu, file=trim(outdir)//"analysis_domain.dat", status='replace', action='write')
      write (iu, '(2F10.2)') settings%xminDOM, settings%xmaxDOM
      write (iu, '(2F10.2)') settings%yminDOM, settings%ymaxDOM
      write (iu, '(2F10.2)') settings%zminDOM, settings%zmaxDOM
      close (iu)

   end subroutine
!=========================================================
!=======
!=========================================================
   subroutine write_topography(OBJsettings, x, y, z, n)
      implicit none
      type(MeshSettings), INTENT(IN) :: OBJsettings
      real(dp), intent(in)    :: x(:), y(:)
      integer, intent(in)     :: n
      real(dp), INTENT(INOUT) :: z(:)
      integer                 :: i, iu

      z = z/1000.0d0

      open (newunit=iu, file=trim(outdir)//OBJsettings%topography_file, status='replace', action='write')
      do i = 1, n
         if (.not. OBJsettings%has_sea) then
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
   subroutine write_bathymetry(OBJsettings, x, y, z, n)

      implicit none
      TYPE(MeshSettings), INTENT(IN) :: OBJsettings
      real(dp), intent(in):: x(:), y(:)
      real(dp), INTENT(INOUT) :: z(:)
      integer, intent(in):: n
      integer:: j, iu

      z = z/1000.0d0

      open (newunit=iu, file=trim(outdir)//OBJsettings%bathymetry_file, status='replace')
      do j = 1, n
         if (.not. OBJsettings%has_sea) then
            ! Caso SIN mar: no hay batimetría
            write (iu, '(3F15.5)') x(j), y(j), -1.0_dp/1000.0d0
         else
            ! Caso CON mar
            if (z(i) < OBJsettings%sea_level) then
               write (iu, '(3F15.5)') x(j), y(j), -z(j)
            else
               write (iu, '(3F15.5)') x(j), y(j), -1.0_dp/1000.0d0
            end if
         end if
      end do
      close (iu)
   end subroutine
!=========================================================
!=======
!=========================================================
   subroutine write_coast_line(OBJsettings)

      implicit none
      type(MeshSettings), INTENT(IN) :: OBJsettings
      integer:: iu

      if (.not. OBJsettings%has_sea) then
         open (newunit=iu, file=trim(outdir)//OBJsettings%coastLine_file, status='replace')
         write (iu, *) 1
         write (iu, '(2(F21.15, 1X), 2I2)') OBJsettings%xminDOM - 5.0d0, OBJsettings%yminDOM - 5.0d0, 0, 0
         write (iu, '(2(F21.15, 1X), 2I2)') OBJsettings%xmaxDOM + 5.0d0, OBJsettings%yminDOM - 5.0d0, 0, 0
         write (iu, '(2(F21.15, 1X), 2I2)') OBJsettings%xmaxDOM + 5.0d0, OBJsettings%ymaxDOM + 5.0d0, 0, 0
         write (iu, '(2(F21.15, 1X), 2I2)') OBJsettings%xminDOM - 5.0d0, OBJsettings%ymaxDOM + 5.0d0, 1, 0
         close (iu)
      else
         stop 'Sea case not implemented yet'
      end if
   end subroutine
!=========================================================
!=======
!=========================================================
   subroutine run_makeTetraMesh_and_assign_regions(OBJmodRegions)
      implicit none

      type(ModelRegion), INTENT(IN) :: OBJmodRegions
      ! integer, intent(in) :: N_regions, ID_regions(N_regions)
      ! real(8), intent(in) :: coord_regions(3, N_regions)

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

      call execute_command_line('cd preprocessing/buildMesh && cp ../geometry/* .', wait=.true., exitstat=stat)
      if (stat /= 0) stop 'ERROR: copying geometry failed'
      call execute_command_line('echo " "')

      ! -----------------------------
      ! 3. Ejecutar steps 1–4
      ! -----------------------------
      call execute_command_line('echo "step 1 --> Defining Land/Sea Boundary"')
      call execute_command_line('cd preprocessing/buildMesh && makeTetraMesh -stp 1', wait=.true., exitstat=stat)
      if (stat /= 0) then
         error stop 'ERROR: makeTetraMesh step 1 failed'
      end if
      call execute_command_line('echo "done..." && sleep 1')

      call execute_command_line('echo "step 2 --> Building 2D mesh..."')
      call execute_command_line('cd preprocessing/buildMesh && makeTetraMesh -stp 2', wait=.true., exitstat=stat)
      if (stat /= 0) then
         error stop 'ERROR: makeTetraMesh step 2 failed'
      end if
      call execute_command_line('echo "done..." && sleep 1')

      call execute_command_line('echo "step 3 --> Interpolating altitudes"')
      call execute_command_line('cd preprocessing/buildMesh && makeTetraMesh -stp 3', wait=.true., exitstat=stat)
      if (stat /= 0) then
         error stop 'ERROR: makeTetraMesh step 3 failed'
      end if
      call execute_command_line('echo "done..." && sleep 1')

      call execute_command_line('echo "step 4 --> Making surface mesh"')
      call execute_command_line('cd preprocessing/buildMesh && makeTetraMesh -stp 4', wait=.true., exitstat=stat)
      if (stat /= 0) then
         error stop 'ERROR: makeTetraMesh step 4 failed'
      end if
      call execute_command_line('echo "done..." && sleep 1')

      ! -----------------------------
      ! 4. Verificar output.poly
      ! -----------------------------
      inquire (file='preprocessing/buildMesh/output.poly', exist=ex)
      if (.not. ex) stop 'ERROR: output.poly was not generated'

      ! -----------------------------
      ! 5. Parchear regiones
      ! -----------------------------

      ! after_part4 = .false.
      found_part4 = .false.

      call execute_command_line('echo " " ')
      call execute_command_line('echo "Assigning regions in output.poly file"')
      open (newunit=iu_in, file='preprocessing/buildMesh/output.poly', status='old')
      open (newunit=iu_out, file='preprocessing/buildMesh/output.poly.tmp', status='replace')

      !
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
            write (iu_out, '(I0)') OBJmodRegions%Nregions

            do k = 1, OBJmodRegions%Nregions
               write (iu_out, '(I3,3F10.3,I4,1PE12.4)') k, OBJmodRegions%coord(:, k), OBJmodRegions%ID(k), OBJmodRegions%rho(k)
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

      call execute_command_line('mv preprocessing/buildMesh/output.poly.tmp preprocessing/buildMesh/output.poly')

      write (*, *) 'Ok: ✅ makeTetraMesh steps 1–4 done and regions added to output.poly'

   end subroutine run_makeTetraMesh_and_assign_regions
!=========================================================
!=======
!=========================================================
   pure real(dp) function skin_depth_km(rho_ohmm, f_hz) result(delta)
      real(dp), intent(in) :: rho_ohmm, f_hz

      ! Aproximación común: δ[km] ≈ 500 * sqrt(rho/f)
      if (rho_ohmm <= 0.0_dp .or. f_hz <= 0.0_dp) then
         delta = 0.0_dp
      else
         delta = 500.0_dp*sqrt(rho_ohmm/f_hz)
         delta = delta/1000.0d0
      end if

   end function skin_depth_km
!=========================================================
!=======
!=========================================================
   subroutine write_makeMtr_param(xcenter, ycenter, site_x, site_y, Nsites, OBJrefiParam, cfg)

      implicit none

      type(GlobalRefinement), intent(in)  :: cfg
      type(ParamRefinement), INTENT(IN)   :: OBJrefiParam
      integer, intent(in)                 :: Nsites
      real(dp), intent(in)                :: xcenter, ycenter
      real(dp), intent(in)                :: site_x(Nsites), site_y(Nsites)

      integer :: iu, i, ios
      real(dp) :: z0
      real(dp) :: dist, max_r
      real(dp) :: a1, ai
      real(dp) :: max_elem_size
      real(dp) :: depth_min, depth_last
      real(dp) :: len_i
      real(dp) :: fh, fvp, fvm, fh_final, fvp_final, depth_ratio_core, depth_ratio_trans
      real(dp) :: relax_factor, z_core, z_relax_end, n_core, n_trans
      integer :: n_core_i, n_trans_i
      real(dp) :: a_core_end, a_last, a_ratio_core, a_ratio_far, step
      real(dp) :: len_mult(5)
      real(dp) :: fh_tab(5), fvp_tab(5), fvm_tab(5)
      character(len=512) :: fname

      z0 = 0.0_dp

      if (cfg%levels < 2) stop "ERROR: mesh.volume.levels debe ser >= 2"

      ! ----------------------------------------------------------
      ! 1) distancia maxima de un site al centro del dominio
      ! ----------------------------------------------------------
      max_r = -1.0_dp
      do i = 1, Nsites
         dist = sqrt(site_x(i)**2 + site_y(i)**2)
         if (dist > max_r) max_r = dist
      end do

      ! ----------------------------------------------------------
      ! 2) a1 (radio core horizontal)
      ! ----------------------------------------------------------
      a1 = max_r + cfg%horizontal_padding

      ! Evitar geometría imposible: necesitamos ai > depth_i
      ! (más adelante ajustamos ai si hiciera falta, pero a1 debe ser razonable)
      a1 = max(a1, 1.2_dp*(0.3_dp*cfg%target_depth))  ! usando depth_min interno

      ! ----------------------------------------------------------
      ! 3) cap físico para len (por input file)
      ! ----------------------------------------------------------
      max_elem_size = cfg%farfield_resolution

      ! 4) PROFUNDIDADES POR NIVEL (V1 SIMPLE, INTERNO)

      ! ----------------------------------------------------------
      ! 4) cap físico para len (por input file)
      ! ----------------------------------------------------------
      relax_factor = 2.0_dp   ! transición hasta 2*z_core (puedes ajustar interno)

      ! --- 4.1) Profundidad core controlada por HRL
      depth_min = cfg%depth_min_factor*cfg%target_depth   ! mínimo interno (km)
      z_core = max(0.5_dp, depth_min)

      ! --- 4.2) Extensión total del dominio
      depth_last = cfg%vertical_padding_factor*cfg%target_depth

      ! --- 4.3) Fin de transición vertical
      ! z_relax_end = depth_last
      z_relax_end = min(depth_last, max(relax_factor*z_core, 0.7*depth_last))

      ! --- 4.4) Reparto de niveles
      n_core = real(cfg%levels/2, dp)          ! mitad niveles core
      n_trans = real(cfg%levels - n_core, dp)   ! mitad transición

      if (n_core < 1.0_dp) n_core = 1.0_dp
      if (n_trans < 1.0_dp) n_trans = 1.0_dp

      ! --- 5) Ratios geométricos independientes
      depth_ratio_core = (z_core/0.5_dp)**(1.0_dp/max(n_core - 1.0_dp, 1.0_dp))
      depth_ratio_trans = (depth_last/z_core)**(1.0_dp/max(n_trans, 1.0_dp))

      ! ----------------------------------------------------------
      ! 6) fh/fvp: relajación lineal a 0
      ! ----------------------------------------------------------
      fh_final = 0.0_dp
      fvp_final = 0.0_dp

      ! --- Core: 5 niveles como el autor (o menos si levels es chico)
      n_core_i = min(5, cfg%levels - 2)   ! deja al menos 2 niveles para trans+far
      if (n_core_i < 1) n_core_i = 1

      ! --- Transición: 1 nivel (simple y efectivo)
      n_trans_i = 1
      if (cfg%levels - n_core_i <= 1) n_trans_i = 0   ! si no hay espacio

      ! Tablas tipo autor (se recortan si n_core_i < 5)
      len_mult = (/0.35_dp, 1.15_dp, 1.8_dp, 3.0_dp, 5.0_dp/)

      ! fh_tab = (/0.5_dp, 0.5_dp, 0.5_dp, 0.3_dp, 0.2_dp/)    !----> refinamiento algo estrecho en y, y cargado en hacia direccion x
      fh_tab = (/0.3_dp, 0.2_dp, 0.1_dp, 0.05_dp, 0.0_dp/)      !----> caso de refinamiento casi circular en xy plane
      ! fh_tab = (/0.8_dp, 0.7_dp, 0.5_dp, 0.5_dp, 0.4_dp/)      !----> caso de refinamiento estecho en y y marcado en x

      !tierra
      fvp_tab = (/0.3_dp, 0.3_dp, 0.2_dp, 0.1_dp, 0.1_dp/)

      !aire
      fvm_tab = (/0.7_dp, 0.6_dp, 0.5_dp, 0.4_dp, 0.3_dp/)

      ! Core A suave: usa ratio fijo suave (parecido al autor)
      a_ratio_core = 1.15_dp   ! puedes dejarlo fijo interno por ahora

      ! Donde termina el core en radio:
      a_core_end = a1*a_ratio_core**real(n_core_i - 1, dp)

      ! A_last: el último radio del archivo (aquí sí usa tu crecimiento global)
      a_last = a1*cfg%a_ratio**real(cfg%levels - 1, dp)

      ! Asegura estrictamente creciente (muy importante)
      if (a_last <= a_core_end) a_last = 1.05_dp*a_core_end

      ! Ratio farfield para llegar a a_last con los niveles restantes
      if (cfg%levels > n_core_i + n_trans_i) then
         a_ratio_far = (a_last/a_core_end)**(1.0_dp/real(cfg%levels - (n_core_i + n_trans_i), dp))
      else
         a_ratio_far = 1.0_dp
      end if

      ! ----------------------------------------------------------
      ! 7) Escribir makeMtr.param
      ! ----------------------------------------------------------
      fname = trim(outdir)//'/makeMtr.param'
      open (newunit=iu, file=fname, status='replace', action='write', iostat=ios)
      if (ios /= 0) stop "ERROR: no se pudo abrir makeMtr.param"

      write (iu, '(3(F4.2,1x))') xcenter, ycenter, z0
      write (iu, '(F4.2)') OBJrefiParam%rotation
      write (iu, '(I0)') cfg%levels

      do i = 1, cfg%levels

         if (i <= n_core_i) then
            ! -------- CORE (patrón autor) --------
            ai = a1*a_ratio_core**real(i - 1, dp)

            ! growth = cfg%vol_len_factor * len_mult(i)
            ! len_i = min(cfg%core_resolution*growth, max_elem_size)

            len_i = min(cfg%core_resolution*len_mult(i), max_elem_size)

            fh = fh_tab(i)
            fvp = fvp_tab(i)*1.6
            fvm = fvm_tab(i)

         else if (i <= n_core_i + n_trans_i) then
            ! -------- TRANSICIÓN (apaga anisotropía) --------
            ! ai = a_core_end * a_ratio_far   ! primer paso fuera del core (simple)
            step = i - n_core_i
            ai = a_core_end*a_ratio_far**step

            ! aquí puedes escoger: o un len intermedio, o ya dejar que crezca
            len_i = min(max(len_i, cfg%core_resolution*10.0_dp), max_elem_size)

            !    fh = (1.0_dp - t)*cfg%fh_core_vol + t*fh_final
            !    fh = max(0.0_dp, min(fh, 0.99_dp))
            fh = 0.0_dp
            fvp = 0.0_dp
            fvm = 0.0_dp

         else
            ! -------- FARFIELD (todo 0) --------
            ai = a_core_end*a_ratio_far**real(i - (n_core_i + n_trans_i) + 1.25, dp)

            ! crecimiento suave de len hasta cap
            len_i = min(len_i*cfg%len_growth, max_elem_size)
            print'(A,F15.5)', 'growth', len_i*cfg%len_growth
            print'(A,F15.5)', 'max_elem_size o FAR_FIEL_RESOL', max_elem_size

            !    fh = (1.0_dp - t)*cfg%fh_core_vol + t*fh_final
            !    fh = max(0.0_dp, min(fh, 0.99_dp))
            fh = 0.0_dp
            fvp = 0.0_dp
            fvm = 0.0_dp
         end if

         ! clamp seguridad
         fh = max(0.0_dp, min(fh, 0.99_dp))
         fvp = max(0.0_dp, min(fvp, 0.99_dp))
         fvm = max(0.0_dp, min(fvm, 0.99_dp))

         write (iu, 10) ai, len_i, fh, fvp, fvm
      end do

      close (iu)
10    format(F8.2, 2X, F6.2, 2X, F5.3, 1X, F6.3, 1X, F6.3)

   end subroutine write_makeMtr_param
!=========================================================
!=======
!=========================================================
   subroutine write_obs_sites(site_coorX, site_coorY, site_coorZ, nSites, OBJparamRefi)

      implicit none
      integer, intent(in) :: nSites
      real(dp), intent(in) :: site_coorX(nSites), site_coorY(nSites), site_coorZ(nSites)
      type(ParamRefinement), intent(in) :: OBJparamRefi

      integer :: iu, ij, k
      real(dp) :: r_ratio, e_ratio
      real(dp) :: r_k, L_k
      real(dp) :: t, obl_k, oblateness_core

      character(len=512) :: fname

      oblateness_core = 0.3_dp

      fname = trim(outdir)//'/obs_site.dat'
      open (newunit=iu, file=fname, status='replace', action='write')

      write (iu, '(I0)') nSites

      print *, ' Prints en obs_site'
      ! ---- Growth ratios ----
      if (OBJparamRefi%Nsph > 1) then
         print *, OBJparamRefi%maxRAD
         r_ratio = (OBJparamRefi%maxRAD/OBJparamRefi%minRAD)** &
                   (1.0_dp/real(OBJparamRefi%Nsph - 1, dp))

         e_ratio = (OBJparamRefi%coreResol/OBJparamRefi%minEDGE)** &
                   (1.0_dp/real(OBJparamRefi%Nsph - 1, dp))
         print *, 'entra al if'

         print'(A,f15.5)', 'esto es r_ratio en if', r_ratio
         print'(A,f15.5)', 'esto es e_ratio en if', e_ratio
      else
         r_ratio = 1.0_dp
         e_ratio = 1.0_dp
      end if
      print *, r_ratio
      print *, e_ratio

      do ij = 1, nSites

         write (iu, '(3(F14.9,1x))') site_coorX(ij), site_coorY(ij), site_coorZ(ij)
         write (iu, '(I0)') OBJparamRefi%Nsph

         r_k = OBJparamRefi%minRAD
         L_k = OBJparamRefi%minEDGE

         do k = 1, OBJparamRefi%Nsph

            t = real(k - 1, dp)/real(OBJparamRefi%Nsph - 1, dp)

            obl_k = (1.0_dp - t)*oblateness_core
            obl_k = max(0.0_dp, min(obl_k, 0.99_dp))

            write (iu, '(3(F5.2,1x))') r_k, L_k, obl_k

            r_k = r_k*r_ratio
            L_k = L_k*e_ratio

         end do

      end do

      close (iu)

   end subroutine write_obs_sites

   ! subroutine setGlobalMeshRefinement(xcenter, ycenter, Nglob_ellipses, corXsite, corYsite, nSites, OBJglobRefi)
   !
   !    implicit none
   !    type(GlobalRefinement), INTENT(INOUT) :: OBJglobRefi
   !    integer, intent(in)   :: Nglob_ellipses
   !    real(dp), intent(in)  :: ycenter, xcenter!, rot_glob
   !    character(len=512)    :: fname
   !    real(dp)              :: z0, oblatness_on_ZXplane
   !    integer               :: iu, i, nSites, k
   !    ! Parámetros hardcodeados (pueden ser arreglos en el futuro)
   !    real(dp) :: a(Nglob_ellipses), lengths(Nglob_ellipses), fh(Nglob_ellipses), fvp(Nglob_ellipses), fvm(Nglob_ellipses)
   !    ! real(dp), intent(in) :: len_alon_x_axis(n_site_ellipses), max_edge_len_within_ellipse(n_site_ellipses)
   !    real(dp), intent(in) :: corXsite(nSites), corYsite(nSites)
   !
   !    z0 = 0.0d0
   !
   !    ! Datos del ejemplo del autor
   !    a = (/40.0, 45.0, 50.0, 60.0, 80.0, 100.0, 200.0, 300.0, 400.0, 500.0/)
   !    lengths = (/1.0, 1.5, 3.0, 5.0, 8.0, 10.0, 15.0, 20.0, 30.0, 45.0/)
   !    fh = (/0.5, 0.5, 0.5, 0.3, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0/)
   !    fvp = (/0.7, 0.5, 0.4, 0.3, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0/)
   !    fvm = (/0.9, 0.7, 0.7, 0.5, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0/)
   !
   !    a = a*1.0d0
   !
   !    fname = trim(outdir)//'makeMtr.param'
   !    open (newunit=iu, file=fname, status='replace', action='write')
   !
   !    ! 1. Coordenadas del centro (Y, X, Z)
   !    write (iu, '(3F12.5)') xcenter, ycenter, z0
   !
   !    ! 2. Ángulo de rotación
   !    write (iu, '(F12.2)') OBJglobRefi%rotation
   !
   !    ! 3. Número de elipsoides
   !    ! write (iu, '(I5)') Nglob_ellipses
   !    write (iu, '(I5)') OBJglobRefi%n_ellipses_site
   !
   !    ! 4. Bloque de elipsoides
   !    do i = 1, min(10, 10)
   !       write (iu, '(F6.1, F6.1, 1x, 3F6.2)') a(i), lengths(i), fh(i), fvp(i), fvm(i)
   !    end do
   !
   !    close (iu)
   !    print *, 'Archivo makeMtr.param creado exitosamente en ', trim(outdir)
   !
   !    oblatness_on_ZXplane = 0.3
   !
   !    fname = trim(outdir)//'/obs_site.dat'
   !    open (newunit=iu, file=fname, status='replace', action='write', form='formatted')
   !
   !    ! ======================
   !    ! Write number of sites
   !    ! ======================
   !    write (iu, '(I4)') nSites
   !
   !    ! ======================
   !    ! Write each observation site
   !    ! ======================
   !    do i = 1, nSites
   !       ! X  Y  z0
   !       write (iu, '(3F15.6 )') corXsite(i), corYsite(i), z0
   !       write (iu, '(I1)') OBJglobRefi%n_ellipses_site
   !
   !       ! (Ri, Li) pairs
   !       do k = 1, OBJglobRefi%n_ellipses_site
   !          write (iu, '(3F5.2 )') &
   !             OBJglobRefi%lenEllipseSite(k), OBJglobRefi%maxSiteEdge(k), oblatness_on_zxplane
   !       end do
   !    end do
   !
   !    close (iu)
   !    print *, 'Archivo obs_site.dat creado exitosamente en ', trim(outdir)
   !
   ! end subroutine setGlobalMeshRefinement
!=========================================================
!=======
!=========================================================
   subroutine PARAM_ELLIPSOIDS(Nsites, siteX, siteY, OBJsettings, obj_modReg, OBJglobRefi, A, LEN, FH, FV)
      ! target_depth_km, fmin_hz, rho_ref,
      ! Ne, A, LEN, FH, FV)

      implicit none

      type(MeshSettings), INTENT(IN) :: OBJsettings
      type(ModelRegion), intent(in) :: obj_modreg
      type(GlobalRefinement), INTENT(IN) :: OBJglobRefi
      integer, intent(in) :: Nsites
      real(dp), intent(in) :: siteX(Nsites), siteY(Nsites)
      real(dp), allocatable, intent(out) :: A(:), LEN(:), FH(:), FV(:)

      ! --------------------------
      ! Local variables
      ! --------------------------
      integer :: k, Ne
      real(dp) :: max_r, dist
      real(dp) :: xmin, xmax, ymin, ymax
      real(dp) :: domain_r, rho_ref, fmin_hz
      real(dp) :: a1, a_last, a_ratio
      real(dp) :: len_core_param, len_far_param
      real(dp) :: g_param
      real(dp) :: fh_core, fv_core
      real(dp) :: fh_final, fv_final
      real(dp) :: delta_max, frac_delta
      real(dp) :: t, relax_t
      real(dp) :: relax_factor_param
      real(dp) :: depth_relax_start, depth_relax_end
      real(dp) :: len_cap_phys
      real(dp) :: depth_eff

      ! --------------------------
      ! Inicializaciones internas
      ! --------------------------

      Ne = obj_modreg%NpEllipses
      g_param = 1.8_dp
      fh_core = 0.4_dp
      fv_core = 0.7_dp
      fh_final = 0.0_dp
      fv_final = 0.0_dp
      frac_delta = OBJglobRefi%frac_skin_depth
      relax_factor_param = 2.0_dp

      allocate (A(Ne), LEN(Ne), FH(Ne), FV(Ne))

      ! ----------------------------------------------------------
      ! 1) Radio máximo del survey
      ! ----------------------------------------------------------
      max_r = 0.0_dp
      do k = 1, Nsites
         dist = sqrt(siteX(k)**2 + siteY(k)**2)
         if (dist > max_r) max_r = dist
      end do

      ! ----------------------------------------------------------
      ! 2) Radio del dominio desde LOCATION
      ! coord(:, :) viene de LOCATION xyz xyz
      ! ----------------------------------------------------------
      xmin = OBJsettings%xminDOM
      ymin = OBJsettings%xmaxDOM
      xmax = OBJsettings%yminDOM
      ymax = OBJsettings%ymaxDOM

      domain_r = max(sqrt(xmin**2 + ymin**2), &
                     sqrt(xmin**2 + ymax**2), &
                     sqrt(xmax**2 + ymin**2), &
                     sqrt(xmax**2 + ymax**2))

      ! ----------------------------------------------------------
      ! 3) Radios elipsoidales primera y ultima topada al 95% del maxDOMAIN
      ! ----------------------------------------------------------
      a1 = max_r + OBJglobRefi%horizontal_padding
      a_last = 0.95_dp*domain_r

      if (a_last <= a1) a_last = 1.05_dp*a1

      if (Ne > 1) then
         a_ratio = max((a_last/a1)**(1.0_dp/real(Ne - 1, dp)), 1.35_dp)
      else
         a_ratio = 1.0_dp
      end if

      rho_ref = OBJglobRefi%background_rho
      fmin_hz = OBJglobRefi%fmin_hz
      ! ----------------------------------------------------------
      ! 4) Skin depth con fmin
      ! δ(km) ≈ 503 * sqrt(rho/f) / 1000
      ! ----------------------------------------------------------
      delta_max = skin_depth_km(rho_ref, fmin_hz)

      print'(A,f15.5)', '1sto es delta_max', delta_max
      len_cap_phys = frac_delta*delta_max

      print'(A,f15.5)', 'Esto es len_cap_sys', len_cap_phys

      ! ----------------------------------------------------------
      ! 5) Tamaños de parámetro
      ! ----------------------------------------------------------
      len_core_param = 2.0_dp
      print'(A,f15.5)', 'Esto es domain_r', domain_r

      len_far_param = min(1.00_dp*domain_r, len_cap_phys)
      ! len_far_param  = min(len_far_param, len_cap_phys)
      !
      if (len_far_param < len_core_param) len_far_param = len_core_param

      ! ----------------------------------------------------------
      ! 6) Construccion de radios
      ! ----------------------------------------------------------
      depth_relax_start = OBJglobRefi%target_depth
      depth_relax_end = relax_factor_param*OBJglobRefi%target_depth

      do k = 1, Ne

         ! Radios
         if (k == 1) then
            A(k) = a1
            LEN(k) = len_core_param
         else
            A(k) = A(k - 1)*a_ratio
            LEN(k) = min(LEN(k - 1)*g_param, len_far_param)
         end if

         ! Interpolación horizontal
         if (Ne > 1) then
            t = real(k - 1, dp)/real(Ne - 1, dp)
         else
            t = 0.0_dp
         end if

         FH(k) = (1.0_dp - t)*fh_core + t*fh_final
         FH(k) = max(0.0_dp, min(FH(k), 0.99_dp))

         ! Relajación vertical basada en profundidad proxy
         depth_eff = OBJglobRefi%target_depth*sqrt(A(k)/a1)

         if (depth_eff <= depth_relax_start) then
            FV(k) = fv_core
         else if (depth_eff >= depth_relax_end) then
            FV(k) = 0.0_dp
         else
            relax_t = (depth_eff - depth_relax_start)/ &
                      (depth_relax_end - depth_relax_start)
            FV(k) = (1.0_dp - relax_t)*fv_core
         end if

         FV(k) = max(0.0_dp, min(FV(k), 0.99_dp))

      end do

   end subroutine PARAM_ELLIPSOIDS
!=========================================================
!=======
!=========================================================
   subroutine resistivitty_attribute(Nsites, coorXsite, coorYsite, coorZsite, OBJsettings, OBJglobRefi, OBJmodReg)

      implicit none

      type(ModelRegion), INTENT(IN) :: OBJmodReg
      type(MeshSettings), INTENT(IN) :: OBJsettings
      type(GlobalRefinement), INTENT(IN) :: OBJglobRefi
      ! ======================
      ! Inputs
      ! ======================
      integer, intent(in) :: Nsites!, N_regions, N_paramEsfer
      real(dp), intent(in) :: coorXsite(Nsites), coorYsite(Nsites), coorZsite(Nsites)
      ! real(dp), intent(in) :: esferEdges(N_paramEsfer), radiusEsfer(N_paramEsfer)
      ! real(dp), intent(in) :: rhoRegions(N_regions)
      ! integer, intent(in) :: repRegion(N_regions), fixed(N_regions), ID_regions(N_regions)

      ! ======================
      ! Local variables
      ! ======================
      integer :: iu, i, j, Ne
      character(len=512) :: fname

      ! Parámetros para la Parte 2 (Elipsoides Regionales)
      real(dp), DIMENSION(:), allocatable :: a_reg, len_reg, fh_reg, fv_reg

      call param_ellipsoids(Nsites, coorXsite, coorYsite, OBJsettings, OBJmodReg, OBJglobRefi, a_reg, len_reg, fh_reg, fv_reg)

      Ne = OBJmodReg%NpEllipses

      fname = trim(outdir)//'resistivity_attr.dat'
      open (newunit=iu, file=fname, status='replace', action='write')

      ! -----------------------------------------------------------
      ! PARTE 1: Atributos de Región (Nreg = 3: Aire, Mar, Subsuelo)
      ! Formato: ID  Resistividad  RepeatNumber  FixFlag
      ! -----------------------------------------------------------
      write (iu, '(I0)') OBJmodReg%Nregions
      do j = 1, OBJmodReg%Nregions
         !write (iu, '(I0, 2X, ES10.3, 1X, I3, 1X, I3)') ID_regions(j), rhoRegions(j), repRegion(j), fixed(j)
        write (iu, '(I0,2X,ES10.3,2(1X,I3))') OBJmodReg%ID(j), OBJmodReg%rho(j), OBJmodReg%repeatPartition(j), OBJmodReg%isRHOfix(j)
      end do

      ! -----------------------------------------------------------
      ! PARTE 2: Elipsoides de Control Regional
      ! -----------------------------------------------------------
      write (iu, '(3(F9.6,1x))') 0.0, 0.0, 0.0  ! Centro (Y, X, Z) en km
      write (iu, '(F4.1)') 0.0            ! Rotación (grados)
      write (iu, '(I0)') Ne                 ! Número de elipsoides
      do i = 1, Ne
         write (iu, '(4(F6.2,1x))') a_reg(i), len_reg(i), fh_reg(i), fv_reg(i)
      end do

      ! -----------------------------------------------------------
      ! PARTE 3: Refinamiento Local por Puntos (Estaciones)
      ! -----------------------------------------------------------
      write (iu, '(I0)') Nsites           ! Número de estaciones
      do i = 1, Nsites
         ! Coordenadas (X, Y, Z) en km. Nota: Z viene en metros en tu código.
         write (iu, '(3(F15.8,2x))') coorXsite(i), coorYsite(i), coorZsite(i)

         ! Número de esferas por sitio (Ejemplo: 2 esferas)
         write (iu, '(I0)') OBJmodReg%NparamEsfer
         ! Formato: Radio_km  len_km
         do j = 1, OBJmodReg%NparamEsfer
            write (iu, '(2(F3.1, 1X))') OBJmodReg%radiusForEsfer(j), OBJmodReg%edgesForEsfer(j)
         end do

      end do

      close (iu)
      print *, 'file resistivity_attr.dat escrito exitosamente en ', trim(outdir)

   end subroutine resistivitty_attribute
!=========================================================
!=======
!=========================================================
   subroutine mesh_convertion(OBJsetting)
      implicit none

      type(MeshSettings), intent(in) :: OBJsetting
      character(len=20) :: check
      integer :: stat

      check = OBJsetting%mesh_nature
      check = trim(check)

      select case (check)
      case ("native")
         print *, "No mesh convertion, building native mesh"

      case ("external")

         call EXECUTE_COMMAND_LINE('echo " " ')
         call execute_command_line('echo "step1 --> running meshio to convert the mesh " ')
         call EXECUTE_COMMAND_LINE('cd preprocessing/buildMesh && cp ../geometry/external.nas . ')
         call EXECUTE_COMMAND_LINE('meshio convert external.nas output.ele', wait=.true., exitstat=stat)
         if (stat /= 0) then
            error stop 'ERROR: tetgen execution failed'
         end if
         call execute_command_line('echo "done..." && sleep 1')

         call EXECUTE_COMMAND_LINE('echo " " ')
         call execute_command_line('echo "step 2 --> reindexing converted files to 1 based enumeration" ')
         call EXECUTE_COMMAND_LINE('cd preprocessing/buildMesh && cp ../geometry/output.* .')
         call EXECUTE_COMMAND_LINE('cd preprocessing/buildMesh && ./reindex_tetgen_1based.py output', wait=.true., exitstat=stat)
         if (stat /= 0) then
            error stop 'ERROR: tetgen execution failed'
         end if

         call execute_command_line('echo "step 3 --> running tetgen to build face, neigh and edge files" ')
         call execute_command_line('cd preprocessing/buildMesh && tetgen -n output_1b.ele', wait=.true., exitstat=stat)
         if (stat /= 0) then
            error stop 'ERROR: tetgen execution failed'
         end if
      end select

   end subroutine mesh_convertion
!=========================================================
!=======
!=========================================================
   ! subroutine reindex_tetgen_files(basename)
   !
   !    implicit none
   !    character(len=*), intent(in) :: basename
   !
   !    character(len=256) :: node_in, node_out
   !    character(len=256) :: ele_in, ele_out
   !    character(len=512) :: line
   !    integer :: iu_in, iu_out
   !    integer :: ios
   !    logical :: header_done
   !    integer :: id, n1, n2, n3, n4
   !    integer :: attr
   !    integer :: min_id
   !    integer :: first_data_id
   !    logical :: need_reindex
   !
   !    node_in = trim(basename)//'.node'
   !    node_out = trim(basename)//'_1b.node'
   !    ele_in = trim(basename)//'.ele'
   !    ele_out = trim(basename)//'_1b.ele'
   !
   !    ! ----------------------------------------
   !    ! Detect if node file is 0-based
   !    ! ----------------------------------------
   !
   !    need_reindex = .false.
   !    header_done = .false.
   !
   !    open (newunit=iu_in, file=node_in, status='old', action='read')
   !
   !    do
   !       read (iu_in, '(A)', iostat=ios) line
   !       if (ios /= 0) exit
   !
   !       if (line(1:1) == '#' .or. len_trim(line) == 0) cycle
   !
   !       if (.not. header_done) then
   !          header_done = .true.
   !          cycle
   !       end if
   !
   !       read (line, *) first_data_id
   !       if (first_data_id == 0) need_reindex = .true.
   !       exit
   !    end do
   !
   !    close (iu_in)
   !
   !    if (.not. need_reindex) then
   !       print *, 'Mesh already 1-based. No reindex needed.'
   !       return
   !    end if
   !
   !    print *, 'Reindexing mesh from 0-based to 1-based...'
   !
   !    ! ----------------------------------------
   !    ! Reindex NODE file
   !    ! ----------------------------------------
   !
   !    header_done = .false.
   !
   !    open (newunit=iu_in, file=node_in, status='old', action='read')
   !    open (newunit=iu_out, file=node_out, status='replace', action='write')
   !
   !    do
   !       read (iu_in, '(A)', iostat=ios) line
   !       if (ios /= 0) exit
   !
   !       if (line(1:1) == '#' .or. len_trim(line) == 0) then
   !          write (iu_out, '(A)') trim(line)
   !          cycle
   !       end if
   !
   !       if (.not. header_done) then
   !          write (iu_out, '(A)') trim(line)
   !          header_done = .true.
   !          cycle
   !       end if
   !
   !       read (line, *) id
   !       write (iu_out, '(I0,1X,A)') id + 1, adjustl(line(index(line, ' '):))
   !
   !    end do
   !
   !    close (iu_in)
   !    close (iu_out)
   !
   !    ! ----------------------------------------
   !    ! Reindex ELE file
   !    ! ----------------------------------------
   !
   !    header_done = .false.
   !
   !    open (newunit=iu_in, file=ele_in, status='old', action='read')
   !    open (newunit=iu_out, file=ele_out, status='replace', action='write')
   !
   !    do
   !       read (iu_in, '(A)', iostat=ios) line
   !       if (ios /= 0) exit
   !
   !       if (line(1:1) == '#' .or. len_trim(line) == 0) then
   !          write (iu_out, '(A)') trim(line)
   !          cycle
   !       end if
   !
   !       if (.not. header_done) then
   !          write (iu_out, '(A)') trim(line)
   !          header_done = .true.
   !          cycle
   !       end if
   !
   !       ! Leer primeros 5 enteros (id + 4 nodos)
   !       read (line, *) id, n1, n2, n3, n4
   !
   !       ! Escribir reindexado y conservar resto de línea (atributos)
   !       write (iu_out, '(I0,1X,I0,1X,I0,1X,I0,1X,I0,1X,A)') &
   !          id + 1, n1 + 1, n2 + 1, n3 + 1, n4 + 1, &
   !          adjustl(line(scan(line, ' ', .false., 5):))
   !
   !    end do
   !
   !    close (iu_in)
   !    close (iu_out)
   !
   !    print *, 'Reindex completed.'
   !    print *, 'Created: ', trim(node_out)
   !    print *, 'Created: ', trim(ele_out)
   !
   ! end subroutine reindex_tetgen_files
!=========================================================
!=======
!=========================================================
   subroutine run_TETGEN_and_refine_mesh(OBJmodReg)
      implicit none

      type(GlobalRefinement), intent(in) :: OBJmodReg
      integer :: r, stat
      character(len=512) :: cmd

      ! -----------------------------
      ! Runing tetgen on output.poly file
      ! -----------------------------
      call execute_command_line('echo " "')
      call execute_command_line('echo "running tetgen --> Building 2D Mesh including topography"')
      call execute_command_line('cd preprocessing/buildMesh && tetgen -nVpYAakq3.0/0 output.poly > meshtranTetGen.log 2>&1', wait=.true., exitstat=stat)
      if (stat /= 0) then
         error stop 'ERROR: tetgen execution failed'
      end if
      call execute_command_line('echo "done..." && sleep 3')
      call execute_command_line('echo " " ')
      call execute_command_line('echo " Preparing mesh refinement " && cd preprocessing/buildMesh && mkdir -p refinement', &
                                wait=.true., exitstat=stat)
      if (stat /= 0) error stop 'ERROR: could not create refinement directory'

      ! call execute_command_line('cd preprocessing/buildMesh && cp ../geometry/* .', &
      !                           wait=.true., exitstat=stat)
      ! if (stat /= 0) stop 'ERROR: copying geometry failed'
      ! call execute_command_line('echo " "')

      call execute_command_line('cd preprocessing/buildMesh && cp output.1* refinement')
     call execute_command_line('cd preprocessing/buildMesh && cp makeMtr.param obs_site.dat refinement', wait=.true., exitstat=stat)
      if (stat /= 0) stop 'ERROR: there is no files content refinement parameters'

      ! --------------------------------------------------
      ! Iterative refinement
      ! for i in 1 2 ... OBJmodReg%n_iterative_refi
      ! --------------------------------------------------
      call execute_command_line('echo " "')
      call execute_command_line('echo " Refinement step => Iterative tetgen2femtic execution"')
      r = 1
      do while (r .le. OBJmodReg%n_iterative_refi)
         write (*, '(A,I0)') '  🔁 Refinement iteration: ', r
         ! makeMtr output.$r
     write (cmd, '(A,I0, A)') 'cd preprocessing/buildMesh/refinement && makeMtr output.', r, ' >> meshtranRefinementTetGen.log 2>&1'
         call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
         if (stat /= 0) error stop 'ERROR in makeMtr'
         ! tetgen ... output.$r
         write (cmd, '(A,I0, A)') 'cd preprocessing/buildMesh/refinement && tetgen -nmpYVrAakq3.0/0 output.', r, ' >> meshtranRefinementTetGen.log 2>&1'
         call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
         if (stat /= 0) error stop 'ERROR in tetgen refinement'

         r = r + 1
      end do

      call execute_command_line('echo "Refinement done..." && sleep 1')

   end subroutine run_TETGEN_and_refine_mesh
!=========================================================
!=======
!=========================================================
   subroutine run_TetGen2Femtic(OBJsettings, OBJglobRefi)
      implicit none

      type(GlobalRefinement), intent(inout) :: OBJglobRefi
      type(MeshSettings), intent(in) :: OBJsettings
      character(len=256) :: cmd, fname
      integer :: stat, last_refinement
      logical :: ex1, ex2, ex3

      last_refinement = OBJglobRefi%n_iterative_refi + 1

      call execute_command_line('echo " "')
      call execute_command_line('echo " Last step => tetgen2femtic execution"')
      call execute_command_line('cd preprocessing/buildMesh && mkdir -p tetgenTOfemtic', wait=.true., exitstat=stat)
      if (stat /= 0) then
         error stop 'ERROR: could not be created tetgenTOfemtic directory'
      end if

      select case (OBJsettings%mesh_nature)
      case ('native')
         ! Copiar el último refinement
         write (cmd, '(A,I0,A)') 'cd preprocessing/buildMesh && cp refinement/output.', last_refinement, '* tetgenTOfemtic'
         call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
         if (stat /= 0) error stop 'ERROR copying final refinement to tetgenTOfemtic folder'

         call execute_command_line('cd preprocessing/buildMesh && cp resistivity_attr.dat tetgenTOfemtic')

         ! Ejecutar TetGen2Femtic con el último número
         write (cmd, '(A,I0,A)') 'cd preprocessing/buildMesh/tetgenTOfemtic && TetGen2Femtic output.', last_refinement, &
            ' > meshtratTetGen2Femtic.log 2>&1'
         call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
         if (stat /= 0) error stop 'ERROR running TetGen2Femtic'

         call execute_command_line('echo "done..." && sleep 1')

         ! write (cmd, '(A,I0,A)') 'cp ../preprocessing/buildMesh/tetgenTOfemtic/{mesh.dat,resistivity_block_iter0.dat,output.', &
         !    OBJglobRefi%n_iterative_refi, '.femtic.vtk} ../preprocessing/inv'
         ! call execute_command_line(trim(cmd), wait=.true., exitstat=stat)

      case ('external')

         OBJglobRefi%n_iterative_refi = 1

         write (cmd, '(A,I0,A)') 'cd preprocessing/buildMesh && cp output_1b.1* tetgenTOfemtic'
         call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
         if (stat /= 0) error stop 'ERROR copying converted final files to tetgenTOfemtic folder'

         call execute_command_line('cd preprocessing/buildMesh && cp resistivity_attr.dat tetgenTOfemtic')

         ! Ejecutar TetGen2Femtic con el último número
         write (cmd, '(A,I0)') 'cd preprocessing/buildMesh/tetgenTOfemtic && TetGen2Femtic output_1b.', last_refinement
         call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
         if (stat /= 0) error stop 'ERROR running TetGen2Femtic'

      end select

      write (cmd, '(A,I0,A)') 'cp preprocessing/buildMesh/tetgenTOfemtic/{mesh.dat,resistivity_block_iter0.dat,output.', &
         last_refinement, '.femtic.vtk} computing'
      call execute_command_line(trim(cmd), wait=.true., exitstat=stat)

      ! Verificación
      inquire (file='computing/mesh.dat', exist=ex1)
      inquire (file='computing/resistivity_block_iter0.dat', exist=ex2)
      write (fname, '(A,I0,A)') 'computing/output.', last_refinement, '.femtic.vtk'
      inquire (file=fname, exist=ex3)

      if (ex1 .and. ex2 .and. ex3) then
         write (*, *) 'Input files for running Femtic are successfully copied to computing/'
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
