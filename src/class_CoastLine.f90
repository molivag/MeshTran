module class_CoastLine
   !
   ! This class module define the object CoastLine
   ! Implements the generate procedure which is an orchestation of the entire
   ! workflow of extract, simplify, close and visualization of the coast line
   !
   ! The final product is a closed-simplified coast line in UTM coordinates.
   ! whit as same as the DEM, x pointing to the north and y to the east.

   use mesh_config
   implicit none

   private

   public :: Coast_Line

   type :: Coast_Line
      character(len=5) :: has_sea
      integer :: npoints = 0
      real(dp) :: sea_level

      real(dp), allocatable :: x(:)
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: lat(:)
      real(dp), allocatable :: long(:)

      logical :: is_closed = .false.

   contains
      procedure :: generate
      procedure :: extract
      ! procedure :: simplify
      procedure :: closing

   end type Coast_Line

contains
   subroutine generate(this, OBJsettings, dem_x, dem_y, dem_z, nDEM)!, closed_Clx, closed_CLy, npoinClosedCL)

      use geo_utils, only: lat_long_to_UTM_km
      implicit none

      class(Coast_Line), intent(inout) :: this

      type(MeshSettings), intent(in)                                 :: OBJsettings
      integer, intent(in)                                            :: nDEM
      real(dp), intent(in)                                           :: dem_x(nDEM), dem_y(nDEM), dem_z(nDEM)
      ! real(dp), dimension(:), allocatable, intent(out)               :: closed_CLx, closed_Cly
      ! integer,                             intent(out)               :: npoinClosedCL
      integer                             :: npoincoastLine
      ! real(dp), dimension(:), allocatable :: cordX_cl, cordY_cl
      real(dp), DIMENSION(:), allocatable :: east_cl_km, north_cl_km

      !step 1 we get the coast line from DEM in lat-long
      call this%extract(OBJsettings, npoincoastLine)!, cordX_cl, cordY_cl)

      ! step 2 we convert coast line coordinate lat-long to km UTM
      call lat_long_to_UTM_km(this%lat, this%long, this%npoints, east_cl_km, north_cl_km)
      ! Despues de que mande llamar la rutina me entregara las coordenadas norte y este de la linea de costa

      !at this point the previous routine lat_long... gives east and north as is common in the world
      ! y=north and x=east

      allocate (this%x(this%npoints), this%y(this%npoints))
      this%x = north_cl_km
      this%y = east_cl_km

      !however for MT convention where x point to north, we feed the following routine with x=north
      call this%closing( dem_x, dem_y, dem_z, nDEM)!, closed_CLx, closed_Cly, npoinClosedCL)

   end subroutine generate
   !=========================================================
   !=======
   !=========================================================

   subroutine extract(this, OBJsettings, pointsCoastLine)!, cl_x, cl_y)
      implicit none

      class(Coast_Line), intent(inout) :: this
      !
      character(len=*), parameter:: tmpfile = "curves_list.tmp"
      type(MeshSettings), intent(in)   :: OBJsettings
      ! real(dp), dimension(:), allocatable, INTENT(OUT)  :: cl_x, cl_y
      integer, intent(out)  :: pointsCoastLine
      character(len=250)               :: DEM_file, folder, full_path
      character(len=350)               :: cmd
      integer                          :: stat, best_id, iu, i
      real(dp) :: longg, latt

      print *, 'Extracting MAIN coast line from DEM 🗾'
      DEM_file = OBJsettings%dem_file
      DEM_file = trim(DEM_file)
      folder = DEM_file(:len_trim(DEM_file) - 4)
      print *, folder
      full_path = "preprocessing/DEM/coast_line_"//trim(folder)

      call EXECUTE_COMMAND_LINE("mkdir -p preprocessing/DEM/coast_line_"//trim(folder), exitstat=stat)
      if (stat /= 0) error stop 'Error creating coast_line'

      !Crea el directorio de trabajo gis
      cmd = "cd preprocessing/DEM/coast_line_"//trim(folder)//" && cp ../"//trim(DEM_file)//" ."
      call EXECUTE_COMMAND_LINE(trim(cmd), exitstat=stat)
      if (stat /= 0) error stop 'ERROR extracting coast line from DEM'

      ! ------ Starting The Geospatial Data Abstraction Library commands ------ !

      !step -- 1
      !ejecuta el primer comando, extraccion de la linea de costa
      write (cmd, '(A)') &
         "cd preprocessing/DEM/coast_line_"//trim(folder)//" && gdal_contour -fl 0 "//trim(DEM_file)// &
         " coastline.shp -a elevation"
      call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
      if (stat /= 0) error stop 'ERROR extracting coast line from DEM'

      !step -- condition_2
      write (cmd, '(A)') "cd preprocessing/DEM/coast_line_"//trim(folder)// &
         " && ogrinfo -al -geom=SUMMARY coastline.shp > "//tmpfile
      call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
      if (stat /= 0) error stop 'ERROR ccreating coast line step 2 summary of z=0 curves'

      !step -- 3
      !Buscar en el archivo temporal con las cotas z=0, la de mayor numero de puntos y extraer su id
      call get_main_coastline_id(full_path, tmpfile, best_id, pointsCoastLine)
      print *, best_id

      ! Paso 4: con el id obtenido, extraer solo esa linea
      write (cmd, '(A,I0,A)') &
         'cd preprocessing/DEM/coast_line_'//trim(folder)// &
         ' && ogr2ogr -f "ESRI Shapefile" coastline_main.shp '//'coastline.shp -where "ID = ', best_id, '" '
      print '(A,I0,A)', cmd

      call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
      if (stat /= 0) error stop 'ERROR in extracting main coastline step 4 in extract_coastLine'

      !STEP 5 -
      write (cmd, '(A)') "cd preprocessing/DEM/coast_line_"//trim(folder)// &
         "&&  ogr2ogr -f GeoJSON coast.geojson coastline_main.shp"
      call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
      if (stat /= 0) error stop 'ERROR in coinvertig *.shp file into *.GEOjson step 5 in extract_coastLine'

      !step 6 - Converting to UTM
      write (cmd, '(A)') ' cd preprocessing/DEM/coast_line_'//trim(folder)// &
         '&& jq -r ''.features[0].geometry.coordinates[] | '//'"\(.[0]) \(.[1])"'' coast.geojson > coastline_raw.csv'
      call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
      if (stat /= 0) error stop 'ERROR in coinvertig *.GEOjson into *.csv step 6 in extract_coastLine'

      !The very lst step is to read again the csv just generated and store x an y coast line coordinate
      allocate (this%long(pointsCoastLine), this%lat(pointsCoastLine))

      open (newunit=iu, file=trim(full_path)//'/coastline_raw.csv', status='old')
      do i = 1, pointsCoastLine, 1
         read (iu, *) longg , latt 
         this%long(i) = longg
         this%lat(i) = latt
      end do
      close (iu)
      this%npoints = pointsCoastLine
      !
   end subroutine extract
   !=========================================================
   !=======
   !=========================================================
   subroutine get_main_coastline_id(coast_line_path, tmpfile, best_id, max_pts)
      implicit none
      character(len=*), intent(in)  :: coast_line_path, tmpfile
      integer, intent(out) :: best_id, max_pts
      CHARACTER(len=100) ::  fname
      integer :: unit

      character(len=256) :: line
      integer :: ios, current_id, npts
      character(len=50) :: dummy

      ! Paso 1: ejecutar ogrinfo y guardar en fichero temporal
      ! cmd = 'ogrinfo -al -geom=SUMMARY ' // trim(shpfile) // &
      !       ' > /tmp/coast_summary.txt 2>&1'
      ! call execute_command_line(trim(cmd))

      ! Paso 2: parsear el fichero buscando el LINESTRING mas grande
      max_pts = 0
      best_id = 0
      current_id = 0

      fname = trim(coast_line_path)//"/"//tmpfile
      print *, 'Esto es el folder en get_main_coastline ', trim(fname)
      open (newunit=unit, file=trim(fname), status='old', iostat=ios, action='read')
      if (ios /= 0) then
         write (*, *) "ERROR: no se pudo abrir ", trim(fname)
         return
      end if

      do
         read (unit, '(A)', iostat=ios) line
         if (ios /= 0) exit

         line = adjustl(line)  ! quitar espacios al inicio

         ! Buscar linea tipo:  ID (Integer) = 13
         if (index(line, 'ID (Integer) =') > 0) then
            read (line, *, iostat=ios) dummy, dummy, dummy, current_id
         end if

         ! Buscar linea tipo:  LINESTRING : 1349 points
         if (index(line, 'LINESTRING :') > 0) then
            read (line, *, iostat=ios) dummy, dummy, npts
            if (npts > max_pts) then
               max_pts = npts
               best_id = current_id
            end if
         end if

      end do

      close (unit)

      write (*, '(A,I0,A,I0,A)') &
         ' >> Costa principal: ID=', best_id, &
         ' con ', max_pts, ' puntos'

   end subroutine get_main_coastline_id
   !!=========================================================
   !!=======
   !!=========================================================
   subroutine closing(this, dem_x, dem_y, dem_z, nDEM)!, closed_x, closed_y, n_closed)
      !   
      !   
      !   
      ! This routine closes the raw coastline extracted from the DEM.
      !
      ! The coastline obtained from GDAL is an open polyline, not a closed polygon.
      ! A closed polygon is required later to determine whether DEM points belong
      ! to land or sea using point-in-polygon algorithms (Ray Casting).
      !
      ! The routine analyses the DEM boundaries (west/east/north/south) and
      ! determines where the sea is predominantly located. The coastline is then
      ! automatically closed toward the sea side by adding artificial corner points
      ! along the domain boundary.
      !
      ! This avoids hardcoding assumptions such as:
      !   - sea always at west
      !   - sea always at north
      ! etc.
      !
      ! The routine works for general coastal DEMs as long as the sea is mainly
      ! located along one side of the domain.
      !
      ! Input:
      !   - raw coastline coordinates
      !   - DEM coordinates and elevations
      !   - domain limits
      !
      ! Output:
      !   - closed coastline polygon
      !
      ! Notes:
      !   - DEM NoData values must be cleaned before using this routine.
      !   - The coastline points must already be ordered.
      !   - This routine does NOT classify land/sea yet.
      ! -------------------------------------------------------------------------



      implicit none

      class(Coast_Line), intent(inout) :: this

      integer, intent(in) :: nDEM

      real(dp), intent(in) :: dem_x(nDEM)
      real(dp), intent(in) :: dem_y(nDEM)
      real(dp), intent(in) :: dem_z(nDEM)

      ! ============================================================
      ! LOCAL
      ! ============================================================

      integer :: i, nsea, n_coast, ntotal, n_closed

      real(dp) :: frac, xmin, xmax, ymin, ymax, area,  tol
      real(dp), DIMENSION(:), allocatable :: coast_x, coast_y, closed_x, closed_y

      logical :: sea_west, sea_xmin, sea_xmax, sea_ymin, sea_ymax
      logical :: sea_east
      logical :: sea_north
      logical :: sea_south

      character(len=10) :: close_to


      xmin = minval(dem_x)       !min North
      xmax = maxval(dem_x)        !max North  
      ymin = minval(dem_y)       !min East
      ymax = maxval(dem_y)        !max East    

      n_coast =  this%npoints 
      allocate( coast_x(n_coast), coast_y(n_coast) )

     coast_x = this%x 
     coast_y = this%y



      ! ============================================================
      ! INITIALIZE
      ! ============================================================

      sea_xmin = .false. 
      sea_xmax = .false.
      sea_ymin = .false.
      sea_ymax = .false.

      sea_west = .false.
      sea_east = .false.
      sea_north = .false.
      sea_south = .false.

      ! ============================================================
      ! WEST EDGE is ymin
      ! ============================================================

      nsea = 0
      tol = 1.d-3
      ntotal = 0

      do i = 1, nDEM

         if (abs(dem_y(i) - ymin) < tol) then
            if (dem_z(i) < 0.d0) then
               nsea = nsea + 1
            end if
            ! este punto pertenece al borde oeste
            ntotal = ntotal + 1
         end if

      end do
      frac = real(nsea, 8)/ntotal

      if (frac > 0.7d0) then
         sea_ymin = .true.
      end if

         sea_west = sea_ymin


      ! ============================================================
      ! EAST EDGE is ymax
      ! ============================================================
     
      nsea = 0
      tol = 1.d-3
      ntotal = 0
     
      do i = 1, nDEM

         if (abs(dem_y(i) - ymax) < tol) then
            if (dem_z(i) < 0.d0) then
               nsea = nsea + 1
            end if
            ! este punto pertenece al borde oeste
            ntotal = ntotal + 1
         end if

      end do
      frac = real(nsea, 8)/ntotal

      if (frac > 0.7d0) then
         sea_ymax = .true.
      end if

         sea_east = sea_ymax
     
      ! ============================================================
      ! SOUTH EDGE is xmin
      ! ============================================================
     

      nsea = 0
      tol = 1.d-3
      ntotal = 0
     
      do i = 1, nDEM

         if (abs(dem_x(i) - xmin) < tol) then
            if (dem_z(i) < 0.d0) then
               nsea = nsea + 1
            end if
            ! este punto pertenece al borde oeste
            ntotal = ntotal + 1
         end if

      end do
      frac = real(nsea, 8)/ntotal

      if (frac > 0.7d0) then
         sea_xmin = .true.
      end if

         sea_south = sea_xmin

     
      ! ============================================================
      ! NORTH EDGE is xmax
      ! ============================================================
     
      nsea = 0
      tol = 1.d-3
      ntotal = 0
     
      do i = 1, nDEM

         if (abs(dem_x(i) - xmax) < tol) then
            if (dem_z(i) < 0.d0) then
               nsea = nsea + 1
            end if
            ! este punto pertenece al borde oeste
            ntotal = ntotal + 1
         end if

      end do
      frac = real(nsea, 8)/ntotal

      if (frac > 0.7d0) then
         sea_xmax = .true.
      end if

         sea_north = sea_xmax



      ! ============================================================
      ! PRINT INFO
      ! ============================================================
     
      print *, "Sea west  :", sea_west
      print *, "Sea east  :", sea_east
      print *, "Sea south :", sea_south
      print *, "Sea north :", sea_north
     
      ! ============================================================
      ! DETERMINE CLOSURE SIDE
      ! ============================================================
     
      close_to = "none"
     
      if (sea_east) then
         close_to = "west"
      elseif (sea_west) then
         close_to = "east"
      elseif (sea_north) then
         close_to = "south"
      elseif (sea_south) then
         close_to = "north"
      end if
     
      print *, "Closing coastline toward :", trim(close_to)
     
      ! ============================================================
      ! CHECK ORIENTATION (signed area)
      ! ============================================================
     
      area = 0.d0
     
      do i = 1, n_coast - 1
     
         area = area + &
                (coast_x(i)*coast_y(i + 1) - &
                 coast_x(i + 1)*coast_y(i))
     
      end do
     
      area = 0.5d0*area
     
      print *, "Signed area =", area
     
      if (area > 0.d0) then
         print *, "Orientation : CCW"
      else
         print *, "Orientation : CW"
      end if
     
      ! ============================================================
      ! ALLOCATE CLOSED POLYGON
      ! ============================================================
     
      n_closed = n_coast + 4
     
      allocate (closed_x(n_closed))
      allocate (closed_y(n_closed))
     
      ! ============================================================
      ! COPY ORIGINAL COASTLINE
      ! ============================================================
     
      closed_x(1:n_coast) = coast_x(:)
      closed_y(1:n_coast) = coast_y(:)
     
      ! ============================================================
      ! CLOSE POLYGON
      ! ============================================================

      ! Here    west/east move y
      ! and   north/south move x
     
      select case (trim(close_to))
     
         ! ============================================================
         ! EAST
         ! ============================================================
     
      case ("east")
     
         closed_x(n_coast + 1) = coast_x(n_coast)
         closed_y(n_coast + 1) = ymax
     
         closed_x(n_coast + 2) = xmax
         closed_y(n_coast + 2) = ymax
     
         closed_x(n_coast + 3) = coast_x(1)
         closed_y(n_coast + 3) = ymax
     
         closed_x(n_coast + 4) = coast_x(1)
         closed_y(n_coast + 4) = coast_y(1)
     
         ! ============================================================
         ! WEST
         ! ============================================================
     
      case ("west")
         !horizontal proyection to the west boundary
     
         closed_x(n_coast + 1) = coast_x(n_coast)
         closed_y(n_coast + 1) = ymin
     
         closed_x(n_coast + 2) = xmax
         closed_y(n_coast + 2) = ymin
     
         closed_x(n_coast + 3) = coast_x(1)
         closed_y(n_coast + 3) = ymin
     
         closed_x(n_coast + 4) = coast_x(1)
         closed_y(n_coast + 4) = coast_y(1)
     
         ! ============================================================
         ! NORTH
         ! ============================================================
     
      case ("north")
     
         closed_x(n_coast + 1) = xmax
         closed_y(n_coast + 1) = coast_x(n_coast)
     
         closed_x(n_coast + 2) = xmax
         closed_y(n_coast + 2) = ymax
     
         closed_x(n_coast + 3) = xmax
         closed_y(n_coast + 3) = coast_x(1)
     
         closed_x(n_coast + 4) = coast_x(1)
         closed_y(n_coast + 4) = coast_y(1)
     
         ! ============================================================
         ! SOUTH
         ! ============================================================
     
      case ("south")
     
         closed_x(n_coast + 1) = xmin
         closed_y(n_coast + 1) = coast_y(n_coast)
     
         closed_x(n_coast + 2) = xmin
         closed_y(n_coast + 2) = ymin
     
         closed_x(n_coast + 3) = xmin
         closed_y(n_coast + 3) = coast_y(1)
     
         closed_x(n_coast + 4) = coast_x(1)
         closed_y(n_coast + 4) = coast_y(1)
     
         ! ============================================================
         ! DEFAULT
         ! ============================================================
     
      case default
     
         print *, "ERROR : could not determine closure side"
         stop
     
      end select

      this%npoints = n_closed
      this%x = closed_x
      this%y = closed_y

   end subroutine closing
!
end module class_CoastLine
