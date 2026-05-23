module class_CoastLine

   use mesh_config
   implicit none

   private 

   public :: CoastLine



   type :: Coast_Line
      integer :: npoints = 0

      real(dp), allocatable :: x(:)
      real(dp), allocatable :: y(:)

      logical :: is_closed = .false.
      

   contains
      procedure :: generate
      procedure :: extract
      procedure :: closing

   end type Coast_Line

contains
   subroutine generate(this, OBJsettings, dem_x, dem_y, dem_z,nDEM)

      use geo_utils, only: lat_long_to_UTM_km
      implicit none

      class(Coast_Line), intent(inout) :: this

      type(MeshSettings), intent(in)      :: OBJsettings
      integer, intent(in)                 :: nDEM
      real(dp), intent(in)                :: dem_x(nDEM)
      real(dp), intent(in)                :: dem_y(nDEM)
      real(dp), intent(in)                :: dem_z(nDEM)
      character(len=250)                  :: fname, DEM_file, folder, full_path
      character(len=350)                  :: cmd
      integer                             :: stat, best_id, npoincoastLine
      real(dp), dimension(:), allocatable :: cordX_cl, cordY_cl
      real(dp), DIMENSION(:), allocatable :: east_cl_km, north_cl_km

      !step 1 we get the coast line from DEM in lat-long
      call this%extract(OBJsettings, npoincoastLine, cordX_cl, cordY_cl)

      ! step 2 we convert coast line coordinate lat-long to km UTM
      call lat_long_to_UTM_km(cordX_cl, cordY_cl, npoincoastLine, east_cl_km, north_cl_km)
      ! Despues de que mande llamar la rutina me entregara las coordenadas norte y este de la linea de costa

      !at this point the previous routine lat_long... gives east and north as is common in the world
      ! y=north and x=east

      !however for MT convention where x point to north, we feed the following routine with x=north
      call this%closing(OBJsettings, npoincoastLine, north_cl_km, east_cl_km, dem_x, dem_y, dem_z,nDEM )

   end subroutine generate
   !=========================================================
   !=======
   !=========================================================

   subroutine extract(this, OBJsettings, pointsCoastLine, cl_x, cl_y)
      implicit none

      class(Coast_Line), intent(inout) :: this
      !
      character(len=*), parameter:: tmpfile = "curves_list.tmp"
      type(MeshSettings), intent(in)   :: OBJsettings
      real(dp), dimension(:), allocatable, INTENT(OUT)  :: cl_x, cl_y
      integer, intent(out)  :: pointsCoastLine
      character(len=250)               :: fname, DEM_file, folder, full_path
      character(len=350)               :: cmd
      integer                          :: stat, best_id, iu, i
      real(dp) :: cl_coordX, cl_coordY

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
      allocate (cl_x(pointsCoastLine), cl_y(pointsCoastLine))
      open (newunit=iu, file=trim(full_path)//'coastline_raw.csv', status='old')
      do i = 1, pointsCoastLine, 1
         read (iu, *) cl_coordX, cl_coordY
         cl_x(i) = cl_coordX
         cl_y(i) = cl_coordY
      end do
      close (iu)
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

     character(len=256) :: cmd, line
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
     pause
     open (newunit=unit, file=trim(fname), status='old', iostat=ios, action='read')
     if (ios /= 0) then
        write (*, *) "ERROR: no se pudo abrir ", trim(fname)
        return
     end if
     pause

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
   subroutine closing(this, OBJsettings, pointsCoastLine, coast_x, coast_y, dem_x, dem_y, dem_z, nDEM)!, closed_x, closed_y, n_closed)
!    !
     implicit none

     class(Coast_Line), intent(inout) :: this

     type(MeshSettings), INTENT(IN) :: OBJsettings
     integer, intent(in) :: pointsCoastLine, nDEM
     real(dp), intent(in) :: coast_x(pointsCoastLine)
     real(dp), intent(in) :: coast_y(pointsCoastLine)


     real(dp), intent(in) :: dem_x(nDEM)
     real(dp), intent(in) :: dem_y(nDEM)
     real(dp), intent(in) :: dem_z(nDEM)

     ! ============================================================
     ! OUTPUT
     ! ============================================================

     ! integer, intent(out) :: n_closed
     ! real(dp), allocatable, intent(out) :: closed_x(:)
     ! real(dp), allocatable, intent(out) :: closed_y(:)

     ! ============================================================
     ! LOCAL
     ! ============================================================

     integer :: i, j
     integer :: nsea

     real(dp) :: frac, xmin, xmax, ymin, ymax

     logical :: sea_west
     logical :: sea_east
     logical :: sea_north
     logical :: sea_south

     character(len=10) :: close_to

     real(dp) :: area, n_coast


     xmin = OBJsettings%xminDOM
     xmax = OBJsettings%xmaxDOM
     ymin = OBJsettings%yminDOM
     ymax = OBJsettings%ymaxDOM
     n_coast = pointsCoastLine
!    !
!    !
!    !
!    !   ! ============================================================
!    !   ! INITIALIZE
!    !   ! ============================================================
!    !
!    !   sea_west = .false.
!    !   sea_east = .false.
!    !   sea_north = .false.
!    !   sea_south = .false.
!    !
!    !   ! ============================================================
!    !   ! WEST EDGE
!    !   ! ============================================================
!    !
!    !   nsea = 0
!    !
!    !   do j = 1, nDEM
!    !      if (dem_z(j) < 0.d0) then
!    !         nsea = nsea + 1
!    !      end if
!    !   end do
!    !
!    !   frac = real(nsea, 8)/real(nDEM, 8)
!    !
!    !   if (frac > 0.7d0) sea_west = .true.
!    !   !
!    !   ! ! ============================================================
!    !   ! ! EAST EDGE
!    !   ! ! ============================================================
!    !   !
!    !   ! nsea = 0
!    !   !
!    !   ! do j = 1, ny
!    !   !    if (dem_z(nx, j) < 0.d0) then
!    !   !       nsea = nsea + 1
!    !   !    end if
!    !   ! end do
!    !   !
!    !   ! frac = real(nsea, 8)/real(ny, 8)
!    !   !
!    !   ! if (frac > 0.7d0) sea_east = .true.
!    !   !
!    !   ! ! ============================================================
!    !   ! ! SOUTH EDGE
!    !   ! ! ============================================================
!    !   !
!    !   ! nsea = 0
!    !   !
!    !   ! do i = 1, nx
!    !   !    if (dem_z(i, 1) < 0.d0) then
!    !   !       nsea = nsea + 1
!    !   !    end if
!    !   ! end do
!    !   !
!    !   ! frac = real(nsea, 8)/real(nx, 8)
!    !   !
!    !   ! if (frac > 0.7d0) sea_south = .true.
!    !   !
!    !   ! ! ============================================================
!    !   ! ! NORTH EDGE
!    !   ! ! ============================================================
!    !   !
!    !   ! nsea = 0
!    !   !
!    !   ! do i = 1, nx
!    !   !    if (dem_z(i, ny) < 0.d0) then
!    !   !       nsea = nsea + 1
!    !   !    end if
!    !   ! end do
!    !   !
!    !   ! frac = real(nsea, 8)/real(nx, 8)
!    !   !
!    !   ! if (frac > 0.7d0) sea_north = .true.
!    !   !
!    !   ! ! ============================================================
!    !   ! ! PRINT INFO
!    !   ! ! ============================================================
!    !   !
!    !   ! print *, "Sea west  :", sea_west
!    !   ! print *, "Sea east  :", sea_east
!    !   ! print *, "Sea south :", sea_south
!    !   ! print *, "Sea north :", sea_north
!    !   !
!    !   ! ! ============================================================
!    !   ! ! DETERMINE CLOSURE SIDE
!    !   ! ! ============================================================
!    !   !
!    !   ! close_to = "none"
!    !   !
!    !   ! if (sea_east) then
!    !   !    close_to = "east"
!    !   ! elseif (sea_west) then
!    !   !    close_to = "west"
!    !   ! elseif (sea_north) then
!    !   !    close_to = "north"
!    !   ! elseif (sea_south) then
!    !   !    close_to = "south"
!    !   ! end if
!    !   !
!    !   ! print *, "Closing coastline toward :", trim(close_to)
!    !   !
!    !   ! ! ============================================================
!    !   ! ! CHECK ORIENTATION (signed area)
!    !   ! ! ============================================================
!    !   !
!    !   ! area = 0.d0
!    !   !
!    !   ! do i = 1, n_coast - 1
!    !   !
!    !   !    area = area + &
!    !   !           (coast_x(i)*coast_y(i + 1) - &
!    !   !            coast_x(i + 1)*coast_y(i))
!    !   !
!    !   ! end do
!    !   !
!    !   ! area = 0.5d0*area
!    !   !
!    !   ! print *, "Signed area =", area
!    !   !
!    !   ! if (area > 0.d0) then
!    !   !    print *, "Orientation : CCW"
!    !   ! else
!    !   !    print *, "Orientation : CW"
!    !   ! end if
!    !   !
!    !   ! ! ============================================================
!    !   ! ! ALLOCATE CLOSED POLYGON
!    !   ! ! ============================================================
!    !   !
!    !   n_closed = n_coast + 4
!    !   !
!    !   ! allocate (closed_x(n_closed))
!    !   ! allocate (closed_y(n_closed))
!    !   !
!    !   ! ! ============================================================
!    !   ! ! COPY ORIGINAL COASTLINE
!    !   ! ! ============================================================
!    !   !
!    !   ! closed_x(1:n_coast) = coast_x(:)
!    !   ! closed_y(1:n_coast) = coast_y(:)
!    !   !
!    !   ! ! ============================================================
!    !   ! ! CLOSE POLYGON
!    !   ! ! ============================================================
!    !   !
!    !   ! select case (trim(close_to))
!    !   !
!    !   !    ! ============================================================
!    !   !    ! EAST
!    !   !    ! ============================================================
!    !   !
!    !   ! case ("east")
!    !   !
!    !   !    closed_x(n_coast + 1) = xmax
!    !   !    closed_y(n_coast + 1) = coast_y(n_coast)
!    !   !
!    !   !    closed_x(n_coast + 2) = xmax
!    !   !    closed_y(n_coast + 2) = ymax
!    !   !
!    !   !    closed_x(n_coast + 3) = xmax
!    !   !    closed_y(n_coast + 3) = coast_y(1)
!    !   !
!    !   !    closed_x(n_coast + 4) = coast_x(1)
!    !   !    closed_y(n_coast + 4) = coast_y(1)
!    !   !
!    !   !    ! ============================================================
!    !   !    ! WEST
!    !   !    ! ============================================================
!    !   !
!    !   ! case ("west")
!    !   !
!    !   !    closed_x(n_coast + 1) = xmin
!    !   !    closed_y(n_coast + 1) = coast_y(n_coast)
!    !   !
!    !   !    closed_x(n_coast + 2) = xmin
!    !   !    closed_y(n_coast + 2) = ymax
!    !   !
!    !   !    closed_x(n_coast + 3) = xmin
!    !   !    closed_y(n_coast + 3) = coast_y(1)
!    !   !
!    !   !    closed_x(n_coast + 4) = coast_x(1)
!    !   !    closed_y(n_coast + 4) = coast_y(1)
!    !   !
!    !   !    ! ============================================================
!    !   !    ! NORTH
!    !   !    ! ============================================================
!    !   !
!    !   ! case ("north")
!    !   !
!    !   !    closed_x(n_coast + 1) = coast_x(n_coast)
!    !   !    closed_y(n_coast + 1) = ymax
!    !   !
!    !   !    closed_x(n_coast + 2) = xmax
!    !   !    closed_y(n_coast + 2) = ymax
!    !   !
!    !   !    closed_x(n_coast + 3) = coast_x(1)
!    !   !    closed_y(n_coast + 3) = ymax
!    !   !
!    !   !    closed_x(n_coast + 4) = coast_x(1)
!    !   !    closed_y(n_coast + 4) = coast_y(1)
!    !   !
!    !   !    ! ============================================================
!    !   !    ! SOUTH
!    !   !    ! ============================================================
!    !   !
!    !   ! case ("south")
!    !   !
!    !   !    closed_x(n_coast + 1) = coast_x(n_coast)
!    !   !    closed_y(n_coast + 1) = ymin
!    !   !
!    !   !    closed_x(n_coast + 2) = xmax
!    !   !    closed_y(n_coast + 2) = ymin
!    !   !
!    !   !    closed_x(n_coast + 3) = coast_x(1)
!    !   !    closed_y(n_coast + 3) = ymin
!    !   !
!    !   !    closed_x(n_coast + 4) = coast_x(1)
!    !   !    closed_y(n_coast + 4) = coast_y(1)
!    !   !
!    !   !    ! ============================================================
!    !   !    ! DEFAULT
!    !   !    ! ============================================================
!    !   !
!    !   ! case default
!    !   !
!    !   !    print *, "ERROR : could not determine closure side"
!    !   !    stop
!    !   !
!    !   ! end select
!    !
   end subroutine closing
!
end module class_CoastLine
