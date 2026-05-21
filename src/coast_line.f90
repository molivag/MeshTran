module coast_line

   use mesh_config
   implicit none

   contains
      subroutine extract_coastLine(OBJsettings)
         implicit none

         character(len=*), parameter:: tmpfile = "curves_list.tmp"
         type(MeshSettings), intent(in)   :: OBJsettings
         character(len=250)               :: fname,DEM_file, folder, full_path
         character(len=350)               :: cmd
         integer                          :: stat, best_id

         print*, 'Extracting coast line from DEM 🗾'
         DEM_file = OBJsettings%dem_file
         DEM_file = trim(DEM_file)
         folder = DEM_file(:len_trim(DEM_file)-4)
         print*, folder
         full_path = "preprocessing/DEM/coast_line_" // trim(folder) 

         call EXECUTE_COMMAND_LINE("mkdir -p preprocessing/DEM/coast_line_"//trim(folder), exitstat = stat)
         if (stat /= 0) error stop 'Error creating coast_line'


         !Crea el directorio de trabajo gis
         cmd = "cd preprocessing/DEM/coast_line_" // trim(folder) // " && cp ../"//trim(DEM_file)//" ."
         call EXECUTE_COMMAND_LINE(trim(cmd),exitstat=stat)
         if (stat /= 0) error stop 'ERROR extracting coast line from DEM'

         ! ------ Starting The Geospatial Data Abstraction Library commands ------ !

         !step -- 1
         !ejecuta el primer comando, extraccion de la linea de costa
         write (cmd, '(A)')&
         "cd preprocessing/DEM/coast_line_" // trim(folder)// " && gdal_contour -fl 0 " // trim(DEM_file) //" coastline.shp -a elevation"
         call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
         if (stat /= 0) error stop 'ERROR extracting coast line from DEM'

         !step -- condition_2
         write(cmd, '(A)') "cd preprocessing/DEM/coast_line_" // trim(folder)// " && ogrinfo -al -geom=SUMMARY coastline.shp > "// tmpfile
         call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
         if (stat /= 0) error stop 'ERROR ccreating coast line step 2 summary of z=0 curves'

         !step -- 3
         !Buscar en el archivo temporal con las cotas z=0, la de mayor numero de puntos y extraer su id
         call get_main_coastline_id(full_path, tmpfile, best_id)
         print*, best_id

         pause


         ! Paso 4: con el id obtenido, extraer solo esa linea
         write(cmd, '(A,I0,A)') &
      'cd preprocessing/DEM/coast_line_' // trim(folder)// ' && ogr2ogr -f "ESRI Shapefile" coastline_main.shp ' // 'coastline.shp -where "ID = ',best_id,'" '
         print '(A,I0,A)', cmd



         call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
         if (stat /= 0) error stop 'ERROR in extracting main coastline step 4 in extract_coastLine'

         !STEP 5 - 
         write(cmd, '(A)') "cd preprocessing/DEM/coast_line_" // trim(folder)//"&&  ogr2ogr -f GeoJSON coast.geojson coastline_main.shp" 
         call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
         if (stat /= 0) error stop 'ERROR in coinvertig *.shp file into *.GEOjson step 5 in extract_coastLine'

         !step 6 - Converting to UTM
         write(cmd, '(A)') ' cd preprocessing/DEM/coast_line_'//trim(folder)//'&& jq -r ''.features[0].geometry.coordinates[] | ' //  '"\(.[0]) \(.[1])"'' coast.geojson > coastline_raw.csv'
         call execute_command_line(trim(cmd), wait=.true., exitstat=stat)
         if (stat /= 0) error stop 'ERROR in coinvertig *.GEOjson into *.csv step 6 in extract_coastLine'






         
      end subroutine extract_coastLine
!=========================================================
!=======
!=========================================================
      subroutine get_main_coastline_id(coast_line_path, tmpfile,  best_id)
        implicit none
        character(len=*), intent(in)  :: coast_line_path, tmpfile
        integer,          intent(out) :: best_id
      CHARACTER(len=100) ::  fname
      integer :: unit
      
        character(len=256) :: cmd, line
        integer :: ios, current_id, npts, max_pts
        character(len=50) :: dummy
      
        ! Paso 1: ejecutar ogrinfo y guardar en fichero temporal
        ! cmd = 'ogrinfo -al -geom=SUMMARY ' // trim(shpfile) // &
        !       ' > /tmp/coast_summary.txt 2>&1'
        ! call execute_command_line(trim(cmd))
      
        ! Paso 2: parsear el fichero buscando el LINESTRING mas grande
        max_pts    = 0
        best_id    = 0
        current_id = 0
      
      fname = trim(coast_line_path)//"/"//tmpfile
      print*, 'Esto es el folder en get_main_coastline ', trim(fname)
      pause
        open(newunit=unit, file=trim(fname), status='old', iostat=ios, action='read')
        if (ios /= 0) then
          write(*,*) "ERROR: no se pudo abrir ",trim(fname)
          return
        end if
         pause
      
        do
          read(unit, '(A)', iostat=ios) line
          if (ios /= 0) exit
      
          line = adjustl(line)  ! quitar espacios al inicio
      
          ! Buscar linea tipo:  ID (Integer) = 13
          if (index(line, 'ID (Integer) =') > 0) then
            read(line, *, iostat=ios) dummy, dummy, dummy, current_id
          end if
      
          ! Buscar linea tipo:  LINESTRING : 1349 points
          if (index(line, 'LINESTRING :') > 0) then
            read(line, *, iostat=ios) dummy, dummy, npts
            if (npts > max_pts) then
              max_pts = npts
              best_id = current_id
            end if
          end if
      
        end do
      
        close(unit)
      
        write(*,'(A,I0,A,I0,A)') &
          ' >> Costa principal: ID=', best_id, &
          ' con ', max_pts, ' puntos'
      
      end subroutine get_main_coastline_id



end module coast_line
