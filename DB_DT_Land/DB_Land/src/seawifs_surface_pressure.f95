! TODO Add function get_surface_pressure(lat,lon,res)
module seawifs_surface_pressure
 
  implicit none
  
! -- restrict access to only module methods listed as public.
  private
  public  ::  load_surfterr_table, get_elevation, get_surface_pressure
!  public  ::  ps2
  
!  real, dimension(360,720)    ::  ps        ! 0.5x0.5 degree
 	real, dimension(2160, 4320) ::  ps2       ! 0.09x0.09 degree
  real, dimension(2160, 4320) ::  elev      ! elevation in meters
  
  contains
  
  integer function load_surfterr_table(surfterr_file)
    implicit none
    
    include 'hdf.inc'
    include 'dffunc.inc'
  
    character(len=*), intent(in)  ::  surfterr_file

		character(len=255)    ::  sds_name
	  integer					      ::  sd_id, sds_id, sds_index
		integer, dimension(2) ::  start2, stride2, dims2
		
 		integer  				      ::  status
	  
	  load_surfterr_table = 0
    
    sd_id = sfstart(surfterr_file, DFACC_READ)
    if (sd_id == FAIL ) then
    	print *,"ERROR: failed to start SDS interface on look up table file."
    	load_surfterr_table = sd_id
    	return
    end if
  
! 		sds_name = 'ps'
! 		sds_index = sfn2index(sd_id, sds_name)
!		if (sds_index == FAIL ) then
!			print *,"ERROR: failed to convert ps SDS name to index."
!			load_surfterr_table = sds_index
!			return
!		end if
!							
!		sds_id = sfselect(sd_id, sds_index)
!		if ( sds_id == FAIL ) then
!			print *,"ERROR: failed to select ps SDS."
!			load_surfterr_table = sds_id
!			return
!		end if
!       	
!    start2 = (/0,0/)
!    stride2 = (/1,1/)
!    dims2 = (/360,720/)
!		status = sfrdata(sds_id, start2, stride2, dims2, ps)
!		if ( status == FAIL ) then
!			print *,"ERROR: failed to read ps SDS: ", status
!			load_surfterr_table = status
!			return
!		end if
!	  
!		status = sfendacc(sds_id)
!		if ( status == FAIL ) then
!			print *,"ERROR: failed to close ps SDS."
!			load_surfterr_table = status
!			return
!		end if
! 		
      sds_name = 'surface_pressure'
    sds_index = sfn2index(sd_id, sds_name)
    if (sds_index == FAIL ) then
      print *,"ERROR: failed to convert surface pressure SDS name to index."
      load_surfterr_table = sds_index
      return
    end if
            
    sds_id = sfselect(sd_id, sds_index)
    if ( sds_id == FAIL ) then
      print *,"ERROR: failed to select surface pressure SDS."
      load_surfterr_table = sds_id
      return
    end if
    
    start2  = (/0,0/)
    stride2 = (/1,1/)
    dims2 = (/2160, 4320/)
    status = sfrdata(sds_id, start2, stride2, dims2, ps2)
    if ( status == FAIL ) then
      print *,"ERROR: failed to read surface pressure SDS."
      load_surfterr_table = status
      return
    end if
  
    status = sfendacc(sds_id)
    if ( status == FAIL ) then
      print *,"ERROR: failed to close surface pressure SDS."
      load_surfterr_table = status
      return
    end if
 	  
 	  sds_name = 'surface_elevation'
	  sds_index = sfn2index(sd_id, sds_name)
		if (sds_index == FAIL ) then
			print *,"ERROR: failed to convert surface elevation SDS name to index."
			load_surfterr_table = sds_index
			return
		end if
							
		sds_id = sfselect(sd_id, sds_index)
		if ( sds_id == FAIL ) then
			print *,"ERROR: failed to select surface elevation SDS."
			load_surfterr_table = sds_id
			return
		end if
			
		start2  = (/0,0/)
		stride2 = (/1,1/)
		dims2 = (/2160, 4320/)
		status = sfrdata(sds_id, start2, stride2, dims2, elev)
		if ( status == FAIL ) then
			print *,"ERROR: failed to read surface elevation SDS."
			load_surfterr_table = status
			return
		end if
		
		status = sfendacc(sds_id)
		if ( status == FAIL ) then
			print *,"ERROR: failed to close surface elevation SDS."
			load_surfterr_table = status
			return
		end if
 	
 		status = sfend(sd_id)
    if ( status == FAIL ) then
  	  print *,"ERROR: failed to close surface terrain HDF file."
  	  load_surfterr_table = status
  	  return
    end if
    
  end function load_surfterr_table

!----------------------------------------------------------------------
  subroutine load_surface_pressure(lutable_file)
    implicit none
		
		include 'hdf.inc'
	  include 'dffunc.inc'
	  
    character(len=*)  ::  lutable_file
 		real, dimension(2160, 4320) :: sfc_pressure

    integer					      ::  sd_id, data_type, n_attrs
		integer, dimension(2) ::  start2, stride2, dims2
		character*255 	      ::  sds_name
 		integer  				      ::  status, rank, sds_id, sds_index
 		
     sd_id = sfstart(lutable_file, DFACC_READ)
    if (sd_id == FAIL ) then
    	print *,"ERROR: failed to start SDS interface on surface pressure file."
      stop
    end if
    
    sds_name = 'surface_pressure'
	  sds_index = sfn2index(sd_id, sds_name)
		if (sds_index == FAIL ) then
			print *,"ERROR: failed to convert surface pressure SDS name to index."
			stop
		end if
							
		sds_id = sfselect(sd_id, sds_index)
		if ( sds_id == FAIL ) then
			print *,"ERROR: failed to select surface pressure SDS."
			stop
		end if
	
		status = sfginfo(sds_id, sds_name, rank, dims2, data_type, n_attrs)
		if (status == FAIL) then
		  print *, "ERROR: Unable to fetch surface pressure SDS info: ", status
		  stop
		end if
		
		start2  = (/0,0/)
		stride2 = (/1,1/)
		status = sfrdata(sds_id, start2, stride2, dims2, sfc_pressure)
		if ( status == FAIL ) then
			print *,"ERROR: failed to read surface pressure SDS."
			stop
		end if
		
		status = sfendacc(sds_id)
		if ( status == FAIL ) then
			print *,"ERROR: failed to close surface pressure SDS."
			stop
		end if
	
		status = sfend(sd_id)
    if ( status == FAIL ) then
  	  print *,"ERROR: failed to close surface pressure HDF file."
    end if
    
    ps2 = sfc_pressure
    
  end subroutine load_surface_pressure

  real function get_elevation(lat, lon, status, debug) result(elevation)
    implicit none
    character(len=255)    :: func_name = "get_elevation"
    
    real, intent(in)                ::  lat
    real, intent(in)                ::  lon
    integer, intent(inout)          ::  status
    logical, intent(in), optional   ::  debug
    
    integer               ::  ilat, ilon
    logical               ::  dflag
    
    dflag = .false.
    if (present(debug)) then
      dflag = debug
    end if
    
    status = -1
    
    ilat = lat2index(lat)
    ilon = lon2index(lon)
    
    elevation = elev(ilat,ilon)
    
    if (dflag) then
      print *, trim(func_name), ": lat, lon, ilat, ilon, elev: ", lat, lon, ilat, ilon, elevation
    end if
    
    status = 0
    return
    
    contains
    
    integer function lat2index(lat) result(index)
      implicit none
      
      real, intent(in)  ::  lat
      index = floor((90.0-lat) * 12.0) + 1
      if (index > 2160) index = 2160
      if (index < 1)   index = 1
      
      return
      
    end function lat2index
    
    integer function lon2index(lon) result(index)
      implicit none
      
      real, intent(in)    ::  lon
      index = floor((lon+180.0) * 12.0) + 1
		  if (index > 4320) index = 4320
		  if (index < 1)   index = 1
		  
		  return
		end function lon2index
		
  end function get_elevation
  
  real function get_surface_pressure(lat, lon, status, debug) result(sp)
    implicit none
    character(len=255)    :: func_name = "get_surface_pressure"
    
    real, intent(in)                ::  lat
    real, intent(in)                ::  lon
    integer, intent(inout)          ::  status
    logical, intent(in), optional   ::  debug
    
    integer               ::  ilat, ilon
    logical               ::  dflag
    
    dflag = .false.
    if (present(debug)) then
      dflag = debug
    end if
    
    status = 0
    
    ilat = lat2index(lat)
    ilon = lon2index(lon)
    
    sp = ps2(ilat,ilon)/1013.25
    
    if (dflag) then
      print *, trim(func_name), ": lat, lon, ilat, ilon, pressure: ", lat, lon, ilat, ilon, sp
    end if
    
    return
    
    contains
    
    integer function lat2index(lat) result(index)
      implicit none
      
      real, intent(in)  ::  lat
      index = floor((90.0-lat) * 12.0) + 1
      if (index > 2160) index = 2160
      if (index < 1)   index = 1
      
      return
      
    end function lat2index
    
    integer function lon2index(lon) result(index)
      implicit none
      
      real, intent(in)    ::  lon
      index = floor((lon+180.0) * 12.0) + 1
		  if (index > 4320) index = 4320
		  if (index < 1)   index = 1
		  
		  return
		end function lon2index
	
	end function get_surface_pressure
end module seawifs_surface_pressure
