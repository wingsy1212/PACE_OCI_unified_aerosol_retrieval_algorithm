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
    
    use netcdf
    USE OCIUAAER_Config_Module
  
    character(len=*), intent(in)  ::  surfterr_file

		character(len=255)    ::  sds_name
	  integer					      ::  sd_id, sds_id, sds_index
		integer, dimension(2) ::  start2, stride2, dims2
		
 		integer  				      ::  status
	  
	  load_surfterr_table = 0
    
    character(len=255)    ::  sds_name
    character(len=255)    ::  dset_name
    character(len=255)    ::  attr_name
    character(len=255)    ::  group_name

    integer               ::  nc_id
    integer               ::  dim_id
    integer               ::  dset_id
    integer               ::  grp_id

    start  = (/ 1,1 /)
    edge   = (/ 2160, 4320 /)
    stride = (/ 1,1 /)

    dbdt_file = cfg%db_nc4
    status = nf90_open(dbdt_file, nf90_nowrite, nc_id)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to open deepblue lut_nc4 file: ", status
        return
    end if

    group_name = 'surface_pressure'
    status = nf90_inq_ncid(nc_id, group_name, grp_id)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to get ID of group "//trim(group_name)//": ", status
        return
    end if

    dset_name = 'surface_pressure'
    status = nf90_inq_varid(grp_id, dset_name, dset_id)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
        return
    end if
    status = nf90_get_var(grp_id, dset_id, ps2, start=start, &
                          stride=stride, count=edge)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
        return
    end if

    dset_name = 'surface_elevation'
    status = nf90_inq_varid(grp_id, dset_name, dset_id)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
        return
    end if
    status = nf90_get_var(grp_id, dset_id, elev, start=start, &
                          stride=stride, count=edge)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
        return
    end if

    status = nf90_close(nc_id)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to close lut_nc4 file: ", status
        return
    end if
    
  end function load_surfterr_table

!----------------------------------------------------------------------
  subroutine load_surface_pressure(lutable_file)
    implicit none
		
    use netcdf
    USE OCIUAAER_Config_Module
	  
    character(len=*)  ::  lutable_file
 		real, dimension(2160, 4320) :: sfc_pressure

    integer					      ::  sd_id, data_type, n_attrs
		integer, dimension(2) ::  start2, stride2, dims2
		character*255 	      ::  sds_name
 		integer  				      ::  status, rank, sds_id, sds_index

    character(len=255)    ::  sds_name
    character(len=255)    ::  dset_name
    character(len=255)    ::  attr_name
    character(len=255)    ::  group_name

    integer               ::  nc_id
    integer               ::  dim_id
    integer               ::  dset_id
    integer               ::  grp_id

    start  = (/ 1,1 /)
    edge   = (/ 2160, 4320 /)
    stride = (/ 1,1 /)

    dbdt_file = cfg%db_nc4
    status = nf90_open(dbdt_file, nf90_nowrite, nc_id)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to open deepblue lut_nc4 file: ", status
        return
    end if

    group_name = 'dbdt_regions'
    status = nf90_inq_ncid(nc_id, group_name, grp_id)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to get ID of group "//trim(group_name)//": ", status
        return
    end if

    dset_name = 'surface_pressure'
    status = nf90_inq_varid(grp_id, dset_name, dset_id)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
        return
    end if
    status = nf90_get_var(grp_id, dset_id, sfc_pressure, start=start, &
                          stride=stride, count=edge)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
        return
    end if

    status = nf90_close(nc_id)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to close lut_nc4 file: ", status
        return
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
