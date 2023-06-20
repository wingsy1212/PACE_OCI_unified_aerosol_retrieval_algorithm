module landcover

  implicit none
  
! -- restrict access to only module methods listed as public.
  private
  public  ::  load_landcover, unload_landcover, get_landcover


  integer, dimension(:,:), allocatable  ::  lcdata
  
  contains
  
  integer function load_landcover(lc_file) result(status)
    implicit none
    
    include 'hdf.f90'
	  include 'dffunc.f90'
  
    character(len=*), intent(in)  ::  lc_file

    character(len=255)    ::  sds_name
	  integer					      ::  sd_id, sds_id, sds_index
		integer, dimension(2) ::  start2, stride2, dims2
		integer               ::  rank, dtype, nattrs

    sd_id = sfstart(lc_file, DFACC_READ)
    if (sd_id == -1 ) then
    	print *,"ERROR: failed to start SDS interface on look up table file."
    	status = sd_id
    	return
    end if
  
 		sds_name = 'Land_Vegetation_Type'
 		sds_index = sfn2index(sd_id, sds_name)
		if (sds_index == -1 ) then
			print *,"ERROR: failed to convert land cover SDS name to index."
			status = sds_index
			return
		end if
							
		sds_id = sfselect(sd_id, sds_index)
		if ( sds_id == -1) then
			print *,"ERROR: failed to select ps SDS."
			status = sds_id
			return
		end if
		
		status = sfginfo(sds_id, sds_name, rank, dims2, dtype, nattrs)
		if (status /= 0) then
		  print *, "ERROR: Unable to get info on land cover SDS: ", status
		  return
		end if
    
    allocate(lcdata(dims2(1), dims2(2)), stat=status)
    if (status /= 0) then
      print *, "ERROR: Unable to allocate lcdata array: ", status
      return
    end if
    
    start2 = (/0,0/)
    stride2 = (/1,1/)
		status = sfrdata(sds_id, start2, stride2, dims2, lcdata)
		if ( status == FAIL ) then
			print *,"ERROR: failed to read land cover SDS: ", status
			return
		end if
	  
		status = sfendacc(sds_id)
		if ( status == FAIL ) then
			print *,"ERROR: failed to close land cover SDS."
			return
		end if
		
		status = sfend(sd_id)
    if ( status == FAIL ) then
  	  print *,"ERROR: failed to close land cover file."
  	  return
    end if
  end function load_landcover
  
  subroutine unload_landcover(status)
    implicit none
    
    integer, intent(inout) :: status
    
    deallocate(lcdata, stat=status)
    if (status /= 0) then
      print *, "ERROR: Unable to deallocate land cover array: ", status
      return
    end if
  
  end subroutine unload_landcover
  
  integer function get_landcover(lat,lon,status) result(lc)
    implicit none
    
    real, intent(in)        ::  lat
    real, intent(in)        ::  lon
    integer, intent(inout)  ::  status
    
    integer                 ::  ilat,ilon
    
    status = 0
    
    ilat = floor(lat * 20.0) + 1800 + 1
    ilon = floor(lon * 20.0) + 3600 + 1
    if (lon == 180.0 .AND. ilon == 3601) ilon = 1
    if (ilat > size(lcdata,1) .OR. ilon > size(lcdata,2)) then
      print *, "ERROR: Indices out of bounds for land cover array: ", ilat, ilon
      status = -1
      return
    end if
    lc = lcdata(ilat,ilon)
    
    return
  end function get_landcover
  
    
end module landcover