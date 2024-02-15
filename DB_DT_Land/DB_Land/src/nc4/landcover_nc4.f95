module landcover

  implicit none
  
! -- restrict access to only module methods listed as public.
  private
  public  ::  load_landcover, unload_landcover, get_landcover


  integer, dimension(:,:), allocatable  ::  lcdata
  
  contains
  
  integer function load_landcover(lc_file) result(status)
    
!   include 'hdf.f90'
!   include 'dffunc.f90'
    use netcdf
    USE OCIUAAER_Config_Module

    implicit none
  
    character(len=255), intent(in)  ::  lc_file

    character(len=255)    ::  sds_name
    integer               ::  sd_id, sds_id, sds_index
    integer, dimension(2) ::  start2, stride2, dims2
    integer, dimension(32)::  dimids
    integer               ::  rank, dtype, nattrs

    integer               ::  nc_id
    integer               ::  dim_id
    integer               ::  dset_id
    integer               ::  grp_id

    character(len=255)    ::  dset_name
    character(len=255)    ::  attr_name
    character(len=255)    ::  group_name
    character(len=255)    ::  err_msg

    lc_file = cfg%db_nc4
    status = nf90_open(lc_file, nf90_nowrite, nc_id)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to open deepblue lut_nc4 file: ", status
        return
    end if

    group_name = 'land_cover'
    status = nf90_inq_ncid(nc_id, group_name, grp_id)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to get ID of group "//trim(group_name)//": ", status
        return
    end if

    dset_name = 'Land_Vegetation_Type'
    status = nf90_inq_varid(grp_id, dset_name, dset_id)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
        return
    end if
    status = nf90_inquire_variable(grp_id, dset_id, dimids=dimids)
    status = nf90_inquire_dimension(grp_id, dimids(1), len = dims2(1))
    status = nf90_inquire_dimension(grp_id, dimids(2), len = dims2(2))
    allocate(lcdata(dims2(1), dims2(2)), stat=status)
    if (status /= 0) then
      print *, "ERROR: Unable to allocate lcdata array: ", status
      return
    end if

    start2=(/1,1/)
    stride2=(/1,1/)
    status = nf90_get_var(grp_id, dset_id, lcdata, start=start2, &
                            stride=stride2, count=dims2)
    err_msg = NF90_STRERROR(status)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
        return
    end if

    status = nf90_close(nc_id)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to close lut_nc4 file: ", status
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
