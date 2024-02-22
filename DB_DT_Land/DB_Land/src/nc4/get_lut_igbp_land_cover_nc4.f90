subroutine get_lut_igbp_land_cover( lut_file, status)
    !
    !f90
    !description:
    !   This subroutine reads core SDSs in the IGBP Land Cover
    !   LUT hdf file.
    !   file= 'MCD12C1.A2004001.005.Global_IGBP_Land_Cover_0.10deg.hdf'
    !   includes:       IGBP_Land_Cover & Region_Index
    !   array dimension [3600, 1800]

    !   ver. 1.0  written by CES (jul 2010)
    !   ver. 1.1  modified by MJ (aug 2011)
    !
    !   use GeneralAuxType
    !   use lut_arrays
   
    use netcdf
    USE OCIUAAER_Config_Module

    implicit none
    include '../newaottbl90.inc'
   
    character(*),           intent (in)   :: lut_file
    integer,                intent (out)  :: status
    integer, dimension (2)                :: start, edge, stride

    integer               :: number_type, nattrs
    integer               :: sds_id,sds_index,attr_index, hdfid
    character(len=255)    ::  sds_name
    character(len=255)    ::  dset_name
    character(len=255)    ::  attr_name
    character(len=255)    ::  group_name

    integer               ::  nc_id
    integer               ::  dim_id
    integer               ::  dset_id
    integer               ::  grp_id

    start  = (/ 1,1 /)
    edge   = (/ 3600,1800 /)
    stride = (/ 1,1 /)
    
    status = nf90_open(cfg%db_nc4, nf90_nowrite, nc_id)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to open deepblue lut_nc4 file: ", status
        return
    end if

    group_name = 'veg_land_cover'
    status = nf90_inq_ncid(nc_id, group_name, grp_id)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to get ID of group "//trim(group_name)//": ", status
        return
    end if

    dset_name = 'IGBP_Land_Cover'
    status = nf90_inq_varid(grp_id, dset_name, dset_id)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
        return
    end if
    status = nf90_get_var(grp_id, dset_id, xlcvr_2, start=start, &
                          stride=stride, count=edge)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
        return
    end if

    dset_name = 'Region_Index'
    status = nf90_inq_varid(grp_id, dset_name, dset_id)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
        return
    end if
    status = nf90_get_var(grp_id, dset_id, regid_2, start=start, &
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

end subroutine get_lut_igbp_land_cover

