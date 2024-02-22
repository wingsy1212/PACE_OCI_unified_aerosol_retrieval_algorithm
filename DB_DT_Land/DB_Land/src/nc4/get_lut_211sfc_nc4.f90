subroutine get_lut_211sfc( lut_file, status)
    !
    !f90
    !description:
    !   This subroutine reads LUTs for 2.1um sfc refl. calculation
    !   file= 'nvalx21um4sfc.hdf'
    !   array dimension [10,46,30,4]

    !   ver. 1.0  modified by MJ (aug 2011)
    !
    !   use GeneralAuxType
    !   use lut_arrays

    use netcdf
    USE OCIUAAER_Config_Module

    implicit none
    include '../sfc21tbl90.inc'
   
    character(*),           intent (in)      :: lut_file
    integer,                intent (out)     :: status
    integer, dimension (4)                  :: start, edges, stride

    integer                ::  number_type, nattrs
    integer                ::  sds_id, sds_index, attr_index, hdfid
    character(len=255)     ::  sds_name
    character(len=255)    ::  dset_name
    character(len=255)    ::  attr_name
    character(len=255)    ::  group_name

    integer               ::  nc_id
    integer               ::  dim_id
    integer               ::  dset_id
    integer               ::  grp_id

    start  = (/ 1,1,1,1 /)
    edges   = (/ 10,46,30,4 /)
    stride = (/ 1,1,1,1 /)

    status = nf90_open(cfg%db_nc4, nf90_nowrite, nc_id)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to open deepblue lut_nc4 file: ", status
        return
    end if

    group_name = 'veg_21sfc'
    status = nf90_inq_ncid(nc_id, group_name, grp_id)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to get ID of group "//trim(group_name)//": ", status
        return
    end if

    dset_name = 'NVALX21_SFC'
    status = nf90_inq_varid(grp_id, dset_name, dset_id)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
        return
    end if
    status = nf90_get_var(grp_id, dset_id, nvalx21, start=start, &
                          stride=stride, count=edges)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
        return
    end if

    dset_name = 'R0X21_SFC'
    status = nf90_inq_varid(grp_id, dset_name, dset_id)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
        return
    end if
    status = nf90_get_var(grp_id, dset_id, r0x_21, start=start, &
                          stride=stride, count=edges)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
        return
    end if

    dset_name = 'SX21_SFC'
    status = nf90_inq_varid(grp_id, dset_name, dset_id)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
        return
    end if
    status = nf90_get_var(grp_id, dset_id, sx_21, start=start, &
                          stride=stride, count=edges)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
        return
    end if

    dset_name = 'TX21_SFC'
    status = nf90_inq_varid(grp_id, dset_name, dset_id)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
        return
    end if
    status = nf90_get_var(grp_id, dset_id, tx_21, start=start, &
                          stride=stride, count=edges)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
        return
    end if

    dset_name = 'NVALX672_SFC'
    status = nf90_inq_varid(grp_id, dset_name, dset_id)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
        return
    end if
    status = nf90_get_var(grp_id, dset_id, nvalx672, start=start, &
                          stride=stride, count=edges)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
        return
    end if

    dset_name = 'R0X672_SFC'
    status = nf90_inq_varid(grp_id, dset_name, dset_id)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
        return
    end if
    status = nf90_get_var(grp_id, dset_id, r0x_672, start=start, &
                          stride=stride, count=edges)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
        return
    end if

    dset_name = 'SX672_SFC'
    status = nf90_inq_varid(grp_id, dset_name, dset_id)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
        return
    end if
    status = nf90_get_var(grp_id, dset_id, sx_672, start=start, &
                          stride=stride, count=edges)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
        return
    end if

    dset_name = 'TX672_SFC'
    status = nf90_inq_varid(grp_id, dset_name, dset_id)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
        return
    end if
    status = nf90_get_var(grp_id, dset_id, tx_672, start=start, &
                          stride=stride, count=edges)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
        return
    end if

    dset_name = 'NVALX865_SFC'
    status = nf90_inq_varid(grp_id, dset_name, dset_id)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
        return
    end if
    status = nf90_get_var(grp_id, dset_id, nvalx865, start=start, &
                          stride=stride, count=edges)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
        return
    end if

    dset_name = 'R0X865_SFC'
    status = nf90_inq_varid(grp_id, dset_name, dset_id)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
        return
    end if
    status = nf90_get_var(grp_id, dset_id, r0x_865, start=start, &
                          stride=stride, count=edges)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
        return
    end if

    dset_name = 'SX865_SFC'
    status = nf90_inq_varid(grp_id, dset_name, dset_id)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
        return
    end if
    status = nf90_get_var(grp_id, dset_id, sx_865, start=start, &
                          stride=stride, count=edges)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
        return
    end if

    dset_name = 'TX865_SFC'
    status = nf90_inq_varid(grp_id, dset_name, dset_id)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
        return
    end if
    status = nf90_get_var(grp_id, dset_id, tx_865, start=start, &
                          stride=stride, count=edges)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
        return
    end if

    status = nf90_close(nc_id)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to close lut_nc4 file: ", status
        return
    end if

end subroutine get_lut_211sfc

