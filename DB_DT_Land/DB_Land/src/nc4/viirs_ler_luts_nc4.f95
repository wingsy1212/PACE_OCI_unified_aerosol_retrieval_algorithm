module viirs_ler_luts

!private

public    ::  load_viirs_ler_luts
!public    ::  unload_viirs_ler_luts

real, dimension(20800)    ::  logi0
real, dimension(20800)    ::  z1i0
real, dimension(20800)    ::  z2i0
real, dimension(20800)    ::  ti0
real, dimension(260)      ::  sb
real, dimension(160)      ::  li0r
real, dimension(160)      ::  z1i0r
real, dimension(160)      ::  z2i0r
real, dimension(160)      ::  ti0r
real, dimension(2)        ::  sbr

real, dimension(10)       ::  xxzlog = (/0.0000000, 0.00977964, 0.0395086, 0.0904221, &
&                                        0.164818, 0.266515, 0.401776, 0.581261,      &
&                                        0.824689, 1.17436/)
!     -- these numbers correspond to satellite zenith angle
!     -- node points of 0.,16.,30.,40.,48.,54.,58.,60. degrees
real, dimension(8)       ::  xxlog = (/0.000000,0.0395086,0.143841,0.266515,0.401776, &
&                                       0.531394,0.635031,0.693147/)

real, dimension(10)       ::  xzlog
real, dimension(8)        ::  xlog
real, dimension(4,7)      ::  densol
real, dimension(4,5)      ::  denscn
real, dimension(4)        ::  cthet0
real, dimension(4)        ::  ctheta
real, dimension(16)       ::  cofs
integer                   ::  indsol
integer                   ::  indscn
integer                   ::  iofset
real                      ::  p1
real                      ::  pr
real, dimension(10)       ::  dum10

common /lpoly/ xzlog, xlog, densol, denscn, cthet0, ctheta, cofs, &
&              indsol, indscn, iofset, p1, pr, dum10
              
real :: xlat
real :: xlong
real :: sza
real :: xthet
real :: xphi
real :: cphi
real :: c2phi
real :: pteran
real :: xnvalm
integer ::  isnow
integer ::  partial
real :: so2ind
real :: resn
real :: sens
real :: rsens
real :: ozbst
real :: ref
real :: estozn
real :: ozcld
real :: pcloud
real :: prfrac
real :: clfrac
real :: rayval
real :: r412
real :: r470
real :: sfref412
real :: sfref470
real :: sfref650
real :: qdif412
real :: qdif470
real :: qdif650
real :: stdv

common /sample/ xlat,xlong,sza,xthet,xphi,cphi,c2phi,pteran,    &
& xnvalm(6),isnow,so2ind,resn(5),sens(5),                  &
& rsens,ozbst,ref,estozn,ozcld,pcloud,prfrac,clfrac,partial, &
& rayval(5),r412,r470,sfref412,sfref470,sfref650,qdif412,       &
& qdif470,qdif650, stdv
        
real, dimension(26) :: realbuf
common /bufout/ realbuf
                 
common /table/  logi0, z1i0, z2i0, ti0, sb, li0r, z1i0r, z2i0r, ti0r, sbr

contains

subroutine load_viirs_ler_luts(lut_ler_file, status)
    use netcdf
    USE OCIUAAER_Config_Module

    implicit none

    character(len=255), intent(in)  ::  lut_ler_file

    real        ::  xdenom
    integer     ::  i, j, k
    integer     ::  status

    character(len=255)    ::  sds_name
    integer               ::  sd_id, sds_id, sds_index
    integer, dimension(1) ::  start1, stride1, dims1
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

    status = nf90_open(trim(cfg%db_nc4), nf90_nowrite, nc_id)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to open deepblue lut_nc4 file: ", status
        return
    end if

    group_name = 'viirs_ler'
    status = nf90_inq_ncid(nc_id, group_name, grp_id)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to get ID of group "//trim(group_name)//": ", status
        return
    end if

    start1=(/1/)
    stride1=(/1/)

    dset_name = 'logi0'
    status = nf90_inq_varid(grp_id, dset_name, dset_id)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
        return
    end if
    status = nf90_inquire_variable(grp_id, dset_id, dimids=dimids)
    status = nf90_inquire_dimension(grp_id, dimids(1), len = dims1(1))
    status = nf90_get_var(grp_id, dset_id, logi0, start=start1, &
                            stride=stride1, count=dims1)
    err_msg = NF90_STRERROR(status)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
        return
    end if

    dset_name = 'z1i0'
    status = nf90_inq_varid(grp_id, dset_name, dset_id)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
        return
    end if
    status = nf90_inquire_variable(grp_id, dset_id, dimids=dimids)
    status = nf90_inquire_dimension(grp_id, dimids(1), len = dims1(1))
    status = nf90_get_var(grp_id, dset_id, z1i0, start=start1, &
                            stride=stride1, count=dims1)
    err_msg = NF90_STRERROR(status)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
        return
    end if

    dset_name = 'z2i0'
    status = nf90_inq_varid(grp_id, dset_name, dset_id)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
        return
    end if
    status = nf90_inquire_variable(grp_id, dset_id, dimids=dimids)
    status = nf90_inquire_dimension(grp_id, dimids(1), len = dims1(1))
    status = nf90_get_var(grp_id, dset_id, z2i0, start=start1, &
                            stride=stride1, count=dims1)
    err_msg = NF90_STRERROR(status)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
        return
    end if

    dset_name = 'ti0'
    status = nf90_inq_varid(grp_id, dset_name, dset_id)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
        return
    end if
    status = nf90_inquire_variable(grp_id, dset_id, dimids=dimids)
    status = nf90_inquire_dimension(grp_id, dimids(1), len = dims1(1))
    status = nf90_get_var(grp_id, dset_id, ti0, start=start1, &
                            stride=stride1, count=dims1)
    err_msg = NF90_STRERROR(status)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
        return
    end if

    dset_name = 'sb'
    status = nf90_inq_varid(grp_id, dset_name, dset_id)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
        return
    end if
    status = nf90_inquire_variable(grp_id, dset_id, dimids=dimids)
    status = nf90_inquire_dimension(grp_id, dimids(1), len = dims1(1))
    status = nf90_get_var(grp_id, dset_id, sb, start=start1, &
                            stride=stride1, count=dims1)
    err_msg = NF90_STRERROR(status)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
        return
    end if

    dset_name = 'logi0r'
    status = nf90_inq_varid(grp_id, dset_name, dset_id)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
        return
    end if
    status = nf90_inquire_variable(grp_id, dset_id, dimids=dimids)
    status = nf90_inquire_dimension(grp_id, dimids(1), len = dims1(1))
    status = nf90_get_var(grp_id, dset_id, li0r, start=start1, &
                            stride=stride1, count=dims1)
    err_msg = NF90_STRERROR(status)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
        return
    end if

    dset_name = 'z1i0r'
    status = nf90_inq_varid(grp_id, dset_name, dset_id)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
        return
    end if
    status = nf90_inquire_variable(grp_id, dset_id, dimids=dimids)
    status = nf90_inquire_dimension(grp_id, dimids(1), len = dims1(1))
    status = nf90_get_var(grp_id, dset_id, z1i0r, start=start1, &
                            stride=stride1, count=dims1)
    err_msg = NF90_STRERROR(status)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
        return
    end if

    dset_name = 'z2i0r'
    status = nf90_inq_varid(grp_id, dset_name, dset_id)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
        return
    end if
    status = nf90_inquire_variable(grp_id, dset_id, dimids=dimids)
    status = nf90_inquire_dimension(grp_id, dimids(1), len = dims1(1))
    status = nf90_get_var(grp_id, dset_id, z2i0r, start=start1, &
                            stride=stride1, count=dims1)
    err_msg = NF90_STRERROR(status)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
        return
    end if

    dset_name = 'ti0r'
    status = nf90_inq_varid(grp_id, dset_name, dset_id)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
        return
    end if
    status = nf90_inquire_variable(grp_id, dset_id, dimids=dimids)
    status = nf90_inquire_dimension(grp_id, dimids(1), len = dims1(1))
    status = nf90_get_var(grp_id, dset_id, ti0r, start=start1, &
                            stride=stride1, count=dims1)
    err_msg = NF90_STRERROR(status)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
        return
    end if

    dset_name = 'sbr'
    status = nf90_inq_varid(grp_id, dset_name, dset_id)
    if (status /= NF90_NOERR) then
        print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
        return
    end if
    status = nf90_inquire_variable(grp_id, dset_id, dimids=dimids)
    status = nf90_inquire_dimension(grp_id, dimids(1), len = dims1(1))
    status = nf90_get_var(grp_id, dset_id, sbr, start=start1, &
                            stride=stride1, count=dims1)
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

  xzlog(:) = xxzlog(:)
  xlog(:) = xxlog(:)
  
! -- calculate values needed in table interpolation
  do j = 1, 7
    do k = j, j + 3
      xdenom = 1.0
      do i = j, j + 3
        if (i .ne. k) xdenom = xdenom * (xzlog(k) - xzlog(i))
      end do
      densol(k - j + 1, j) = xdenom
    end do
  end do
  do j = 1, 5
    do k = j, j + 3
      xdenom = 1.0
      do i = j, j + 3
        if (i .ne. k) xdenom = xdenom * (xlog(k) - xlog(i))
      end do
      denscn(k - j + 1, j) = xdenom
    end do
  end do
      
  return
end subroutine load_viirs_ler_luts

integer function wrap_total() result(status)
  implicit none
  
  status = 0
  return
  
end function wrap_total

end module viirs_ler_luts
