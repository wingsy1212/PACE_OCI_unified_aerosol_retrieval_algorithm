
! Draft netCDF reader for interm file
! Simplified to fewer variables for testing
! Yingxi
! Created 2022-08-08

 module read_intrim_forUVtau
 
! Declare use of Fortran 90 netCDF package
! If this part doesn't compile: module load sips/gcc

use netcdf

implicit none

contains

! ---------------------------------------------------------------------
!
! Begin subroutine make_nc(ALL VARIABLES GO HERE,XL,YL)
!
! ---------------------------------------------------------------------

      subroutine read_Output_forUVtau(Ret_Xtrack,Ret_Lines,uvdbdtaod)
                   
      include 'output_Variables.inc' 
!     Declare input variables

       real(8) :: scale_factor, add_offset, cprime
       character(len=200) :: SDS_name
       integer ::  i, j, zl3, k ,length
       integer(2) :: fill_val
       integer(2), parameter :: fillvalue_int2 = -9999
       integer(2), dimension(:,:,:), allocatable :: temp
       REAL,DIMENSION(:,:,:),ALLOCATABLE :: uvdbdtaod
       integer :: file_id, grp_id,var_id, err_code, start(2), stride(2), edge(2), attr_id
       integer :: xtype, ndims,nAtts,deflate_level,endianness,quantize_mode,nsd
       integer(4), dimension(:),allocatable :: dimids
!       logical :: contiguous,shuffle,fletcher32
!       integer :: nvars
!       integer(4), dimension(:),allocatable :: var_ids

       character (len=*), parameter :: nc_name = 'Interm_file.nc'
      call check(nf90_open(nc_name,NF90_NOWRITE,file_id))
      SDS_name='Optical_Depth_Land'
!     Open the netCDF file already created via ncgen
      print *,nf90_noerr
      ZL3 = 5
!      stride = 1
!      edge(1) = Ret_Xtrack
!      edge(2) = Ret_Lines
      allocate(temp(zl3,Ret_Xtrack,Ret_Lines))
!      allocate(uvdbdtaod(zl3,Ret_Xtrack,Ret_Lines))
      err_code = nf90_inq_ncid(file_id,'geophysical_data',grp_id)
      !print *,grp_id
      !err_code = nf90_inquire(grp_id, nvariables=nvars)
      !allocate(var_ids(nvars))
      !err_code = nf90_inq_varids(grp_id, nvars, var_ids)
      err_code = nf90_inq_varid(grp_id, trim(SDS_name), var_id)
      err_code = nf90_inquire_variable(grp_id,var_id,ndims=ndims)
      allocate(dimids(ndims))
      err_code = nf90_inquire_variable(grp_id,var_id,dimids=dimids)
!      print *,dimids
      do i = 1, ndims
          err_code = nf90_inquire_dimension(file_id,dimids(i),len=length)
!          print *,length
      enddo
      err_code = nf90_get_var(grp_id, var_id, temp)
      call handle_err(err_code)
      
      err_code = nf90_get_att(grp_id, var_id, "scale_factor", scale_factor)
      err_code = nf90_get_att(grp_id, var_id, "add_offset", add_offset)
      err_code = nf90_get_att(grp_id, var_id, "_FillValue", fill_val)
      call handle_err(err_code)
      if (err_code /= 0) fill_val = fillvalue_int2
!      print *,'fill: ',fill_val
!      err_code = nf90_get_att(file_id, var_id, "cprime", cprime)
      do j=1,Ret_Lines
         do i=1,Ret_Xtrack
            do k=1,zl3     
                if (temp(k,i,j) /= fill_val .AND. temp(k,i,j) /= 28008 ) then 
!                    if (temp(k,i,j) .gt. 0) then
!                        print *,temp(k,i,j)* scale_factor + add_offset
!                    endif
                !print *, '',temp(k,i,j),scale_factor,add_offset
                    uvdbdtaod(i,j,k) = (temp(k,i,j) * scale_factor + add_offset)
                else            
                    uvdbdtaod(i,j,k) = -999.
                endif
            end do
         end do
      end do
      
      deallocate(temp)
      deallocate(dimids)
          
      call check( nf90_close(file_id) )
      print *, 'Read variables from ',nc_name
      return
      end subroutine read_Output_forUVtau

      subroutine check(status)

      integer, intent ( in) :: status

      if(status /= nf90_noerr) then
         print *, trim(nf90_strerror(status))
         stop "Stopped"
      end if

      end subroutine check

      subroutine handle_err(status)
          integer, intent ( in) :: status
          if(status /= nf90_noerr) then
          print *, trim(nf90_strerror(status))
          stop "Stopped"
          end if
      end subroutine handle_err     
end module read_intrim_forUVtau






