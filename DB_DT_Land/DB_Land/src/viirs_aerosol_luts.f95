module viirs_aerosol_luts

private

public    ::  load_viirs_aerosol_luts
public    ::  unload_viirs_aerosol_luts

public    ::  aero_412, aero_470, aero_650
public    ::  aero_412_abs, aero_470_abs
public    ::  aero_412_dust, aero_470_dust, aero_650_dust
public    ::  aero_412_abs_dust, aero_470_abs_dust
public    ::  new_intep, polint

integer, parameter		:: 	SGL = SELECTED_REAL_KIND(p=3)

real, dimension(30)             ::  lut_raa
real, dimension(46)             ::  lut_vza

real, dimension(12,46,30,10,8,20) ::  nvalx412
real, dimension(12,46,30,10,4,24) ::  nvalx470
real, dimension(12,46,30,10,24)   ::  nvalx650
    
type  ::  viirs_aerosol_lut

  integer                                   ::  nsza
  integer                                   ::  nvza
  integer                                   ::  nraa
  integer                                   ::  naot
  integer                                   ::  nssa
  integer                                   ::  nsfc
  
  real, dimension(:), allocatable           ::  sza
  real, dimension(:), allocatable           ::  vza
  real, dimension(:), allocatable           ::  raa
  real, dimension(:), allocatable           ::  aot
  real, dimension(:), allocatable           ::  ssa
  real, dimension(:), allocatable           ::  sfc
  
  real, dimension(:,:,:,:,:,:), allocatable ::  nvalx
end type viirs_aerosol_lut

type(viirs_aerosol_lut)    ::  default_lut412, default_lut488, default_lut672 
type(viirs_aerosol_lut)    ::  dust_lut412, dust_lut488, dust_lut672

common /aottbl/ nvalx412, nvalx(12,46,30,10),             &
&               theta0(12), theta(46),                                      &
&               phi(30), sfc_ref412(20), tau(10), w0(8), w0_470(4),         &
&               nvalx650, sfc_ref650(24), sfc_ref470(24),   &
&               sfcprs(360,720), nvalx470 

contains

integer function load_viirs_aerosol_luts(lut_land_file, lut_dust_file) result(status)
  implicit none

  character(len=*), intent(in)    ::  lut_land_file
  character(len=*), intent(in)    ::  lut_dust_file
        
!  print *, 'Reading land LUT...'
  status = read_aerosol_lut_file(lut_land_file, default_lut412, default_lut488, default_lut672)
  if (status /= 0) then 
    print *, 'ERROR: Failed to read in default land aerosol model LUT: ', status
    print *, 'File: ', trim(lut_land_file)
    return
  end if
!  print *, 'done.'
  
!  print *, 'Reading dust LUT...'
  status = read_aerosol_lut_file(lut_dust_file, dust_lut412, dust_lut488, dust_lut672)
  if (status /= 0) then 
    print *, 'ERROR: Failed to read in dust aerosol model LUT: ', status
!    print *, 'File: ', trim(lut_dust_file)
    return
  end if
!  print *, 'done.'

! -- translate new LUT node arrays to old variabled in common block  
! -- copy over the data to the deprecated common block 
! -- still used in the vegetated retrieval -- see find_v_vegset.f
  nvalx412 = default_lut412%nvalx(:,:,:,:,:,:)
  nvalx470 = default_lut488%nvalx(:,:,:,:,:,:)
  nvalx650 = default_lut672%nvalx(:,:,:,:,1,:)
  
  tau         = default_lut412%aot
  theta0      = default_lut412%sza
  phi         = default_lut412%raa
  theta       = default_lut412%vza
  w0          = default_lut412%ssa
  
  w0_470      = default_lut488%ssa 
  sfc_ref412  = default_lut412%sfc 
  sfc_ref470  = default_lut488%sfc 
  sfc_ref650  = default_lut672%sfc
  
  return
  
end function load_viirs_aerosol_luts
  
subroutine unload_viirs_aerosol_luts(status)
  implicit none
  
  integer,intent(inout) ::  status
  
  status = 0 
!  deallocate(nvalx412, nvalx470, nvalx650, stat=status)
!  if (status /= 0) then
!    print *, "ERROR: Failed to deallocate aerosol LUT arrays: ", status
!    return
!  end if

end subroutine unload_viirs_aerosol_luts

integer function read_aerosol_lut_file(lut_file, lut412, lut488, lut672) result(status)
  implicit none
  
  include 'hdf.inc'
  include 'dffunc.inc'
  
  character(len=*), intent(in)            ::  lut_file 
  type(viirs_aerosol_lut), intent(inout)  ::  lut412
  type(viirs_aerosol_lut), intent(inout)  ::  lut488
  type(viirs_aerosol_lut), intent(inout)  ::  lut672
  
  character(len=255)  ::  sds_name
  integer :: sd_id, sds_index, sds_id, rank, n_attrs, data_type
  integer, dimension(1) ::  start1, stride1, edges1, dim_sizes1
  integer, dimension(5) ::  edges5, start5, stride5, dim_sizes5
  integer, dimension(6) ::  edges6, start6, stride6, dim_sizes6
  integer, dimension(32)::  dim_sizes
  
  integer                         ::  i

  sd_id = sfstart(lut_file, DFACC_READ)
  if (sd_id == FAIL) then
    print *, "ERROR: Unable to start SD interface on land aerosol model table: ", sd_id
    status = -1
    return
  end if
  
! -- 412nm
  sds_name = 'SZA412_Nodes'
  sds_index = sfn2index(sd_id, sds_name)
  if (sds_index == FAIL) then
    print *, "ERROR: Unable to find index of "//trim(sds_name)//": ", sds_index
    status = -1
    return
  end if
  
  sds_id = sfselect(sd_id, sds_index)
  if (sds_id == FAIL) then
    print *, "ERROR: Unable to select SDS "//trim(sds_name)//": ", sds_id
    status = -1
    return
  end if
  
  status = sfginfo(sds_id, sds_name, rank, dim_sizes, data_type, n_attrs)
  if (status == FAIL) then
    print *, "ERROR: Unable to get info on SDS "//trim(sds_name)//": ", status
    status = -1
    return
  end if
  
  lut412%nsza = dim_sizes(1)
  allocate(lut412%sza(lut412%nsza), stat=status)
  if (status /= 0) then
    print *, "ERROR: Unable to allocate SZA412 data array: ", status
    return
  end if
  
  start1  = (/0/)
  stride1 = (/1/)
  edges1  = (/lut412%nsza/)
  status = sfrdata(sds_id, start1, stride1, edges1, lut412%sza)
  if (status == FAIL) then
    print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
    return
  end if
  
  status = sfendacc(sds_id)
  if (status == FAIL) then
    print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
    return
  end if
  
  sds_name = 'VZA412_Nodes'
  sds_index = sfn2index(sd_id, sds_name)
  if (sds_index == FAIL) then
    print *, "ERROR: Unable to find index of "//trim(sds_name)//": ", sds_index
    status = -1
    return
  end if
  
  sds_id = sfselect(sd_id, sds_index)
  if (sds_id == FAIL) then
    print *, "ERROR: Unable to select SDS "//trim(sds_name)//": ", sds_id
    status = -1
    return
  end if
  
  status = sfginfo(sds_id, sds_name, rank, dim_sizes, data_type, n_attrs)
  if (status == FAIL) then
    print *, "ERROR: Unable to get info on SDS "//trim(sds_name)//": ", status
    status = -1
    return
  end if
  
  lut412%nvza = dim_sizes(1)
  allocate(lut412%vza(lut412%nvza), stat=status)
  if (status /= 0) then
    print *, "ERROR: Unable to allocate VZA412 data array: ", status
    return
  end if
  
  start1  = (/0/)
  stride1 = (/1/)
  edges1  = (/lut412%nvza/)
  status = sfrdata(sds_id, start1, stride1, edges1, lut412%vza)
  if (status == FAIL) then
    print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
    return
  end if
  
  status = sfendacc(sds_id)
  if (status == FAIL) then
    print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
    return
  end if
  
  sds_name = 'RAA412_Nodes'
  sds_index = sfn2index(sd_id, sds_name)
  if (sds_index == FAIL) then
    print *, "ERROR: Unable to find index of "//trim(sds_name)//": ", sds_index
    status = -1
    return
  end if
  
  sds_id = sfselect(sd_id, sds_index)
  if (sds_id == FAIL) then
    print *, "ERROR: Unable to select SDS "//trim(sds_name)//": ", sds_id
    status = -1
    return
  end if
  
  status = sfginfo(sds_id, sds_name, rank, dim_sizes, data_type, n_attrs)
  if (status == FAIL) then
    print *, "ERROR: Unable to get info on SDS "//trim(sds_name)//": ", status
    status = -1
    return
  end if
  
  lut412%nraa = dim_sizes(1)
  allocate(lut412%raa(lut412%nraa), stat=status)
  if (status /= 0) then
    print *, "ERROR: Unable to allocate RAA412 data array: ", status
    return
  end if
  
  start1  = (/0/)
  stride1 = (/1/)
  edges1  = (/lut412%nraa/)
  status = sfrdata(sds_id, start1, stride1, edges1, lut412%raa)
  if (status == FAIL) then
    print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
    return
  end if
  
  status = sfendacc(sds_id)
  if (status == FAIL) then
    print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
    return
  end if
  
  sds_name = 'AOT412_Nodes'
  sds_index = sfn2index(sd_id, sds_name)
  if (sds_index == FAIL) then
    print *, "ERROR: Unable to find index of "//trim(sds_name)//": ", sds_index
    status = -1
    return
  end if
  
  sds_id = sfselect(sd_id, sds_index)
  if (sds_id == FAIL) then
    print *, "ERROR: Unable to select SDS "//trim(sds_name)//": ", sds_id
    status = -1
    return
  end if
  
  status = sfginfo(sds_id, sds_name, rank, dim_sizes, data_type, n_attrs)
  if (status == FAIL) then
    print *, "ERROR: Unable to get info on SDS "//trim(sds_name)//": ", status
    status = -1
    return
  end if
  
  lut412%naot = dim_sizes(1)
  allocate(lut412%aot(lut412%naot), stat=status)
  if (status /= 0) then
    print *, "ERROR: Unable to allocate AOT412 data array: ", status
    return
  end if
  
  start1  = (/0/)
  stride1 = (/1/)
  edges1  = (/lut412%naot/)
  status = sfrdata(sds_id, start1, stride1, edges1, lut412%aot)
  if (status == FAIL) then
    print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
    return
  end if
  
  status = sfendacc(sds_id)
  if (status == FAIL) then
    print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
    return
  end if
  
  sds_name = 'SSA412_Nodes'
  sds_index = sfn2index(sd_id, sds_name)
  if (sds_index == FAIL) then
    print *, "ERROR: Unable to find index of "//trim(sds_name)//": ", sds_index
    status = -1
    return
  end if
  
  sds_id = sfselect(sd_id, sds_index)
  if (sds_id == FAIL) then
    print *, "ERROR: Unable to select SDS "//trim(sds_name)//": ", sds_id
    status = -1
    return
  end if
  
  status = sfginfo(sds_id, sds_name, rank, dim_sizes, data_type, n_attrs)
  if (status == FAIL) then
    print *, "ERROR: Unable to get info on SDS "//trim(sds_name)//": ", status
    status = -1
    return
  end if
  
  lut412%nssa = dim_sizes(1)
  allocate(lut412%ssa(lut412%nssa), stat=status)
  if (status /= 0) then
    print *, "ERROR: Unable to allocate SSA412 data array: ", status
    return
  end if
  
  start1  = (/0/)
  stride1 = (/1/)
  edges1  = (/lut412%nssa/)
  status = sfrdata(sds_id, start1, stride1, edges1, lut412%ssa)
  if (status == FAIL) then
    print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
    return
  end if
  
  status = sfendacc(sds_id)
  if (status == FAIL) then
    print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
    return
  end if
  
  sds_name = 'SR412_Nodes'
  sds_index = sfn2index(sd_id, sds_name)
  if (sds_index == FAIL) then
    print *, "ERROR: Unable to find index of "//trim(sds_name)//": ", sds_index
    status = -1
    return
  end if
  
  sds_id = sfselect(sd_id, sds_index)
  if (sds_id == FAIL) then
    print *, "ERROR: Unable to select SDS "//trim(sds_name)//": ", sds_id
    status = -1
    return
  end if
  
  status = sfginfo(sds_id, sds_name, rank, dim_sizes, data_type, n_attrs)
  if (status == FAIL) then
    print *, "ERROR: Unable to get info on SDS "//trim(sds_name)//": ", status
    status = -1
    return
  end if
  
  lut412%nsfc = dim_sizes(1)
  allocate(lut412%sfc(lut412%nsfc), stat=status)
  if (status /= 0) then
    print *, "ERROR: Unable to allocate SR412 data array: ", status
    return
  end if
  
  start1  = (/0/)
  stride1 = (/1/)
  edges1  = (/lut412%nsfc/)
  status = sfrdata(sds_id, start1, stride1, edges1, lut412%sfc)
  if (status == FAIL) then
    print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
    return
  end if
  
  status = sfendacc(sds_id)
  if (status == FAIL) then
    print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
    return
  end if
  
  sds_name = 'nvalx412'
  allocate(lut412%nvalx(lut412%nsza,lut412%nvza,lut412%nraa,lut412%naot,lut412%nssa, &
  & lut412%nsfc), stat=status)
  if (status /= 0) then
    print *, "ERROR: Unable to allocate nvalx412 data array: ", status
    return
  end if
  
  sds_index = sfn2index(sd_id, sds_name)
  if (sds_index == FAIL) then
    print *, "ERROR: Unable to find index of SDS "//trim(sds_name)//": ", sds_index
    status = -1
    return
  end if
  
  sds_id = sfselect(sd_id, sds_index)
  if (sds_id == FAIL) then
    print *, "ERROR: Unable to select SDS "//trim(sds_name)//": ", sds_id
    status = -1
    return
  end if
  
  start6  = (/0,0,0,0,0,0/)
  stride6 = (/1,1,1,1,1,1/)
  edges6  = shape(lut412%nvalx)
  status = sfrdata(sds_id, start6, stride6, edges6, lut412%nvalx)
  if (status == FAIL) then
    print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
    return
  end if
  
  status = sfendacc(sds_id)
  if (status == FAIL) then
    print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
    return
  end if
  
! -- 488nm
  sds_name = 'SZA488_Nodes'
  sds_index = sfn2index(sd_id, sds_name)
  if (sds_index == FAIL) then
    print *, "ERROR: Unable to find index of "//trim(sds_name)//": ", sds_index
    status = -1
    return
  end if
  
  sds_id = sfselect(sd_id, sds_index)
  if (sds_id == FAIL) then
    print *, "ERROR: Unable to select SDS "//trim(sds_name)//": ", sds_id
    status = -1
    return
  end if
  
  status = sfginfo(sds_id, sds_name, rank, dim_sizes, data_type, n_attrs)
  if (status == FAIL) then
    print *, "ERROR: Unable to get info on SDS "//trim(sds_name)//": ", status
    status = -1
    return
  end if
  
  lut488%nsza = dim_sizes(1)
  allocate(lut488%sza(lut488%nsza), stat=status)
  if (status /= 0) then
    print *, "ERROR: Unable to allocate SZA488 data array: ", status
    return
  end if
  
  start1  = (/0/)
  stride1 = (/1/)
  edges1  = (/lut488%nsza/)
  status = sfrdata(sds_id, start1, stride1, edges1, lut488%sza)
  if (status == FAIL) then
    print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
    return
  end if
  
  status = sfendacc(sds_id)
  if (status == FAIL) then
    print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
    return
  end if
  
  sds_name = 'VZA488_Nodes'
  sds_index = sfn2index(sd_id, sds_name)
  if (sds_index == FAIL) then
    print *, "ERROR: Unable to find index of "//trim(sds_name)//": ", sds_index
    status = -1
    return
  end if
  
  sds_id = sfselect(sd_id, sds_index)
  if (sds_id == FAIL) then
    print *, "ERROR: Unable to select SDS "//trim(sds_name)//": ", sds_id
    status = -1
    return
  end if
  
  status = sfginfo(sds_id, sds_name, rank, dim_sizes, data_type, n_attrs)
  if (status == FAIL) then
    print *, "ERROR: Unable to get info on SDS "//trim(sds_name)//": ", status
    status = -1
    return
  end if
  
  lut488%nvza = dim_sizes(1)
  allocate(lut488%vza(lut488%nvza), stat=status)
  if (status /= 0) then
    print *, "ERROR: Unable to allocate VZA488 data array: ", status
    return
  end if
  
  start1  = (/0/)
  stride1 = (/1/)
  edges1  = (/lut488%nvza/)
  status = sfrdata(sds_id, start1, stride1, edges1, lut488%vza)
  if (status == FAIL) then
    print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
    return
  end if
  
  status = sfendacc(sds_id)
  if (status == FAIL) then
    print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
    return
  end if
  
  sds_name = 'RAA488_Nodes'
  sds_index = sfn2index(sd_id, sds_name)
  if (sds_index == FAIL) then
    print *, "ERROR: Unable to find index of "//trim(sds_name)//": ", sds_index
    status = -1
    return
  end if
  
  sds_id = sfselect(sd_id, sds_index)
  if (sds_id == FAIL) then
    print *, "ERROR: Unable to select SDS "//trim(sds_name)//": ", sds_id
    status = -1
    return
  end if
  
  status = sfginfo(sds_id, sds_name, rank, dim_sizes, data_type, n_attrs)
  if (status == FAIL) then
    print *, "ERROR: Unable to get info on SDS "//trim(sds_name)//": ", status
    status = -1
    return
  end if
  
  lut488%nraa = dim_sizes(1)
  allocate(lut488%raa(lut488%nraa), stat=status)
  if (status /= 0) then
    print *, "ERROR: Unable to allocate RAA488 data array: ", status
    return
  end if
  
  start1  = (/0/)
  stride1 = (/1/)
  edges1  = (/lut488%nraa/)
  status = sfrdata(sds_id, start1, stride1, edges1, lut488%raa)
  if (status == FAIL) then
    print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
    return
  end if
  
  status = sfendacc(sds_id)
  if (status == FAIL) then
    print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
    return
  end if
  
  sds_name = 'AOT488_Nodes'
  sds_index = sfn2index(sd_id, sds_name)
  if (sds_index == FAIL) then
    print *, "ERROR: Unable to find index of "//trim(sds_name)//": ", sds_index
    status = -1
    return
  end if
  
  sds_id = sfselect(sd_id, sds_index)
  if (sds_id == FAIL) then
    print *, "ERROR: Unable to select SDS "//trim(sds_name)//": ", sds_id
    status = -1
    return
  end if
  
  status = sfginfo(sds_id, sds_name, rank, dim_sizes, data_type, n_attrs)
  if (status == FAIL) then
    print *, "ERROR: Unable to get info on SDS "//trim(sds_name)//": ", status
    status = -1
    return
  end if
  
  lut488%naot = dim_sizes(1)
  allocate(lut488%aot(lut488%naot), stat=status)
  if (status /= 0) then
    print *, "ERROR: Unable to allocate AOT488 data array: ", status
    return
  end if
  
  start1  = (/0/)
  stride1 = (/1/)
  edges1  = (/lut488%naot/)
  status = sfrdata(sds_id, start1, stride1, edges1, lut488%aot)
  if (status == FAIL) then
    print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
    return
  end if
  
  status = sfendacc(sds_id)
  if (status == FAIL) then
    print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
    return
  end if
  
  sds_name = 'SSA488_Nodes'
  sds_index = sfn2index(sd_id, sds_name)
  if (sds_index == FAIL) then
    print *, "ERROR: Unable to find index of "//trim(sds_name)//": ", sds_index
    status = -1
    return
  end if
  
  sds_id = sfselect(sd_id, sds_index)
  if (sds_id == FAIL) then
    print *, "ERROR: Unable to select SDS "//trim(sds_name)//": ", sds_id
    status = -1
    return
  end if
  
  status = sfginfo(sds_id, sds_name, rank, dim_sizes, data_type, n_attrs)
  if (status == FAIL) then
    print *, "ERROR: Unable to get info on SDS "//trim(sds_name)//": ", status
    status = -1
    return
  end if
  
  lut488%nssa = dim_sizes(1)
  allocate(lut488%ssa(lut488%nssa), stat=status)
  if (status /= 0) then
    print *, "ERROR: Unable to allocate SSA488 data array: ", status
    return
  end if
  
  start1  = (/0/)
  stride1 = (/1/)
  edges1  = (/lut488%nssa/)
  status = sfrdata(sds_id, start1, stride1, edges1, lut488%ssa)
  if (status == FAIL) then
    print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
    return
  end if
  
  status = sfendacc(sds_id)
  if (status == FAIL) then
    print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
    return
  end if
  
  sds_name = 'SR488_Nodes'
  sds_index = sfn2index(sd_id, sds_name)
  if (sds_index == FAIL) then
    print *, "ERROR: Unable to find index of "//trim(sds_name)//": ", sds_index
    status = -1
    return
  end if
  
  sds_id = sfselect(sd_id, sds_index)
  if (sds_id == FAIL) then
    print *, "ERROR: Unable to select SDS "//trim(sds_name)//": ", sds_id
    status = -1
    return
  end if
  
  status = sfginfo(sds_id, sds_name, rank, dim_sizes, data_type, n_attrs)
  if (status == FAIL) then
    print *, "ERROR: Unable to get info on SDS "//trim(sds_name)//": ", status
    status = -1
    return
  end if
  
  lut488%nsfc = dim_sizes(1)
  allocate(lut488%sfc(lut488%nsfc), stat=status)
  if (status /= 0) then
    print *, "ERROR: Unable to allocate SR488 data array: ", status
    return
  end if
  
  start1  = (/0/)
  stride1 = (/1/)
  edges1  = (/lut488%nsfc/)
  status = sfrdata(sds_id, start1, stride1, edges1, lut488%sfc)
  if (status == FAIL) then
    print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
    return
  end if
  
  status = sfendacc(sds_id)
  if (status == FAIL) then
    print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
    return
  end if
  
  sds_name = 'nvalx488'
  allocate(lut488%nvalx(lut488%nsza,lut488%nvza,lut488%nraa,lut488%naot,lut488%nssa, &
  & lut488%nsfc), stat=status)
  if (status /= 0) then
    print *, "ERROR: Unable to allocate nvalx488 data array: ", status
    return
  end if
  
  sds_index = sfn2index(sd_id, sds_name)
  if (sds_index == FAIL) then
    print *, "ERROR: Unable to find index of SDS "//trim(sds_name)//": ", sds_index
    status = -1
    return
  end if
  
  sds_id = sfselect(sd_id, sds_index)
  if (sds_id == FAIL) then
    print *, "ERROR: Unable to select SDS "//trim(sds_name)//": ", sds_id
    status = -1
    return
  end if
  
  start6  = (/0,0,0,0,0,0/)
  stride6 = (/1,1,1,1,1,1/)
  edges6  = shape(lut488%nvalx)
  status = sfrdata(sds_id, start6, stride6, edges6, lut488%nvalx)
  if (status == FAIL) then
    print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
    return
  end if
  
  status = sfendacc(sds_id)
  if (status == FAIL) then
    print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
    return
  end if
  
  sds_name = 'SZA672_Nodes'
  sds_index = sfn2index(sd_id, sds_name)
  if (sds_index == FAIL) then
    print *, "ERROR: Unable to find index of "//trim(sds_name)//": ", sds_index
    status = -1
    return
  end if
  
  sds_id = sfselect(sd_id, sds_index)
  if (sds_id == FAIL) then
    print *, "ERROR: Unable to select SDS "//trim(sds_name)//": ", sds_id
    status = -1
    return
  end if
  
  status = sfginfo(sds_id, sds_name, rank, dim_sizes, data_type, n_attrs)
  if (status == FAIL) then
    print *, "ERROR: Unable to get info on SDS "//trim(sds_name)//": ", status
    status = -1
    return
  end if
  
! -- 672 nm
  lut672%nsza = dim_sizes(1)
  allocate(lut672%sza(lut672%nsza), stat=status)
  if (status /= 0) then
    print *, "ERROR: Unable to allocate SZA670 data array: ", status
    return
  end if
  
  start1  = (/0/)
  stride1 = (/1/)
  edges1  = (/lut672%nsza/)
  status = sfrdata(sds_id, start1, stride1, edges1, lut672%sza)
  if (status == FAIL) then
    print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
    return
  end if
  
  status = sfendacc(sds_id)
  if (status == FAIL) then
    print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
    return
  end if
  
  sds_name = 'VZA672_Nodes'
  sds_index = sfn2index(sd_id, sds_name)
  if (sds_index == FAIL) then
    print *, "ERROR: Unable to find index of "//trim(sds_name)//": ", sds_index
    status = -1
    return
  end if
  
  sds_id = sfselect(sd_id, sds_index)
  if (sds_id == FAIL) then
    print *, "ERROR: Unable to select SDS "//trim(sds_name)//": ", sds_id
    status = -1
    return
  end if
  
  status = sfginfo(sds_id, sds_name, rank, dim_sizes, data_type, n_attrs)
  if (status == FAIL) then
    print *, "ERROR: Unable to get info on SDS "//trim(sds_name)//": ", status
    status = -1
    return
  end if
  
  lut672%nvza = dim_sizes(1)
  allocate(lut672%vza(lut672%nvza), stat=status)
  if (status /= 0) then
    print *, "ERROR: Unable to allocate VZA670 data array: ", status
    return
  end if
  
  start1  = (/0/)
  stride1 = (/1/)
  edges1  = (/lut672%nvza/)
  status = sfrdata(sds_id, start1, stride1, edges1, lut672%vza)
  if (status == FAIL) then
    print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
    return
  end if
  
  status = sfendacc(sds_id)
  if (status == FAIL) then
    print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
    return
  end if
  
  sds_name = 'RAA672_Nodes'
  sds_index = sfn2index(sd_id, sds_name)
  if (sds_index == FAIL) then
    print *, "ERROR: Unable to find index of "//trim(sds_name)//": ", sds_index
    status = -1
    return
  end if
  
  sds_id = sfselect(sd_id, sds_index)
  if (sds_id == FAIL) then
    print *, "ERROR: Unable to select SDS "//trim(sds_name)//": ", sds_id
    status = -1
    return
  end if
  
  status = sfginfo(sds_id, sds_name, rank, dim_sizes, data_type, n_attrs)
  if (status == FAIL) then
    print *, "ERROR: Unable to get info on SDS "//trim(sds_name)//": ", status
    status = -1
    return
  end if
  
  lut672%nraa = dim_sizes(1)
  allocate(lut672%raa(lut672%nraa), stat=status)
  if (status /= 0) then
    print *, "ERROR: Unable to allocate RAA670 data array: ", status
    return
  end if
  
  start1  = (/0/)
  stride1 = (/1/)
  edges1  = (/lut672%nraa/)
  status = sfrdata(sds_id, start1, stride1, edges1, lut672%raa)
  if (status == FAIL) then
    print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
    return
  end if
  
  status = sfendacc(sds_id)
  if (status == FAIL) then
    print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
    return
  end if
  
  sds_name = 'AOT672_Nodes'
  sds_index = sfn2index(sd_id, sds_name)
  if (sds_index == FAIL) then
    print *, "ERROR: Unable to find index of "//trim(sds_name)//": ", sds_index
    status = -1
    return
  end if
  
  sds_id = sfselect(sd_id, sds_index)
  if (sds_id == FAIL) then
    print *, "ERROR: Unable to select SDS "//trim(sds_name)//": ", sds_id
    status = -1
    return
  end if
  
  status = sfginfo(sds_id, sds_name, rank, dim_sizes, data_type, n_attrs)
  if (status == FAIL) then
    print *, "ERROR: Unable to get info on SDS "//trim(sds_name)//": ", status
    status = -1
    return
  end if
  
  lut672%naot = dim_sizes(1)
  allocate(lut672%aot(lut672%naot), stat=status)
  if (status /= 0) then
    print *, "ERROR: Unable to allocate AOT670 data array: ", status
    return
  end if
  
  start1  = (/0/)
  stride1 = (/1/)
  edges1  = (/lut672%naot/)
  status = sfrdata(sds_id, start1, stride1, edges1, lut672%aot)
  if (status == FAIL) then
    print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
    return
  end if
  
  status = sfendacc(sds_id)
  if (status == FAIL) then
    print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
    return
  end if
  
  sds_name = 'SSA672_Nodes'
  sds_index = sfn2index(sd_id, sds_name)
  if (sds_index == FAIL) then
    print *, "ERROR: Unable to find index of "//trim(sds_name)//": ", sds_index
    status = -1
    return
  end if
  
  sds_id = sfselect(sd_id, sds_index)
  if (sds_id == FAIL) then
    print *, "ERROR: Unable to select SDS "//trim(sds_name)//": ", sds_id
    status = -1
    return
  end if
  
  status = sfginfo(sds_id, sds_name, rank, dim_sizes, data_type, n_attrs)
  if (status == FAIL) then
    print *, "ERROR: Unable to get info on SDS "//trim(sds_name)//": ", status
    status = -1
    return
  end if
  
  lut672%nssa = dim_sizes(1)
  allocate(lut672%ssa(lut672%nssa), stat=status)
  if (status /= 0) then
    print *, "ERROR: Unable to allocate SSA670 data array: ", status
    return
  end if
  
  start1  = (/0/)
  stride1 = (/1/)
  edges1  = (/lut672%nssa/)
  status = sfrdata(sds_id, start1, stride1, edges1, lut672%ssa)
  if (status == FAIL) then
    print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
    return
  end if
  
  status = sfendacc(sds_id)
  if (status == FAIL) then
    print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
    return
  end if
  
  
  sds_name = 'SR672_Nodes'
  sds_index = sfn2index(sd_id, sds_name)
  if (sds_index == FAIL) then
    print *, "ERROR: Unable to find index of "//trim(sds_name)//": ", sds_index
    status = -1
    return
  end if
  
  sds_id = sfselect(sd_id, sds_index)
  if (sds_id == FAIL) then
    print *, "ERROR: Unable to select SDS "//trim(sds_name)//": ", sds_id
    status = -1
    return
  end if
  
  status = sfginfo(sds_id, sds_name, rank, dim_sizes, data_type, n_attrs)
  if (status == FAIL) then
    print *, "ERROR: Unable to get info on SDS "//trim(sds_name)//": ", status
    status = -1
    return
  end if
  
  lut672%nsfc = dim_sizes(1)
  allocate(lut672%sfc(lut672%nsfc), stat=status)
  if (status /= 0) then
    print *, "ERROR: Unable to allocate SR670 data array: ", status
    return
  end if
  
  start1  = (/0/)
  stride1 = (/1/)
  edges1  = (/lut672%nsfc/)
  status = sfrdata(sds_id, start1, stride1, edges1, lut672%sfc)
  if (status == FAIL) then
    print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
    return
  end if
  
  status = sfendacc(sds_id)
  if (status == FAIL) then
    print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
    return
  end if

  sds_name = 'nvalx672'
  allocate(lut672%nvalx(lut672%nsza,lut672%nvza,lut672%nraa,lut672%naot, &
  & lut672%nssa, lut672%nsfc), stat=status)
  if (status /= 0) then
    print *, "ERROR: Unable to allocate nvalx670 data array: ", status
    return
  end if
  
  sds_index = sfn2index(sd_id, sds_name)
  if (sds_index == FAIL) then
    print *, "ERROR: Unable to find index of SDS "//trim(sds_name)//": ", sds_index
    status = -1
    return
    return
  end if
  
  sds_id = sfselect(sd_id, sds_index)
  if (sds_id == FAIL) then
    print *, "ERROR: Unable to find index of SDS "//trim(sds_name)//": ", sds_index
    status = -1
    return
  end if
  
  start6  = (/0,0,0,0,0,0/)
  stride6 = (/1,1,1,1,1,1/)
  edges6  = shape(lut672%nvalx)
  status = sfrdata(sds_id, start6, stride6, edges6, lut672%nvalx)
  if (status == FAIL) then
    print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
    return
  end if
  
  status = sfendacc(sds_id)
  if (status == FAIL) then
    print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
    return
  end if
    
  status = sfend(sd_id)
  if (status /= 0) then
    print *, "ERROR: Unable to close land aerosol model file: ", status
    return
  end if
      
  return


end function read_aerosol_lut_file


!--------------------------------------------------------
!-- NOTE: if model_frac = 0.0, the aerosol model = imod
!-- if model_frac > 0.0, the aerosol model will be interpolated between
!--  imod and imod + 1, using the value of model_frac
!
subroutine aero_470(dflag, refl, x1, x2, x3, mm, nn, ll, ma, imod,  &
&   r470, tau_x470, tau_x470_flag, trflg, model_frac, debug)
     
  implicit none
        
  logical, intent(inout)              ::  dflag
  real, intent(in)                    ::  refl
  real, intent(in)                    ::  x1
  real, intent(in)                    ::  x2
  real, intent(in)                    ::  x3
  integer, intent(in)                 ::  mm
  integer, intent(in)                 ::  nn
  integer, intent(in)                 ::  ll
  integer, intent(in)                 ::  ma
  integer, intent(in)                 ::  imod
  real, intent(in)                    ::  r470
  real, intent(inout)                 ::  tau_x470
  integer, intent(inout)              ::  tau_x470_flag
  real, intent(inout)                 ::  trflg
  real, intent(in)                    ::  model_frac
  logical, intent(in)                 ::  debug
              
  real, dimension(:), allocatable     ::  yy
  real, dimension(8)                  ::  yy2

  integer                             ::  ii
  real                                ::  frac
  
  integer                             ::  status
          
  tau_x470      = -999.0
  tau_x470_flag = -999
  
  allocate(yy(default_lut488%naot), stat=status)
  if (status /= 0) then
    print *, "ERROR: Failed to allocate array for reduced AOT table: ", status
    return
  end if
        
  status = create_reduced_lut_aot(default_lut488, refl, x1,x2,x3, imod,  &
  &   r470, model_frac, yy, debug)
  if (status /= 0) then
    dflag = .true.
    return
  end if

  if (refl <= yy(1)) then
    tau_x470 = 0.06 
    if (trflg > 0.0) tau_x470 = 0.02
    tau_x470_flag = -10
    if (debug) print *, 'aero_470, hit low bound: ', refl, yy(1)
    return
  end if
 
! Reflc off the charts!  Set AOT to max and set flag.
  if (refl >= yy(size(yy))) then
    tau_x470 = extrap(refl, yy, default_lut488%aot(1:default_lut488%naot), status)
    if (status /= 0) then
      if (status == 1) then
        tau_x470 = 5.0
      else
!         print *, 'ERROR: Failed to extrapolate AOT: ', status
        return
      end if
    end if
    if (tau_x470 > 5.0) tau_x470 = 5.0
    tau_x470_flag = 1
    if (debug) print *, 'aero_470, hit hi bound: ', refl, yy(10)
    return
  endif

!
!     Check if the reflectance increase with AOT
!
 
  if (yy(1) < yy(2)) go to 650

  if (refl < yy(4)) return

  yy2(:) = -999.0
  yy2(1:7) = yy(4:10)  ! last cell is empty, how it was originally written!

  if (yy2(2) < yy2(1)) return
  
  ii = search(refl, yy2, status, frac=frac)
  if (status /= 0) then
    dflag = .true.
    return
  end if

  tau_x470 = frac*default_lut488%aot(ii+1+3) + (1.-frac)*default_lut488%aot(ii+3)
  tau_x470_flag = 0
!  if (debug) print *, 'aero_470, exit 2358, aot: ', tau_x470
  return
 
650     continue

!
!     Pass the monotonic order check
!
  ii = search(refl, yy, status, frac=frac)
  if (status /= 0) then
    dflag = .true.
    return
  end if
  
  tau_x470 = frac*default_lut488%aot(ii+1) + (1.-frac)*default_lut488%aot(ii) 
  tau_x470_flag = 0
!  if (debug) print *, 'aero_470, exit 2371, aot: ', tau_x470
!      print *,'tau_x470 =', tau_x470
  
  deallocate(yy, stat=status)
  if (status /= 0) then
    print *, "WARNING: Failed to deallocate reduced AOT table: ", status
  end if
  
  return

end subroutine aero_470


!--------------------------------------------------------
! TODO: check/use/add status flag      
subroutine aero_650(dflag,refl,x1,x2,x3,mm,nn,ll,ma,r650,tau_x650,     &
&   tau_x650_flag,tau_x470_flag,tau_x412,tau_x470,tau_x412_flag_91,trflg)
  implicit none
        
  logical, intent(inout)              ::  dflag
  real, intent(in)                    ::  refl
  real, intent(in)                    ::  x1
  real, intent(in)                    ::  x2
  real, intent(in)                    ::  x3
  integer, intent(in)                 ::  mm
  integer, intent(in)                 ::  nn
  integer, intent(in)                 ::  ll
  integer, intent(in)                 ::  ma
  real, intent(in)                    ::  r650
  real, intent(inout)                 ::  tau_x650
  integer, intent(inout)              ::  tau_x650_flag
  integer, intent(in)                 ::  tau_x470_flag
  real, intent(in)                    ::  tau_x412
  real, intent(in)                    ::  tau_x470
  integer, intent(in)                 ::  tau_x412_flag_91
  real, intent(in)                    ::  trflg
  !logical, intent(in), optional      ::  debug
        
  real, dimension(:), allocatable     ::  yy
  real, dimension(8)                  ::  yy2
  real, dimension(4)                  ::  yy3
  real, dimension(6)                  ::  yy5
     
  real                                ::  frac
  real                                ::  w0_x
  integer                             ::  ii
        
  integer                             :: status
  logical                             :: debug
        
  tau_x650_flag = -999
  tau_x650      = -999.0    
  debug         = .false.    
  
  allocate(yy(default_lut672%naot), stat=status)
  if (status /= 0) then
    print *, "ERROR: Failed to allocate array for reduced AOT table: ", status
    return
  end if
  
  status = create_reduced_lut_aot(default_lut672, refl, x1,x2,x3, 1,  &
  &   r650, 1.0, yy, debug)
  if (status /= 0) then
    dflag = .true.
    return
  end if

  if (refl <= yy(1) .and. yy(1) < yy(2)) then
    tau_x650 = 0.002
    if (trflg > 0.0) tau_x650 = 0.02
    tau_x650_flag = -10
    return
  end if
 
! Reflc off the charts!  Set AOT to max and set flag.
  if (refl >= yy(10)) then
    tau_x650  = extrap(refl, yy, default_lut672%aot(1:default_lut672%naot), status)
    if (status /= 0) then
      if (status == 1) then
        tau_x650 = 5.0
      else
!         print *, 'ERROR: Failed to extrapolate AOT: ', status
        return
      end if
    end if
    if (tau_x650 > 5.0) tau_x650 = 5.0
    w0_x      = -999.
    tau_x650_flag = 1 
    return
  end if

!     -- check for absorbing dust over bright surface
  if (refl >= yy(7)) then
    yy3(:) = yy(7:10)
 
    if (yy3(2) < yy3(1)) return
        
    ii = search(refl, yy3, status, frac=frac)
    if (status /= 0) then
      dflag = .true.
      return
    end if
        
    tau_x650 = frac*default_lut672%aot(ii+1+6) + (1.-frac)*default_lut672%aot(ii+6)
    tau_x650_flag = 0
  
    return
  
  end if
      
!     -- check if 470 or 412 retrievals were off the charts (refl > yy(10))
!      if (tau_x470_flag > 0) go to 670
!      if (tau_x412_flag_91 > 0) go to 680

!
!     Check if the reflectance increase with AOT
!
  if (yy(1) < yy(2)) go to 650
 
  if (refl < yy(4)) return
      
  yy2(1:7) = yy(4:10)
 
  if (yy2(2) < yy2(1)) return
    
  ii = search(refl, yy2, status, frac=frac)
  if (status /= 0) then
    dflag = .true.
    return
  end if

  tau_x650 = frac*default_lut672%aot(ii+1+3) + (1.-frac)*default_lut672%aot(ii+3)
  tau_x650_flag = 0
  return

670   continue
!      if (refl < yy(8)) return
!
!      do i = 1, 3
!        yy3(i) = yy(i+7)
!      end do
!
!      if (yy3(2) < yy3(1)) return
!      call search2(dflag,refl,yy3,3,ii,frac)
!      if (dflag) return
!      tau_x650 = frac*tau(ii+1+7) + (1.-frac)*tau(ii+7)
!      tau_x650_flag = 0
!      return

680   continue
  if (refl < yy(5)) return

  yy5(:) = yy(5:10)

  if (yy5(2) < yy5(1)) return

  ii = search(refl, yy5, status, frac=frac)
  if (status /= 0) then
    dflag = .true.
    return
  end if

! if (yy3(1) > yy3(2)) print *,'yy=',refl,(yy(i),I=1,10)
  tau_x650 = frac*default_lut672%aot(ii+1+4) + (1.-frac)*default_lut672%aot(ii+4)
  tau_x650_flag = 0

  return
 
650   continue
!
!     Pass the monotonic order check
!
  ii = search(refl, yy, status, frac=frac)
  if (status /= 0) then
    dflag = .true.
    return
  end if

! print *,'after 2nd search'
  tau_x650 = frac*default_lut672%aot(ii+1) + (1.-frac)*default_lut672%aot(ii)
  tau_x650_flag = 0
      
  deallocate(yy, stat=status)
  if (status /= 0) then
    print *, "WARNING: Failed to deallocate reduced AOT table: ", status
  end if
  
  return
      
end subroutine aero_650 
        
! -- NOTE: if model_frac = 0.0, the aerosol model = imod
! -- if model_frac > 0.0, the aerosol model will be interpolated between
! --  imod and imod + 1, using the value of model_frac
!
!--------------------------------------------------------
  subroutine aero_412(dflag,refl,x1,x2,x3,mm,nn,ll,ma,imod,r412,        &
  &   tau_x412,tau_x412_flag,trflg,model_frac,debug)
        
  implicit none
        

  logical, intent(inout)            ::  dflag
  real, intent(in)                  ::  refl
  real, intent(in)                  ::  x1
  real, intent(in)                  ::  x2
  real, intent(in)                  ::  x3
  integer, intent(in)               ::  mm
  integer, intent(in)               ::  nn
  integer, intent(in)               ::  ll
  integer, intent(in)               ::  ma
  integer, intent(in)               ::  imod
  real, intent(in)                  ::  r412
  real, intent(inout)               ::  tau_x412
  integer, intent(inout)            ::  tau_x412_flag
  real, intent(in)                  ::  trflg
  real, intent(in)                  ::  model_frac
  logical, intent(in)               ::  debug

  real, dimension(:), allocatable   ::  yy
  real, dimension(8)                ::  yy2
              
  integer                           ::  ii
  real                              ::  frac
        
  integer                           :: status
        
  tau_x412      = -999.0
  tau_x412_flag = -999
             
  allocate(yy(default_lut412%naot), stat=status)
  if (status /= 0) then
    print *, "ERROR: Failed to allocate array for reduced AOT table: ", status
    return
  end if
  
  status = create_reduced_lut_aot(default_lut412, refl, x1,x2,x3, imod,  &
  &   r412, model_frac, yy, debug)
  if (status /= 0) then
    dflag = .true.
    return
  end if
  
  if (refl <= yy(1)) then
    tau_x412 = 0.06 
    if (trflg > 0.0) tau_x412 = 0.02
    tau_x412_flag = -10       
    if (debug) print *, 'aero_412, hit low bound: ', refl, yy(1)
    return
  end if

! Reflc off the charts!  Set AOT to max and set flag.
  if (refl >= yy(10)) then
    tau_x412 = extrap(refl, yy, default_lut412%aot(1:default_lut412%naot), status)
    if (status /= 0) then
      if (status == 1) then
        tau_x412 = 5.0
      else
!         print *, 'ERROR: Failed to extrapolate AOT: ', status
        return
      end if
    end if
    if (tau_x412 > 5.0) tau_x412 = 5.0
    
    tau_x412_flag = 1
    if (debug) print *, 'aero_412, hit hi bound: ', refl, yy(10)
    return
  end if
!
!     Check if the reflectance increase with AOT
!
  if (yy(1) < yy(2)) go to 650

  if (refl < yy(4)) return
  
  yy2(:) = -999.0
  yy2(1:7) = yy(4:10)
 
  if (yy2(2) < yy2(1)) return
  
  ii = search(refl, yy2, status, frac=frac)
  if (status /= 0) then
    dflag = .true.
    return
  end if

  tau_x412 = frac*default_lut412%aot(ii+1+3) + (1.0-frac)*default_lut412%aot(ii+3)
  tau_x412_flag = 0
!  if (debug) print *, 'aero_412, exit 2355, aot: ', tau_x412
  return
 
650 continue
!
!     Pass the monotonic order check
!       
  ii = search(refl, yy, status, frac=frac)
  if (status /= 0) then
    dflag = .true.
    return
  end if

  tau_x412 = frac*default_lut412%aot(ii+1) + (1.-frac)*default_lut412%aot(ii) 
  tau_x412_flag = 0
!  if (debug) print *, 'aero_412, exit 2367, aot: ', tau_x412, refl
  
  deallocate(yy, stat=status)
  if (status /= 0) then
    print *, "WARNING: Failed to deallocate reduced AOT table: ", status
  end if
  
  return
  
end subroutine aero_412
            
!--------------------------------------------------------
subroutine aero_412_abs(dflag,refl,x1,x2,x3,mm,nn,ll,r412,tau_x,w0_x)      
  implicit none
  
  logical, intent(inout)        ::  dflag
  real, intent(in)              ::  refl
  real, intent(in)              ::  x1
  real, intent(in)              ::  x2
  real, intent(in)              ::  x3
  integer, intent(in)           ::  mm
  integer, intent(in)           ::  nn
  integer, intent(in)           ::  ll
  real, intent(in)              ::  r412
  real, intent(in)              ::  tau_x
  real, intent(inout)           ::  w0_x
  
  integer                       ::  index_ii, index_ia
  real                          ::  frac, frac_ia
  integer                       ::  status
  
  real, dimension(:,:,:,:), allocatable   ::  nnvalxw
  real, dimension(:), allocatable         ::  yyw
               
  dflag = .false.
        
  allocate(nnvalxw(4,4,2,default_lut412%nssa), yyw(default_lut412%nssa), stat=status)
  if (status /= 0) then 
    print *, "ERROR: Failed to allocate arrays: ", status
    dflag = .false.
    return
  end if
  
  index_ia = search(tau_x, default_lut412%aot, status, frac=frac_ia)
  if (status /= 0) then
    dflag = .false.
    return
  end if
  
  status = create_reduced_lut_ssa(default_lut412, refl, x1, x2, x3, index_ia, &
  &   r412, frac_ia, yyw)
  if (status /= 0) then
    dflag = .true.
    print *, "ERROR: create_reduced_lut_ssa: aero_412_abs", refl,x1, x2, x3, index_ia, r412  
    return
  end if
  
  if (refl.le.yyw(1)) then
    w0_x = 0.82
    return
  endif

  if (refl.ge.yyw(8)) then
    w0_x = 1.0
    return
  endif
      
  index_ii = search(refl, yyw, status, frac=frac)
  if (status /= 0) then
    dflag = .false.
    return
  end if      
  w0_x = frac*default_lut412%ssa(index_ii+1) + (1.-frac)*default_lut412%ssa(index_ii)
      
  deallocate(nnvalxw, yyw, stat=status)
  if (status /= 0) then 
    print *, "WARNING: Failed to deallocate arrays: ", status
    return
  end if
        
  return
  
end subroutine aero_412_abs
    
!--------------------------------------------------------
subroutine aero_470_abs(dflag2,refl,x1,x2,x3,mm,nn,ll,r470,tau_x,w0_x)
  implicit none
  
  logical, intent(inout)        ::  dflag2
  real, intent(in)              ::  refl
  real, intent(in)              ::  x1
  real, intent(in)              ::  x2
  real, intent(in)              ::  x3
  integer, intent(in)           ::  mm
  integer, intent(in)           ::  nn
  integer, intent(in)           ::  ll
  real, intent(in)              ::  r470
  real, intent(in)              ::  tau_x
  real, intent(inout)           ::  w0_x
  
  integer                       ::  index_ii, index_ia
  real                          ::  frac, frac_ia
  integer                       ::  status
  
  real, dimension(:,:,:,:), allocatable   ::  nnvalxw
  real, dimension(:), allocatable         ::  yyw        
                  
  dflag2 = .false.
      
  allocate(nnvalxw(4,4,2,default_lut488%nssa), yyw(default_lut488%nssa), stat=status)
  if (status /= 0) then 
    print *, "ERROR: Failed to allocate arrays: ", status
    dflag2 = .false.
    return
  end if
  
  index_ia = search(tau_x, default_lut488%aot, status, frac=frac_ia)
  if (status /= 0) then
    dflag2 = .false.
    return
  end if
      
  status = create_reduced_lut_ssa(default_lut488, refl, x1, x2, x3, index_ia, &
  &     r470, frac_ia, yyw)
  if (status /= 0) then
    dflag2 = .true.
    print *, "ERROR: create_reduced_lut_ssa: aero_470_abs", refl,x1, x2, x3, index_ia, r470      
    return
  end if
       
  if (refl.le.yyw(1)) then
    w0_x = -999.
    return
  endif

  if (refl.ge.yyw(4)) then
    w0_x = 1.0
    return
  endif
    
  index_ii = search(refl, yyw, status, frac=frac)
  if (status /= 0) then
    dflag2 = .false.
    return
  end if      
  w0_x = frac*default_lut488%ssa(index_ii+1) + (1.-frac)*default_lut488%ssa(index_ii)

  deallocate(nnvalxw, yyw, stat=status)
  if (status /= 0) then 
    print *, "WARNING: Failed to deallocate arrays: ", status
    return
  end if
 
  return

end subroutine aero_470_abs

!--------------------------------------------------------
!-- NOTE: if model_frac = 0.0, the aerosol model = imod
!-- if model_frac > 0.0, the aerosol model will be interpolated between
!--  imod and imod + 1, using the value of model_frac
!
subroutine aero_470_dust(dflag, refl, x1, x2, x3, mm, nn, ll, ma, imod,  &
&   r470, tau_x470, tau_x470_flag, trflg, model_frac, debug)

  implicit none

  logical, intent(inout)              ::  dflag
  real, intent(in)                    ::  refl
  real, intent(in)                    ::  x1
  real, intent(in)                    ::  x2
  real, intent(in)                    ::  x3
  integer, intent(in)                 ::  mm
  integer, intent(in)                 ::  nn
  integer, intent(in)                 ::  ll
  integer, intent(in)                 ::  ma
  integer, intent(in)                 ::  imod
  real, intent(in)                    ::  r470
  real, intent(inout)                 ::  tau_x470
  integer, intent(inout)              ::  tau_x470_flag
  real, intent(inout)                 ::  trflg
  real, intent(in)                    ::  model_frac
  logical, intent(in)                 ::  debug

  real, dimension(:), allocatable     ::  yy
  real, dimension(8)                  ::  yy2

  integer                             ::  ii
  real                                ::  frac

  integer                             ::  status

  tau_x470      = -999.0
  tau_x470_flag = -999

  allocate(yy(dust_lut488%naot), stat=status)
  if (status /= 0) then
    print *, "ERROR: Failed to allocate array for reduced AOT table: ", status
    return
  end if
  
  status = create_reduced_lut_aot(dust_lut488, refl, x1,x2,x3, imod,  &
  &   r470, model_frac, yy, debug)
  if (status /= 0) then
    dflag = .true.
    return
  end if

  if (refl <= yy(1)) then
    tau_x470 = 0.06
    if (trflg > 0.0) tau_x470 = 0.02
    tau_x470_flag = -10
    if (debug) print *, 'aero_470_dust, hit low bound: ', refl, yy(1)
    return
  end if

! Reflc off the charts!  Set AOT to max and set flag.
  if (refl >= yy(size(yy))) then
    tau_x470 = extrap(refl, yy, dust_lut488%aot(1:dust_lut488%naot), status)
    
    if (status /= 0) then
      if (status == 1) then
        tau_x470 = 5.0
      else
!         print *, 'ERROR: Failed to extrapolate AOT: ', status
        return
      end if
    end if
    if (tau_x470 > 5.0) tau_x470 = 5.0
    tau_x470_flag = 1
    if (debug) print *, 'aero_470_dust, hit hi bound: ', refl, yy(10)
    return
  endif

!
!     Check if the reflectance increase with AOT
!

  if (yy(1) < yy(2)) go to 650

  if (refl < yy(4)) return

  yy2(:) = -999.0
  yy2(1:7) = yy(4:10)  ! last cell is empty, how it was originally written!

  if (yy2(2) < yy2(1)) return

  ii = search(refl, yy2, status, frac=frac)
  if (status /= 0) then
    dflag = .true.
    return
  end if

  tau_x470 = frac*dust_lut488%aot(ii+1+3) + (1.-frac)*dust_lut488%aot(ii+3)
  tau_x470_flag = 0
  if (debug) print *, 'aero_470_dust, exit 2358, aot: ', tau_x470
  return

650     continue

!
!     Pass the monotonic order check
!
  ii = search(refl, yy, status, frac=frac)
  if (status /= 0) then
    dflag = .true.
    return
  end if

  tau_x470 = frac*dust_lut488%aot(ii+1) + (1.-frac)*dust_lut488%aot(ii)
  tau_x470_flag = 0
  if (debug) print *, 'aero_470_dust, exit 2371, aot: ', tau_x470
!      print *,'tau_x470 =', tau_x470

  deallocate(yy, stat=status)
  if (status /= 0) then
    print *, "WARNING: Failed to deallocate reduced AOT table: ", status
  end if

  return

end subroutine aero_470_dust

!--------------------------------------------------------
! TODO: check/use/add status flag
subroutine aero_650_dust(dflag,refl,x1,x2,x3,mm,nn,ll,ma,r650,tau_x650,     &
&   tau_x650_flag,tau_x470_flag,tau_x412,tau_x470,tau_x412_flag_91,trflg)
  implicit none

  logical, intent(inout)              ::  dflag
  real, intent(in)                    ::  refl
  real, intent(in)                    ::  x1
  real, intent(in)                    ::  x2
  real, intent(in)                    ::  x3
  integer, intent(in)                 ::  mm
  integer, intent(in)                 ::  nn
  integer, intent(in)                 ::  ll
  integer, intent(in)                 ::  ma
  real, intent(in)                    ::  r650
  real, intent(inout)                 ::  tau_x650
  integer, intent(inout)              ::  tau_x650_flag
  integer, intent(in)                 ::  tau_x470_flag
  real, intent(in)                    ::  tau_x412
  real, intent(in)                    ::  tau_x470
  integer, intent(in)                 ::  tau_x412_flag_91
  real, intent(in)                    ::  trflg
  !logical, intent(in), optional      ::  debug

  real, dimension(:), allocatable     ::  yy
  real, dimension(8)                  ::  yy2
  real, dimension(4)                  ::  yy3
  real, dimension(6)                  ::  yy5

  real                                ::  frac
  real                                ::  w0_x
  integer                             ::  ii

  integer                             :: status
  logical                             :: debug

  tau_x650_flag = -999
  tau_x650      = -999.0
  debug         = .false.

  allocate(yy(dust_lut672%naot), stat=status)
  if (status /= 0) then
    print *, "ERROR: Failed to allocate array for reduced AOT table: ", status
    return
  end if

  status = create_reduced_lut_aot(dust_lut672, refl, x1,x2,x3, 1,  &
  &   r650, 1.0, yy, debug)
  if (status /= 0) then
    dflag = .true.
    return
  end if

  if (refl <= yy(1) .and. yy(1) < yy(2)) then
    tau_x650 = 0.002
    if (trflg > 0.0) tau_x650 = 0.02
    tau_x650_flag = -10
    return
  end if

! Reflc off the charts!  Set AOT to max and set flag.
  if (refl >= yy(10)) then
    tau_x650  = extrap(refl, yy, dust_lut672%aot(1:dust_lut672%naot), status)
    if (status /= 0) then
      if (status == 1) then
        tau_x650 = 5.0
      else
!         print *, 'ERROR: Failed to extrapolate AOT: ', status
        return
      end if
    end if
    if (tau_x650 > 5.0) tau_x650 = 5.0
    w0_x      = -999.
    tau_x650_flag = 1
    return
  end if

!     -- check for absorbing dust over bright surface
  if (refl >= yy(7)) then
    yy3(:) = yy(7:10)

    if (yy3(2) < yy3(1)) return

    ii = search(refl, yy3, status, frac=frac)
    if (status /= 0) then
      dflag = .true.
      return
    end if

    tau_x650 = frac*dust_lut672%aot(ii+1+6) + (1.-frac)*dust_lut672%aot(ii+6)
    tau_x650_flag = 0

    return

  end if

!     -- check if 470 or 412 retrievals were off the charts (refl > yy(10))
!      if (tau_x470_flag > 0) go to 670
!      if (tau_x412_flag_91 > 0) go to 680

!
!     Check if the reflectance increase with AOT
!
  if (yy(1) < yy(2)) go to 650

  if (refl < yy(4)) return

  yy2(1:7) = yy(4:10)

  if (yy2(2) < yy2(1)) return

  ii = search(refl, yy2, status, frac=frac)
  if (status /= 0) then
    dflag = .true.
    return
  end if

  tau_x650 = frac*dust_lut672%aot(ii+1+3) + (1.-frac)*dust_lut672%aot(ii+3)
  tau_x650_flag = 0
  return

670   continue
!      if (refl < yy(8)) return
!
!      do i = 1, 3
!        yy3(i) = yy(i+7)
!      end do
!
!      if (yy3(2) < yy3(1)) return
!      call search2(dflag,refl,yy3,3,ii,frac)
!      if (dflag) return
!      tau_x650 = frac*tau(ii+1+7) + (1.-frac)*tau(ii+7)
!      tau_x650_flag = 0
!      return

680   continue
  if (refl < yy(5)) return

  yy5(:) = yy(5:10)

  if (yy5(2) < yy5(1)) return

  ii = search(refl, yy5, status, frac=frac)
  if (status /= 0) then
    dflag = .true.
    return
  end if

! if (yy3(1) > yy3(2)) print *,'yy=',refl,(yy(i),I=1,10)
  tau_x650 = frac*dust_lut672%aot(ii+1+4) + (1.-frac)*dust_lut672%aot(ii+4)
  tau_x650_flag = 0

  return

650   continue
!
!     Pass the monotonic order check
!
  ii = search(refl, yy, status, frac=frac)
  if (status /= 0) then
    dflag = .true.
    return
  end if

! print *,'after 2nd search'
  tau_x650 = frac*dust_lut672%aot(ii+1) + (1.-frac)*dust_lut672%aot(ii)
  tau_x650_flag = 0

  deallocate(yy, stat=status)
  if (status /= 0) then
    print *, "WARNING: Failed to deallocate reduced AOT table: ", status
  end if

  return

end subroutine aero_650_dust

! -- NOTE: if model_frac = 0.0, the aerosol model = imod
! -- if model_frac > 0.0, the aerosol model will be interpolated between
! --  imod and imod + 1, using the value of model_frac
!
!--------------------------------------------------------
  subroutine aero_412_dust(dflag,refl,x1,x2,x3,mm,nn,ll,ma,imod,r412,        &
  &   tau_x412,tau_x412_flag,trflg,model_frac,debug)

  implicit none


  logical, intent(inout)            ::  dflag
  real, intent(in)                  ::  refl
  real, intent(in)                  ::  x1
  real, intent(in)                  ::  x2
  real, intent(in)                  ::  x3
  integer, intent(in)               ::  mm
  integer, intent(in)               ::  nn
  integer, intent(in)               ::  ll
  integer, intent(in)               ::  ma
  integer, intent(in)               ::  imod
  real, intent(in)                  ::  r412
  real, intent(inout)               ::  tau_x412
  integer, intent(inout)            ::  tau_x412_flag
  real, intent(in)                  ::  trflg
  real, intent(in)                  ::  model_frac
  logical, intent(in)               ::  debug

  real, dimension(:), allocatable   ::  yy
  real, dimension(8)                ::  yy2

  integer                           ::  ii
  real                              ::  frac

  integer                           :: status

  tau_x412      = -999.0
  tau_x412_flag = -999

  allocate(yy(dust_lut412%naot), stat=status)
  if (status /= 0) then
    print *, "ERROR: Failed to allocate array for reduced AOT table: ", status
    return
  end if

  status = create_reduced_lut_aot(dust_lut412, refl, x1,x2,x3, imod,  &
  &   r412, model_frac, yy, debug)
  if (status /= 0) then
    dflag = .true.
    return
  end if

  if (refl <= yy(1)) then
    tau_x412 = 0.06
    if (trflg > 0.0) tau_x412 = 0.02
    tau_x412_flag = -10
    if (debug) print *, 'aero_412_dust, hit low bound: ', refl, yy(1)
    return
  end if

! Reflc off the charts!  Set AOT to max and set flag.
  if (refl >= yy(10)) then
    tau_x412 = extrap(refl, yy, dust_lut412%aot(1:dust_lut412%naot), status)
    if (status /= 0) then
      if (status == 1) then
        tau_x412 = 5.0
      else
!         print *, 'ERROR: Failed to extrapolate AOT: ', status
        return
      end if
    end if
    if (tau_x412 > 5.0) tau_x412 = 5.0

    tau_x412_flag = 1
    if (debug) print *, 'aero_412_dust, hit hi bound: ', refl, yy(10)
    return
  end if
!
!     Check if the reflectance increase with AOT
!
  if (yy(1) < yy(2)) go to 650

  if (refl < yy(4)) return

  yy2(:) = -999.0
  yy2(1:7) = yy(4:10)

  if (yy2(2) < yy2(1)) return

  ii = search(refl, yy2, status, frac=frac)
  if (status /= 0) then
    dflag = .true.
    return
  end if

  tau_x412 = frac*dust_lut412%aot(ii+1+3) + (1.0-frac)*dust_lut412%aot(ii+3)
  tau_x412_flag = 0
  if (debug) print *, 'aero_412_dust, exit 2355, aot: ', tau_x412
  return

650 continue
!
!     Pass the monotonic order check
!
  ii = search(refl, yy, status, frac=frac)
  if (status /= 0) then
    dflag = .true.
    return
  end if

  tau_x412 = frac*dust_lut412%aot(ii+1) + (1.-frac)*dust_lut412%aot(ii)
  tau_x412_flag = 0
  if (debug) print *, 'aero_412_dust, exit 2367, aot: ', tau_x412, refl

  deallocate(yy, stat=status)
  if (status /= 0) then
    print *, "WARNING: Failed to deallocate reduced AOT table: ", status
  end if

  return

end subroutine aero_412_dust

!--------------------------------------------------------
subroutine aero_412_abs_dust(dflag,refl,x1,x2,x3,mm,nn,ll,r412,tau_x,w0_x)
  implicit none

  logical, intent(inout)        ::  dflag
  real, intent(in)              ::  refl
  real, intent(in)              ::  x1
  real, intent(in)              ::  x2
  real, intent(in)              ::  x3
  integer, intent(in)           ::  mm
  integer, intent(in)           ::  nn
  integer, intent(in)           ::  ll
  real, intent(in)              ::  r412
  real, intent(in)              ::  tau_x
  real, intent(inout)           ::  w0_x

  integer                       ::  index_ii, index_ia
  real                          ::  frac, frac_ia
  integer                       ::  status

  real, dimension(:,:,:,:), allocatable   ::  nnvalxw
  real, dimension(:), allocatable         ::  yyw

  dflag = .false.

  allocate(nnvalxw(4,4,2,dust_lut412%nssa), yyw(dust_lut412%nssa), stat=status)
  if (status /= 0) then
    print *, "ERROR: Failed to allocate arrays: ", status
    dflag = .false.
    return
  end if

  index_ia = search(tau_x, dust_lut412%aot, status, frac=frac_ia)
  if (status /= 0) then
    dflag = .false.
    return
  end if

  status = create_reduced_lut_ssa(dust_lut412, refl, x1, x2, x3, index_ia, &
  &   r412, frac_ia, yyw)
  if (status /= 0) then
    dflag = .true.
    print *, "ERROR: create_reduced_lut_ssa: aero_412_abs_dust", refl,x1, x2, x3, index_ia, r412
    return
  end if

  if (refl.le.yyw(1)) then
    w0_x = 0.82
    return
  endif

  if (refl.ge.yyw(8)) then
    w0_x = 1.0
    return
  endif

  index_ii = search(refl, yyw, status, frac=frac)
  if (status /= 0) then
    dflag = .false.
    return
  end if
  w0_x = frac*dust_lut412%ssa(index_ii+1) + (1.-frac)*dust_lut412%ssa(index_ii)

  deallocate(nnvalxw, yyw, stat=status)
  if (status /= 0) then
    print *, "WARNING: Failed to deallocate arrays: ", status
    return
  end if

  return

end subroutine aero_412_abs_dust

!--------------------------------------------------------
subroutine aero_470_abs_dust(dflag2,refl,x1,x2,x3,mm,nn,ll,r470,tau_x,w0_x)
  implicit none

  logical, intent(inout)        ::  dflag2
  real, intent(in)              ::  refl
  real, intent(in)              ::  x1
  real, intent(in)              ::  x2
  real, intent(in)              ::  x3
  integer, intent(in)           ::  mm
  integer, intent(in)           ::  nn
  integer, intent(in)           ::  ll
  real, intent(in)              ::  r470
  real, intent(in)              ::  tau_x
  real, intent(inout)           ::  w0_x

  integer                       ::  index_ii, index_ia
  real                          ::  frac, frac_ia
  integer                       ::  status

  real, dimension(:,:,:,:), allocatable   ::  nnvalxw
  real, dimension(:), allocatable         ::  yyw

  dflag2 = .false.

  allocate(nnvalxw(4,4,2,dust_lut488%nssa), yyw(dust_lut488%nssa), stat=status)
  if (status /= 0) then
    print *, "ERROR: Failed to allocate arrays: ", status
    dflag2 = .false.
    return
  end if

  index_ia = search(tau_x, dust_lut488%aot, status, frac=frac_ia)
  if (status /= 0) then
    dflag2 = .false.
    return
  end if

  status = create_reduced_lut_ssa(dust_lut488, refl, x1, x2, x3, index_ia, &
  &     r470, frac_ia, yyw)
  if (status /= 0) then
    dflag2 = .true.
    print *, "ERROR: create_reduced_lut_ssa: aero_470_abs_dust", refl,x1, x2, x3, index_ia, r470
    return
  end if

  if (refl.le.yyw(1)) then
    w0_x = -999.
    return
  endif

  if (refl.ge.yyw(4)) then
    w0_x = 1.0
    return
  endif

  index_ii = search(refl, yyw, status, frac=frac)
  if (status /= 0) then
    dflag2 = .false.
    return
  end if
  w0_x = frac*dust_lut488%ssa(index_ii+1) + (1.-frac)*dust_lut488%ssa(index_ii)

  deallocate(nnvalxw, yyw, stat=status)
  if (status /= 0) then
    print *, "WARNING: Failed to deallocate arrays: ", status
    return
  end if

  return

end subroutine aero_470_abs_dust

!--------------------------------------------------------
integer function create_reduced_lut_aot(lut, refl, sza, vza, raa, imod,  &
&   rXXX, model_frac, yy, debug) result(status)
      
  implicit none
  
  type(viirs_aerosol_lut)           ::  lut
  real, intent(in)                  ::  refl
  real, intent(in)                  ::  sza
  real, intent(in)                  ::  vza
  real, intent(in)                  ::  raa
  integer, intent(in)               ::  imod
  real, intent(in)                  ::  rXXX
  real, intent(in)                  ::  model_frac
  real, intent(inout),dimension(:)  ::  yy
  logical, intent(in),optional      ::  debug

  real, dimension(:,:,:,:), allocatable  ::  nnvalx, nnvalx1, nnvalx2
  
  integer                     ::  index_ii, ii, jj, mbeg, nbeg
  integer                     ::  ia, j, i
  real                        ::  frac, xfrac, y, dy

  real, parameter             ::  pi = 3.14159
              
!  if (debug) print *, 'aero_XXX, in: ', refl, sza, vza, raa, rXXX, imod, model_frac
        
  status = -1
        
  if (rXXX < 0.0) return
      
  index_ii = search(rXXX, lut%sfc, status, frac=frac) 
  if (status /= 0) then 
!     print *, 'ERROR: Specified SFC not within table: '
    return
  endif
!  if (debug) print *, 'aero_XXX, sfc indx: ', rXXX, index_ii
  
  ii = search(raa, lut%raa, status, frac=xfrac)
  if (status /= 0) then
    print *, 'ERROR: Specified RAA not within table: ', raa, lut%raa(1), lut%raa(lut%nraa)
    return
  end if
!  if (debug) print *, 'aero_XXX, raa: ', raa, ii, lut%raa(ii), lut%raa(ii+1)
  
  mbeg = search(sza, lut%sza, status) 
  if (status /= 0) then 
    print *, 'ERROR: Specified SZA not within table: ', sza, lut%sza(1), lut%sza(lut%nsza)
    return
  end if
  mbeg = max(0, mbeg-2)
  if (mbeg > lut%nsza-4) mbeg = lut%nsza - 4
  
  nbeg = search(vza, lut%vza, status)    
  if (status /= 0) then 
    print *, 'ERROR: Specified VZA not within table: ', vza, lut%vza(1), lut%vza(lut%nvza)
    return
  end if  
  nbeg = max(0, nbeg-2)
  if (nbeg > lut%nvza-4) nbeg = lut%nvza - 4
!  if (debug) print *, 'mbeg, nbeg: ', mbeg, nbeg

! -- interpolate between models if requested.
! if (present(model_frac)) then
  allocate(nnvalx(4,4,2,lut%naot), nnvalx1(4,4,2,lut%naot), &
  &   nnvalx2(4,4,2,lut%naot), stat=status)
  if (status /= 0) then
    print *, "ERROR: Failed to allocate reduced LUT arrays: ", status
    return
  end if
        
  if (imod < lut%nssa) then              ! temp check to avoid out-of-bounds issue with imod=nssa
    nnvalx1(:,:,:,:) = -999.0
    nnvalx2(:,:,:,:) = -999.0
    do ia = 1, lut%naot
      do i = 1, 4
        do j = 1, 4
          nnvalx1(i,j,1,ia) = lut%nvalx(mbeg+i,nbeg+j,ii,ia,imod,index_ii)*     &
  &            (1.0-frac) + lut%nvalx(mbeg+i,nbeg+j,ii,ia,imod,index_ii+1)*frac
          nnvalx1(i,j,2,ia) = lut%nvalx(mbeg+i,nbeg+j,ii+1,ia,imod,index_ii)*   &
  &             (1.0-frac) + lut%nvalx(mbeg+i,nbeg+j,ii+1,ia,imod,index_ii+1)*frac

          nnvalx2(i,j,1,ia) = lut%nvalx(mbeg+i,nbeg+j,ii,ia,imod+1,index_ii)*   &
  &             (1.0-frac) + lut%nvalx(mbeg+i,nbeg+j,ii,ia,imod+1,index_ii+1)*frac
          nnvalx2(i,j,2,ia) = lut%nvalx(mbeg+i,nbeg+j,ii+1,ia,imod+1,index_ii)* &
  &             (1.0-frac) + lut%nvalx(mbeg+i,nbeg+j,ii+1,ia,imod+1,index_ii+1)*frac

          nnvalx(i,j,1,ia) = (1.0-model_frac) * nnvalx1(i,j,1,ia) +   &  
  &                                 model_frac * nnvalx2(i,j,1,ia)
          nnvalx(i,j,2,ia) = (1.0-model_frac) * nnvalx1(i,j,2,ia) +   &
  &                              model_frac * nnvalx2(i,j,2,ia)
        end do
      end do
    end do
          
! -- just use imod, no interpolation
  else      
    nnvalx(:,:,:,:) = -999.0
    do ia = 1, lut%naot
      do i = 1, 4
        do j = 1, 4
          nnvalx(i,j,1,ia) = lut%nvalx(mbeg+i,nbeg+j,ii,ia,imod,index_ii)*      &     
  &              (1.-frac) + lut%nvalx(mbeg+i,nbeg+j,ii,ia,imod,index_ii+1)*frac   
          nnvalx(i,j,2,ia) = lut%nvalx(mbeg+i,nbeg+j,ii+1,ia,imod,index_ii)*    &
  &              (1.-frac) + lut%nvalx(mbeg+i,nbeg+j,ii+1,ia,imod,index_ii+1)*frac   
        end do
      end do
    end do
  end if
        
!---     interpolating AOT tables
  yy(:) = -999.0
  do ia = 1, lut%naot
    call new_intep(lut%sza, lut%vza, lut%raa, nnvalx, lut%nsza, lut%nvza, lut%nraa, ia,   &    
    &   sza,vza,raa,y,dy,mbeg,nbeg,xfrac)

    yy(ia) = y/pi
    !print *,'tau, i/f=', tau(ia), y/pi, dy
  end do
        
  deallocate(nnvalx, nnvalx1, nnvalx2, stat=status)
  if (status /= 0) then
    print *, "WARNING: Failed to deallocate reduced LUT arrays: ", status
  end if

  return
        
end function create_reduced_lut_aot

!--------------------------------------------------------
integer function create_reduced_lut_ssa(lut, refl, sza, vza, raa, index_ia,  &
&   rXXX, aot_frac, yy, debug) result(status)
  implicit none

  type(viirs_aerosol_lut)     ::  lut
  real, intent(in)            ::  refl
  real, intent(in)            ::  sza
  real, intent(in)            ::  vza
  real, intent(in)            ::  raa
  integer, intent(in)         ::  index_ia
  real, intent(in)            ::  rXXX
  real, intent(in)            ::  aot_frac
  real, intent(inout),dimension(:)  ::  yy
  logical, intent(in),optional::  debug

  real, dimension(:,:,:,:), allocatable ::  nnvalx

  real, parameter             ::  pi = 3.14159

  integer                     ::  index_ii, ii, mbeg, nbeg
  integer                     ::  iw, j, i
  real                        ::  frac, xfrac, y, dy
  real                        ::  dd1, dd2

! -- calculate the interpolation indices for surf. reflectance, RAA, SZA, VZA.
  index_ii = search(rXXX, lut%sfc, status, frac=frac)
  if (status /= 0) then
    print *, "ERROR: Failed to get interpolation index for surface: ", rXXX, status
    return
  end if

  ii = search(raa, lut%raa, status, frac=xfrac)
  if (status /= 0) then
    print *, "ERROR: Failed to get interpolation index for RAA: ", raa, status
    return
  end if

  mbeg = search(sza, lut%sza, status)
  if (status /= 0) then
    print *, "ERROR: Failed to get interpolation index for SZA: ", sza, status
    return
  end if
  mbeg = max(0, mbeg-2)
  if (mbeg > lut%nsza-4) mbeg = lut%nsza - 4

  nbeg = search(vza, lut%vza, status)
  if (status /= 0) then
    print *, "ERROR: Failed to get interpolation index for VZA: ", vza, status
    return
  end if
  nbeg = max(0, nbeg-2)
  if (nbeg > lut%nvza-4) nbeg = lut%nvza - 4

! -- perform interpolation of LUT as function of SSA.
  allocate(nnvalx(4,4,2,lut%nssa), stat=status)
  if (status /= 0) then
    print *, "ERROR: Failed to allocate nnvalx or yy array: ", status
    return
  end if

  do iw = 1, lut%nssa
    do i = 1, 4
      do j = 1, 4
      dd1 = lut%nvalx(mbeg+i,nbeg+j,ii,index_ia,iw,index_ii)*  &
  &      (1.-frac) +                                             &
  &  lut%nvalx(mbeg+i,nbeg+j,ii,index_ia,iw,index_ii+1)*frac
    dd2= lut%nvalx(mbeg+i,nbeg+j,ii,index_ia+1,iw,index_ii)* &
  &  (1.-frac)+                                              &
  &   lut%nvalx(mbeg+i,nbeg+j,ii,index_ia+1,iw,index_ii+1)    &
  &   *frac

     nnvalx(i,j,1,iw) = dd1* (1.-aot_frac) + dd2*aot_frac

    dd1 = lut%nvalx(mbeg+i,nbeg+j,ii+1,index_ia,iw,index_ii)*&
  &  (1.-frac) +                                             &
  &   lut%nvalx(mbeg+i,nbeg+j,ii+1,index_ia,iw,index_ii+1)*frac
    dd2= lut%nvalx(mbeg+i,nbeg+j,ii+1,index_ia+1,iw,index_ii)*&
  &  (1.-frac)+                                              &
  &   lut%nvalx(mbeg+i,nbeg+j,ii+1,index_ia+1,iw,index_ii+1)  &
  &   *frac

    nnvalx(i,j,2,iw) = dd1* (1.-aot_frac) + dd2*aot_frac

      end do
    end do
  end do

!-- interpolating W0 tables
  do iw = 1, lut%nssa
    call new_intep(lut%sza, lut%vza, lut%raa, nnvalx,  &
    &   lut%nsza, lut%nvza, lut%nraa, iw, sza,vza,raa,y,   &
    &   dy,mbeg,nbeg,xfrac)

    yy(iw) = y/pi
  end do

end function create_reduced_lut_ssa

!-------------------------------------------------------- 
  integer function search(xbar,x,status,frac) result(i)
!     
!     purpose
!       locate position in table of point at which interpolation is
!       required
!
!     usage
!       i = search (xbar, x, status, frac=frac)
!
!     description of parameters
!       xbar   - point at which interpolation is required
!       x      - array containing independent variable
!       status - integer indicates success (0) or failure (/=0)
!       frac   - optional parameter for fraction.
!
    implicit none
      
    real, intent(in)                  ::  xbar
    real, dimension(:), intent(in)    ::  x
    integer, intent(inout)            ::  status
    real, intent(inout), optional     ::  frac
      
    integer                           ::  n
      
    integer                           ::  m, k
    integer                           ::  icnt
    real, parameter                   ::  b = 0.69314718  ! (ln(2))

    status = -1
    
    n = size(x)
    icnt = 0
    if (n < 2) then
      print *, "Search n is less than 2."
      return
    end if
      
    if (x(1) > x(2)) then
      print *, "Search table is not in increasing order."
      status = -1
      return
    end if
    
    m = int((log(float(n)))/b)
    i = 2**m
    if (i >= n) i = n-1
    k = i
    do 
      k=k/2
      
      if (k == 0) icnt = icnt + 1
      if (icnt >= 2) exit     ! starting an infinite loop, bail out!

!     -- check if we found our indices and exit the loop if so.
!     -- otherwise, reset i and loop again.
      if (xbar >= x(i) .AND. xbar < x(i+1)) then
        status = 0
        exit
      else
        if (xbar > x(i)) then
          i = i + k
          if (i >= n) then
            i = n - 1
          end if
        else
          i = i - k
        end if
      end if

    end do

!   -- detect failed search above via icnt. if >= 2, do sequential search,
!   -- else above search was successful.    
    if (icnt >= 2 .OR. status /= 0) then
      do i = 1, n-1
        if (xbar >= x(i) .AND. xbar <= x(i+1)) then
          status = 0
          exit
         end if
      end do
    end if
    
!   -- calculate fraction if desired and search was successful.
    if (status == 0 .AND. present(frac)) then
      if (i < n) then
        frac = (xbar-x(i))/ (x(i+1)-x(i))
      else
        frac = 1.0
      end if
    end if
    
    return

  end function search
       
! -- perform a binary search of array xx for j such that x lies between xx(j) and xx(j+1).
! -- see Numerical Recipes in Fortran, Second Edition, p.110
! -- returns  values:
! --     0 if xx(1) > x and size(xx) if xx(size(xx)) < x
! --    -1, 1 if x < xx(1) or x > xx(size(xx)) respectively.
! --
! -- xx must be sorted.
  integer function locate(x,xx,status,frac) result(j)
    implicit none
    
    real, dimension(:), intent(in)    ::  xx
    real, intent(in)                  ::  x
    integer, intent(inout)            ::  status
    real, intent(inout), optional     ::  frac
    
    integer                           ::  jl, jm, ju
    integer                           ::  i, n
    
    n = size(xx)
    status = 0
    
!   -- start binary search.
    jl = 0
    ju = n + 1
    do
      if (ju-jl > 1) then
          jm = (ju + jl) / 2
        if ((xx(n) >= xx(1)) .EQV. (x >= xx(jm))) then
          jl = jm
        else
          ju = jm
        end if
      else
        exit
      endif
    end do

!   -- check endpoint equality, otherwise use jl from above.
    if (x == xx(1)) then 
      j = 1
    else if (x == xx(n)) then
      j = n-1
    else
      j = jl
    end if
    
!   -- set status, indicate success or appropriate failure condition.
    status = 0
    if (j >= n) status = 1
    if (j < 1) status = -1
    
    if (present(frac) .AND. status == 0) then
      if (j < n) then
        frac = (x-xx(j))/ (xx(j+1)-xx(j))
      else
        frac = 1.0
      end if
    end if
    
    return
     
  end function locate
  
! -- https://en.wikipedia.org/wiki/Extrapolation
  real function extrap(x, xx, yy, status) result(res)
    implicit none
    
    real, intent(in)                ::  x
    real, dimension(:), intent(in)  ::  xx
    real, dimension(:), intent(in)  ::  yy
    integer, intent(out)            ::  status
    
    real, dimension(2)              ::  r
    integer                         ::  n, i
    
    res = -999.0
    status = -1
    
    n = size(xx, 1) 
    status = linfit(xx(n-1:n), yy(n-1:n), r)
    if (status /= 0) then 
!       print *, "ERROR: linfit failed, skipping: ", status
      return
    end if
    
    res = r(1) + (x)*r(2)

!   -- if extrapolation produces a negative AOT, exclude the last node point
!   -- and try again.
    if (res < 3.5) then
      i = n-1
      do while (res < 3.5 .AND. i >= 2)
        status = linfit(xx(i-1:i), yy(i-1:i), r)
        if (status /= 0) then 
!           print *, "ERROR: linfit failed, skipping: ", status
          return
        end if
    
        res = r(1) + (x)*r(2)
        i = i-1 
      end do
      if (i < 2) then 
        status = 1
        return
      end if
    end if

    status = 0
    return    
    
  end function extrap
  
  integer function linfit(x, y, r) result(status)
    implicit none
    
    real, dimension(:), intent(in)    ::  x
    real, dimension(:), intent(in)    ::  y
    real, dimension(2), intent(inout) ::  r
    
    real                              ::  sx, sy
    real                              ::  sxx, syy, sxy
    
    integer                           ::  i, n
    
    n   = size(x)
    sx  = 0.0
    sy  = 0.0
    sxy = 0.0
    sxx = 0.0
    do i = 1, n
      sx = sx + x(i)
      sy = sy + y(i)
      sxx = sxx + (x(i) * x(i))
      sxy = sxy + (x(i) * y(i))
      syy = syy + (y(i) * y(i))
    end do

    if (abs((n*sxy) - (sx*sy)) < 1.0e-10 .or. abs((n*sxx)-(sx*sx))<1.0e-10 .or. n==0.) then 
     status = -1
    else 
     r(2) = ((n*sxy) - (sx*sy))/((n*sxx)-(sx*sx))
     r(1) = (sy/n)-(r(2)*sx/n)
     status = 0
    end if
    
    return
    
  end function linfit
  
!-------------------------------------------------------------
      subroutine new_intep(x1a,x2a,x3a,ya,m,n,l,ia,x1,x2,x3,y,dy, &
     &                     mbeg,nbeg,frac)
        implicit none
        
        real, dimension(:), intent(in)    ::  x1a
        real, dimension(:), intent(in)    ::  x2a
        real, dimension(:), intent(in)    ::  x3a
        real, dimension(:,:,:,:), intent(in)  ::  ya
        integer, intent(in)               ::  m, n, l
        integer, intent(in)               ::  ia
        real, intent(in)                  ::  x1, x2, x3
        real, intent(inout)               ::  y, dy
        integer, intent(in)               ::  mbeg, nbeg
        real, intent(in)                  ::  frac

        real, dimension(4)                ::  xx2a, xx1a          
        real, dimension(4)                ::  yntmp, ymtmp
        real, dimension(2)                ::  yltmp
        
        integer                           ::  j, k
        
      !dimension x1a(m),x2a(n),x3a(l),ya(4,4,2,10)
      !dimension xx2a(4), xx1a(4)
      !dimension yntmp(4),ymtmp(4),yltmp(2)
 
      do j=1,4
        do k=1,4
          yltmp(1)=ya(j,k,1,ia)
          yltmp(2)=ya(j,k,2,ia)
          yntmp(k) = yltmp(1)*(1.-frac) + yltmp(2)*frac
          xx2a(k) = x2a(k+nbeg)
        end do
        call polint(xx2a,yntmp,4,x2,ymtmp(j),dy)
        xx1a(j) = x1a(j+mbeg)
        
      end do

      call polint(xx1a,ymtmp,4,x1,y,dy)

      return
      end subroutine new_intep
      

      subroutine polint(xa,ya,n,x,y,dy)
      implicit none
      real, dimension(:), intent(in)    ::  xa
      real, dimension(:), intent(in)    ::  ya
      integer, intent(in)               ::  n
      real, intent(in)                  ::  x
      real, intent(inout)               ::  y
      real, intent(inout)               ::  dy
      
      integer, parameter    ::  nmax = 50
      real, dimension(nmax) ::  c
      real, dimension(nmax) ::  d
       
      integer   ::  i, m, ns
      real      ::  ho, hp, dif, dift
      real      ::  w, den
      
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n 
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
!          if(den.eq.0.) print *,'den=',den
!          if(den.eq.0.)pause
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue

      return

      end subroutine polint
      
    
end module viirs_aerosol_luts
