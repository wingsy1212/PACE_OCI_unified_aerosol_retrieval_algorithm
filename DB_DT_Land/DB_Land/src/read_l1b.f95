module read_l1b

private

! -- types
public  ::  viirs_db_svm
! -- functions
public  ::  load_viirs_db_data
public  ::  load_viirs_db_data_nasa
public  ::  load_h8_data
public  ::  load_goes_data
public  ::  load_pace_data,load_pace_test_data,write_pace_test_data

public  ::  latitude, longitude

integer, parameter  ::  i8  = selected_int_kind(3)
integer, parameter  ::  i64 = selected_int_kind(18)
integer, parameter	:: 	r64 = selected_real_kind(p=10)

type  ::  viirs_db_svm
  integer                             ::  scan
  integer                             ::  xscan
  integer                             ::  yr
  integer                             ::  mo
  integer                             ::  dy
  integer                             ::  doy
  integer                             ::  hr
  integer                             ::  min
  
  real(kind=r64), dimension(:), allocatable  ::  scan_time
   
  real, dimension(:,:), allocatable   ::  lat
  real, dimension(:,:), allocatable   ::  lon
  real, dimension(:,:), allocatable   ::  sza
  real, dimension(:,:), allocatable   ::  saa
  real, dimension(:,:), allocatable   ::  vza
  real, dimension(:,:), allocatable   ::  vaa
  real, dimension(:,:), allocatable   ::  raa
  real, dimension(:,:), allocatable   ::  sca
  real, dimension(:,:), allocatable   ::  gla
  real, dimension(:,:), allocatable   ::  amf
  
  real, dimension(:,:), allocatable   ::  uv01_refl
  real, dimension(:,:), allocatable   ::  uv02_refl
  real, dimension(:,:), allocatable   ::  m01_refl
  real, dimension(:,:), allocatable   ::  m02_refl
  real, dimension(:,:), allocatable   ::  m03_refl
  real, dimension(:,:), allocatable   ::  m04_refl
  real, dimension(:,:), allocatable   ::  m05_refl
  real, dimension(:,:), allocatable   ::  m07_refl
  real, dimension(:,:), allocatable   ::  m08_refl
  real, dimension(:,:), allocatable   ::  m09_refl
  real, dimension(:,:), allocatable   ::  m10_refl
  real, dimension(:,:), allocatable   ::  m11_refl

  real, dimension(:,:), allocatable   ::  m14_rad
  real, dimension(:,:), allocatable   ::  m14_bt
  real, dimension(:,:), allocatable   ::  m15_rad
  real, dimension(:,:), allocatable   ::  m15_bt
  real, dimension(:,:), allocatable   ::  m16_rad
  real, dimension(:,:), allocatable   ::  m16_bt
  real, dimension(:,:), allocatable   ::  m1120_bt
  real, dimension(:,:), allocatable   ::  m13_bt
  
  real, dimension(:,:), allocatable   ::  btd4
  real, dimension(:,:), allocatable   ::  btd8
  real, dimension(:,:), allocatable   ::  btd11
  real, dimension(:,:), allocatable   ::  dstar
  real, dimension(:,:), allocatable   ::  ndvi
  
  integer, dimension(:,:), allocatable  ::  land_mask
  real, dimension(:,:), allocatable     ::  elev
  real, dimension(:,:), allocatable     ::  ps
  integer, dimension(:), allocatable    ::  mirror_side
!  real(kind=r64), dimension(:,:), allocatable       ::  mnorm
  
end type viirs_db_svm

real, dimension(:,:), allocatable     ::  latitude, longitude, lat_copy, lon_copy

contains





subroutine load_pace_test_data(geo_file,lat,lon,sza,saz,vza,vaz,elev,&
                        UVtoSWIR_Radiances,Year,Month,Day,doy, &
                        grn_lwmask,UVAI,nXTrack,nLines,DTspec, DTAOD,dtlat, dtlon,&
                        Ret_CldMsk_500_Land,Ret_land_Quality_Flag,dtfmf,ncXTrack,ncLines, proxy, synth)
  use calendars, only: gdatetime, gregorian_from_doy,doy_from_gregorian

  implicit none
  
  include 'hdf.inc'
  include 'dffunc.inc'
  

                  
  character(len=*), intent(in)          ::  geo_file, proxy, synth
  integer, intent(inout)                  ::  year, month, day, doy

  real, dimension(:,:), allocatable , intent(out)   ::  lat,lon,sza, saz,vza,vaz,elev,UVAI
  real, dimension(:,:), intent(inout)   ::  dtlat, dtlon,Ret_land_Quality_Flag,dtfmf
  integer, dimension(:,:), intent(inout)::  Ret_CldMsk_500_Land
  
  real, dimension(:,:,:), allocatable , intent(out)   ::  UVtoSWIR_Radiances, DTspec,DTAOD
  integer, dimension(:,:), allocatable, intent(out)   ::  grn_lwmask
  real, dimension(:,:), allocatable     ::  dum_flag
  integer, intent(inout)                :: nXTrack,nLines 
  integer                               :: ncXTrack,ncLines,dtncLine,dtnLine

  integer                               ::  nc_id
  integer                               ::  dim_id
  integer                               ::  dset_id    
  integer                               ::  grp_id     
  
  character(len=255)                    ::  dset_name
  character(len=255)                    ::  attr_name
  character(len=255)                    ::  group_name
  character(len=5)                      ::  day_or_night
  character(len=24)                     ::  start_datetime
  character(len=10)                     ::  start_date
  character(len=5)                      ::  start_time  
  
  
  real, dimension(:,:), allocatable     ::  dum,cdum,dum2, dtdum
  integer, dimension(:,:), allocatable  ::  i_val
  integer                               ::  i_fillvalue
  integer, dimension(:,:), allocatable  ::  lw_mask
  real                                  ::  scale_factor
  real                                  ::  add_offset
  real, dimension(:), allocatable       ::  rad_lut
  real, dimension(:), allocatable       ::  bt_lut
  integer                               ::  indx
  real                                  ::  frac

  byte, dimension(:), allocatable       ::  flag_masks
  integer(kind=i8), dimension(:,:), allocatable     ::  qa
  integer, dimension(:,:), allocatable  ::  qa_mask
  integer                               ::  dtype
  
  integer                               ::  nvals,status

  integer                               ::  scan, xscan, nlut, nscan
  real, parameter                       ::  pi = 3.1415927, d2r = pi/180.0

  integer                               ::  i, j, k
  real                                  ::  tmp0, tmp1
  
  real(kind=r64), dimension(:), allocatable :: tmp_scan_time    ! tmp buffer
  integer, dimension(2) ::  start2, stride2, edges2, dim_sizes2
  integer :: sd_id, sds_index, sds_id, rank, n_attrs, data_type
  character(len=255)  ::  sds_name
  integer, dimension(32)::  dim_sizes
  type(gdatetime)       ::  gdt1
  real, dimension(14)   ::  UVtoSWIR_IRRadiances



  UVtoSWIR_IRRadiances=(/ 1000.98999, 1080.02991, 1089.37, 1695.31995, 1912.5, 1870.92993, &
  & 1514.62, 1476.18994, 1462.45007, 978.640015, 448.290009, 364.679993, 244.390015, 76.8499985/)
  if (doy /= 0) then 
   gdt1  = gregorian_from_doy(year, doy)
   month=gdt1%month
   day=gdt1%mday
  endif
!   year=2019
!   month=3
!   day=24

  gdt1 = gdatetime(year, month, day,12, 30, 0, 0, 0) 
 
  ncXTrack= nXTrack/7 !181
  ncLines=  nLines/7  !246
  dtnline = nLines
  IF (synth == 'NULL') THEN ;
   if (nxTrack >3000) then 
     ncXTrack= nXTrack/8 !181
     ncLines=  nLines/8  !246  
     if (nxTrack <=3232) then 
       dtncLine  = ncLines+2
       dtnLine = 3248
     endif
   endif
  endif

  IF (proxy == 'NULL') THEN
    dtncLine  = ncLines
  endif
    
  allocate(dum(nXTrack,nLines),lat(nXTrack,nLines),lon(nXTrack,nLines),sza(nXTrack,nLines),&
  & saz(nXTrack,nLines),vza(nXTrack,nLines),vaz(nXTrack,nLines),elev(nXTrack,nLines), &
  & UVtoSWIR_Radiances(14,nXTrack,nLines),dum_flag(nXTrack,nLines),UVAI(nXTrack,nLines),&
  & DTspec(ncXTrack,dtncLine,9), DTAOD(ncXTrack,dtncLine,4),cdum(ncXTrack,dtncLine), &
  & dtdum(nXTrack,dtnLine),stat=status)
!   & DTspec(ncXTrack,500,9), DTAOD(ncXTrack,500,4),cdum(ncXTrack,500),dum2(nXTrack,4000), stat=status) 
  sd_id = sfstart(geo_file, DFACC_READ)
  if (sd_id == FAIL) then
    print *, "ERROR: Unable to start SD interface on land aerosol model table: ", sd_id
    status = -1
    return
  end if

  sds_name = 'lat'
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
  
  start2  = (/0,0/)
  stride2 = (/1,1/)
  edges2  = shape(lat)
  status = sfrdata(sds_id, start2, stride2, edges2, lat)
  if (status == FAIL) then
    print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
    return
  end if
  
  status = sfendacc(sds_id)
  if (status == FAIL) then
    print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
    return
  end if
 
  latitude  = lat    

  sds_name = 'lon'
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
  
  start2  = (/0,0/)
  stride2 = (/1,1/)
  edges2  = shape(lon)
  status = sfrdata(sds_id, start2, stride2, edges2, lon)
  if (status == FAIL) then
    print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
    return
  end if
  
  status = sfendacc(sds_id)
  if (status == FAIL) then
    print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
    return
  end if
  longitude = lon
 
 
 
   sds_name = 'sza'
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
  
  start2  = (/0,0/)
  stride2 = (/1,1/)
  edges2  = shape(lat)
  status = sfrdata(sds_id, start2, stride2, edges2, sza)
  if (status == FAIL) then
    print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
    return
  end if
  
  status = sfendacc(sds_id)
  if (status == FAIL) then
    print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
    return
  end if


   sds_name = 'saz'
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
  
  start2  = (/0,0/)
  stride2 = (/1,1/)
  edges2  = shape(lat)
  status = sfrdata(sds_id, start2, stride2, edges2, saz)
  if (status == FAIL) then
    print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
    return
  end if
  
  status = sfendacc(sds_id)
  if (status == FAIL) then
    print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
    return
  end if     


   sds_name = 'vaz'
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
  
  start2  = (/0,0/)
  stride2 = (/1,1/)
  edges2  = shape(lat)
  status = sfrdata(sds_id, start2, stride2, edges2, vaz)
  if (status == FAIL) then
    print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
    return
  end if
  
  status = sfendacc(sds_id)
  if (status == FAIL) then
    print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
    return
  end if     



   sds_name = 'vza'
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
  
  start2  = (/0,0/)
  stride2 = (/1,1/)
  edges2  = shape(lat)
  status = sfrdata(sds_id, start2, stride2, edges2, vza)
  if (status == FAIL) then
    print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
    return
  end if
  
  status = sfendacc(sds_id)
  if (status == FAIL) then
    print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
    return
  end if     


 
  sds_name = 'Land_Sea_Flag'
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
  
  start2  = (/0,0/)
  stride2 = (/1,1/)
  edges2  = shape(lat)
  status = sfrdata(sds_id, start2, stride2, edges2, dum_flag)
  if (status == FAIL) then
    print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
    return
  end if
  grn_lwmask = int(dum_flag)
  
  
  status = sfendacc(sds_id)
  if (status == FAIL) then
    print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
    return
  end if      
  where (grn_lwmask == 0.0) grn_lwmask = 0. !land. 
  where (grn_lwmask >= 100) grn_lwmask = 3. !ocean. 

  sds_name = 'elev'
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
  
  start2  = (/0,0/)
  stride2 = (/1,1/)
  edges2  = shape(lat)
  status = sfrdata(sds_id, start2, stride2, edges2, elev)
  if (status == FAIL) then
    print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
    return
  end if
  
  status = sfendacc(sds_id)
  if (status == FAIL) then
    print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
    return
  end if       

  if (status /= 0) then
    print *, "ERROR: Failed to allocate reflectance arrays: ", status
    return
  end if

  sds_name = 'UVAI'
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
  
  start2  = (/0,0/)
  stride2 = (/1,1/)
  edges2  = shape(lat)
  status = sfrdata(sds_id, start2, stride2, edges2, UVAI)
  if (status == FAIL) then
    print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
    return
  end if
  
  status = sfendacc(sds_id)
  if (status == FAIL) then
    print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
    return
  end if       

  if (status /= 0) then
    print *, "ERROR: Failed to allocate reflectance arrays: ", status
    return
  end if
!   sds_name = 'rad1'
!   sds_index = sfn2index(sd_id, sds_name)
!   if (sds_index == FAIL) then
!     print *, "ERROR: Unable to find index of "//trim(sds_name)//": ", sds_index
!     status = -1
!     return
!   end if
!   
!   sds_id = sfselect(sd_id, sds_index)
!   if (sds_id == FAIL) then
!     print *, "ERROR: Unable to select SDS "//trim(sds_name)//": ", sds_id
!     status = -1
!     return
!   end if
!   
!   start2  = (/0,0/)
!   stride2 = (/1,1/)
!   edges2  = shape(lat)
!   status = sfrdata(sds_id, start2, stride2, edges2, dum)
!   if (status == FAIL) then
!     print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
!     return
!   end if
!   UVtoSWIR_Radiances(1,:,:)= dum!*cos(sza*d2r)
!   status = sfendacc(sds_id)
!   if (status == FAIL) then
!     print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
!     return
!   end if   
  
  sds_name = 'rad2'
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
  
  start2  = (/0,0/)
  stride2 = (/1,1/)
  edges2  = shape(lat)
  status = sfrdata(sds_id, start2, stride2, edges2, dum)
  if (status == FAIL) then
    print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
    return
  end if
  UVtoSWIR_Radiances(2,:,:)= dum!*cos(sza*d2r)
  
  status = sfendacc(sds_id)
  if (status == FAIL) then
    print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
    return
  end if   
  
  sds_name = 'rad3'
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
  
  start2  = (/0,0/)
  stride2 = (/1,1/)
  edges2  = shape(lat)
  status = sfrdata(sds_id, start2, stride2, edges2, dum)
  if (status == FAIL) then
    print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
    return
  end if
  UVtoSWIR_Radiances(3,:,:)= dum!*cos(sza*d2r)
  status = sfendacc(sds_id)
  if (status == FAIL) then
    print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
    return
  end if  
       
  sds_name = 'rad4'
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
  
  start2  = (/0,0/)
  stride2 = (/1,1/)
  edges2  = shape(lat)
  status = sfrdata(sds_id, start2, stride2, edges2, dum)
  if (status == FAIL) then
    print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
    return
  end if
  UVtoSWIR_Radiances(4,:,:)= dum!*cos(sza*d2r)
  
  status = sfendacc(sds_id)
  if (status == FAIL) then
    print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
    return
  end if     

  sds_name = 'rad5'
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
  
  start2  = (/0,0/)
  stride2 = (/1,1/)
  edges2  = shape(lat)
  status = sfrdata(sds_id, start2, stride2, edges2, dum)
  if (status == FAIL) then
    print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
    return
  end if
  
  status = sfendacc(sds_id)
  if (status == FAIL) then
    print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
    return
  end if    
  UVtoSWIR_Radiances(5,:,:)= dum!*cos(sza*d2r) 
    
  sds_name = 'rad6'
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
  
  start2  = (/0,0/)
  stride2 = (/1,1/)
  edges2  = shape(lat)
  status = sfrdata(sds_id, start2, stride2, edges2, dum)
  if (status == FAIL) then
    print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
    return
  end if
  
  status = sfendacc(sds_id)
  if (status == FAIL) then
    print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
    return
  end if      
  UVtoSWIR_Radiances(6,:,:)= dum!*cos(sza*d2r)
    
  sds_name = 'rad7'
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
  
  start2  = (/0,0/)
  stride2 = (/1,1/)
  edges2  = shape(lat)
  status = sfrdata(sds_id, start2, stride2, edges2, dum)
  if (status == FAIL) then
    print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
    return
  end if
  
  status = sfendacc(sds_id)
  if (status == FAIL) then
    print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
    return
  end if    
  UVtoSWIR_Radiances(7,:,:)= dum!*cos(sza*d2r)  
    
  sds_name = 'rad10'
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
  
  start2  = (/0,0/)
  stride2 = (/1,1/)
  edges2  = shape(lat)
  status = sfrdata(sds_id, start2, stride2, edges2, dum)
  if (status == FAIL) then
    print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
    return
  end if
  
  status = sfendacc(sds_id)
  if (status == FAIL) then
    print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
    return
  end if      
  UVtoSWIR_Radiances(10,:,:)= dum!*cos(sza*d2r)
    
  sds_name = 'rad11'
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
  
  start2  = (/0,0/)
  stride2 = (/1,1/)
  edges2  = shape(lat)
  status = sfrdata(sds_id, start2, stride2, edges2, dum)
  if (status == FAIL) then
    print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
    return
  end if
  
  status = sfendacc(sds_id)
  if (status == FAIL) then
    print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
    return
  end if     
  UVtoSWIR_Radiances(11,:,:)= dum!*cos(sza*d2r) 
    
   sds_name = 'rad12'
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
  
  start2  = (/0,0/)
  stride2 = (/1,1/)
  edges2  = shape(lat)
  status = sfrdata(sds_id, start2, stride2, edges2, dum)
  if (status == FAIL) then
    print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
    return
  end if
  
  status = sfendacc(sds_id)
  if (status == FAIL) then
    print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
    return
  end if   
  UVtoSWIR_Radiances(12,:,:)= dum!*cos(sza*d2r)  
    
  sds_name = 'rad13'
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
  
  start2  = (/0,0/)
  stride2 = (/1,1/)
  edges2  = shape(lat)
  status = sfrdata(sds_id, start2, stride2, edges2, dum)
  if (status == FAIL) then
    print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
    return
  end if
  
  status = sfendacc(sds_id)
  if (status == FAIL) then
    print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
    return
  end if  
  UVtoSWIR_Radiances(13,:,:)= dum!*cos(sza*d2r)    
    
  sds_name = 'rad14'
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
  
  start2  = (/0,0/)
  stride2 = (/1,1/)
  edges2  = shape(lat)
  status = sfrdata(sds_id, start2, stride2, edges2, dum)
  if (status == FAIL) then
    print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
    return
  end if
  
  status = sfendacc(sds_id)
  if (status == FAIL) then
    print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
    return
  end if      
  UVtoSWIR_Radiances(14,:,:)= dum!*cos(sza*d2r)

  sds_name = 'DTcld'
  
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

  start2  = (/0,0/)
  stride2 = (/1,1/)
  edges2  = shape(dum)

  status = sfrdata(sds_id, start2, stride2, edges2, dtdum)
  if (status == FAIL) then
    print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
    return
  end if
  print *, shape(dtdum)
  Ret_CldMsk_500_Land(:,:)= dtdum

  status = sfendacc(sds_id)
  if (status == FAIL) then
    print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
    return
  end if     

  do i = 1, 9
    if (i ==1) sds_name = 'DTspec1'
    if (i ==2) sds_name = 'DTspec2'
    if (i ==3) sds_name = 'DTspec3'
    if (i ==4) sds_name = 'DTspec4'
    if (i ==5) sds_name = 'DTspec5'
    if (i ==6) sds_name = 'DTspec6'
    if (i ==7) sds_name = 'DTspec7'
    if (i ==8) sds_name = 'DTspec8'
    if (i ==9) sds_name = 'DTspec9'
    
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
  
    start2  = (/0,0/)
    stride2 = (/1,1/)
    edges2  = shape(cdum)
    status = sfrdata(sds_id, start2, stride2, edges2, cdum)
    if (status == FAIL) then
      print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
      return
    end if
    DTspec(:,:,i)= cdum


    status = sfendacc(sds_id)
    if (status == FAIL) then
      print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
      return
    end if        
  enddo

  do i = 1, 4
    if (i ==1) sds_name = 'DTAOD1'
    if (i ==2) sds_name = 'DTAOD2'
    if (i ==3) sds_name = 'DTAOD3'
    if (i ==4) sds_name = 'DTAOD4'

    
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
  
    start2  = (/0,0/)
    stride2 = (/1,1/)
    edges2  = shape(cdum)
    status = sfrdata(sds_id, start2, stride2, edges2, cdum)
    if (status == FAIL) then
      print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
      return
    end if
    DTaod(:,:,i)= cdum

    status = sfendacc(sds_id)
    if (status == FAIL) then
      print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
      return
    end if        
  enddo

  sds_name = 'DTlat'
  
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

  start2  = (/0,0/)
  stride2 = (/1,1/)
  edges2  = shape(cdum)
  status = sfrdata(sds_id, start2, stride2, edges2, cdum)
  if (status == FAIL) then
    print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
    return
  end if
  print *, shape(cdum)
  print *, shape(dtlat)
  dtlat(:,:)= cdum

  status = sfendacc(sds_id)
  if (status == FAIL) then
    print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
    return
  end if   

  sds_name = 'DTlon'
  
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

  start2  = (/0,0/)
  stride2 = (/1,1/)
  edges2  = shape(cdum)
  status = sfrdata(sds_id, start2, stride2, edges2, cdum)
  if (status == FAIL) then
    print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
    return
  end if
  dtlon(:,:)= cdum

  status = sfendacc(sds_id)
  if (status == FAIL) then
    print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
    return
  end if   

  sds_name = 'DTQA'
  
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

  start2  = (/0,0/)
  stride2 = (/1,1/)
  edges2  = shape(cdum)
  status = sfrdata(sds_id, start2, stride2, edges2, cdum)
  if (status == FAIL) then
    print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
    return
  end if
  Ret_land_Quality_Flag(:,:)= cdum

  status = sfendacc(sds_id)
  if (status == FAIL) then
    print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
    return
  end if   

  sds_name = 'DTfmf'
  
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

  start2  = (/0,0/)
  stride2 = (/1,1/)
  edges2  = shape(cdum)
  status = sfrdata(sds_id, start2, stride2, edges2, cdum)
  if (status == FAIL) then
    print *, "ERROR: Unable to read SDS "//trim(sds_name)//": ", status
    return
  end if
  dtfmf(:,:)= cdum

  status = sfendacc(sds_id)
  if (status == FAIL) then
    print *, "ERROR: Unable to end access to SDS "//trim(sds_name)//": ", status
    return
  end if   




!    all_output_datasets(39)  = dataset('DTcld',Ret_CldMsk_500_Land)
!    all_output_datasets(40)  = dataset('DTQA',Ret_land_Quality_Flag)  
! Ret_CldMsk_500_Land,Ret_land_Quality_Flag
  status = sfend(sd_id)
  if (status /= 0) then
    print *, "ERROR: Unable to close land aerosol model file: ", status
    return
  end if  
  
end subroutine load_pace_test_data


subroutine write_pace_test_data(geo_file,lat,lon,sza,saz,&
                        vza,vaz,elev,UVtoSWIR_Radiances,Year,Month,Day,&
                        grn_lwmask,UVAI, DTspec,dtaod,dtLat,dtLon, &
                        Ret_CldMsk_500_Land,Ret_land_Quality_Flag,dt_fmf,Ret_Xtrack,Ret_Lines)
  use db_debug 
  use calendars, only: gdatetime, gregorian_from_doy,doy_from_gregorian
                      
  implicit none
  include 'hdf.inc'
  include 'dffunc.inc'
 
                  
  character(len=*), intent(in)          ::  geo_file
  integer, intent(in)                   ::  year, month, day,Ret_Xtrack,Ret_Lines
  integer                               ::  nXTrack,nLines

  real, dimension(:,:), allocatable , intent(in)   :: lat,lon,sza, saz,vza,vaz,elev,UVAI
  real, dimension(:,:), intent(in)   :: dtLat,dtLon,Ret_land_Quality_Flag,dt_fmf
  real, dimension(:,:,:), intent(in)               :: dtspec,dtaod
  integer, dimension(:,:), intent(in)              :: Ret_CldMsk_500_Land
  real, dimension(:,:,:), allocatable , intent(in) :: UVtoSWIR_Radiances
  integer, dimension(:,:), allocatable, intent(in) :: grn_lwmask
  real, dimension(:,:), allocatable                :: dum

  type(dataset), dimension(:), allocatable  :: all_output_datasets
  integer                                   ::  status, i
  type(output_file)   ::  debug_file
  integer, dimension(2) ::  sz



   allocate(all_output_datasets(41), stat=status)
   if (status /= 0) then 
     print *, "ERROR: Failed to allocate output dataset array: ", status
     stop
   end if 
   do i = 1, size(all_output_datasets,1)  
    allocate(all_output_datasets(i)%values(nXTrack,nLines), stat=status)
   end do
  sz  = shape(lat)
  nXTrack=sz(1)
  nLines= sz(2)
  
  allocate(dum(nXTrack,nLines), stat=status)

   all_output_datasets(1)  = dataset('lat', lat)
   all_output_datasets(2)  = dataset('lon', lon)
   all_output_datasets(3)  = dataset('sza', sza)
   all_output_datasets(4)  = dataset('saz', saz) 
   all_output_datasets(5)  = dataset('vza', vza) 
   all_output_datasets(6)  = dataset('vaz', vaz) 
   all_output_datasets(7)  = dataset('elev', elev)  
   all_output_datasets(8)  = dataset('Land_Sea_Flag', real(grn_lwmask)) 

   dum = UVtoSWIR_Radiances(1,:,:)
   all_output_datasets(9)  = dataset('rad1',dum )  
   dum = UVtoSWIR_Radiances(2,:,:)
   all_output_datasets(10)  = dataset('rad2',dum ) 
   dum = UVtoSWIR_Radiances(3,:,:)
   all_output_datasets(11)  = dataset('rad3',dum )  
   dum = UVtoSWIR_Radiances(4,:,:)
   all_output_datasets(12)  = dataset('rad4',dum ) 
   dum = UVtoSWIR_Radiances(5,:,:)
   all_output_datasets(13)  = dataset('rad5',dum )   
   dum = UVtoSWIR_Radiances(6,:,:)
   all_output_datasets(14)  = dataset('rad6',dum ) 
   dum = UVtoSWIR_Radiances(7,:,:)
   all_output_datasets(15)  = dataset('rad7',dum )  
   dum = UVtoSWIR_Radiances(8,:,:)
   all_output_datasets(16)  = dataset('rad8',dum ) 
   dum = UVtoSWIR_Radiances(9,:,:)
   all_output_datasets(17)  = dataset('rad9',dum )    

   dum = UVtoSWIR_Radiances(10,:,:)
   all_output_datasets(18)  = dataset('rad10',dum )   
   dum = UVtoSWIR_Radiances(11,:,:)
   all_output_datasets(19)  = dataset('rad11',dum ) 
   dum = UVtoSWIR_Radiances(12,:,:)
   all_output_datasets(20)  = dataset('rad12',dum )  
   dum = UVtoSWIR_Radiances(13,:,:)
   all_output_datasets(21)  = dataset('rad13',dum ) 
   dum = UVtoSWIR_Radiances(14,:,:)
   all_output_datasets(22)  = dataset('rad14',dum )      
   
   all_output_datasets(23)  = dataset('UVAI',UVAI )      

   all_output_datasets(24)  = dataset('DTspec1',dtspec(:,:,1)) 
   all_output_datasets(25)  = dataset('DTspec2',dtspec(:,:,2))
   all_output_datasets(26)  = dataset('DTspec3',dtspec(:,:,3))
   all_output_datasets(27)  = dataset('DTspec4',dtspec(:,:,4))
   all_output_datasets(28)  = dataset('DTspec5',dtspec(:,:,5))
   all_output_datasets(29)  = dataset('DTspec6',dtspec(:,:,6))
   all_output_datasets(30)  = dataset('DTspec7',dtspec(:,:,7))
   all_output_datasets(31)  = dataset('DTspec8',dtspec(:,:,8))
   all_output_datasets(32)  = dataset('DTspec9',dtspec(:,:,9))
   all_output_datasets(33)  = dataset('DTAOD1',dtaod(:,:,1))
   all_output_datasets(34)  = dataset('DTAOD2',dtaod(:,:,2))
   all_output_datasets(35)  = dataset('DTAOD3',dtaod(:,:,3))
   all_output_datasets(36)  = dataset('DTAOD4',dtaod(:,:,4))

   all_output_datasets(37)  = dataset('DTlon',dtlon)
   all_output_datasets(38)  = dataset('DTlat',dtlat)    
   all_output_datasets(39)  = dataset('DTcld',Ret_CldMsk_500_Land)
   all_output_datasets(40)  = dataset('DTQA',Ret_land_Quality_Flag)    
   all_output_datasets(41)  = dataset('DTfmf',dt_fmf)          
   print *, 'about to call debug...'
   debug_file = output_file(geo_file, all_output_datasets)
   status = db_debug_output(debug_file)
   if (status /= 0) then
     print *, "ERROR: Failed to dump debug output: ", status
     stop
   end if

   do i = 1, size(all_output_datasets,1)
     deallocate(all_output_datasets(i)%values, stat=status)
     if (status /= 0) then
       print *, 'ERROR: Failed to deallocate all_output_datasets(i)%values: ', i, status
       stop
     end if
   end do

   deallocate(all_output_datasets)  
  
end subroutine write_pace_test_data


type (viirs_db_svm) function load_pace_data(lat_in, lon_in, sza, saz, vza, vaz, elev, &
  & Land_Sea_Flag,UVtoSWIR_Reflectances, Year,month,day,doy, status) result (vdbs)


  use calendars, only: gdatetime, gregorian_from_doy,doy_from_gregorian

  implicit none
  

  REAL, dimension(:,:), intent(in)      ::  lat_in, lon_in
  REAL, dimension(:,:), intent(in)      ::  sza, saz, vza, vaz, elev
  REAL, dimension(:,:,:), intent(in)    ::  UVtoSWIR_Reflectances
  Integer, dimension(:,:), intent(in)   ::  Land_Sea_Flag 
                  
!   character(len=*), intent(in)          ::  geo_file
!   character(len=*), intent(in)          ::  mband_file
  integer, intent(in)                   ::  year, month, day,doy
  integer, intent(inout)                ::  status

  integer                               ::  nc_id
  integer                               ::  dim_id
  integer                               ::  dset_id    
  integer                               ::  grp_id     
  
  character(len=255)                    ::  dset_name
  character(len=255)                    ::  attr_name
  character(len=255)                    ::  group_name
  character(len=5)                      ::  day_or_night
  character(len=24)                     ::  start_datetime
  character(len=10)                     ::  start_date
  character(len=5)                      ::  start_time  
  
  real, dimension(:,:), allocatable     ::  saa, vaa
  integer, dimension(:,:), allocatable  ::  i_val
  integer                               ::  i_fillvalue
  integer, dimension(:,:), allocatable  ::  lw_mask
  real                                  ::  scale_factor
  real                                  ::  add_offset
  real, dimension(:), allocatable       ::  rad_lut
  real, dimension(:), allocatable       ::  bt_lut
  integer                               ::  indx
  real                                  ::  frac

  byte, dimension(:), allocatable       ::  flag_masks
  integer(kind=i8), dimension(:,:), allocatable     ::  qa
  integer, dimension(:,:), allocatable  ::  qa_mask
  integer                               ::  dtype
  
  integer                               ::  nvals

  integer                               ::  scan, xscan, nlut, nscan
  real, parameter                       ::  pi = 3.1415927, d2r = pi/180.0

  integer                               ::  i, j, k
  real                                  ::  tmp0, tmp1
  
  real(kind=r64), dimension(:), allocatable :: tmp_scan_time    ! tmp buffer
  integer, dimension(2) ::  start2, stride2, edges2, dim_sizes2
  integer :: sd_id, sds_index, sds_id, rank, n_attrs, data_type
  character(len=255)  ::  sds_name
  integer, dimension(32)::  dim_sizes
  type(gdatetime)       ::  gdt1
  real, dimension(14)   ::  UVtoSWIR_IRRadiances



  UVtoSWIR_IRRadiances=(/ 1000.98999, 1080.02991, 1089.37, 1695.31995, 1912.5, 1870.92993, &
  & 1514.62, 1476.18994, 1462.45007, 978.640015, 448.290009, 364.679993, 244.390015, 76.8499985/)
  vdbs%yr = year
  vdbs%mo = month
  vdbs%dy = day
  vdbs%hr = 12
  vdbs%min = 30
  if (doy==0) then 
    gdt1 = gdatetime(year, month, day,12, 0, 0, 0, 0) 
    vdbs%doy = doy_from_gregorian(gdt1)  
  else
    vdbs%doy = doy
  endif
  
!  print *,'year , doy:', year, doy
  vdbs%scan   = size(lat_in,2)
  vdbs%xscan  = size(lat_in,1)
  
  allocate(vdbs%lat(vdbs%xscan,vdbs%scan), vdbs%lon(vdbs%xscan,vdbs%scan), vdbs%sza(vdbs%xscan,vdbs%scan),  &
  &       vdbs%vza(vdbs%xscan,vdbs%scan), vdbs%raa(vdbs%xscan,vdbs%scan), vdbs%sca(vdbs%xscan,vdbs%scan),   &
  &       vdbs%saa(vdbs%xscan,vdbs%scan), vdbs%vaa(vdbs%xscan,vdbs%scan),   &
  &       vdbs%gla(vdbs%xscan,vdbs%scan), vdbs%amf(vdbs%xscan,vdbs%scan), vdbs%land_mask(vdbs%xscan,vdbs%scan),&
  &       vdbs%elev(vdbs%xscan,vdbs%scan), vdbs%ps(vdbs%xscan,vdbs%scan), stat=status)

  vdbs%lat  = lat_in
  vdbs%lon  = lon_in
  vdbs%sza  = sza
  vdbs%saa  = saz
  vdbs%vaa  = vaz
  vdbs%vza  = vza
  longitude = lon_in
  latitude  = lat_in
  

  vdbs%raa  = -999.0
  vdbs%raa  = (vdbs%vaa - vdbs%saa) - 180.0
  where (vdbs%raa > 180.0)  vdbs%raa = vdbs%raa - 360.0
  where (vdbs%raa < -180.0) vdbs%raa = vdbs%raa + 360.0
  where (vdbs%raa < 0.0)    vdbs%raa = -1.0 * vdbs%raa

  vdbs%sca  = acos(cos(vdbs%sza*d2r)*cos(vdbs%vza*d2r) -       &
                sin(vdbs%sza*d2r)*sin(vdbs%vza*d2r)*cos(vdbs%raa*d2r))
  vdbs%sca  = 180.0 - (vdbs%sca/d2r)

  vdbs%gla  = acos(cos(vdbs%sza*d2r)*cos(vdbs%vza*d2r) +       &
                sin(vdbs%sza*d2r)*sin(vdbs%vza*d2r)*cos(vdbs%raa*d2r))
  vdbs%gla = vdbs%gla/d2r

  vdbs%amf = 1.0/cos(vdbs%sza*d2r)+1.0/cos(vdbs%vza*d2r)

  vdbs%land_mask = Land_Sea_Flag
    
  where (vdbs%land_mask == 0.0) vdbs%land_mask = 0. !land. 
  where (vdbs%land_mask >= 100) vdbs%land_mask = 3. !ocean. 
  vdbs%elev = elev

  allocate(vdbs%m01_refl(vdbs%xscan,vdbs%scan), vdbs%m03_refl(vdbs%xscan,vdbs%scan), vdbs%m04_refl(vdbs%xscan,vdbs%scan),&
  &         vdbs%m05_refl(vdbs%xscan,vdbs%scan),vdbs%m07_refl(vdbs%xscan,vdbs%scan), vdbs%m08_refl(vdbs%xscan,vdbs%scan), &
  &         vdbs%m09_refl(vdbs%xscan,vdbs%scan), vdbs%m10_refl(vdbs%xscan,vdbs%scan), &
  &         vdbs%m11_refl(vdbs%xscan,vdbs%scan),vdbs%uv01_refl(vdbs%xscan,vdbs%scan), &
  &         vdbs%uv02_refl(vdbs%xscan,vdbs%scan) , stat=status)
  if (status /= 0) then
    print *, "ERROR: Failed to allocate reflectance arrays: ", status
    return
  end if
  vdbs%uv01_refl = UVtoSWIR_Reflectances(2,:,:)
  vdbs%uv02_refl = UVtoSWIR_Reflectances(3,:,:)
  vdbs%m01_refl = UVtoSWIR_Reflectances(4,:,:)
  vdbs%m03_refl = UVtoSWIR_Reflectances(5,:,:)
  vdbs%m04_refl = UVtoSWIR_Reflectances(6,:,:)
  vdbs%m05_refl = UVtoSWIR_Reflectances(7,:,:)
  vdbs%m07_refl = UVtoSWIR_Reflectances(10,:,:)
  vdbs%m08_refl = UVtoSWIR_Reflectances(11,:,:)
  vdbs%m09_refl = UVtoSWIR_Reflectances(12,:,:)
  vdbs%m10_refl = UVtoSWIR_Reflectances(13,:,:)
  vdbs%m11_refl = UVtoSWIR_Reflectances(14,:,:)
  
!   vdbs%m01_refl= vdbs%m01_refl*pi/UVtoSWIR_IRRadiances(4)
!   vdbs%m03_refl= vdbs%m03_refl*pi/UVtoSWIR_IRRadiances(5)
!   vdbs%m04_refl= vdbs%m04_refl*pi/UVtoSWIR_IRRadiances(6)
!   vdbs%m05_refl= vdbs%m05_refl*pi/UVtoSWIR_IRRadiances(7)
!   vdbs%m07_refl= vdbs%m07_refl*pi/UVtoSWIR_IRRadiances(10)
!   vdbs%m08_refl= vdbs%m08_refl*pi/UVtoSWIR_IRRadiances(11)
!   vdbs%m09_refl= vdbs%m09_refl*pi/UVtoSWIR_IRRadiances(12)
!   vdbs%m10_refl= vdbs%m10_refl*pi/UVtoSWIR_IRRadiances(13)
!   vdbs%m11_refl= vdbs%m11_refl*pi/UVtoSWIR_IRRadiances(14)


!   vdbs%m01_refl= vdbs%m01_refl/cos(vdbs%sza*d2r)
!   vdbs%m03_refl= vdbs%m03_refl/cos(vdbs%sza*d2r)
!   vdbs%m04_refl= vdbs%m04_refl/cos(vdbs%sza*d2r)
!   vdbs%m05_refl= vdbs%m05_refl/cos(vdbs%sza*d2r)
!   vdbs%m07_refl= vdbs%m07_refl/cos(vdbs%sza*d2r)
!   vdbs%m08_refl= vdbs%m08_refl/cos(vdbs%sza*d2r)
!   vdbs%m09_refl= vdbs%m09_refl/cos(vdbs%sza*d2r)
!   vdbs%m10_refl= vdbs%m10_refl/cos(vdbs%sza*d2r)
!   vdbs%m11_refl= vdbs%m11_refl/cos(vdbs%sza*d2r)

  allocate( vdbs%ndvi(vdbs%xscan,vdbs%scan), stat=status)!vdbs%dstar(vdbs%xscan,vdbs%scan),
  if (status /= 0) then
    print *, "ERROR: Failed to allocate D* arrays: ", status
    return
  end if
  
  vdbs%ndvi(:,:) = -999.0
  where (vdbs%m07_refl > -900.0 .AND.  vdbs%m05_refl > -900.0)
    vdbs%ndvi = (vdbs%m07_refl - vdbs%m05_refl) / (vdbs%m07_refl + vdbs%m05_refl)
  end where    

  
!   print *,'================================================================='
!   print *,'Temporarily reduce L1b size' 
!   print *,'================================================================='
!   where (vdbs%lat > 75. .or. vdbs%lat < 60. .or. vdbs%lon > -90. .or. vdbs%lon < -105.)
!   ! where (vdbs%lat > -22.63 .or. vdbs%lat < -28.03 .or. vdbs%lon > 180 .or. vdbs%lon < -180)
! ! where (vdbs%lat > 90 .or. vdbs%lat < -90 .or. vdbs%lon > -161.63 .or. vdbs%lon < -161.83)
! !   where (vdbs%lat > 36.5 .or. vdbs%lat < 36. .or. vdbs%lon > -75. .or. vdbs%lon < -75.5)
!         vdbs%m01_refl =-999.0
!         vdbs%m03_refl =-999.0
!         vdbs%m04_refl =-999.0
!         vdbs%m05_refl =-999.0
!         vdbs%m07_refl =-999.0
!         vdbs%m08_refl =-999.0
!         vdbs%m09_refl =-999.0 
!         vdbs%m10_refl =-999.0
!         vdbs%m11_refl =-999.0
!         
!   end where
  
  return
  
  
end function load_pace_data




type (viirs_db_svm) function load_h8_data(h8_file, status) result(h8)
  
  use netcdf
  use calendars, only: gdatetime, gregorian_from_doy,doy_from_gregorian
  implicit none
  
  character(len=*), intent(in)    ::  h8_file
  integer, intent(out)            ::  status

  character(len=255)    :: dset_name,str
  character(len=255)    :: attr_name
  character(len=255)    :: grp_name
  character(len=50)     :: fname
  character(len=24)     ::  start_date
  integer               ::  nc_id
  integer               ::  grp_id
  integer               ::  dim_id
  integer               ::  dset_id
  
  integer               ::  scan
  integer               ::  xscan
  integer               ::  ntime,fname_idx
  integer                         ::  i,j,screen_sum,i1,j1,i2,j2
  real, dimension(:,:), allocatable   ::  saa, vaa
  real, parameter                       ::  d2r =  3.14159/180.0  
  integer, dimension(:,:), allocatable  ::  itmp
  real, dimension(:,:), allocatable  ::  tmp
  type(gdatetime)       ::  gdt1
  real    ::  scale_factor, add_offset
  integer ::  fill_value
  integer, dimension(:), allocatable    ::  array_screen
  real, dimension(:), allocatable       ::  sub_array1,sub_array2
  integer, dimension(:,:), allocatable  ::  minsza
  
  status = -1
        
  status = nf90_open(h8_file, nf90_nowrite, nc_id)
  str= "ERROR: Failed to open H8 file: "
  call check(status,str)

  fname_idx = index(h8_file, 'HS_H08_')
  fname = h8_file(fname_idx: fname_idx+50)
  start_date = fname(8:20)
  read(start_date(1:4), fmt='(I4)') h8%yr
  read(start_date(5:6), fmt='(I2)') h8%mo
  read(start_date(7:8), fmt='(I2)') h8%dy
  read(start_date(10:11), fmt='(I2)') h8%hr
  read(start_date(12:13), fmt='(I2)') h8%min
  gdt1 = gdatetime(h8%yr, h8%mo, h8%dy, h8%hr, h8%min, 0, 0, 0) 
  h8%doy = doy_from_gregorian(gdt1)  

! -- 2.0 km bands
! ------------------------------------------------------------------------------
  grp_name = '2.0_km'
  status = nf90_inq_ncid(nc_id, grp_name, grp_id)
  str= "ERROR: Failed to get ID of group "//trim(grp_name)//": "
  call check(status,str)
  
  dset_name = 'line'
  status = nf90_inq_dimid(grp_id, dset_name, dim_id)
  str= "ERROR: Failed to get ID of dimension "//trim(dset_name)//": "
  call check(status,str)
  
  status = nf90_inquire_dimension(grp_id, dim_id, len=scan)
  str= "ERROR: Failed to get size of dimension "//trim(dset_name)//": "
  call check(status,str)

  dset_name = 'column'
  status = nf90_inq_dimid(grp_id, dset_name, dim_id)
  str= "ERROR: Failed to get ID of dimension "//trim(dset_name)//": "
  call check(status,str)
  
  status = nf90_inquire_dimension(grp_id, dim_id, len=xscan)
  str= "ERROR: Failed to size of dimension "//trim(dset_name)//": "
  call check(status,str)

  h8%scan = scan
  h8%xscan  = xscan
  
  print *, 'xscan, scan: ', xscan, scan   
  allocate(h8%lat(xscan,scan), h8%lon(xscan,scan), h8%m10_refl(xscan,scan), h8%m11_refl(xscan,scan), &
  &       h8%m14_bt(xscan,scan),h8%m15_bt(xscan,scan),h8%m16_bt(xscan,scan), &
  &       h8%raa(xscan,scan), h8%sca(xscan,scan), h8%sza(xscan,scan), h8%vza(xscan,scan), &
  &       saa(xscan,scan), vaa(xscan,scan), h8%saa(xscan,scan), h8%vaa(xscan,scan), &
  &       itmp(xscan,scan),h8%land_mask(xscan,scan), h8%elev(xscan,scan), h8%ps(xscan,scan),&
  &       h8%amf(xscan,scan),h8%m1120_bt(xscan,scan),h8%m13_bt(xscan,scan),minsza(xscan,scan), stat=status)
  h8%land_mask(:,:) = -999
  str= "ERROR: Failed to allocate 0.5 km data arrays: "
  call check(status,str)

  allocate(latitude(xscan,scan), longitude(xscan,scan), stat=status)
  str= "ERROR: Failed to allocate latitude and longitude arrays: "
  call check(status,str)
  
  dset_name = 'Latitude'
  status = nf90_inq_varid(grp_id, dset_name, dset_id)
  str= "ERROR: Failed to get ID of dataset "//trim(dset_name)//": "
  call check(status,str)
  
  status = nf90_get_var(grp_id, dset_id, h8%lat)
  str= "ERROR: Failed to read dataset "//trim(dset_name)//": "
  call check(status,str)
!   where (h8%lat < -900.0) h8%lat = -999. 
  latitude = h8%lat
   
  
  dset_name = 'Longitude'
  status = nf90_inq_varid(grp_id, dset_name, dset_id)
  str= "ERROR: Failed to get ID of dataset "//trim(dset_name)//": "
  call check(status,str)
  
  status = nf90_get_var(grp_id, dset_id, h8%lon)
  str= "ERROR: Failed to read dataset "//trim(dset_name)//": "
  call check(status,str)
!   where (h8%lon < -900.0) h8%lon = -999. 
  longitude = h8%lon
  
  dset_name = 'SatelliteAzimuthAngle'
  status = nf90_inq_varid(grp_id, dset_name, dset_id)
  str=  "ERROR: Failed to get ID of dataset "//trim(dset_name)//": "
  call check(status,str)
  
  attr_name = 'scale_factor'
  status = nf90_get_att(grp_id, dset_id, attr_name, scale_factor)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)
  
  attr_name = 'add_offset'
  status = nf90_get_att(grp_id, dset_id, attr_name, add_offset)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)

  attr_name = '_FillValue'
  status = nf90_get_att(grp_id, dset_id, attr_name, fill_value)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)
  
  status = nf90_get_var(grp_id, dset_id, itmp)
  str= "ERROR: Failed to read dataset "//trim(dset_name)//": "
  call check(status,str)
  
  vaa(:,:) = -999.0
  where(itmp /= fill_value)
    vaa = itmp * scale_factor + add_offset
  end where
  h8%vaa = vaa
!   where (h8%vaa < -900.0) h8%vaa = -999. 
  
  dset_name = 'SolarAzimuthAngle'
  status = nf90_inq_varid(grp_id, dset_name, dset_id)
  str= "ERROR: Failed to get ID of dataset "//trim(dset_name)//": "
  call check(status,str)
  
  attr_name = 'scale_factor'
  status = nf90_get_att(grp_id, dset_id, attr_name, scale_factor)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)
  
  attr_name = 'add_offset'
  status = nf90_get_att(grp_id, dset_id, attr_name, add_offset)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)
  
  attr_name = '_FillValue'
  status = nf90_get_att(grp_id, dset_id, attr_name, fill_value)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)
  
  status = nf90_get_var(grp_id, dset_id, itmp)
  str= "ERROR: Failed to read dataset "//trim(dset_name)//": "
  call check(status,str)
  
  saa(:,:) = -999.0
  where(itmp /= fill_value)
    saa = itmp * scale_factor + add_offset
  end where
  
! -- SAA from files is measured from the south while VAA is measured from the north (???).
! -- flip the SAA to the north so they're consistent.
  saa = saa - 180.0
  where (saa < 0.0) saa = saa + 360.0
  h8%saa = saa
!   where (h8%saa < -900.0) h8%saa = -999. 
  
  dset_name = 'SolarZenithAngle'
  status = nf90_inq_varid(grp_id, dset_name, dset_id)
  str= "ERROR: Failed to get ID of dataset "//trim(dset_name)//": "
  call check(status,str)
  
  attr_name = 'scale_factor'
  status = nf90_get_att(grp_id, dset_id, attr_name, scale_factor)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)
  
  attr_name = 'add_offset'
  status = nf90_get_att(grp_id, dset_id, attr_name, add_offset)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)
  
  attr_name = '_FillValue'
  status = nf90_get_att(grp_id, dset_id, attr_name, fill_value)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)
  
  status = nf90_get_var(grp_id, dset_id, itmp)
  str= "ERROR: Failed to read dataset "//trim(dset_name)//": "
  call check(status,str)
  
  h8%sza(:,:) = -999.0
  where(itmp /= fill_value)
    h8%sza = itmp * scale_factor + add_offset
  end where
!   where (h8%sza < -900.0) h8%sza = -999. 
  
  dset_name = 'SatelliteZenithAngle'
  status = nf90_inq_varid(grp_id, dset_name, dset_id)
  str= "ERROR: Failed to get ID of dataset "//trim(dset_name)//": "
  call check(status,str)
  
  attr_name = 'scale_factor'
  status = nf90_get_att(grp_id, dset_id, attr_name, scale_factor)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)
  
  attr_name = 'add_offset'
  status = nf90_get_att(grp_id, dset_id, attr_name, add_offset)
  str=  "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)
  
  attr_name = '_FillValue'
  status = nf90_get_att(grp_id, dset_id, attr_name, fill_value)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)

  
  status = nf90_get_var(grp_id, dset_id, itmp)
  str= "ERROR: Failed to read dataset "//trim(dset_name)//": "
  call check(status,str)

  h8%vza(:,:) = -999.0
  where(itmp /= fill_value)
    h8%vza = itmp * scale_factor + add_offset
  end where
!   where (h8%vza < -900.0) h8%vza = -999. 
  
  dset_name = 'B05_Albedo'
  status = nf90_inq_varid(grp_id, dset_name, dset_id)
  str= "ERROR: Failed to get ID of dataset "//trim(dset_name)//": "
  call check(status,str)
  
  attr_name = 'scale_factor'
  status = nf90_get_att(grp_id, dset_id, attr_name, scale_factor)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)
  
  attr_name = 'add_offset'
  status = nf90_get_att(grp_id, dset_id, attr_name, add_offset)
  str=  "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)
  
  attr_name = '_FillValue'
  status = nf90_get_att(grp_id, dset_id, attr_name, fill_value)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)
  
  status = nf90_get_var(grp_id, dset_id, itmp)
  str= "ERROR: Failed to read dataset "//trim(dset_name)//": "
  call check(status,str)
  
  h8%m10_refl(:,:) = -999.0
  where(itmp /= fill_value)
    h8%m10_refl = itmp * scale_factor + add_offset
  end where
!   where (h8%m10_refl < -900.0) h8%m10_refl = -999. 
  
  dset_name = 'B06_Albedo'
  status = nf90_inq_varid(grp_id, dset_name, dset_id)
  str= "ERROR: Failed to get ID of dataset "//trim(dset_name)//": "
  call check(status,str)
  
  attr_name = 'scale_factor'
  status = nf90_get_att(grp_id, dset_id, attr_name, scale_factor)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)

  attr_name = 'add_offset'
  status = nf90_get_att(grp_id, dset_id, attr_name, add_offset)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)
  
  attr_name = '_FillValue'
  status = nf90_get_att(grp_id, dset_id, attr_name, fill_value)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)
  
  status = nf90_get_var(grp_id, dset_id, itmp)
  str= "ERROR: Failed to read dataset "//trim(dset_name)//": "
  call check(status,str)

  h8%m11_refl(:,:) = -999.0
  where(itmp /= fill_value)
    h8%m11_refl = itmp * scale_factor + add_offset
  end where
  
  dset_name = 'B07_BT'
  status = nf90_inq_varid(grp_id, dset_name, dset_id)
  str= "ERROR: Failed to get ID of dataset "//trim(dset_name)//": "
  call check(status,str)
  
  attr_name = 'scale_factor'
  status = nf90_get_att(grp_id, dset_id, attr_name, scale_factor)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)

  attr_name = 'add_offset'
  status = nf90_get_att(grp_id, dset_id, attr_name, add_offset)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)
  
  attr_name = '_FillValue'
  status = nf90_get_att(grp_id, dset_id, attr_name, fill_value)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)
  
  status = nf90_get_var(grp_id, dset_id, itmp)
  str= "ERROR: Failed to read dataset "//trim(dset_name)//": "
  call check(status,str)

  h8%m13_bt(:,:) = -999.0
  where(itmp /= fill_value)
    h8%m13_bt = itmp * scale_factor + add_offset
  end where
    
!   where (h8%m11_refl < -900.0) h8%m11_refl = -999. 
 
  dset_name = 'B11_BT'
  status = nf90_inq_varid(grp_id, dset_name, dset_id)
  str= "ERROR: Failed to get ID of dataset "//trim(dset_name)//": "
  call check(status,str)
  
  attr_name = 'scale_factor'
  status = nf90_get_att(grp_id, dset_id, attr_name, scale_factor)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)
  
  attr_name = 'add_offset'
  status = nf90_get_att(grp_id, dset_id, attr_name, add_offset)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)
  
  attr_name = '_FillValue'
  status = nf90_get_att(grp_id, dset_id, attr_name, fill_value)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)
  
  status = nf90_get_var(grp_id, dset_id, itmp)
  str= "ERROR: Failed to read dataset "//trim(dset_name)//": "
  call check(status,str)
  
  h8%m14_bt(:,:) = -999.0
  where(itmp /= fill_value)
    h8%m14_bt = itmp * scale_factor + add_offset
  end where
!   where (h8%m14_bt < -900.0) h8%m14_bt = -999. 
  
  dset_name = 'B13_BT'
  status = nf90_inq_varid(grp_id, dset_name, dset_id)
  str= "ERROR: Failed to get ID of dataset "//trim(dset_name)//": "
  call check(status,str)

  attr_name = 'scale_factor'
  status = nf90_get_att(grp_id, dset_id, attr_name, scale_factor)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)
  
  attr_name = 'add_offset'
  status = nf90_get_att(grp_id, dset_id, attr_name, add_offset)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)
  
  attr_name = '_FillValue'
  status = nf90_get_att(grp_id, dset_id, attr_name, fill_value)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)
  
  status = nf90_get_var(grp_id, dset_id, itmp)
  str= "ERROR: Failed to read dataset "//trim(dset_name)//": "
  call check(status,str)
  
  h8%m15_bt(:,:) = -999.0
  where(itmp /= fill_value)
    h8%m15_bt = itmp * scale_factor + add_offset
  end where
!   where (h8%m15_bt < -900.0) h8%m15_bt = -999. 

  dset_name = 'B14_BT'
  status = nf90_inq_varid(grp_id, dset_name, dset_id)
  str= "ERROR: Failed to get ID of dataset "//trim(dset_name)//": "
  call check(status,str)

  attr_name = 'scale_factor'
  status = nf90_get_att(grp_id, dset_id, attr_name, scale_factor)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)
  
  attr_name = 'add_offset'
  status = nf90_get_att(grp_id, dset_id, attr_name, add_offset)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)
  
  attr_name = '_FillValue'
  status = nf90_get_att(grp_id, dset_id, attr_name, fill_value)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)
  
  status = nf90_get_var(grp_id, dset_id, itmp)
  str= "ERROR: Failed to read dataset "//trim(dset_name)//": "
  call check(status,str)
  
  h8%m1120_bt(:,:) = -999.0
  where(itmp /= fill_value)
    h8%m1120_bt = itmp * scale_factor + add_offset
  end where
  
  dset_name = 'B15_BT'
  status = nf90_inq_varid(grp_id, dset_name, dset_id)
  str= "ERROR: Failed to get ID of dataset "//trim(dset_name)//": "
  call check(status,str)
  
  attr_name = 'scale_factor'
  status = nf90_get_att(grp_id, dset_id, attr_name, scale_factor)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)
  
  attr_name = 'add_offset'
  status = nf90_get_att(grp_id, dset_id, attr_name, add_offset)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)
  
  attr_name = '_FillValue'
  status = nf90_get_att(grp_id, dset_id, attr_name, fill_value)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)
  
  status = nf90_get_var(grp_id, dset_id, itmp)
  str=  "ERROR: Failed to read dataset "//trim(dset_name)//": "
  call check(status,str)
  
  h8%m16_bt(:,:) = -999.0
  where(itmp /= fill_value)
    h8%m16_bt = itmp * scale_factor + add_offset
  end where
!   where (h8%m16_bt < -900.0) h8%m16_bt = -999. 
  
  h8%raa  = (vaa - saa) - 180.0
  where (h8%raa > 180.0)  h8%raa = h8%raa - 360.0
  where (h8%raa < -180.0) h8%raa = h8%raa + 360.0
  where (h8%raa < 0.0)    h8%raa = -1.0 * h8%raa

  h8%sca  = acos(cos(h8%sza*d2r)*cos(h8%vza*d2r) -       &
                sin(h8%sza*d2r)*sin(h8%vza*d2r)*cos(h8%raa*d2r))
  h8%sca  = 180.0 - (h8%sca/d2r)
  h8%amf = 1.0/cos(h8%sza*d2r)+1.0/cos(h8%vza*d2r)
  deallocate(itmp, saa, vaa, stat=status)
  str= "ERROR: Failed to deallocate tmp array: "
  call check(status,str)

! -- 1.0 km band
! ------------------------------------------------------------------------------
  grp_name = '1.0_km'
  status = nf90_inq_ncid(nc_id, grp_name, grp_id)
  str= "ERROR: Failed to get ID of group "//trim(grp_name)//": "
  call check(status,str)
  
  dset_name = 'line'
  status = nf90_inq_dimid(grp_id, dset_name, dim_id)
  str= "ERROR: Failed to get ID of dimension "//trim(dset_name)//": "
  call check(status,str)
  
  status = nf90_inquire_dimension(grp_id, dim_id, len=scan)
  str= "ERROR: Failed to get size of dimension "//trim(dset_name)//": "
  call check(status,str)

  dset_name = 'column'
  status = nf90_inq_dimid(grp_id, dset_name, dim_id)
  str= "ERROR: Failed to get ID of dimension "//trim(dset_name)//": "
  call check(status,str)
  
  status = nf90_inquire_dimension(grp_id, dim_id, len=xscan)
  str= "ERROR: Failed to size of dimension "//trim(dset_name)//": "
  call check(status,str)

  print *, 'xscan, scan: ', xscan, scan  
    
  allocate(h8%m03_refl(h8%xscan,h8%scan), h8%m04_refl(h8%xscan,h8%scan), &
  &       h8%m07_refl(h8%xscan,h8%scan), itmp(xscan,scan), tmp(xscan,scan), &
  &       stat=status)
  str= "ERROR: Failed to allocate 1.0 km data arrays: "
  call check(status,str)

  dset_name = 'B01_Albedo'
  status = nf90_inq_varid(grp_id, dset_name, dset_id)
  str= "ERROR: Failed to get ID of dataset "//trim(dset_name)//": "
  call check(status,str)

  
  attr_name = 'scale_factor'
  status = nf90_get_att(grp_id, dset_id, attr_name, scale_factor)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)
  
  attr_name = 'add_offset'
  status = nf90_get_att(grp_id, dset_id, attr_name, add_offset)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)
  
  attr_name = '_FillValue'
  status = nf90_get_att(grp_id, dset_id, attr_name, fill_value)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)

  
  status = nf90_get_var(grp_id, dset_id, itmp)
  str= "ERROR: Failed to read dataset "//trim(dset_name)//": "
  call check(status,str)
  
  h8%m03_refl(:,:) = -999.0
  tmp(:,:) = -999.0
  where(itmp /= fill_value)
    tmp = itmp * scale_factor + add_offset
  end where
  
  status = rescale_array(tmp, h8%m03_refl, 2)
  str= "ERROR: Failed to rescale B1 array: "
  call check(status,str)
!   where (h8%m03_refl < -900.0) h8%m03_refl = -999. 
  
  dset_name = 'B02_Albedo'
  status = nf90_inq_varid(grp_id, dset_name, dset_id)
  str= "ERROR: Failed to get ID of dataset "//trim(dset_name)//": "
  call check(status,str)
  
  attr_name = 'scale_factor'
  status = nf90_get_att(grp_id, dset_id, attr_name, scale_factor)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)
  
  attr_name = 'add_offset'
  status = nf90_get_att(grp_id, dset_id, attr_name, add_offset)
  str=  "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)
  
  attr_name = '_FillValue'
  status = nf90_get_att(grp_id, dset_id, attr_name, fill_value)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)
  
  status = nf90_get_var(grp_id, dset_id, itmp)
  str= "ERROR: Failed to read dataset "//trim(dset_name)//": "
  call check(status,str)
  
  h8%m04_refl(:,:) = -999.0
  tmp(:,:) = -999.0
  where(itmp /= fill_value)
    tmp = itmp * scale_factor + add_offset
  end where
  
  status = rescale_array(tmp, h8%m04_refl, 2)
  str= "ERROR: Failed to rescale B2 array: "
!   where (h8%m04_refl < -900.0) h8%m04_refl = -999. 
  
  dset_name = 'B04_Albedo'
  status = nf90_inq_varid(grp_id, dset_name, dset_id)
  str= "ERROR: Failed to get ID of dataset "//trim(dset_name)//": "
  call check(status,str)
  
  attr_name = 'scale_factor'
  status = nf90_get_att(grp_id, dset_id, attr_name, scale_factor)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)
  
  attr_name = 'add_offset'
  status = nf90_get_att(grp_id, dset_id, attr_name, add_offset)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)
  
  attr_name = '_FillValue'
  status = nf90_get_att(grp_id, dset_id, attr_name, fill_value)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)
  
  status = nf90_get_var(grp_id, dset_id, itmp)
  str= "ERROR: Failed to read dataset "//trim(dset_name)//": "
  call check(status,str)
  
  h8%m07_refl(:,:) = -999.0
  tmp(:,:) = -999.0
  where(itmp /= fill_value)
    tmp = itmp * scale_factor + add_offset
  end where
  
  status = rescale_array(tmp, h8%m07_refl, 2)
  str=  "ERROR: Failed to rescale B4 array: "
  call check(status,str)

!   where (h8%m07_refl < -900.0) h8%m07_refl = -999. 
  
  deallocate(itmp, tmp, stat=status)
  str= "ERROR: Failed to deallocate tmp array: "
  call check(status,str)
  
! -- 0.5 km band
! ------------------------------------------------------------------------------
  grp_name = '0.5_km'
  status = nf90_inq_ncid(nc_id, grp_name, grp_id)
  str= "ERROR: Failed to get ID of group "//trim(grp_name)//": "
  call check(status,str)
  
  dset_name = 'line'
  status = nf90_inq_dimid(grp_id, dset_name, dim_id)
  str= "ERROR: Failed to get ID of dimension "//trim(dset_name)//": "
  call check(status,str)
  
  status = nf90_inquire_dimension(grp_id, dim_id, len=scan)
  str= "ERROR: Failed to get size of dimension "//trim(dset_name)//": "
  call check(status,str)

  dset_name = 'column'
  status = nf90_inq_dimid(grp_id, dset_name, dim_id)
  str= "ERROR: Failed to get ID of dimension "//trim(dset_name)//": "
  call check(status,str)
  
  status = nf90_inquire_dimension(grp_id, dim_id, len=xscan)
  str= "ERROR: Failed to size of dimension "//trim(dset_name)//": "
  call check(status,str)
  print *, 'xscan, scan: ', xscan, scan  
    
  allocate(h8%m05_refl(h8%xscan,h8%scan), itmp(xscan,scan), tmp(xscan,scan), stat=status)
  str= "ERROR: Failed to allocate 0.5 km data arrays: "
  call check(status,str)

  dset_name = 'B03_Albedo'
  status = nf90_inq_varid(grp_id, dset_name, dset_id)
  str= "ERROR: Failed to get ID of dataset "//trim(dset_name)//": "
  call check(status,str)
  
  attr_name = 'scale_factor'
  status = nf90_get_att(grp_id, dset_id, attr_name, scale_factor)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)
  
  attr_name = 'add_offset'
  status = nf90_get_att(grp_id, dset_id, attr_name, add_offset)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)
  
  attr_name = '_FillValue'
  status = nf90_get_att(grp_id, dset_id, attr_name, fill_value)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)
  
  status = nf90_get_var(grp_id, dset_id, itmp)
  str= "ERROR: Failed to read dataset "//trim(dset_name)//": "
  call check(status,str)
  
  h8%m05_refl(:,:) = -999.0
  tmp(:,:) = -999.0
  where(itmp /= fill_value)
    tmp = itmp * scale_factor + add_offset
  end where
  status = rescale_array(tmp, h8%m05_refl, 4)
  str= "ERROR: Failed to rescale B3 array: "
  call check(status,str)
!   where (h8%m05_refl < -900.0) h8%m05_refl = -999. 
  
  deallocate(itmp, tmp, stat=status)
  str= "ERROR: Failed to deallocate tmp array: "
  call check(status,str)  

! -- by default, the reflectances are in units I*pi / F0. In order to match VIIRS and MODIS,
! -- convert to I*pi / u0 * F0 by dividing by the cos(sza).
 where (h8%m05_refl > -900.0) h8%m05_refl = h8%m05_refl / cos(h8%sza*d2r)
 where (h8%m04_refl > -900.0) h8%m04_refl = h8%m04_refl / cos(h8%sza*d2r)
 where (h8%m03_refl > -900.0) h8%m03_refl = h8%m03_refl / cos(h8%sza*d2r)
 where (h8%m07_refl > -900.0) h8%m07_refl = h8%m07_refl / cos(h8%sza*d2r)
 where (h8%m10_refl > -900.0)h8%m10_refl = h8%m10_refl / cos(h8%sza*d2r)
 where (h8%m11_refl > -900.0)h8%m11_refl = h8%m11_refl / cos(h8%sza*d2r)
  
! -- secondary, derived products  
  allocate(h8%ndvi(h8%xscan,h8%scan), h8%btd8(h8%xscan,h8%scan),  &
  &       h8%btd11(h8%xscan,h8%scan), h8%dstar(h8%xscan,h8%scan), h8%btd4(h8%xscan,h8%scan), stat=status)
  str= "ERROR: Failed to allocate NDVI array: "
  call check(status,str)  
  
  h8%ndvi(:,:) = -999.0
  where (h8%m07_refl > -900.0 .AND. h8%m05_refl > -900.0)
    h8%ndvi = (h8%m07_refl - h8%m05_refl) / (h8%m07_refl + h8%m05_refl)
  end where
  
  h8%btd8   = h8%m14_bt - h8%m15_bt
  h8%btd11  = h8%m15_bt - h8%m16_bt
  h8%btd4   = h8%m13_bt - h8%m1120_bt
  status = calc_dstar(h8%btd8, h8%btd11, h8%dstar)

! remove lon/lat grid value when reflectance is nan value
  where (h8%m03_refl < -900. .and. h8%m05_refl < -900.)
       h8%lat =-999.0
       h8%lon =-999.0
  end where  
  latitude          = h8%lat
  longitude         = h8%lon  

!  print *,'================================================================='
!  print *,'Temporarily reduce L1b size' 
!  print *,'================================================================='
!  where (h8%lat > 80. .or. h8%lat < -80. .or. h8%lon > 180 .or. h8%lon < 120.)
!        h8%m10_refl =-999.0
!        h8%m11_refl =-999.0
!        h8%m14_bt   =-999.0
!        h8%m15_bt   =-999.0
!        h8%m16_bt   =-999.0
!        h8%m03_refl =-999.0
!        h8%m04_refl =-999.0
!        h8%m05_refl =-999.0
!        h8%m07_refl =-999.0
!        h8%ndvi     =-999.0
!        h8%btd8     =-999.0
!        h8%btd11    =-999.0
!        h8%dstar    =-999.0
!  end where

  minsza  = h8%sza
  where (minsza < -900.)
       minsza =9999.
  end where  
  
  if (minval(minsza) > 85.) then
    print *, "ERROR: Granule is a night granule. Skipping.", minval(minsza)
    status = -1
    return
  end if

end function load_h8_data


type (viirs_db_svm) function load_goes_data(goes_file,goes_f1,goes_f2,goes_f3,goes_f4,&
        & goes_f5,goes_f6,goes_f7,goes_f11,goes_f13,goes_f14,goes_f15, segment, status) result(goes_seg)
  
  use netcdf
  use calendars, only: gdatetime, gregorian_from_doy,doy_from_gregorian
  implicit none
  
  type(viirs_db_svm)       ::  goes
  character(len=*), intent(in)    ::  goes_file,goes_f1,goes_f2,goes_f3,goes_f4,goes_f5,goes_f6
  character(len=*), intent(in)    ::  goes_f7,goes_f11,goes_f13,goes_f14,goes_f15
  character(len=*), intent(in)    ::  segment
  integer, intent(out)            ::  status

  character(len=255)    :: dset_name, prefix, surfix, goes_file_read
  character(len=255)    :: attr_name,goes_time,str
  character(len=255)    :: grp_name
  character(len=50)     :: fname
  character(len=24)     :: start_date
  character(len=2)      :: band_number, goes_ver
  
  integer               ::  nc_id,do_bt,ridx
  integer               ::  grp_id
  integer               ::  dim_id
  integer               ::  dset_id,seg_int,seg_sz,st_indx,fi_indx
  
  integer               ::  scan
  integer               ::  xscan, fdoy
  integer               ::  ntime,fname_idx
  integer                         ::  i,j,screen_sum,i1,j1,i2,j2
  real, dimension(:,:), allocatable   ::  dummy
  real, parameter                       ::  d2r =  3.14159/180.0  
  integer, dimension(:,:), allocatable  ::  itmp
  real, dimension(:,:), allocatable  ::  tmp
  type(gdatetime)       ::  gdt1
  real                  ::  scale_factor, add_offset, pi, dist, esun,glon, glat
  real*8                ::  time, ftime
  integer ::  fill_value
  integer, dimension(:), allocatable  ::  array_screen
  real, dimension(:), allocatable  ::  sub_array1,sub_array2
  integer, dimension(:,:), allocatable  ::  minsza

  pi =  3.14159 
  status = -1

!   prefix=goes_file(1:32)  
!   surfix=goes_file(35:90)
  ridx  = index(goes_file, '/',.true.)
!   goes_time=goes_file(41:54)
  goes_time=goes_file(ridx+28:ridx+41)
  read(goes_time(1:4), fmt='(I4)') goes_seg%yr
  read(goes_time(5:7), fmt='(I3)') goes_seg%doy 
  read(goes_time(8:9), fmt='(I2)') goes_seg%hr
  read(goes_time(10:11), fmt='(I2)') goes_seg%min

  gdt1  = gregorian_from_doy(goes_seg%yr, goes_seg%doy)
  goes_seg%mo=gdt1%month
  goes_seg%dy=gdt1%mday
  
  ! -- Initialization & Read Band 06 (2km resolution)
  ! -------------------------------------------------------------------------------
  band_number='06'
  goes_file_read = goes_f6  !trim(prefix) // trim(band_number) // trim(surfix)
  print *, trim(goes_file_read)
  goes_ver  = goes_file_read(38:39)
  
  print *, 'read GOES',goes_ver
  status = nf90_open(goes_file_read, nf90_nowrite, nc_id) 
  str= "ERROR: Failed to open GOES file: "
  call check(status,str)

  dset_name = 'y'
  status = nf90_inq_dimid(nc_id, dset_name, dim_id)
  str= "ERROR: Failed to get size of dimension "//trim(dset_name)//": "
  call check(status,str)  
  status = nf90_inquire_dimension(nc_id, dim_id, len = scan)
  str= "ERROR: Failed to size of dimension "//trim(dset_name)//": "
  call check(status,str)
 
  dset_name = 'x'
  status = nf90_inq_dimid(nc_id, dset_name, dim_id)
  str= "ERROR: Failed to get size of dimension "//trim(dset_name)//": "
  call check(status,str)  
  status = nf90_inquire_dimension(nc_id, dim_id, len = xscan)
  str= "ERROR: Failed to size of dimension "//trim(dset_name)//": "
  call check(status,str)    

  dset_name = 't'
  status = nf90_inq_varid(nc_id, dset_name, dset_id)
  str= "ERROR: Failed to get ID of dataset "//trim(dset_name)//": "
  call check(status,str)

  status = nf90_get_var(nc_id, dset_id, time)
  str= "ERROR: Failed to read dataset "//trim(dset_name)//": "
  call check(status,str)  

  dset_name = 'geospatial_lat_lon_extent'
  status = nf90_inq_varid(nc_id, dset_name, dset_id)
  str= "ERROR: Failed to get ID of dataset "//trim(dset_name)//": "
  call check(status,str)
    
  attr_name = 'geospatial_lat_center'
  status = nf90_get_att(nc_id, dset_id, attr_name, glat)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)
    
  attr_name = 'geospatial_lon_center'
  status = nf90_get_att(nc_id, dset_id, attr_name, glon)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)

  goes%scan   = scan
  goes%xscan  = xscan    
  print *, 'xscan, scan: ', xscan, scan
  
  read(segment,'(I2)')  seg_int
  seg_int = seg_int - 1
  seg_sz  = 544   !about 1/10 of scan also multiple of 4
  st_indx = (1 + (seg_sz*seg_int))-4    !4 pixel as buffer
  if (st_indx .lt. 1) st_indx=1
  fi_indx = ((seg_sz*(seg_int+1)))+4  !4 pixels as buffer
  if (fi_indx .gt. scan) fi_indx=scan
  seg_sz  = fi_indx - st_indx +1
  
  allocate(goes%lat(xscan,scan), goes%lon(xscan,scan), goes%m03_refl(xscan,scan),goes%m05_refl(xscan,scan), &
  &       goes%m07_refl(xscan,scan),goes%m09_refl(xscan,scan), goes%m10_refl(xscan,scan), &
  &       goes%m11_refl(xscan,scan), goes%m14_bt(xscan,scan), goes%m15_bt(xscan,scan),  &
  &       goes%m16_bt(xscan,scan),goes%raa(xscan,scan), goes%sca(xscan,scan), goes%sza(xscan,scan), &
  &       goes%vza(xscan,scan), dummy(xscan,scan), goes%saa(xscan,scan), goes%vaa(xscan,scan), &
  &       itmp(xscan,scan),goes%land_mask(xscan,scan), goes%elev(xscan,scan), goes%ps(xscan,scan),&
  &       goes%amf(xscan,scan),goes%m1120_bt(xscan,scan),goes%m13_bt(xscan,scan), stat=status)  

  goes%land_mask(:,:) = -999
  goes%lat(:,:) = -999
  goes%lon(:,:) = -999
  goes%m03_refl(:,:) = -999
  goes%m05_refl(:,:) = -999
  goes%m07_refl(:,:) = -999
  goes%m09_refl(:,:) = -999
  goes%m10_refl(:,:) = -999
  goes%m11_refl(:,:) = -999
  goes%m14_bt(:,:) = -999
  goes%m15_bt(:,:) = -999
  goes%m16_bt(:,:) = -999
  goes%raa(:,:) = -999
  goes%sca(:,:) = -999
  goes%sza(:,:) = -999
  goes%vza(:,:) = -999
  goes%saa(:,:) = -999
  goes%vaa(:,:) = -999
  goes%elev(:,:) = -999
  goes%ps(:,:) = -999
  goes%amf(:,:) = -999
  goes%m1120_bt(:,:) = -999
  goes%m13_bt(:,:) = -999
  
  str= "ERROR: Failed to allocate 2 km data arrays: "
  call check(status,str)    
  
  allocate(latitude(xscan,scan), longitude(xscan,scan),lat_copy(xscan,scan), lon_copy(xscan,scan), stat=status)
  str= "ERROR: Failed to allocate latitude and longitude arrays: "
  call check(status,str)  
  
  !define 2km lon lat grid
  call goes_geo(nc_id,xscan,scan,longitude,latitude, 2, goes_ver)
  goes%lat  = latitude 
  goes%lon  = longitude
  lat_copy  = latitude 
  lon_copy  = longitude  
  dset_name = 'Rad'
  do_bt     = 0
  call goes_read_rad(nc_id, dset_name,itmp, dummy,do_bt ) 
  goes%m11_refl = dummy 
  
  !define time, sat angle, solar angle
  fdoy  = goes_seg%doy
  ftime = ((mod(time,86400.))/3600.)+12.
 
  if (ftime>24) then
    fdoy  = fdoy+1
    ftime = ftime-24
  end if 
  call POSSOL(fdoy, ftime, longitude,latitude, goes%sza, goes%saa)
  call COMPUTE_SATELLITE_ANGLES(glon, glat, longitude,latitude, goes%vza, goes%vaa)
  
!   goes%raa  = (goes%vaa - goes%saa) - 180.0
!   where (goes%raa > 180.0)  goes%raa = goes%raa - 360.0
!   where (goes%raa < -180.0) goes%raa = goes%raa + 360.0
!   where (goes%raa < 0.0)    goes%raa = -1.0 * goes%raa
  goes%raa  = (goes%saa - goes%vaa)
  where (goes%raa < 0.0)    goes%raa = -1.0 * goes%raa
  where (goes%raa > 180.0)  goes%raa = 360.0 -goes%raa
  
  goes%sca  = acos(cos(goes%sza*d2r)*cos(goes%vza*d2r) -       &
                sin(goes%sza*d2r)*sin(goes%vza*d2r)*cos(goes%raa*d2r))
  goes%sca  = 180.0 - (goes%sca/d2r)
  goes%amf = 1.0/cos(goes%sza*d2r)+1.0/cos(goes%vza*d2r)
  
  ! -- Read Band 04 (2km resolution)
  ! ------------------------------------------------------------------------------ 
  band_number='04'
  goes_file_read = goes_f4  !trim(prefix) // trim(band_number) // trim(surfix)
  print *, trim(goes_file_read)
  status = nf90_open(goes_file_read, nf90_nowrite, nc_id) 
  str= "ERROR: Failed to open GOES file: "
  call check(status,str)
    
  dset_name = 'Rad'
  do_bt     = 0
  call goes_read_rad(nc_id, dset_name,itmp, dummy,do_bt ) 
  goes%m09_refl = dummy
  str= "ERROR: Failed to rescale m04_refl array: "  
  
  ! -- Read Band 7 BT (2km resolution)
  ! ------------------------------------------------------------------------------ 
  band_number='07'
  goes_file_read = goes_f7  !trim(prefix) // trim(band_number) // trim(surfix)
  print *, trim(goes_file_read)
  status = nf90_open(goes_file_read, nf90_nowrite, nc_id) 
  str= "ERROR: Failed to open GOES file: "
  call check(status,str)
    
  dset_name = 'Rad'
  do_bt     = 1
  call goes_read_rad(nc_id, dset_name,itmp, dummy,do_bt ) 
  goes%m13_bt = dummy 
    
  ! -- Read Band 11 BT (2km resolution)
  ! ------------------------------------------------------------------------------ 
  band_number='11'
  goes_file_read = goes_f11 !trim(prefix) // trim(band_number) // trim(surfix)
  print *, trim(goes_file_read)
  status = nf90_open(goes_file_read, nf90_nowrite, nc_id) 
  str= "ERROR: Failed to open GOES file: "
  call check(status,str)
    
  dset_name = 'Rad'
  do_bt     = 1
  call goes_read_rad(nc_id, dset_name,itmp, dummy,do_bt ) 
  goes%m14_bt = dummy 

  ! -- Read Band 13 BT (2km resolution)
  ! ------------------------------------------------------------------------------ 
  band_number='13'
  goes_file_read = goes_f13 !trim(prefix) // trim(band_number) // trim(surfix)
  print *, trim(goes_file_read)
  status = nf90_open(goes_file_read, nf90_nowrite, nc_id) 
  str= "ERROR: Failed to open GOES file: "
  call check(status,str)
    
  dset_name = 'Rad'
  do_bt     = 1
  call goes_read_rad(nc_id, dset_name,itmp, dummy,do_bt ) 
  goes%m15_bt = dummy 

  ! -- Read Band 15 BT (2km resolution)
  ! ------------------------------------------------------------------------------ 
  band_number='15'
  goes_file_read = goes_f15 !trim(prefix) // trim(band_number) // trim(surfix)
  print *, trim(goes_file_read)
  status = nf90_open(goes_file_read, nf90_nowrite, nc_id) 
  str= "ERROR: Failed to open GOES file: "
  call check(status,str)
    
  dset_name = 'Rad'
  do_bt     = 1
  call goes_read_rad(nc_id, dset_name,itmp, dummy,do_bt ) 
  goes%m16_bt = dummy  
  
  ! -- Read Band 14 BT (2km resolution)
  ! ------------------------------------------------------------------------------ 
  band_number='14'
  goes_file_read = goes_f14 !trim(prefix) // trim(band_number) // trim(surfix)
  print *, trim(goes_file_read)
  status = nf90_open(goes_file_read, nf90_nowrite, nc_id) 
  str= "ERROR: Failed to open GOES file: "
  call check(status,str)
    
  dset_name = 'Rad'
  do_bt     = 1
  call goes_read_rad(nc_id, dset_name,itmp, dummy,do_bt ) 
  goes%m1120_bt = dummy  

  ! -- Initialization & Read Band 01 (1km resolution)
  ! -------------------------------------------------------------------------------
  band_number='01'
  goes_file_read = goes_f1  !trim(prefix) // trim(band_number) // trim(surfix)
  print *, trim(goes_file_read)
  status = nf90_open(goes_file_read, nf90_nowrite, nc_id) 
  str= "ERROR: Failed to open GOES file: "
  call check(status,str)

  dset_name = 'y'
  status = nf90_inq_dimid(nc_id, dset_name, dim_id)
  str= "ERROR: Failed to get size of dimension "//trim(dset_name)//": "
  call check(status,str)  
  status = nf90_inquire_dimension(nc_id, dim_id, len = scan)
  str= "ERROR: Failed to size of dimension "//trim(dset_name)//": "
  call check(status,str)
 
  dset_name = 'x'
  status = nf90_inq_dimid(nc_id, dset_name, dim_id)
  str= "ERROR: Failed to get size of dimension "//trim(dset_name)//": "
  call check(status,str)  
  status = nf90_inquire_dimension(nc_id, dim_id, len = xscan)
  str= "ERROR: Failed to size of dimension "//trim(dset_name)//": "
  call check(status,str)     

  deallocate(dummy, itmp, stat=status)
  str= "ERROR: Failed to deallocate tmp array: "
  call check(status,str)    

  allocate(dummy(xscan,scan),itmp(xscan,scan), stat=status)    
  str= "ERROR: Failed to allocate 1 km data arrays: "
  call check(status,str)    
   
  dset_name = 'Rad'
  do_bt     = 0
  call goes_read_rad(nc_id, dset_name,itmp, dummy,do_bt ) 
  
  status = rescale_array(dummy, goes%m03_refl, 2)
  str= "ERROR: Failed to rescale m03_refl array: "

  ! -- Read Band 03 (1km resolution)
  ! ------------------------------------------------------------------------------ 
  band_number='03'
  goes_file_read = goes_f3  !trim(prefix) // trim(band_number) // trim(surfix)
  print *, trim(goes_file_read)
  status = nf90_open(goes_file_read, nf90_nowrite, nc_id) 
  str= "ERROR: Failed to open GOES file: "
  call check(status,str)
    
  dset_name = 'Rad'
  do_bt     = 0
  call goes_read_rad(nc_id, dset_name,itmp, dummy,do_bt ) 
  status = rescale_array(dummy, goes%m07_refl, 2)
  str= "ERROR: Failed to rescale m07_refl array: "

  ! -- Read Band 05 (1km resolution)
  ! ------------------------------------------------------------------------------ 
  band_number='05'
  goes_file_read = goes_f5  !trim(prefix) // trim(band_number) // trim(surfix)
  print *, trim(goes_file_read)
  status = nf90_open(goes_file_read, nf90_nowrite, nc_id) 
  str= "ERROR: Failed to open GOES file: "
  call check(status,str)
    
  dset_name = 'Rad'
  do_bt     = 0
  call goes_read_rad(nc_id, dset_name,itmp, dummy,do_bt ) 
  status = rescale_array(dummy, goes%m10_refl, 2)
  str= "ERROR: Failed to rescale m10_refl array: "

  ! -- Initialization & Read Band 02 (0.5km resolution)
  ! -------------------------------------------------------------------------------
  band_number='02'
  goes_file_read = goes_f2  !trim(prefix) // trim(band_number) // trim(surfix)
  print *, trim(goes_file_read)
  status = nf90_open(goes_file_read, nf90_nowrite, nc_id) 
  str= "ERROR: Failed to open GOES file: "
  call check(status,str)

  dset_name = 'y'
  status = nf90_inq_dimid(nc_id, dset_name, dim_id)
  str= "ERROR: Failed to get size of dimension "//trim(dset_name)//": "
  call check(status,str)  
  status = nf90_inquire_dimension(nc_id, dim_id, len = scan)
  str= "ERROR: Failed to size of dimension "//trim(dset_name)//": "
  call check(status,str)
 
  dset_name = 'x'
  status = nf90_inq_dimid(nc_id, dset_name, dim_id)
  str= "ERROR: Failed to get size of dimension "//trim(dset_name)//": "
  call check(status,str)  
  status = nf90_inquire_dimension(nc_id, dim_id, len = xscan)
  str= "ERROR: Failed to size of dimension "//trim(dset_name)//": "
  call check(status,str)
       
  deallocate(dummy, itmp, stat=status)
  str= "ERROR: Failed to deallocate tmp array: "
  call check(status,str)    
  allocate(dummy(xscan,scan),itmp(xscan,scan), stat=status)    
  str= "ERROR: Failed to allocate 0.5 km data arrays: "
  call check(status,str)    
   
  dset_name = 'Rad'
  do_bt     = 0
  call goes_read_rad(nc_id, dset_name,itmp, dummy,do_bt ) 
  
  status = rescale_array(dummy, goes%m05_refl, 4)
  str= "ERROR: Failed to rescale m03_refl array: "  

! -- secondary, derived products  
  allocate(goes%ndvi(goes%xscan,goes%scan), goes%btd8(goes%xscan,goes%scan),  &
  &       goes%btd11(goes%xscan,goes%scan), goes%dstar(goes%xscan,goes%scan), &
  &       goes%btd4(goes%xscan,goes%scan), stat=status)
  str= "ERROR: Failed to allocate NDVI array: "
  call check(status,str)  
  
  goes%ndvi(:,:) = -999.0
  goes%btd8(:,:) = -999.0
  goes%btd11(:,:) = -999.0
  goes%dstar(:,:) = -999.0
  goes%btd4(:,:) = -999.0
  
  where (goes%m07_refl > -900.0 .AND. goes%m05_refl > -900.0)
    goes%ndvi = (goes%m07_refl - goes%m05_refl) / (goes%m07_refl + goes%m05_refl)
  end where
  
  goes%btd8  = goes%m14_bt - goes%m15_bt
  goes%btd11 = goes%m15_bt - goes%m16_bt
  goes%btd4  = goes%m13_bt - goes%m1120_bt
  status = calc_dstar(goes%btd8, goes%btd11, goes%dstar)

! -- reduce to segment size
  deallocate(latitude, longitude, stat=status)
  allocate(goes_seg%lat(goes%xscan,seg_sz), goes_seg%lon(goes%xscan,seg_sz), &
  &       goes_seg%m03_refl(goes%xscan,seg_sz),goes_seg%m05_refl(goes%xscan,seg_sz), &
  &       goes_seg%m07_refl(goes%xscan,seg_sz),goes_seg%m09_refl(goes%xscan,seg_sz), &
  &       goes_seg%m10_refl(goes%xscan,seg_sz),goes_seg%m11_refl(goes%xscan,seg_sz), &
  &       goes_seg%m14_bt(goes%xscan,seg_sz), goes_seg%m15_bt(goes%xscan,seg_sz),  &
  &       goes_seg%m16_bt(goes%xscan,seg_sz),goes_seg%raa(goes%xscan,seg_sz), &
  &       goes_seg%sca(goes%xscan,seg_sz), goes_seg%sza(goes%xscan,seg_sz), &
  &       goes_seg%vza(goes%xscan,seg_sz),goes_seg%saa(goes%xscan,seg_sz),  &
  &       goes_seg%vaa(goes%xscan,seg_sz),goes_seg%land_mask(goes%xscan,seg_sz),&
  &       goes_seg%elev(goes%xscan,seg_sz), goes_seg%ps(goes%xscan,seg_sz),&
  &       goes_seg%amf(goes%xscan,seg_sz),goes_seg%m1120_bt(goes%xscan,seg_sz),&
  &       goes_seg%m13_bt(goes%xscan,seg_sz), goes_seg%ndvi(goes%xscan,seg_sz), &
  &       goes_seg%btd8(goes%xscan,seg_sz), goes_seg%btd11(goes%xscan,seg_sz),&
  &       goes_seg%dstar(goes%xscan,seg_sz), goes_seg%btd4(goes%xscan,seg_sz),&
  &       latitude(goes%xscan,seg_sz), longitude(goes%xscan,seg_sz),minsza(goes%xscan,seg_sz),stat=status) 
!   latitude          =lat_copy(:,st_indx:fi_indx)
!   longitude         =lon_copy(:,st_indx:fi_indx)
  goes_seg%xscan    =goes%xscan
  goes_seg%scan     =seg_sz
  goes_seg%lat      =goes%lat(:,st_indx:fi_indx)
  goes_seg%lon      =goes%lon(:,st_indx:fi_indx)
  goes_seg%m03_refl =goes%m03_refl(:,st_indx:fi_indx)
  goes_seg%m05_refl =goes%m05_refl(:,st_indx:fi_indx)
  goes_seg%m07_refl =goes%m07_refl(:,st_indx:fi_indx)
  goes_seg%m09_refl =goes%m09_refl(:,st_indx:fi_indx)
  goes_seg%m10_refl =goes%m10_refl(:,st_indx:fi_indx)
  goes_seg%m11_refl =goes%m11_refl(:,st_indx:fi_indx)
  goes_seg%m14_bt   =goes%m14_bt(:,st_indx:fi_indx)
  goes_seg%m15_bt   =goes%m15_bt(:,st_indx:fi_indx)
  goes_seg%m16_bt   =goes%m16_bt(:,st_indx:fi_indx)
  goes_seg%raa      =goes%raa(:,st_indx:fi_indx)
  goes_seg%sca      =goes%sca(:,st_indx:fi_indx)
  goes_seg%sza      =goes%sza(:,st_indx:fi_indx)
  goes_seg%vza      =goes%vza(:,st_indx:fi_indx)
  goes_seg%saa      =goes%saa(:,st_indx:fi_indx)
  goes_seg%vaa      =goes%vaa(:,st_indx:fi_indx)
  goes_seg%land_mask=goes%land_mask(:,st_indx:fi_indx)
  goes_seg%elev     =goes%elev(:,st_indx:fi_indx)
  goes_seg%ps       =goes%ps(:,st_indx:fi_indx)
  goes_seg%amf      =goes%amf(:,st_indx:fi_indx)
  goes_seg%m1120_bt =goes%m1120_bt(:,st_indx:fi_indx)
  goes_seg%m13_bt   =goes%m13_bt(:,st_indx:fi_indx)
  goes_seg%ndvi     =goes%ndvi(:,st_indx:fi_indx)
  goes_seg%btd8     =goes%btd8(:,st_indx:fi_indx)
  goes_seg%btd11    =goes%btd11(:,st_indx:fi_indx)
  goes_seg%dstar    =goes%dstar(:,st_indx:fi_indx)
  goes_seg%btd4     =goes%btd4(:,st_indx:fi_indx)

! remove lon/lat grid value when reflectance is nan value
 where (goes_seg%m03_refl < -900. .and. goes_seg%m05_refl < -900.)
       goes_seg%lat =-999.0
       goes_seg%lon =-999.0
 end where  
  latitude          = goes_seg%lat
  longitude         = goes_seg%lon

!   print *,'================================================================='
!   print *,'Temporarily reduce L1b size' 
!   print *,'================================================================='
! !   where (goes%lat > 61. .or. goes%lat < 59. .or. goes%lon > -89.0 .or. goes%lon < -91.)
!   where (goes%lat > 80. .or. goes%lat < 0. .or. goes%lon > -60. .or. goes%lon < -150.)
!         goes%m10_refl =-999.0
!         goes%m11_refl =-999.0
!         goes%m14_bt   =-999.0
!         goes%m15_bt   =-999.0
!         goes%m16_bt   =-999.0
!         goes%m03_refl =-999.0
!         goes%m05_refl =-999.0
!         goes%m07_refl =-999.0
!         goes%ndvi     =-999.0
!         goes%btd8     =-999.0
!         goes%btd11    =-999.0
!         goes%dstar    =-999.0
!   end where
  
  minsza  = goes_seg%sza
  where (minsza < -900.)
       minsza =9999.
  end where  
  
  if (minval(minsza) > 85.) then
    print *, "ERROR: Granule is a night granule. Skipping."
    status = -1
    return
  end if
end function load_goes_data



 subroutine COMPUTE_SATELLITE_ANGLES(glon, glat, xlon, xlat, zenith, azimuth)
!--------------------------------------------------------------
! Subroutine to make geostationary satellite azimuth field
!   
!     xlon = longitude of the location (positive for western hemisphere) 
!     xlat = latitude of the location  (positive for northern hemisphere) 
!
!     zenith  = satellite zenith view angle 
!     azimuth = satellite azimuth angle clockwise from north
!--------------------------------------------------------------
!  arguments
 !  integer,                         intent(in)  :: isat
   real,               intent(in)  :: glon, glat
   real, dimension(:,:), intent(in)  :: xlon, xlat
   real, dimension(:,:), intent(out) :: zenith, azimuth

!  Local variables
   real, parameter :: pi = 3.14159265, dtor = pi / 180.0
   real :: satlon, satlat, lat, lon, beta, sin_beta
   integer :: i, j, nx, ny

   nx = size(xlon,dim=1)
   ny = size(xlon,dim=2)
   
   zenith = -999.
   azimuth = -999.

   satlon = glon
   satlat = glat


   do j = 1, ny
     do i = 1, nx
     
      if (xlon(i,j) < -900) cycle

      lon = (xlon(i,j) - satlon) * dtor   ! in radians
      lat = (xlat(i,j) - satlat) * dtor   ! in radians

      beta = acos( cos(lat) * cos(lon) )
      sin_beta = sin(beta)

!     zenith angle      
      zenith(i,j) = asin(max(-1.0, min(1.0, &
              42164.0* sin_beta/ sqrt(1.808e09 - 5.3725e08*cos(beta)))))
      zenith(i,j) = zenith(i,j) / dtor

!     azimuth angle
      azimuth(i,j) = sin(lon) / sin_beta
      azimuth(i,j) = min(1.0, max(-1.0,azimuth(i,j)))
      azimuth(i,j) = asin(azimuth(i,j))
      azimuth(i,j) = azimuth(i,j) / dtor
      if (lat < 0.0) then
         azimuth(i,j) = 180.0 - azimuth(i,j)
      endif
      if (azimuth(i,j) < 0.0) then
         azimuth(i,j) = azimuth(i,j) + 360.0
      endif

     enddo
   enddo

 end subroutine COMPUTE_SATELLITE_ANGLES
 

subroutine POSSOL(jday, tu, xlon, xlat, asol, phis)
!---------------------------------------------------------------------
! Compute solar angles
!
! input:
!         jday = julian day
!         tu   = time of day - fractional hours
!         lon  = latitude in degrees
!         lat  = latitude in degrees
!      
!  output:
!         asol = solar zenith angle in degrees
!         phis = solar azimuth angle in degrees
!
!----------------------------------------------------------------------
    implicit none
!   Arguments
    integer,                        intent(in)  :: jday 
    real*8,                         intent(in)  :: tu
    real, dimension(:, :), intent(in)           :: xlon, xlat
    real, dimension(:, :), intent(out)          :: asol, phis

!   Local variables
    real, parameter :: pi = 3.14159265, dtor = pi / 180.0
    real :: tsm, xlo, xla, xj, a1, a2, et, tsv, ah, a3, delta, amuzero, elev, &
         az, caz, azim, pi2, cos_delta, sin_delta, cos_ah, sin_xla, cos_xla, cos_elev

    integer :: i, j, nx, ny

    nx = size(xlon,dim=1)
    ny = size(xlon,dim=2)
    
    asol = -999.
    phis = -999.

    do j = 1, ny
    
      do i = 1, nx
      
      if (xlon(i,j) < -900) cycle

! jday is the number of the day in the month
!
!
!      mean solar time
       tsm = tu + xlon(i,j)/15.0
       xlo = xlon(i,j)*dtor
       xla = xlat(i,j)*dtor
       xj = real(jday)

!      time equation (mn.dec)
       a1 = (1.00554*xj - 6.28306)  * dtor
       a2 = (1.93946*xj + 23.35089) * dtor
       et = -7.67825*sin(a1) - 10.09176*sin(a2)

!      true solar time
       tsv = tsm + et/60.0
       tsv = tsv - 12.0

!      hour angle
       ah = tsv*15.0*dtor

!      solar declination (in radian)
       a3 = (0.9683*xj - 78.00878) * dtor
       delta = 23.4856*sin(a3)*dtor

!     elevation, azimuth
      cos_delta = cos(delta)
      sin_delta = sin(delta)
      cos_ah = cos(ah)
      sin_xla = sin(xla)
      cos_xla = cos(xla)
      
      amuzero = sin_xla*sin_delta + cos_xla*cos_delta*cos_ah
      elev = asin(amuzero)
      cos_elev = cos(elev)
      az = cos_delta*sin(ah)/cos_elev
      caz = (-cos_xla*sin_delta + sin_xla*cos_delta*cos_ah) / cos_elev

      if (az >= 1.0) then
         azim = asin(1.0)
      elseif (az <= -1.0) then
         azim = asin(-1.0)
      else
         azim = asin(az)
      endif

      if (caz <= 0.0) then
         azim = pi - azim
      endif

      if ((caz > 0.0) .and. (az <= 0.0)) then
         azim = 2 * pi + azim
      endif
      azim = azim + pi
      pi2 = 2 * pi
      if (azim > pi2) then
         azim = azim - pi2
      endif

!     conversion in degrees
      elev = elev / dtor
      asol(i,j) = 90.0 - elev
!     phis(i,j) = azim / dtor - 180.0
      phis(i,j) = azim / dtor  !akh - try to get 0 - 360

     enddo
   enddo
  return
end subroutine POSSOL 






subroutine goes_geo(nc_id,xscan,scan, lon, lat, goes_res, goes_ver)
    
  use h5lt
  use hdf5
  use h5ds
  use netcdf
  implicit none
  integer                                 ::  nc_id
  character(len=255)                      ::  dset_name,attr_name, str
  character(len=2), intent(in)            ::  goes_ver
  integer                                 ::  status, dset_id, nx, ny,pixels_per_dim
  integer, intent(in)                     ::  scan, xscan, goes_res
  integer, dimension(2)                   ::  sz
  real*8                                  ::  s_d, s_n,c1, c2, f, fp,H_GOESR, R_EQ,R_pol, d, s_1, s_2, s_3,s_xy
  real*8                                  ::  pi, dtor,sub_lon_radians,sub_lon_degrees,geographic_lon, geographic_lat
  real*8                                  ::  lon_degrees, lat_degrees, rtod,angle_step,max_angle
  real*8                                  ::  lambda0, c, a,b,rs,sx,sy,szg
  real, dimension(:,:), intent(inout)     ::  lon, lat
  real*8, dimension(:), allocatable         ::  x, y

  real    ::  scale_factor, add_offset

  status = -1
  allocate(x(xscan), y(scan), stat=status)  
  pi    =  3.14159265359
  dtor  = (pi/180.0)
  rtod  = (180.0/PI)
  dset_name = 'y'
  status = nf90_inq_varid(nc_id, dset_name, dset_id)
  str= "ERROR: Failed to get ID of dataset "//trim(dset_name)//": "
  call check(status,str)
    
  attr_name = 'scale_factor'
  status = nf90_get_att(nc_id, dset_id, attr_name, scale_factor)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)
    
  attr_name = 'add_offset'
  status = nf90_get_att(nc_id, dset_id, attr_name, add_offset)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)

  status = nf90_get_var(nc_id, dset_id, y)
  str= "ERROR: Failed to read dataset "//trim(dset_name)//": "
  call check(status,str)      
  y = y * scale_factor + add_offset  

  dset_name = 'x'
  status = nf90_inq_varid(nc_id, dset_name, dset_id)
  str= "ERROR: Failed to get ID of dataset "//trim(dset_name)//": "
  call check(status,str)
    
  attr_name = 'scale_factor'
  status = nf90_get_att(nc_id, dset_id, attr_name, scale_factor)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)
    
  attr_name = 'add_offset'
  status = nf90_get_att(nc_id, dset_id, attr_name, add_offset)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)

  status = nf90_get_var(nc_id, dset_id, x)
  str= "ERROR: Failed to read dataset "//trim(dset_name)//": "
  call check(status,str)   

  dset_name = 'goes_imager_projection'
  status = nf90_inq_varid(nc_id, dset_name, dset_id)
  str= "ERROR: Failed to get ID of dataset "//trim(dset_name)//": "
  call check(status,str)
    
  attr_name = 'longitude_of_projection_origin'
  status = nf90_get_att(nc_id, dset_id, attr_name, sub_lon_degrees)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)

  x = x * scale_factor + add_offset       !lamda = x  !theta = y
!   x = atan( tan(x)/cos(y) )
!   y = asin( sin(y)*cos(x) ) 

! longitude(xscan,scan)
  H_GOESR = 42164160      !orbit radius
  R_EQ    = 6378137.         !semi-major axis (equatorial radius km) */
  R_pol   = 6356752.31414
  if (goes_ver == '16') lambda0 = -1.308996939
  if (goes_ver == '17') lambda0 = -2.391101075  
  F       = (1.0/298.257222101) !flattening */
  lon(:,:)= -999.0
  lat(:,:)= -999.0
  sz  = shape(lon)
  c   = H_GOESR**2. - R_EQ**2.
    
  FP      = (1.0/((1.0-F)*(1.0-F)))
  d       = (H_GOESR*H_GOESR - R_EQ*R_EQ)  
!   sub_lon_radians = sub_lon_degrees * dtor
!   yadd  = 0.151844  
!   yscal = -0.000056

  do nx  = 1, sz(1)
  do ny  = 1, sz(2)
!     c1      = (H_GOESR * cos(x(nx)) * cos(y(ny))) * (H_GOESR * cos(x(nx)) * cos(y(ny)))  
!     c2      = (cos(y(ny)) * cos(y(ny)) + FP * sin(y(ny)) * sin(y(ny))) * d  
!     if (c1>c2) then
    a   = (dsin(x(nx))**2.)+(dcos(x(nx))**2.)*((dcos(y(ny))**2.)+ (R_EQ/R_pol)**2.*(dsin(y(ny))**2.)  )
    b   = -2.* H_GOESR*dcos(x(nx))*dcos(y(ny))
    rs  = (-b-sqrt(b**2.-4.*a*c)) / (2.*a)
    sx  = rs*dcos(x(nx))*dcos(y(ny))
    sy  = -rs*dsin(x(nx))
    szg  = rs*dcos(x(nx))*dsin(y(ny))
    lat(nx, ny)=atan(  ((R_EQ/R_pol)**2.) *  szg  /  sqrt((H_GOESR -sx)**2. + sy**2. )    )*rtod 
    lon(nx, ny)=(lambda0  - datan(sy/(H_GOESR -sx)))*rtod
    if (lon(nx, ny) > 180.) lon(nx, ny)= lon(nx, ny)-360.
    if (lon(nx, ny) < -180.) lon(nx, ny)= lon(nx, ny)+360.
    if (isnan(lat(nx, ny))) then
      lat(nx, ny)=-999.0
      lon(nx, ny)=-999.0
    endif 

!     
!       s_d = sqrt(c1 - c2)
!       s_n = (H_GOESR * cos(x(nx)) * cos(y(ny)) - s_d) / (cos(y(ny)) * cos(y(ny)) + FP * sin(y(ny)) * sin(y(ny)))
!       s_1 = H_GOESR - s_n * cos(x(nx)) * cos(y(ny))
!       s_2 = s_n * sin(x(nx)) * cos(y(ny))
!       s_3 = (-1.) * s_n * sin(y(ny))
!       s_xy = sqrt(s_1*s_1 + s_2*s_2)    
!       geographic_lon = atan(s_2/s_1) + sub_lon_radians
!       geographic_lat = atan((-1)*FP*(s_3/s_xy))
!       lon_degrees = (rtod*geographic_lon)
!       lat_degrees = rtod*geographic_lat
!       
!       if (lon_degrees < -180.0) lon_degrees = lon_degrees+360. 
!       if (lon_degrees > 180.0)  lon_degrees = lon_degrees-360. 
!       lon(nx, ny)=lon_degrees
!       lat(nx, ny)=lat_degrees
!     end if   
        
  end do !ny
  end do !nx
  

  return
end subroutine goes_geo


subroutine goes_read_rad(nc_id, dset_in,itmp, dummy,do_bt)
    
  use h5lt
  use hdf5
  use h5ds
  use netcdf
  implicit none
  integer                                 ::  nc_id
  character(len=255) , intent(in)         ::  dset_in
  character(len=255)                      ::  dset_name,attr_name,str
  real, dimension(:,:), intent(inout)     ::  dummy
  integer, dimension(:,:), intent(inout)  ::  itmp
  integer, intent(in)                     ::  do_bt 
  integer                                 ::  status, dset_id
  real    ::  scale_factor, add_offset, pi, dist, esun, bc1, bc2, fk1, fk2
  integer ::  fill_value  

  pi =  3.14159 
  status = -1
    
  status = nf90_inq_varid(nc_id, dset_in, dset_id)
  str= "ERROR: Failed to get ID of dataset "//trim(dset_name)//": "
  call check(status,str)
    
  attr_name = 'scale_factor'
  status = nf90_get_att(nc_id, dset_id, attr_name, scale_factor)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)
    
  attr_name = 'add_offset'
  status = nf90_get_att(nc_id, dset_id, attr_name, add_offset)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)

  attr_name = '_FillValue'
  status = nf90_get_att(nc_id, dset_id, attr_name, fill_value)
  str= "ERROR: Failed to read attribute "//trim(attr_name)//": "
  call check(status,str)

  status = nf90_get_var(nc_id, dset_id, itmp)
  str= "ERROR: Failed to read dataset "//trim(dset_name)//": "
  call check(status,str)  

  dset_name = 'esun'
  status = nf90_inq_varid(nc_id, dset_name, dset_id)
  str= "ERROR: Failed to get ID of dataset "//trim(dset_name)//": "
  call check(status,str)  
  
  status = nf90_get_var(nc_id, dset_id, esun)
  str= "ERROR: Failed to read dataset "//trim(dset_name)//": "
  call check(status,str)  

  dset_name = 'earth_sun_distance_anomaly_in_AU'
  status = nf90_inq_varid(nc_id, dset_name, dset_id)
  str= "ERROR: Failed to get ID of dataset "//trim(dset_name)//": "
  call check(status,str)  

  status = nf90_get_var(nc_id, dset_id, dist)
  str= "ERROR: Failed to read dataset "//trim(dset_name)//": "
  call check(status,str) 

  dset_name = 'planck_bc1'
  status = nf90_inq_varid(nc_id, dset_name, dset_id)
  str= "ERROR: Failed to get ID of dataset "//trim(dset_name)//": "
  call check(status,str)  
  
  status = nf90_get_var(nc_id, dset_id, bc1)
  str= "ERROR: Failed to read dataset "//trim(dset_name)//": "
  call check(status,str)  

  dset_name = 'planck_bc2'
  status = nf90_inq_varid(nc_id, dset_name, dset_id)
  str= "ERROR: Failed to get ID of dataset "//trim(dset_name)//": "
  call check(status,str)  
  
  status = nf90_get_var(nc_id, dset_id, bc2)
  str= "ERROR: Failed to read dataset "//trim(dset_name)//": "
  call check(status,str)    

  dset_name = 'planck_fk1'
  status = nf90_inq_varid(nc_id, dset_name, dset_id)
  str= "ERROR: Failed to get ID of dataset "//trim(dset_name)//": "
  call check(status,str)  
  
  status = nf90_get_var(nc_id, dset_id, fk1)
  str= "ERROR: Failed to read dataset "//trim(dset_name)//": "
  call check(status,str)     

  dset_name = 'planck_fk2'
  status = nf90_inq_varid(nc_id, dset_name, dset_id)
  str= "ERROR: Failed to get ID of dataset "//trim(dset_name)//": "
  call check(status,str)  
  
  status = nf90_get_var(nc_id, dset_id, fk2)
  str= "ERROR: Failed to read dataset "//trim(dset_name)//": "
  call check(status,str)    

  dummy(:,:) = -999.0
  if (do_bt == 0) then 
    where(itmp /= fill_value)
      dummy = ((itmp * scale_factor + add_offset)*(dist**2)*pi)/esun
    end where
  end if
  if (do_bt == 1) then 
    where(itmp /= fill_value)     
      dummy = (fk2/log((fk1/((itmp * scale_factor + add_offset)*dist**2.)+1))-bc1)/bc2
    end where
  end if
  return
end subroutine goes_read_rad




! type (viirs_db_svm) function load_viirs_db_data(gmtco_file, svm01_file, svm02_file, &
! &         svm03_file, svm04_file, svm05_file, svm07_file, svm08_file, svm09_file,   &
! &         svm10_file, svm11_file, svm14_file, svm15_file, svm16_file, iicmo_file,   &
! &         status) result (vdbs)
!   
!   use h5lt
!   use hdf5
!   
!   use viirs_obpg_corrections
!   
!   use calendars, only:  gdatetime,         &
!                         doy_from_gregorian
!   
!   implicit none
! 
!   character(len=*), intent(in)          ::  gmtco_file
!   character(len=*), intent(in)          ::  svm01_file
!   character(len=*), intent(in)          ::  svm02_file
!   character(len=*), intent(in)          ::  svm03_file
!   character(len=*), intent(in)          ::  svm04_file
!   character(len=*), intent(in)          ::  svm05_file
!   character(len=*), intent(in)          ::  svm07_file
!   character(len=*), intent(in)          ::  svm08_file
!   character(len=*), intent(in)          ::  svm09_file
!   character(len=*), intent(in)          ::  svm10_file
!   character(len=*), intent(in)          ::  svm11_file
!   character(len=*), intent(in)          ::  svm14_file
!   character(len=*), intent(in)          ::  svm15_file
!   character(len=*), intent(in)          ::  svm16_file
!   character(len=*), intent(in)          ::  iicmo_file
!   integer, intent(inout)                ::  status
!   
!   integer(hid_t)                        ::  h5f_id          ! HDF5 file id
!   integer(hid_t)                        ::  dset_id         ! HDF5 dataset id
!   integer(hid_t)                        ::  dspace_id       ! HDF5 dataspace id
!   integer(hid_t)                        ::  attr_id         ! HDF5 attribute id
!   integer(hsize_t), dimension(1)        ::  dims1
!   integer(hsize_t), dimension(2)        ::  dims2
!   integer(hsize_t), dimension(3)        ::  max_dims
!   integer                               ::  dclass
!   integer                               ::  dtype!integer(hid_t)->integer,wvkim
!   integer(hid_t)                        ::  mtype
!   integer(size_t)                       ::  tsize
!   
!   character(len=255)                    ::  dset_name
!   character(len=255)                    ::  attr_name
!   character(len=5)                      ::  day_or_night
!   character(len=8)                      ::  start_date, end_date
!   character(len=14)                     ::  start_time
!   character(len=255)                    ::  day
!   character(len=255)                    ::  month
!   character(len=255)                    ::  year
!   character(len=255)                    ::  doy
! 
!   
!   real, dimension(:,:), allocatable     ::  vaa
!   integer, dimension(:,:), allocatable  ::  i_refl
!   integer, dimension(:,:), allocatable  ::  i_rad
!   integer, dimension(:,:), allocatable  ::  i_bt
!   integer, dimension(:,:), allocatable  ::  cld_mask
!   integer, dimension(:,:), allocatable  ::  cld_mask_qa
!   integer, dimension(:,:), allocatable  ::  lw_mask
!   real, dimension(2)                    ::  refl_factors
!   real, dimension(2)                    ::  rad_factors
!   real, dimension(2)                    ::  bt_factors
!   integer                               ::  scan, xscan
!   integer, dimension(1)                 ::  tmp_nmirror
!   integer                               ::  nmirror
!   byte                                  ::  bit_mask
!   
!   real, dimension(:,:), allocatable     ::  sc_position
!   real, dimension(:,:), allocatable     ::  sc_velocity
!   real, dimension(:,:), allocatable     ::  sc_attitude
!   real, dimension(:,:), allocatable     ::  sen_mat
!   real, dimension(:), allocatable       ::  coeffs
!   
!   type(gdatetime)                       ::  gdt1
!   real, parameter                       ::  d2r =  3.14159/180.0
!   
!   integer                               ::  i, j
!   real                                  ::  tmp0, tmp1
!   
!   type(c_ptr)                           :: f_ptr
!   integer(kind=i64), dimension(:), allocatable, target :: tmp_scan_time    ! tmp buffer
! 
!   
!   status = -1
!   
!   call h5open_f(status)
!   if (status < 0) then
!     print *, "ERROR: Unable to start HDF interface: ", status
!     return
!   end if
!   
! ! -- open and read geolocation and geometry information
!   call h5fopen_f(trim(gmtco_file), H5F_ACC_RDONLY_F, h5f_id, status)
!   if (status < 0) then 
!     print *, "ERROR: Unable to open VIIRS geolocation file: ", status
!     print *, "File: ", trim(gmtco_file)
!     return
!   end if
!   
! ! -- check that granule is a day granule
!   dset_name = "/Data_Products/VIIRS-MOD-GEO-TC/VIIRS-MOD-GEO-TC_Gran_0"
!   attr_name = "N_Day_Night_Flag"
!   call h5ltget_attribute_string_f(h5f_id, dset_name, attr_name, &
!                                       day_or_night, status)
!   if (status /= 0) then
!     print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
!     return
!   end if
!                                   
!   if (day_or_night == "Night") then
!     print *, "ERROR: Granule is a night granule. Skipping."
!     status = -1
!     return
!   end if
! 
! ! -- get start date attribute
!   dset_name = "/Data_Products/VIIRS-MOD-GEO-TC/VIIRS-MOD-GEO-TC_Gran_0"
!   attr_name = "Beginning_Date"
!   call h5aopen_by_name_f(h5f_id, dset_name, attr_name, attr_id, status)
!   if (status < 0) then
!     print *, "ERROR: Unable to access "//trim(attr_name)//" attribute: ", status
!     return
!   end if
!   
!   call h5aget_type_f(attr_id, mtype, status)
!   if (status /= 0) then
!     print *, "ERROR: Failed to get type for attribute "//trim(attr_name)//": ", status
!     return
!   end if
!   
!   dims2 = (/0,0/)   ! ignored in call to h5aread_f below.
!   call h5aread_f(attr_id, mtype, start_date, dims2, status)
!   if (status /= 0) then
!     print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
!     return
!   end if
! 
!   call h5aclose_f(attr_id, status)
!   if (status /= 0) then
!     print *, "ERROR: Failed to close "//trim(attr_name)//":", status 
!     return
!   end if
!   
!   read(start_date(1:4), fmt='(I4)') vdbs%yr
!   read(start_date(5:6), fmt='(I4)') vdbs%mo
!   read(start_date(7:8), fmt='(I4)') vdbs%dy
!   
!   gdt1 = gdatetime(vdbs%yr, vdbs%mo, vdbs%dy, 0, 0, 0, 0, 0) 
!   vdbs%doy = doy_from_gregorian(gdt1)
!   
!   attr_name = "Beginning_Time"
!   call h5aopen_by_name_f(h5f_id, dset_name, attr_name, attr_id, status)
!   if (status < 0) then
!     print *, "ERROR: Unable to access "//trim(attr_name)//" attribute: ", status
!     return
!   end if
! 
!   call h5aget_type_f(attr_id, mtype, status)
!   if (status /= 0) then
!     print *, "ERROR: Failed to get type for attribute "//trim(attr_name)//": ", status
!     return
!   end if
!   
!   dims2 = (/0,0/)   ! ignored in call to h5aread_f below.
!   call h5aread_f(attr_id, mtype, start_time, dims2, status) 
!   if (status /= 0) then
!     print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
!     return
!   end if
!   
!   call h5aclose_f(attr_id, status)
!   if (status /= 0) then
!     print *, "ERROR: Failed to close "//trim(attr_name)//":", status 
!     return
!   end if
!   
!   read(start_time(1:2), fmt='(I2)') vdbs%hr
!   read(start_time(3:4), fmt='(I2)') vdbs%min
!   
! ! -- get dimensions of latitude data set, save, and allocate other data arrays.
!   dset_name = "/All_Data/VIIRS-MOD-GEO-TC_All/Latitude"
!   call h5ltget_dataset_info_f(h5f_id, trim(dset_name), dims2, dtype, tsize, status)
!   if (status < 0) then
!     print *, "ERROR: Failed to get dimensions of data set: ", status
!     print *, "Data set: ", trim(dset_name)
!     return
!   end if
! 
!    xscan       = dims2(1)
!    scan        = dims2(2)
!    vdbs%scan   = scan
!    vdbs%xscan  = xscan
! 
!   allocate(vdbs%scan_time(scan), stat=status)
!   if (status /= 0) then
!     print *, "ERROR: Failed to allocate scan time arrays: ", status
!     return
!   end if
! 
!   allocate(vdbs%lat(xscan,scan), vdbs%lon(xscan,scan), vdbs%sza(xscan,scan),  &
!   &       vdbs%saa(xscan,scan), vdbs%vza(xscan,scan), vdbs%vaa(xscan,scan),   &
!   &       vdbs%raa(xscan,scan), vdbs%sca(xscan,scan), vdbs%gla(xscan,scan), stat=status)
!   if (status /= 0) then
!     print *, "ERROR: Failed to allocate geolocation and/or geometry data arrays: ", status
!     return
!   end if
! 
!   allocate(vdbs%m01_refl(xscan,scan), vdbs%m02_refl(xscan,scan), &
!   &         vdbs%m03_refl(xscan,scan), vdbs%m04_refl(xscan,scan),&
!   &         vdbs%m05_refl(xscan,scan), i_refl(xscan,scan), stat=status)
!   if (status /= 0) then
!     print *, "ERROR: Failed to allocate reflectance arrays: ", status
!     return
!   end if
!   
!   allocate(vdbs%m07_refl(xscan,scan), vdbs%m08_refl(xscan,scan),  &
!   &         vdbs%m09_refl(xscan,scan), vdbs%m10_refl(xscan,scan), &
!   &         vdbs%m11_refl(xscan,scan), stat=status)
!   if (status /= 0) then
!     print *, "ERROR: Failed to allocate reflectance arrays: ", status
!     return
!   end if
! 
!   allocate(vdbs%m14_rad(xscan,scan), vdbs%m15_rad(xscan,scan),  &
!   &         vdbs%m16_rad(xscan,scan), i_rad(xscan,scan), stat=status)
!   if (status /= 0) then
!     print *, "ERROR: Failed to allocate radiance arrays: ", status
!     return
!   end if
!   
!   allocate(vdbs%m14_bt(xscan, scan), vdbs%m15_bt(xscan, scan),  &
!   & vdbs%m16_bt(xscan, scan), i_bt(xscan,scan), stat=status)
!   if (status /= 0) then
!     print *, "ERROR: Failed to allocate brightness temperature arrays: ", status
!     return
!   end if
! 
!   allocate(vdbs%land_mask(xscan,scan), vdbs%elev(xscan,scan), vdbs%ps(xscan,scan), &
!   &   stat=status)
!   if (status /= 0) then
!     print *, "ERROR: Failed to allocate mask arrays: ", status
!     return
!   end if
!   
!   allocate(vdbs%btd8(xscan,scan), vdbs%btd11(xscan,scan), stat=status)
!   if (status /= 0) then
!     print *, "ERROR: Failed to allocate BTD arrays: ", status
!     return
!   end if
!   
!   allocate(vdbs%dstar(xscan,scan), vdbs%ndvi(xscan, scan), stat=status)
!   if (status /= 0) then
!     print *, "ERROR: Failed to allocate D* arrays: ", status
!     return
!   end if
!   
!   allocate(cld_mask_qa(xscan,scan), &
!   &  cld_mask(xscan,scan), lw_mask(xscan,scan), stat=status)
!   if (status /= 0) then 
!     print *, "ERROR: Unable to allocate azimuth angle arrays: ", status
!     return
!   end if
! 
!   allocate(latitude(xscan,scan), longitude(xscan,scan), stat=status)
!   if (status /= 0) then
!     print *, "ERROR: Failed to allocate shared geolocation arrays: ", status
!     return
!   end if
!    
! ! -- read in quality flag dataset and pluck out mirror side information
! ! --  see documentation here, p. 151:
! ! --  http://npp.gsfc.nasa.gov/sciencedocs/2015-08/474-00448-02-06_JPSS-DD-Vol-II-Part-6_0200D.pdf
!   dset_name = "/All_Data/VIIRS-MOD-GEO-TC_All/NumberOfScans"
!   dims1(1)  = 1
!   call h5ltread_dataset_int_f(h5f_id, trim(dset_name), tmp_nmirror, dims1, status)
!   if (status < 0) then 
!     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
!     return
!   end if
!   nmirror = tmp_nmirror(1)
!   nmirror = 48 
!   
!   allocate(vdbs%mirror_side(nmirror), stat=status)
!   if (status /= 0) then
!     print *, "ERROR: Failed to allocate mirror side arrays: ", status
!     return
!   end if
!   
!   dset_name = "/All_Data/VIIRS-MOD-GEO-TC_All/QF1_SCAN_VIIRSSDRGEO"
!   call h5ltread_dataset_int_f(h5f_id, trim(dset_name), vdbs%mirror_side, dims1, status)
!   if (status < 0) then
!     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
!     return
!   end if
!   
! ! -- shift bits by 7 to the right and 'AND' with 00000001.
! ! -- mirror side info is last bit in byte.
!   bit_mask = b'1000000'
!   vdbs%mirror_side = iand(vdbs%mirror_side, bit_mask)  
! 
! ! -- RELATED TO POLARIZATION CORRECTION, CURRENTLY NOT NEEDED.
! !  allocate(vdbs%mnorm(3, nscan), stat=status)
! !  if (status /= 0) then
! !    print *, "ERROR: Failed to allocate mirror normal vector: ", status
! !    return
! !  end if
!  
! ! allocate(sc_position(3,nscan), sc_attitude(3,nscan), &
! !   & sc_velocity(3,nscan), stat=status)
! !   if (status /= 0) then
! !     print *, "ERROR: Failed to allocate spacecraft orientation arrays: ", status
! !     return
! !   end if
! ! 
! !   dset_name = "/All_Data/VIIRS-MOD-GEO-TC_All/SCPosition"
! !   dims2     = (/3,nscan/)
! !   call h5ltread_dataset_float_f(h5f_id, trim(dset_name), sc_position, dims2, status)
! !   if (status < 0) then 
! !     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
! !     return
! !   end if
! !   sc_position = sc_position / 1000.0  ! convert from m to km.
!   
! !   dset_name = "/All_Data/VIIRS-MOD-GEO-TC_All/SCAttitude"
! !   dims2     = (/3,48/)
! !   call h5ltread_dataset_float_f(h5f_id, trim(dset_name), sc_attitude, dims2, status)
! !   if (status < 0) then 
! !     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
! !     return
! !   end if
! !   sc_attitude = sc_attitude / 3600.0  ! convert arcsec to degrees.  
! !   
! !   dset_name = "/All_Data/VIIRS-MOD-GEO-TC_All/SCVelocity"
! !   dims2     = (/3,48/)
! !   call h5ltread_dataset_float_f(h5f_id, trim(dset_name), sc_velocity, dims2, status)
! !   if (status < 0) then
! !     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
! !     return
! !   end if
! !   sc_velocity = sc_velocity / 1000.0  ! convert from m to km.
! 
! ! -- read scan times. HDF5 high-level API doesn't have a function to read
! ! -- long ints (or anything larger than ints), so we have to use the full library.
!   dset_name = "/All_Data/VIIRS-MOD-GEO-TC_All/MidTime"
!   call h5dopen_f (h5f_id, dset_name, dset_id, status)
!   if (status /= 0) then
!     print *, "ERROR: Failed to open dataset "//trim(dset_name)//": ", status
!     return
!   end if
!   
!   call h5dget_space_f(dset_id, dspace_id, status) 
!   if (status /= 0) then
!     print *, "ERROR: Failed to get simple dataspace for "//trim(dset_name)//": ", status
!     return
!   end if
!   
!   call h5sget_simple_extent_dims_f(dspace_id, dims1, max_dims, status)
!   if (status < 0) then    ! NOTE: status contains rank of dspace_id on exit.
!     print *, "ERROR: Failed to get dataspace dims for "//trim(dset_name)//": ", status
!     return
!   end if
!   
!   allocate(tmp_scan_time(dims1(1)), stat=status)
!   if (status /= 0) then
!     print *, "ERROR: Failed to allocate tmp array for scan times: ", status
!     return
!   end if
!   
!   f_ptr = c_loc(tmp_scan_time(1))   
!   call h5dread_f(dset_id, h5kind_to_type(i64,H5_INTEGER_KIND), f_ptr, status)
!   if (status /= 0) then
!     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
!     return
!   end if
!   
!   call h5sclose_f(dspace_id, status)
!   if (status /= 0) then
!     print *, "ERROR: Failed to close dataspace on "//trim(dset_name)//": ", status
!     return
!   end if
!   
!   call h5dclose_f(dset_id, status)
!   if (status /= 0) then
!     print *, "ERROR: Failed to close dataset "//trim(dset_name)//": ", status
!     return
!   end if
! 
! ! -- each "scan" technically is 16 pixels wide. Expand times
! ! -- to cover full granule.
!   do i = 1, scan
!     j = (i-1)/16 + 1    ! integer division
!     vdbs%scan_time(i) = tmp_scan_time(j)
!   end do
!     
!   deallocate(tmp_scan_time, stat=status)
!   if (status /= 0) then
!     print *, "ERROR: Failed to deallocate tmp array for time data: ", status
!     return
!   end if
!     
! ! -- read all geolocation and geometry data sets
!   dset_name = "/All_Data/VIIRS-MOD-GEO-TC_All/Latitude"
!   call h5ltread_dataset_float_f(h5f_id, trim(dset_name), vdbs%lat, dims2, status)
!   if (status < 0) then 
!     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
!     return
!   end if
! 
!   latitude = vdbs%lat
!   
!   dset_name = "/All_Data/VIIRS-MOD-GEO-TC_All/Longitude"
!   call h5ltread_dataset_float_f(h5f_id, trim(dset_name), vdbs%lon, dims2, status)
!   if (status < 0) then 
!     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
!     return
!   end if
!   
!   longitude = vdbs%lon
!   
!   dset_name = "/All_Data/VIIRS-MOD-GEO-TC_All/SolarZenithAngle"
!   call h5ltread_dataset_float_f(h5f_id, trim(dset_name), vdbs%sza, dims2, status)
!   if (status < 0) then 
!     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
!     return
!   end if
!   
!   dset_name = "/All_Data/VIIRS-MOD-GEO-TC_All/SatelliteZenithAngle"
!   call h5ltread_dataset_float_f(h5f_id, trim(dset_name), vdbs%vza, dims2, status)
!   if (status < 0) then 
!     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
!     return
!   end if
!   
!   dset_name = "/All_Data/VIIRS-MOD-GEO-TC_All/SolarZenithAngle"
!   call h5ltread_dataset_float_f(h5f_id, trim(dset_name), vdbs%sza, dims2, status)
!   if (status < 0) then 
!     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
!     return
!   end if
! 
! ! -- read azimuth angle data and calculate relative azimuth angle and scattering angle.
!   dset_name = "/All_Data/VIIRS-MOD-GEO-TC_All/SolarAzimuthAngle"
!   call h5ltread_dataset_float_f(h5f_id, trim(dset_name), vdbs%saa, dims2, status)
!   if (status < 0) then
!     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
!     return
!   end if
!   
!   dset_name = "/All_Data/VIIRS-MOD-GEO-TC_All/SatelliteAzimuthAngle"
!   call h5ltread_dataset_float_f(h5f_id, trim(dset_name), vdbs%vaa, dims2, status)
!   if (status < 0) then 
!     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
!     return
!   end if
! 
!   vdbs%raa  = (vdbs%vaa - vdbs%saa) - 180.0
!   where (vdbs%raa > 180.0)  vdbs%raa = vdbs%raa - 360.0
!   where (vdbs%raa < -180.0) vdbs%raa = vdbs%raa + 360.0
!   where (vdbs%raa < 0.0)    vdbs%raa = -1.0 * vdbs%raa
!     
!   vdbs%sca  = acos(cos(vdbs%sza*d2r)*cos(vdbs%vza*d2r) -       &
!                 sin(vdbs%sza*d2r)*sin(vdbs%vza*d2r)*cos(vdbs%raa*d2r))
!   vdbs%sca  = 180.0 - (vdbs%sca/d2r)
!   
!   vdbs%gla  = acos(cos(vdbs%sza*d2r)*cos(vdbs%vza*d2r) +       &
!                 sin(vdbs%sza*d2r)*sin(vdbs%vza*d2r)*cos(vdbs%raa*d2r))
!   vdbs%gla = vdbs%gla/d2r
! 
!   dset_name = "/All_Data/VIIRS-MOD-GEO-TC_All/Height"
!   call h5ltread_dataset_float_f(h5f_id, trim(dset_name), vdbs%elev, dims2, status)
!   if (status < 0) then 
!     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
!     return
!   end if
! 
! ! -- convert surface elevation to surface pressure
! ! -- output of atmosphere is fraction of mean sea-level pressure.
!   do j = 1, scan
!     do i = 1, xscan
!       call atmosphere(vdbs%elev(i,j)/1000.0, tmp0, vdbs%ps(i,j), tmp1)
!       vdbs%ps(i,j) = vdbs%ps(i,j) * 1013.25
!     end do
!   end do
!   
!   call h5fclose_f(h5f_id, status)
!   if (status < 0) then
!     print *, "ERROR: Failed to close GMTCO file: ", status
!     return
!   end if
!     
!   call h5fopen_f(trim(svm01_file), H5F_ACC_RDONLY_F, h5f_id, status)
!   if (status < 0) then 
!     print *, "ERROR: Unable to open VIIRS SMV01 file: ", status
!     print *, "File: ", trim(svm01_file)
!     return
!   end if
!   
!   dset_name = "/All_Data/VIIRS-M1-SDR_All/ReflectanceFactors"
!   dims1(1) = 2
!   call h5ltread_dataset_float_f(h5f_id, trim(dset_name), refl_factors, dims1, status)
!   if (status < 0) then 
!     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
!     return
!   end if
!   
!   dset_name = "/All_Data/VIIRS-M1-SDR_All/Reflectance"
!   call h5ltread_dataset_int_f(h5f_id, trim(dset_name), i_refl, dims2, status)
!   if (status < 0) then
!     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
!     return
!   end if
!   
!   vdbs%m01_refl(:,:) = -999.0
!   where(i_refl < 65528) vdbs%m01_refl = i_refl * refl_factors(1)
!     
!   call h5fopen_f(trim(svm02_file), H5F_ACC_RDONLY_F, h5f_id, status)
!   if (status < 0) then 
!     print *, "ERROR: Unable to open VIIRS SMV02 file: ", status
!     print *, "File: ", trim(svm02_file)
!     return
!   end if
!   
!   dset_name = "/All_Data/VIIRS-M2-SDR_All/ReflectanceFactors"
!   dims1(1) = 2
!   call h5ltread_dataset_float_f(h5f_id, trim(dset_name), refl_factors, dims1, status)
!   if (status < 0) then 
!     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
!     return
!   end if
!   
!   dset_name = "/All_Data/VIIRS-M2-SDR_All/Reflectance"
!   call h5ltread_dataset_int_f(h5f_id, trim(dset_name), i_refl, dims2, status)
!   if (status < 0) then 
!     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
!     return
!   end if
!   
!   vdbs%m02_refl(:,:) = -999.0
!   where(i_refl < 65528) vdbs%m02_refl = i_refl * refl_factors(1)
! 
!   call h5fclose_f(h5f_id, status)
!   if (status < 0) then
!     print *, "ERROR: Failed to close SVM02 file: ", status
!     return
!   end if
! 
!   call h5fopen_f(trim(svm03_file), H5F_ACC_RDONLY_F, h5f_id, status)
!   if (status < 0) then 
!     print *, "ERROR: Unable to open VIIRS SMV03 file: ", status
!     print *, "File: ", trim(svm03_file)
!     return
!   end if
!   
!   dset_name = "/All_Data/VIIRS-M3-SDR_All/ReflectanceFactors"
!   dims1(1) = 2
!   refl_factors(:) = 0.0
!   call h5ltread_dataset_float_f(h5f_id, trim(dset_name), refl_factors, dims1, status)
!   if (status < 0) then 
!     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
!     return
!   end if
!   
!   dset_name = "/All_Data/VIIRS-M3-SDR_All/Reflectance"
!   i_refl(:,:) = -9999
!   call h5ltread_dataset_int_f(h5f_id, trim(dset_name), i_refl, dims2, status)
!   if (status < 0) then 
!     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
!     return
!   end if
!   
!   vdbs%m03_refl(:,:) = -999.0
!   where(i_refl < 65528) vdbs%m03_refl = i_refl * refl_factors(1)
! 
!   call h5fclose_f(h5f_id, status)
!   if (status < 0) then
!     print *, "ERROR: Failed to close SVM03 file: ", status
!     return
!   end if
!   
!   call h5fopen_f(trim(svm04_file), H5F_ACC_RDONLY_F, h5f_id, status)
!   if (status < 0) then 
!     print *, "ERROR: Unable to open VIIRS SMV04 file: ", status
!     print *, "File: ", trim(svm04_file)
!     return
!   end if
!   
!   dset_name = "/All_Data/VIIRS-M4-SDR_All/ReflectanceFactors"
!   dims1(1) = 2
!   refl_factors(:) = 0.0
!   call h5ltread_dataset_float_f(h5f_id, trim(dset_name), refl_factors, dims1, status)
!   if (status < 0) then
!     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
!     return
!   end if
!   
!   dset_name = "/All_Data/VIIRS-M4-SDR_All/Reflectance"
!   i_refl(:,:) = -9999
!   call h5ltread_dataset_int_f(h5f_id, trim(dset_name), i_refl, dims2, status)
!   if (status < 0) then 
!     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
!     return
!   end if
!   
!   vdbs%m04_refl(:,:) = -999.0
!   where(i_refl < 65528) vdbs%m04_refl = i_refl * refl_factors(1)
!     
!   call h5fclose_f(h5f_id, status)
!   if (status < 0) then
!     print *, "ERROR: Failed to close SVM04 file: ", status
!     return
!   end if
!   
!   call h5fopen_f(trim(svm05_file), H5F_ACC_RDONLY_F, h5f_id, status)
!   if (status < 0) then
!     print *, "ERROR: Unable to open VIIRS SMV05 file: ", status
!     print *, "File: ", trim(svm05_file)
!     return
!   end if
!   
!   dset_name = "/All_Data/VIIRS-M5-SDR_All/ReflectanceFactors"
!   dims1(1) = 2
!   refl_factors(:) = 0.0
!   call h5ltread_dataset_float_f(h5f_id, trim(dset_name), refl_factors, dims1, status)
!   if (status < 0) then 
!     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
!     return
!   end if
!   
!   dset_name = "/All_Data/VIIRS-M5-SDR_All/Reflectance"
!   i_refl(:,:) = -9999
!   call h5ltread_dataset_int_f(h5f_id, trim(dset_name), i_refl, dims2, status)
!   if (status < 0) then 
!     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
!     return
!   end if
!   
!   vdbs%m05_refl(:,:) = -999.0
!   where(i_refl < 65528) vdbs%m05_refl = i_refl * refl_factors(1)
!   
!   call h5fclose_f(h5f_id, status)
!   if (status < 0) then
!     print *, "ERROR: Failed to close SVM05 file: ", status
!     return
!   end if
!   
!   call h5fopen_f(trim(svm07_file), H5F_ACC_RDONLY_F, h5f_id, status)
!   if (status < 0) then 
!     print *, "ERROR: Unable to open VIIRS SMV07 file: ", status
!     print *, "File: ", trim(svm07_file)
!     return
!   end if
!   
!   dset_name = "/All_Data/VIIRS-M7-SDR_All/ReflectanceFactors"
!   dims1(1) = 2
!   call h5ltread_dataset_float_f(h5f_id, trim(dset_name), refl_factors, dims1, status)
!   if (status < 0) then 
!     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
!     return
!   end if
!   
!   dset_name = "/All_Data/VIIRS-M7-SDR_All/Reflectance"
!   call h5ltread_dataset_int_f(h5f_id, trim(dset_name), i_refl, dims2, status)
!   if (status < 0) then 
!     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
!     return
!   end if
! 
!   vdbs%m07_refl(:,:) = -999.0
!   where(i_refl < 65528) vdbs%m07_refl = i_refl * refl_factors(1)
!     
!   call h5fclose_f(h5f_id, status)
!   if (status < 0) then
!     print *, "ERROR: Failed to close SVM07 file: ", status
!     return
!   end if
!   
!   call h5fopen_f(trim(svm08_file), H5F_ACC_RDONLY_F, h5f_id, status)
!   if (status < 0) then 
!     print *, "ERROR: Unable to open VIIRS SMV01 file: ", status
!     print *, "File: ", trim(svm08_file)
!     return
!   end if
!   
!   dset_name = "/All_Data/VIIRS-M8-SDR_All/ReflectanceFactors"
!   dims1(1) = 2
!   call h5ltread_dataset_float_f(h5f_id, trim(dset_name), refl_factors, dims1, status)
!   if (status < 0) then 
!     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
!     return
!   end if
!   
!   dset_name = "/All_Data/VIIRS-M8-SDR_All/Reflectance"
!   call h5ltread_dataset_int_f(h5f_id, trim(dset_name), i_refl, dims2, status)
!   if (status < 0) then 
!     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
!     return
!   end if
!   
!   vdbs%m08_refl(:,:) = -999.0
!   where(i_refl < 65528) vdbs%m08_refl = i_refl * refl_factors(1)
!     
!   call h5fclose_f(h5f_id, status)
!   if (status < 0) then
!     print *, "ERROR: Failed to close SVM08 file: ", status
!     return
!   end if
!   
!   call h5fopen_f(trim(svm09_file), H5F_ACC_RDONLY_F, h5f_id, status)
!   if (status < 0) then 
!     print *, "ERROR: Unable to open VIIRS SMV10 file: ", status
!     print *, "File: ", trim(svm10_file)
!     return
!   end if
!   
!   dset_name = "/All_Data/VIIRS-M9-SDR_All/ReflectanceFactors"
!   dims1(1) = 2
!   call h5ltread_dataset_float_f(h5f_id, trim(dset_name), refl_factors, dims1, status)
!   if (status < 0) then 
!     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
!     return
!   end if
!   
!   dset_name = "/All_Data/VIIRS-M9-SDR_All/Reflectance"
!   call h5ltread_dataset_int_f(h5f_id, trim(dset_name), i_refl, dims2, status)
!   if (status < 0) then 
!     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
!     return
!   end if
!   
!   vdbs%m09_refl(:,:) = -999.0
!   where(i_refl < 65528) vdbs%m09_refl = i_refl * refl_factors(1)
!     
!   call h5fclose_f(h5f_id, status)
!   if (status < 0) then
!     print *, "ERROR: Failed to close SVM10 file: ", status
!     return
!   end if
!   
!   call h5fopen_f(trim(svm10_file), H5F_ACC_RDONLY_F, h5f_id, status)
!   if (status < 0) then 
!     print *, "ERROR: Unable to open VIIRS SMV10 file: ", status
!     print *, "File: ", trim(svm10_file)
!     return
!   end if
!   
!   dset_name = "/All_Data/VIIRS-M10-SDR_All/ReflectanceFactors"
!   dims1(1) = 2
!   call h5ltread_dataset_float_f(h5f_id, trim(dset_name), refl_factors, dims1, status)
!   if (status < 0) then 
!     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
!     return
!   end if
!   
!   dset_name = "/All_Data/VIIRS-M10-SDR_All/Reflectance"
!   call h5ltread_dataset_int_f(h5f_id, trim(dset_name), i_refl, dims2, status)
!   if (status < 0) then 
!     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
!     return
!   end if
!   
!   vdbs%m10_refl(:,:) = -999.0
!   where(i_refl < 65528) vdbs%m10_refl = i_refl * refl_factors(1)
!     
!   call h5fclose_f(h5f_id, status)
!   if (status < 0) then
!     print *, "ERROR: Failed to close SVM10 file: ", status
!     return
!   end if
!   
!   call h5fopen_f(trim(svm11_file), H5F_ACC_RDONLY_F, h5f_id, status)
!   if (status < 0) then 
!     print *, "ERROR: Unable to open VIIRS SMV11 file: ", status
!     print *, "File: ", trim(svm11_file)
!     return
!   end if
!   
!   dset_name = "/All_Data/VIIRS-M11-SDR_All/ReflectanceFactors"
!   dims1(1) = 2
!   call h5ltread_dataset_float_f(h5f_id, trim(dset_name), refl_factors, dims1, status)
!   if (status < 0) then 
!     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
!     return
!   end if
! 
!   dset_name = "/All_Data/VIIRS-M11-SDR_All/Reflectance"
!   call h5ltread_dataset_int_f(h5f_id, trim(dset_name), i_refl, dims2, status)
!   if (status < 0) then 
!     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
!     return
!   end if
!   
!   vdbs%m11_refl(:,:) = -999.0
!   where(i_refl < 65528) vdbs%m11_refl = i_refl * refl_factors(1)
!     
!   call h5fclose_f(h5f_id, status)
!   if (status < 0) then
!     print *, "ERROR: Failed to close SVM11 file: ", status
!     return
!   end if
!   
!   call h5fopen_f(trim(svm14_file), H5F_ACC_RDONLY_F, h5f_id, status)
!   if (status < 0) then 
!     print *, "ERROR: Unable to open VIIRS SMV05 file: ", status
!     print *, "File: ", trim(svm14_file)
!     return
!   end if
!   
!   dset_name = "/All_Data/VIIRS-M14-SDR_All/RadianceFactors"
!   dims1(1) = 2
!   rad_factors(:) = 0.0
!   call h5ltread_dataset_float_f(h5f_id, trim(dset_name), rad_factors, dims1, status)
!   if (status < 0) then 
!     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
!     return
!   end if
!   
!   dset_name = "/All_Data/VIIRS-M14-SDR_All/Radiance"
!   i_rad(:,:) = -9999
!   call h5ltread_dataset_int_f(h5f_id, trim(dset_name), i_rad, dims2, status)
!   if (status < 0) then 
!     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
!     return
!   end if
!   
!   vdbs%m14_rad(:,:) = -999.0
!   where(i_rad < 65528) vdbs%m14_rad = i_rad * rad_factors(1) + rad_factors(2)
!   
!   dset_name = "/All_Data/VIIRS-M14-SDR_All/BrightnessTemperatureFactors"
!   dims1(1) = 2
!   bt_factors(:) = 0.0
!   call h5ltread_dataset_float_f(h5f_id, trim(dset_name), bt_factors, dims1, status)
!   if (status < 0) then 
!     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
!     return
!   end if
!   
!   dset_name = "/All_Data/VIIRS-M14-SDR_All/BrightnessTemperature"
!   i_bt(:,:) = -9999
!   call h5ltread_dataset_int_f(h5f_id, trim(dset_name), i_bt, dims2, status)
!   if (status < 0) then 
!     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
!     return
!   end if
!   
!   vdbs%m14_bt(:,:) = -999.0
!   where(i_bt < 65528) vdbs%m14_bt = i_bt * bt_factors(1) + bt_factors(2)
! 
!   call h5fclose_f(h5f_id, status)
!   if (status < 0) then
!     print *, "ERROR: Failed to close SVM14 file: ", status
!     return
!   end if
!   
!   call h5fopen_f(trim(svm15_file), H5F_ACC_RDONLY_F, h5f_id, status)
!   if (status < 0) then 
!     print *, "ERROR: Unable to open VIIRS SMV05 file: ", status
!     print *, "File: ", trim(svm15_file)
!     return
!   end if
!   
!   dset_name = "/All_Data/VIIRS-M15-SDR_All/RadianceFactors"
!   dims1(1) = 2
!   rad_factors(:) = 0.0
!   call h5ltread_dataset_float_f(h5f_id, trim(dset_name), rad_factors, dims1, status)
!   if (status < 0) then 
!     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
!     return
!   end if
!   
!   dset_name = "/All_Data/VIIRS-M15-SDR_All/Radiance"
!   i_rad(:,:) = -9999
!   call h5ltread_dataset_int_f(h5f_id, trim(dset_name), i_rad, dims2, status)
!   if (status < 0) then 
!     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
!     return
!   end if
!   
!   vdbs%m15_rad(:,:) = -999.0
!   where(i_rad < 65528) vdbs%m15_rad = i_rad * rad_factors(1) + rad_factors(2)
!   
!   dset_name = "/All_Data/VIIRS-M15-SDR_All/BrightnessTemperatureFactors"
!   dims1(1) = 2
!   bt_factors(:) = 0.0
!   call h5ltread_dataset_float_f(h5f_id, trim(dset_name), bt_factors, dims1, status)
!   if (status < 0) then
!     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
!     return
!   end if
! 
!   dset_name = "/All_Data/VIIRS-M15-SDR_All/BrightnessTemperature"
!   i_bt(:,:) = -9999
!   call h5ltread_dataset_int_f(h5f_id, trim(dset_name), i_bt, dims2, status)
!   if (status < 0) then
!     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
!     return
!   end if
!   
!   vdbs%m15_bt(:,:) = -999.0
!   where(i_bt < 65528) vdbs%m15_bt = i_bt * bt_factors(1) + bt_factors(2)
!   
!   call h5fclose_f(h5f_id, status)
!   if (status < 0) then
!     print *, "ERROR: Failed to close SVM15 file: ", status
!     return
!   end if
!   
!   call h5fopen_f(trim(svm16_file), H5F_ACC_RDONLY_F, h5f_id, status)
!   if (status < 0) then 
!     print *, "ERROR: Unable to open VIIRS SMV16 file: ", status
!     print *, "File: ", trim(svm16_file)
!     return
!   end if
!   
!   dset_name = "/All_Data/VIIRS-M16-SDR_All/RadianceFactors"
!   dims1(1) = 2
!   rad_factors(:) = 0.0
!   call h5ltread_dataset_float_f(h5f_id, trim(dset_name), rad_factors, dims1, status)
!   if (status < 0) then 
!     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
!     return
!   end if
!   
!   dset_name = "/All_Data/VIIRS-M16-SDR_All/Radiance"
!   i_rad(:,:) = -9999
!   call h5ltread_dataset_int_f(h5f_id, trim(dset_name), i_rad, dims2, status)
!   if (status < 0) then 
!     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
!     return
!   end if
!   
!   vdbs%m16_rad(:,:) = -999.0
!   where(i_rad < 65528) vdbs%m16_rad = i_rad * rad_factors(1) + rad_factors(2)
!     
!   dset_name = "/All_Data/VIIRS-M16-SDR_All/BrightnessTemperatureFactors"
!   dims1(1) = 2
!   bt_factors(:) = 0.0
!   call h5ltread_dataset_float_f(h5f_id, trim(dset_name), bt_factors, dims1, status)
!   if (status < 0) then 
!     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
!     return
!   end if
!   
!   dset_name = "/All_Data/VIIRS-M16-SDR_All/BrightnessTemperature"
!   i_bt(:,:) = -9999
!   call h5ltread_dataset_int_f(h5f_id, trim(dset_name), i_bt, dims2, status)
!   if (status < 0) then 
!     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
!     return
!   end if
!   
!   vdbs%m16_bt(:,:) = -999.0
!   where(i_bt < 65528) vdbs%m16_bt = i_bt * bt_factors(1) + bt_factors(2)
!   
!   call h5fclose_f(h5f_id, status)
!   if (status < 0) then
!     print *, "ERROR: Failed to close SVM16 file: ", status
!     return
!   end if
!   
!   call h5fopen_f(trim(iicmo_file), H5F_ACC_RDONLY_F, h5f_id, status)
!   if (status < 0) then 
!     print *, "ERROR: Unable to open VIIRS IICMO file: ", status
!     print *, "File: ", trim(iicmo_file)
!     return
!   end if
!   
!   dset_name = "/All_Data/VIIRS-CM-IP_All/QF1_VIIRSCMIP"
!   call h5ltread_dataset_int_f(h5f_id, trim(dset_name), cld_mask, dims2, status)
!   if (status < 0) then 
!     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
!     return
!   end if
!   
!   bit_mask = b'00000011'
!   cld_mask_qa = iand(cld_mask, bit_mask)
! 
!   dset_name = "/All_Data/VIIRS-CM-IP_All/QF2_VIIRSCMIP"
!   call h5ltread_dataset_int_f(h5f_id, trim(dset_name), lw_mask, dims2, status)
!   if (status < 0) then 
!     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
!     return
!   end if
!   
!   vdbs%land_mask(:,:) = -999
!   bit_mask = b'00000111'
!   where(cld_mask_qa > 0) vdbs%land_mask = iand(lw_mask, bit_mask)
! 
!   call h5fclose_f(h5f_id, status)
!   if (status < 0) then
!     print *, "ERROR: Failed to close IICMO file: ", status
!     return
!   end if
!   
!   call h5close_f(status)
!   if (status /= 0) then 
!     print *, "ERROR: Failed to close HDF interface: ", status
!     return
!   end if
!   
! ! -- calculate brightness temperature differences, D* parameter, and NDVI.
!   vdbs%btd8(:,:) = -999.0
!   where (vdbs%m14_bt > -900.0 .AND. vdbs%m15_bt > -900.0) 
!     vdbs%btd8 = vdbs%m14_bt - vdbs%m15_bt
!   end where
! 
!   vdbs%btd11(:,:) = -999.0
!   where (vdbs%m15_bt > -900.0 .AND. vdbs%m16_bt > -900.0) 
!     vdbs%btd11 = vdbs%m15_bt - vdbs%m16_bt
!   end where
!   
!   vdbs%dstar(:,:) = -999.0
!   status = calc_dstar(vdbs%btd8, vdbs%btd11, vdbs%dstar) 
!   if (status /= 0) then
!     print *, "ERROR: Failed to calculate D* values: ", status
!     return
!   end if
!     
!   vdbs%ndvi(:,:) = -999.0
!   where (vdbs%m07_refl > -900.0 .AND. vdbs%m05_refl > -900.0)
!     vdbs%ndvi = (vdbs%m07_refl - vdbs%m05_refl) / (vdbs%m07_refl + vdbs%m05_refl)
!   end where
!   
! ! -- calculate the mirror normal vector for the polarization correction
! !   allocate(sen_mat(3,3), coeffs(10), stat=status)
! !   if (status /= 0) then
! !     print *, "ERROR: Failed to allocate sensor orientation arrays: ", status
! !     return
! !   end if
! !   
! !   do i = 1, nscan 
! !     sen_mat(:,:) = 0.0
! !     coeffs(:) = 0.0    
! !     call calc_sensor_orientation(sc_position(:,i), sc_velocity(:,i), sc_attitude(:,i), sen_mat, coeffs)
! !     vdbs%mnorm(:,i) = sen_mat(:,1)
! !   end do
!   
! ! -- clean up!
!   deallocate(i_refl, i_rad, i_bt, cld_mask, cld_mask_qa, lw_mask, stat=status)
!   if (status /= 0) then
!     print *, "WARNING: Failed to deallocate integer data arrays: ", status
!   end if
!   return
!   
! end function load_viirs_db_data

type (viirs_db_svm) function load_viirs_db_data_nasa(geo_file, mband_file, status) result (vdbs)
  
  use netcdf
  use calendars, only: gdatetime, doy_from_gregorian
                  
  implicit none

  character(len=*), intent(in)          ::  geo_file
  character(len=*), intent(in)          ::  mband_file
  integer, intent(inout)                ::  status

  integer                               ::  nc_id
  integer                               ::  dim_id
  integer                               ::  dset_id    
  integer                               ::  grp_id     
  
  character(len=255)                    ::  dset_name
  character(len=255)                    ::  attr_name
  character(len=255)                    ::  group_name
  character(len=5)                      ::  day_or_night
  character(len=24)                     ::  start_datetime
  character(len=10)                     ::  start_date
  character(len=5)                      ::  start_time  
  
  real, dimension(:,:), allocatable     ::  saa, vaa
  integer, dimension(:,:), allocatable  ::  i_val
  integer                               ::  i_fillvalue
  integer, dimension(:,:), allocatable  ::  lw_mask
  real                                  ::  scale_factor
  real                                  ::  add_offset
  real, dimension(:), allocatable       ::  rad_lut
  real, dimension(:), allocatable       ::  bt_lut
  integer                               ::  indx
  real                                  ::  frac

  byte, dimension(:), allocatable       ::  flag_masks
  integer(kind=i8), dimension(:,:), allocatable     ::  qa
  integer, dimension(:,:), allocatable  ::  qa_mask
  integer                               ::  dtype
  
  integer                               ::  nvals

  integer                               ::  scan, xscan, nlut, nscan
  real, parameter                       ::  d2r =  3.14159/180.0

  integer                               ::  i, j, k
  real                                  ::  tmp0, tmp1
  
  type(gdatetime)                       ::  gdt1

  real(kind=r64), dimension(:), allocatable :: tmp_scan_time    ! tmp buffer

  
  status = -1
  
  status = nf90_open(geo_file, nf90_nowrite, nc_id)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to open geolocation file: ", status
    return
  end if
  
  attr_name = "DayNightFlag"
  status = nf90_get_att(nc_id, NF90_GLOBAL, attr_name, day_or_night)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if

  if (day_or_night == "Night") then
    print *, "ERROR: Granule is a night granule. Skipping."
    status = -1
    return
  end if
  
  attr_name = "time_coverage_start"
  status = nf90_get_att(nc_id, NF90_GLOBAL, attr_name, start_datetime)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if
  
  start_date = start_datetime(1:10)  
  read(start_date(1:4), fmt='(I4)') vdbs%yr
  read(start_date(6:7), fmt='(I2)') vdbs%mo
  read(start_date(9:10), fmt='(I2)') vdbs%dy
  
  start_time = start_datetime(12:16)
  read(start_time(1:2), fmt='(I2)') vdbs%hr
  read(start_time(4:5), fmt='(I2)') vdbs%min
  
  gdt1 = gdatetime(vdbs%yr, vdbs%mo, vdbs%dy, vdbs%hr, vdbs%min, 0, 0, 0) 
  vdbs%doy = doy_from_gregorian(gdt1)
  
  
  dset_name = 'number_of_lines'
  status = nf90_inq_dimid(nc_id, dset_name, dim_id)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to get size of dimension "//trim(dset_name)//": ", status
    return
  end if
  
  status = nf90_inquire_dimension(nc_id, dim_id, len = scan)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to size of dimension "//trim(dset_name)//": ", status
    return
  end if
  
  dset_name = 'number_of_pixels'
  status = nf90_inq_dimid(nc_id, dset_name, dim_id)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to get size of dimension "//trim(dset_name)//": ", status
    return
  end if
  status = nf90_inquire_dimension(nc_id, dim_id, len = xscan)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to size of dimension "//trim(dset_name)//": ", status
    return
  end if
  
  dset_name = 'number_of_scans'
  status = nf90_inq_dimid(nc_id, dset_name, dim_id)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to get size of dimension "//trim(dset_name)//": ", status
    return
  end if
  status = nf90_inquire_dimension(nc_id, dim_id, len = nscan)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to size of dimension "//trim(dset_name)//": ", status
    return
  end if
  
  vdbs%scan   = scan
  vdbs%xscan  = xscan 
  print *, 'xscan, scan: ', xscan, scan
  
  allocate(vdbs%scan_time(scan), stat=status)
  if (status /= 0) then
    print *, "ERROR: Failed to allocate scan time arrays: ", status
    return
  end if
  
  allocate(vdbs%lat(xscan,scan), vdbs%lon(xscan,scan), vdbs%sza(xscan,scan),  &
  &       vdbs%vza(xscan,scan), vdbs%raa(xscan,scan), vdbs%sca(xscan,scan),   &
  &       vdbs%saa(xscan,scan), vdbs%vaa(xscan,scan),   &
  &       vdbs%gla(xscan,scan), vdbs%amf(xscan,scan), stat=status)
  if (status /= 0) then 
    print *, "ERROR: Failed to allocate geolocation and/or geometry data arrays: ", status
    return
  end if
  
  allocate(vdbs%m01_refl(xscan,scan), vdbs%m02_refl(xscan,scan), &
  &         vdbs%m03_refl(xscan,scan), vdbs%m04_refl(xscan,scan),&
  &         vdbs%m05_refl(xscan,scan), i_val(xscan,scan), stat=status)
  if (status /= 0) then
    print *, "ERROR: Failed to allocate reflectance arrays: ", status
    return
  end if
 
  allocate(vdbs%m07_refl(xscan,scan), vdbs%m08_refl(xscan,scan),  &
  &         vdbs%m09_refl(xscan,scan), vdbs%m10_refl(xscan,scan), &
  &         vdbs%m11_refl(xscan,scan), stat=status)
  if (status /= 0) then
    print *, "ERROR: Failed to allocate reflectance arrays: ", status
    return
  end if
  
  allocate(vdbs%m14_rad(xscan,scan), vdbs%m15_rad(xscan,scan),  &
  &         vdbs%m16_rad(xscan,scan), stat=status)
  if (status /= 0) then
    print *, "ERROR: Failed to allocate radiance arrays: ", status
    return
  end if
  
  allocate(vdbs%m14_bt(xscan, scan), vdbs%m15_bt(xscan, scan),  &
  & vdbs%m16_bt(xscan, scan), stat=status)
  if (status /= 0) then
    print *, "ERROR: Failed to allocate brightness temperature arrays: ", status
    return
  end if

  allocate(vdbs%land_mask(xscan,scan), vdbs%elev(xscan,scan), vdbs%ps(xscan,scan), &
  &   stat=status)
  if (status /= 0) then
    print *, "ERROR: Failed to allocate mask arrays: ", status
    return
  end if
  
  allocate(vdbs%btd8(xscan,scan), vdbs%btd11(xscan,scan), stat=status)
  if (status /= 0) then
    print *, "ERROR: Failed to allocate BTD arrays: ", status
    return
  end if
  
  allocate(vdbs%dstar(xscan,scan), vdbs%ndvi(xscan, scan), stat=status)
  if (status /= 0) then
    print *, "ERROR: Failed to allocate D* arrays: ", status
    return
  end if

  allocate(saa(xscan,scan), vaa(xscan,scan), lw_mask(xscan,scan), stat=status)
  if (status /= 0) then 
    print *, "ERROR: Unable to allocate azimuth angle arrays: ", status
    return
  end if
  
  allocate(latitude(xscan,scan), longitude(xscan,scan), stat=status)
  if (status /= 0) then
    print *, "ERROR: Failed to allocate shared geolocation arrays: ", status
    return
  end if
  
  allocate(qa_mask(xscan,scan), qa(xscan,scan), stat=status)
  if (status /= 0) then
    print *, "ERROR: Failed to allocate L1B quality assurance array: ", status
    return
  end if
  qa_mask(:,:) = b'00000000'
  qa(:,:) = 0
  
!! -- read scan times first. HDF5 high-level API doesn't have a function to read
!! -- long ints (or anything larger than ints), so we have to use the full library.
  group_name = 'scan_line_attributes'
  status = nf90_inq_ncid(nc_id, group_name, grp_id)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to get ID of group "//trim(group_name)//": ", status
    return
  end if
  
  dset_name = "scan_start_time"
  status = nf90_inq_varid(grp_id, dset_name, dset_id)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
    return
  end if
    
  allocate(tmp_scan_time(nscan), stat=status)
  if (status /= 0) then
    print *, "ERROR: Failed to allocate tmp array for scan times: ", status
    return
  end if
  
  status = nf90_get_var(grp_id, dset_id, tmp_scan_time)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
    return
  end if           

! -- each "scan" technically is 16 pixels wide. Expand times
! -- to cover full granule.
  do i = 1, scan
    j = (i-1)/16 + 1    ! integer division
    vdbs%scan_time(i) = tmp_scan_time(j)
  end do
  
  deallocate(tmp_scan_time, stat=status)
  if (status /= 0) then
    print *, "ERROR: Failed to deallocate tmp array for time data: ", status
    return
  end if
  
  dset_name = "HAM_side"
  status = nf90_inq_varid(grp_id, dset_name, dset_id)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
    return
  end if
  
  allocate(vdbs%mirror_side(nscan), stat=status)
  if (status /= 0) then
    print *, "ERROR: Failed to allocate array for mirror side data: ", status
    return
  end if
  
  status = nf90_get_var(grp_id, dset_id, vdbs%mirror_side)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
    return
  end if   
  
  group_name = 'geolocation_data'
  status = nf90_inq_ncid(nc_id, group_name, grp_id)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to get ID of group "//trim(group_name)//": ", status
    return
  end if
  
  dset_name = 'latitude'
  status = nf90_inq_varid(grp_id, dset_name, dset_id)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
    return
  end if
  
  status = nf90_get_var(grp_id, dset_id, latitude)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
    return
  end if
  
  vdbs%lat = latitude

  dset_name = 'longitude'
  status = nf90_inq_varid(grp_id, dset_name, dset_id)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
    return
  end if
  
  status = nf90_get_var(grp_id, dset_id, longitude)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
    return
  end if            
  
  vdbs%lon = longitude
  
  dset_name = 'solar_zenith'
  status = nf90_inq_varid(grp_id, dset_name, dset_id)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
    return
  end if
  
  status = nf90_get_var(grp_id, dset_id, i_val)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
    return
  end if     
  
  attr_name = "scale_factor"
  status = nf90_get_att(grp_id, dset_id, attr_name, scale_factor)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if
  
  attr_name = "add_offset"
  status = nf90_get_att(grp_id, dset_id, attr_name, add_offset)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if

  attr_name = "_FillValue"
  status = nf90_get_att(grp_id, dset_id, attr_name, i_fillvalue)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if

  vdbs%sza(:,:) = -999.0
  where(i_val /= i_fillvalue) vdbs%sza = i_val * scale_factor + add_offset
  print *, 'Count SZA fills: ', count(i_val == i_fillvalue)
  
  dset_name = 'sensor_zenith'
  status = nf90_inq_varid(grp_id, dset_name, dset_id)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
    return
  end if
  
  status = nf90_get_var(grp_id, dset_id, i_val)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
    return
  end if    
  
  attr_name = "scale_factor"
  status = nf90_get_att(grp_id, dset_id, attr_name, scale_factor)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if

  attr_name = "add_offset"
  status = nf90_get_att(grp_id, dset_id, attr_name, add_offset)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if

  attr_name = "_FillValue"
  status = nf90_get_att(grp_id, dset_id, attr_name, i_fillvalue)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if

  vdbs%vza(:,:) = -999.0
  where(i_val /= i_fillvalue) vdbs%vza = i_val * scale_factor + add_offset
  print *, 'Count VZA fills: ', count(i_val == i_fillvalue)
  
  dset_name = 'solar_azimuth'
  status = nf90_inq_varid(grp_id, dset_name, dset_id)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
    return
  end if

  status = nf90_get_var(grp_id, dset_id, i_val)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
    return
  end if  
  
  attr_name = "scale_factor"
  status = nf90_get_att(grp_id, dset_id, attr_name, scale_factor)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if
  
  attr_name = "add_offset"
  status = nf90_get_att(grp_id, dset_id, attr_name, add_offset)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if

  attr_name = "_FillValue"
  status = nf90_get_att(grp_id, dset_id, attr_name, i_fillvalue)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if
  
  saa(:,:) = -999.0
  where(i_val /= i_fillvalue) saa = i_val * scale_factor + add_offset
  
  dset_name = 'sensor_azimuth'
  status = nf90_inq_varid(grp_id, dset_name, dset_id)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
    return
  end if
  
  status = nf90_get_var(grp_id, dset_id, i_val)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
    return
  end if       
  
  attr_name = "scale_factor"
  status = nf90_get_att(grp_id, dset_id, attr_name, scale_factor)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if
  
  attr_name = "add_offset"
  status = nf90_get_att(grp_id, dset_id, attr_name, add_offset)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if

  attr_name = "_FillValue"
  status = nf90_get_att(grp_id, dset_id, attr_name, i_fillvalue)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if

  vaa(:,:) = -999.0
  where(i_val /= i_fillvalue) vaa = i_val * scale_factor + add_offset
  print *, 'Count VAA fills: ', count(i_val == i_fillvalue)
  vdbs%vaa  = vaa
  vdbs%saa  = saa
  vdbs%raa  = -999.0  
  vdbs%raa  = (vaa - saa) - 180.0
  where (vdbs%raa > 180.0)  vdbs%raa = vdbs%raa - 360.0
  where (vdbs%raa < -180.0) vdbs%raa = vdbs%raa + 360.0
  where (vdbs%raa < 0.0)    vdbs%raa = -1.0 * vdbs%raa
  where(i_val == i_fillvalue) vdbs%raa = -999.0
    
  vdbs%sca  = acos(cos(vdbs%sza*d2r)*cos(vdbs%vza*d2r) -       &
                sin(vdbs%sza*d2r)*sin(vdbs%vza*d2r)*cos(vdbs%raa*d2r))
  vdbs%sca  = 180.0 - (vdbs%sca/d2r)
  
  vdbs%gla  = acos(cos(vdbs%sza*d2r)*cos(vdbs%vza*d2r) +       &
                sin(vdbs%sza*d2r)*sin(vdbs%vza*d2r)*cos(vdbs%raa*d2r))
  vdbs%gla = vdbs%gla/d2r

  vdbs%amf = 1.0/cos(vdbs%sza*d2r)+1.0/cos(vdbs%vza*d2r)

  where(i_val == i_fillvalue) 
    vdbs%sca = -999.0
    vdbs%gla = -999.0
    vdbs%amf = -999.0
  endwhere     
 
  deallocate(saa, vaa, stat=status)
  if (status /= 0) then
    print *, "WARNING: Failed to deallocate intermediate angle arrays: ", status
  end if
  
  dset_name = 'height'
  status = nf90_inq_varid(grp_id, dset_name, dset_id)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
    return
  end if
  
  status = nf90_get_var(grp_id, dset_id, i_val)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
    return
  end if

  attr_name = "scale_factor"
  status = nf90_get_att(grp_id, dset_id, attr_name, scale_factor)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if
  
  attr_name = "add_offset"
  status = nf90_get_att(grp_id, dset_id, attr_name, add_offset)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if
 
  attr_name = "_FillValue"
  status = nf90_get_att(grp_id, dset_id, attr_name, i_fillvalue)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if
 
  vdbs%elev(:,:) = -999.0
  where(i_val /= i_fillvalue) vdbs%elev = i_val * scale_factor + add_offset
  
! -- convert surface elevation to surface pressure
! -- output of atmosphere is fraction of mean sea-level pressure.
  do j = 1, scan
    do i = 1, xscan
      call atmosphere(vdbs%elev(i,j)/1000.0, tmp0, vdbs%ps(i,j), tmp1)
      vdbs%ps(i,j) = vdbs%ps(i,j) * 1013.25
    end do
  end do
  
  dset_name = 'land_water_mask'
  status = nf90_inq_varid(grp_id, dset_name, dset_id)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
    return
  end if
  
  status = nf90_get_var(grp_id, dset_id, i_val)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
    return
  end if  

! -- convert land mask values to be consistent with the original NOAA values.  
  vdbs%land_mask(:,:) = -999
  where (i_val == 0 .OR. i_val == 6 .OR. i_val == 7)
    vdbs%land_mask = 3
  end where
  where (i_val == 1)
    vdbs%land_mask = 0
  end where
  where (i_val == 2)
    vdbs%land_mask = 5
  end where
  where (i_val == 3 .OR. i_val == 5)
    vdbs%land_mask = 2
  end where
  where (i_val == 4)
    vdbs%land_mask = 1
  end where
  
!  bit_mask = b'00000111'
!  where(cld_mask_qa > 0) vdbs%land_mask = iand(lw_mask, bit_mask)

  status = nf90_close(nc_id)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to close geolocation file: ", status
    return
  end if
  
  status = nf90_open(mband_file, nf90_nowrite, nc_id)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to open M-band L1B file: ", status
    return
  end if
    
  dset_name = 'number_of_LUT_values'
  status = nf90_inq_dimid(nc_id, dset_name, dim_id)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to get size of dimension "//trim(dset_name)//": ", status
    return
  end if
  
  status = nf90_inquire_dimension(nc_id, dim_id, len = nlut)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to size of dimension "//trim(dset_name)//": ", status
    return
  end if
  
  allocate(rad_lut(nlut), bt_lut(nlut), stat=status)
  if (status /= 0) then
    print *, "ERROR: Failed to allocate arrays for radiance and brightness "// &
    & "temperature look-up tables: ", status
    return
  end if
  rad_lut(:)  = -999.0
  bt_lut(:)   = -999.0
    
  group_name = 'observation_data'
  status = nf90_inq_ncid(nc_id, group_name, grp_id)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to get ID of group "//trim(group_name)//": ", status
    return
  end if
  
  dset_name = 'M01'
  status = nf90_inq_varid(grp_id, dset_name, dset_id)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
    return
  end if

  attr_name = "_FillValue"
  status = nf90_get_att(grp_id, dset_id, attr_name, i_fillvalue)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if
  
  attr_name = "scale_factor"
  status = nf90_get_att(grp_id, dset_id, attr_name, scale_factor)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if
  
  attr_name = "add_offset"
  status = nf90_get_att(grp_id, dset_id, attr_name, add_offset)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if
  
  status = nf90_get_var(grp_id, dset_id, i_val)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
    return
  end if  
  
  vdbs%m01_refl(:,:) = -999.0
  where(i_val /= i_fillvalue) 
    vdbs%m01_refl = i_val * scale_factor + add_offset
    vdbs%m01_refl = vdbs%m01_refl / cos(vdbs%sza*d2r)
  end where
  
  
  dset_name = 'M02'
  status = nf90_inq_varid(grp_id, dset_name, dset_id)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
    return
  end if
  
  attr_name = "_FillValue"
  status = nf90_get_att(grp_id, dset_id, attr_name, i_fillvalue)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if
  
  attr_name = "scale_factor"
  status = nf90_get_att(grp_id, dset_id, attr_name, scale_factor)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if
  
  attr_name = "add_offset"
  status = nf90_get_att(grp_id, dset_id, attr_name, add_offset)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if
  
  status = nf90_get_var(grp_id, dset_id, i_val)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
    return
  end if  
  
  vdbs%m02_refl(:,:) = -999.0
  where(i_val /= i_fillvalue) 
    vdbs%m02_refl = i_val * scale_factor + add_offset
    vdbs%m02_refl = vdbs%m02_refl / cos(vdbs%sza*d2r)
  end where
  
  dset_name = 'M03'
  status = nf90_inq_varid(grp_id, dset_name, dset_id)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
    return
  end if
  
  attr_name = "_FillValue"
  status = nf90_get_att(grp_id, dset_id, attr_name, i_fillvalue)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if
  
  attr_name = "scale_factor"
  status = nf90_get_att(grp_id, dset_id, attr_name, scale_factor)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if
  
  attr_name = "add_offset"
  status = nf90_get_att(grp_id, dset_id, attr_name, add_offset)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if
  
  status = nf90_get_var(grp_id, dset_id, i_val)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
    return
  end if  
  
  vdbs%m03_refl(:,:) = -999.0
  where(i_val /= i_fillvalue) 
    vdbs%m03_refl = i_val * scale_factor + add_offset
    vdbs%m03_refl = vdbs%m03_refl / cos(vdbs%sza*d2r)
  end where
  
  dset_name = 'M04'
  status = nf90_inq_varid(grp_id, dset_name, dset_id)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
    return
  end if
  
  attr_name = "_FillValue"
  status = nf90_get_att(grp_id, dset_id, attr_name, i_fillvalue)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if
  
  attr_name = "scale_factor"
  status = nf90_get_att(grp_id, dset_id, attr_name, scale_factor)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if
  
  attr_name = "add_offset"
  status = nf90_get_att(grp_id, dset_id, attr_name, add_offset)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if
  
  status = nf90_get_var(grp_id, dset_id, i_val)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
    return
  end if  
  
  vdbs%m04_refl(:,:) = -999.0
  where(i_val /= i_fillvalue) 
    vdbs%m04_refl = i_val * scale_factor + add_offset
    vdbs%m04_refl = vdbs%m04_refl / cos(vdbs%sza*d2r)
  end where
  
  dset_name = 'M05'
  status = nf90_inq_varid(grp_id, dset_name, dset_id)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
    return
  end if
  
  attr_name = "_FillValue"
  status = nf90_get_att(grp_id, dset_id, attr_name, i_fillvalue)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if
  
  attr_name = "scale_factor"
  status = nf90_get_att(grp_id, dset_id, attr_name, scale_factor)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if
  
  attr_name = "add_offset"
  status = nf90_get_att(grp_id, dset_id, attr_name, add_offset)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if
  
  status = nf90_get_var(grp_id, dset_id, i_val)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
    return
  end if  
  
  vdbs%m05_refl(:,:) = -999.0
  where(i_val /= i_fillvalue) 
    vdbs%m05_refl = i_val * scale_factor + add_offset
    vdbs%m05_refl = vdbs%m05_refl / cos(vdbs%sza*d2r)
  end where
  
  dset_name = 'M07'
  status = nf90_inq_varid(grp_id, dset_name, dset_id)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
    return
  end if
  
  attr_name = "_FillValue"
  status = nf90_get_att(grp_id, dset_id, attr_name, i_fillvalue)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if
  
  attr_name = "scale_factor"
  status = nf90_get_att(grp_id, dset_id, attr_name, scale_factor)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if
  
  attr_name = "add_offset"
  status = nf90_get_att(grp_id, dset_id, attr_name, add_offset)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if
  
  status = nf90_get_var(grp_id, dset_id, i_val)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
    return
  end if  
  
  vdbs%m07_refl(:,:) = -999.0
  where(i_val /= i_fillvalue)
    vdbs%m07_refl = i_val * scale_factor + add_offset
    vdbs%m07_refl = vdbs%m07_refl / cos(vdbs%sza*d2r)
  end where
  
!   dset_name = 'M07_quality_flags'
!   status = nf90_inq_varid(grp_id, dset_name, dset_id)
!   if (status /= NF90_NOERR) then
!     print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
!     return
!   end if
!   
!   attr_name = "flag_masks"
!   status = nf90_inquire_attribute(grp_id, dset_id, attr_name, len=nvals, xtype=dtype)
!   if (status /= NF90_NOERR) then
!     print *, "ERROR: Failed to inquire "//trim(attr_name)//": ", status
!     return
!   end if  
!   print * ,'flag_masks type: ', dtype
!   
!   allocate(flag_masks(nvals), stat=status)
!   if (status /= 0) then
!     print *, "ERROR: Failed to allocate array for flag masks on dset "//trim(dset_name)//": ", status
!     return
!   end if
!   
!   status = nf90_get_att(grp_id, dset_id, attr_name, flag_masks)
!   if (status /= NF90_NOERR) then
!     print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
!     return
!   end if  
!   
!   print *, 'flag masks: ', flag_masks
!   
!   status = nf90_inquire_variable(grp_id, dset_id, xtype=dtype)
!   print *, 'QA type: ', dtype, NF90_UBYTE
!             
!   status = nf90_get_var(grp_id, dset_id, qa)
!   if (status /= NF90_NOERR) then
!     print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
!     return
!   end if  
!   print *, 'qa(1500,1): ', qa(:,1)
!   stop
!   print *, 'starting QA masking...'  
!   qa_mask(:,:) = 0
!   do k = 1, size(flag_masks,1) 
!     do j = 1, scan
!       do i = 1, xscan
!         if (and(flag_masks(i), qa(i,j)) /= 0) qa_mask(i,j) = 1
!       end do
!     end do
!   end do 
!   print *, 'end qa masking.'
  
  dset_name = 'M08'
  status = nf90_inq_varid(grp_id, dset_name, dset_id)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
    return
  end if
  
  attr_name = "_FillValue"
  status = nf90_get_att(grp_id, dset_id, attr_name, i_fillvalue)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if
  
  attr_name = "scale_factor"
  status = nf90_get_att(grp_id, dset_id, attr_name, scale_factor)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if

  attr_name = "add_offset"
  status = nf90_get_att(grp_id, dset_id, attr_name, add_offset)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if
  
  status = nf90_get_var(grp_id, dset_id, i_val)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
    return
  end if  
  
  vdbs%m08_refl(:,:) = -999.0
  where(i_val /= i_fillvalue) 
    vdbs%m08_refl = i_val * scale_factor + add_offset
    vdbs%m08_refl = vdbs%m08_refl / cos(vdbs%sza*d2r)
  end where
  
  dset_name = 'M09'
  status = nf90_inq_varid(grp_id, dset_name, dset_id)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
    return
  end if
  
  attr_name = "_FillValue"
  status = nf90_get_att(grp_id, dset_id, attr_name, i_fillvalue)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if
  
  attr_name = "scale_factor"
  status = nf90_get_att(grp_id, dset_id, attr_name, scale_factor)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if
  
  attr_name = "add_offset"
  status = nf90_get_att(grp_id, dset_id, attr_name, add_offset)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if
  
  status = nf90_get_var(grp_id, dset_id, i_val)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
    return
  end if  
  
  vdbs%m09_refl(:,:) = -999.0
  where(i_val /= i_fillvalue) 
    vdbs%m09_refl = i_val * scale_factor + add_offset
    vdbs%m09_refl = vdbs%m09_refl / cos(vdbs%sza*d2r)
  end where
  
  dset_name = 'M10'
  status = nf90_inq_varid(grp_id, dset_name, dset_id)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
    return
  end if
  
  attr_name = "_FillValue"
  status = nf90_get_att(grp_id, dset_id, attr_name, i_fillvalue)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if
  
  attr_name = "scale_factor"
  status = nf90_get_att(grp_id, dset_id, attr_name, scale_factor)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if
  
  attr_name = "add_offset"
  status = nf90_get_att(grp_id, dset_id, attr_name, add_offset)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if
  
  status = nf90_get_var(grp_id, dset_id, i_val)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
    return
  end if  
  
  vdbs%m10_refl(:,:) = -999.0
  where(i_val /= i_fillvalue) 
    vdbs%m10_refl = i_val * scale_factor + add_offset
    vdbs%m10_refl = vdbs%m10_refl / cos(vdbs%sza*d2r)
  end where

  dset_name = 'M11'
  status = nf90_inq_varid(grp_id, dset_name, dset_id)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
    return
  end if
  
  attr_name = "_FillValue"
  status = nf90_get_att(grp_id, dset_id, attr_name, i_fillvalue)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if
  
  attr_name = "scale_factor"
  status = nf90_get_att(grp_id, dset_id, attr_name, scale_factor)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if
  
  attr_name = "add_offset"
  status = nf90_get_att(grp_id, dset_id, attr_name, add_offset)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if
  
  status = nf90_get_var(grp_id, dset_id, i_val)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
    return
  end if  
  
  vdbs%m11_refl(:,:) = -999.0
  where(i_val /= i_fillvalue)
    vdbs%m11_refl = i_val * scale_factor + add_offset
    vdbs%m11_refl = vdbs%m11_refl / cos(vdbs%sza*d2r)
  end where
  
  dset_name = 'M14'
  status = nf90_inq_varid(grp_id, dset_name, dset_id)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
    return
  end if
  
  attr_name = "_FillValue"
  status = nf90_get_att(grp_id, dset_id, attr_name, i_fillvalue)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if

  attr_name = "scale_factor"
  status = nf90_get_att(grp_id, dset_id, attr_name, scale_factor)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if
  
  attr_name = "add_offset"
  status = nf90_get_att(grp_id, dset_id, attr_name, add_offset)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if
  
  status = nf90_get_var(grp_id, dset_id, i_val)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
    return
  end if  
  
  vdbs%m14_rad(:,:) = -999.0
  where(i_val /= i_fillvalue) vdbs%m14_rad = i_val * scale_factor + add_offset
  
  dset_name = 'M14_brightness_temperature_lut'
  status = nf90_inq_varid(grp_id, dset_name, dset_id)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
    return
  end if
  
  status = nf90_get_var(grp_id, dset_id, bt_lut)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
    return
  end if  
  
  do j = 1, scan
    do i = 1, xscan
      if (vdbs%m14_rad(i,j) > -900.0) then
        indx = i_val(i,j)
        vdbs%m14_bt(i,j) = bt_lut(indx)
      end if
    end do
  end do
  
  dset_name = 'M15'
  status = nf90_inq_varid(grp_id, dset_name, dset_id)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
    return
  end if
  
  attr_name = "_FillValue"
  status = nf90_get_att(grp_id, dset_id, attr_name, i_fillvalue)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if
  
  attr_name = "scale_factor"
  status = nf90_get_att(grp_id, dset_id, attr_name, scale_factor)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if
  
  attr_name = "add_offset"
  status = nf90_get_att(grp_id, dset_id, attr_name, add_offset)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if
  
  status = nf90_get_var(grp_id, dset_id, i_val)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
    return
  end if  
  
  vdbs%m15_rad(:,:) = -999.0
  where(i_val /= i_fillvalue) vdbs%m15_rad = i_val * scale_factor + add_offset
  
  dset_name = 'M15_brightness_temperature_lut'
  status = nf90_inq_varid(grp_id, dset_name, dset_id)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
    return
  end if
  
  status = nf90_get_var(grp_id, dset_id, bt_lut)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
    return
  end if  
  
  do j = 1, scan
    do i = 1, xscan
      if (vdbs%m15_rad(i,j) > -900.0) then
        indx = i_val(i,j)
        vdbs%m15_bt(i,j) = bt_lut(indx)
      end if
    end do
  end do
  
  dset_name = 'M16'
  status = nf90_inq_varid(grp_id, dset_name, dset_id)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
    return
  end if
  
  attr_name = "_FillValue"
  status = nf90_get_att(grp_id, dset_id, attr_name, i_fillvalue)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if
  
  attr_name = "scale_factor"
  status = nf90_get_att(grp_id, dset_id, attr_name, scale_factor)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if
  
  attr_name = "add_offset"
  status = nf90_get_att(grp_id, dset_id, attr_name, add_offset)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read attribute "//trim(attr_name)//": ", status
    return
  end if
  
  status = nf90_get_var(grp_id, dset_id, i_val)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
    return
  end if  
  
  vdbs%m16_rad(:,:) = -999.0
  where(i_val /= i_fillvalue) vdbs%m16_rad = i_val * scale_factor + add_offset

  dset_name = 'M16_brightness_temperature_lut'
  status = nf90_inq_varid(grp_id, dset_name, dset_id)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
    return
  end if
  
  status = nf90_get_var(grp_id, dset_id, bt_lut)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
    return
  end if
  
  do j = 1, scan
    do i = 1, xscan
      if (vdbs%m16_rad(i,j) > -900.0) then
        indx = i_val(i,j)
        vdbs%m16_bt(i,j) = bt_lut(indx)
      end if
    end do
  end do
  
  status = nf90_close(nc_id)
  if (status /= NF90_NOERR) then
    print *, "ERROR: Failed to close M-band L1B file: ", status
    return
  end if

! -- calculate brightness temperature differences, D* parameter, and NDVI.
  vdbs%btd8(:,:) = -999.0
  where (vdbs%m14_bt > -900.0 .AND. vdbs%m15_bt > -900.0) 
    vdbs%btd8 = vdbs%m14_bt - vdbs%m15_bt
  end where
  
  vdbs%btd11(:,:) = -999.0
  where (vdbs%m15_bt > -900.0 .AND. vdbs%m16_bt > -900.0) 
    vdbs%btd11 = vdbs%m15_bt - vdbs%m16_bt
  end where
  
  vdbs%dstar(:,:) = -999.0
  status = calc_dstar(vdbs%btd8, vdbs%btd11, vdbs%dstar) 
  if (status /= 0) then
    print *, "ERROR: Failed to calculate D* values: ", status
    return
  end if
    
  vdbs%ndvi(:,:) = -999.0
  where (vdbs%m07_refl > -900.0 .AND. vdbs%m05_refl > -900.0)
    vdbs%ndvi = (vdbs%m07_refl - vdbs%m05_refl) / (vdbs%m07_refl + vdbs%m05_refl)
  end where

! NDVI is calculated twice, may be one of them should be erased (deep_blue.f90 line 877)
!     where (vdbs%m01_refl > -900.0)
!        vdbs%ndvi = (vdbs%m07_refl - vdbs%m05_refl) /  &
!       &         (vdbs%m07_refl + vdbs%m05_refl)
!     end where 

! -- clean up!
  deallocate(i_val, lw_mask, stat=status)
  if (status /= 0) then
    print *, "WARNING: Failed to deallocate integer data arrays: ", status
  end if  

!   print *,'================================================================='
!   print *,'Temporarily reduce L1b size' 
!   print *,'================================================================='
!   where (vdbs%lat > 4. .or. vdbs%lat < 0 .or. vdbs%lon > 62. .or. vdbs%lon < 58)
!   ! where (vdbs%lat > -22.63 .or. vdbs%lat < -28.03 .or. vdbs%lon > 180 .or. vdbs%lon < -180)
! ! where (vdbs%lat > 90 .or. vdbs%lat < -90 .or. vdbs%lon > -161.63 .or. vdbs%lon < -161.83)
! !   where (vdbs%lat > 36.5 .or. vdbs%lat < 36. .or. vdbs%lon > -75. .or. vdbs%lon < -75.5)
!         vdbs%m01_refl =-999.0
!         vdbs%m02_refl =-999.0
!         vdbs%m03_refl =-999.0
!         vdbs%m04_refl =-999.0
!         vdbs%m05_refl =-999.0
!         vdbs%m07_refl =-999.0
!         vdbs%m08_refl =-999.0
!         vdbs%m09_refl =-999.0 
!         vdbs%m10_refl =-999.0
!         vdbs%m11_refl =-999.0
!         vdbs%m14_rad =-999.0
!         vdbs%m15_rad =-999.0
!         vdbs%m16_rad =-999.0
!         
!         vdbs%m14_bt   =-999.0
!         vdbs%m15_bt   =-999.0
!         vdbs%m16_bt   =-999.0
! 
!         vdbs%ndvi     =-999.0
!         vdbs%btd8     =-999.0
!         vdbs%btd11    =-999.0
!         vdbs%dstar    =-999.0
!   end where
!   
!   return
  
end function load_viirs_db_data_nasa


subroutine check(status, str)
  integer, intent ( in) :: status
  character(len=255)    :: str
  
  if(status /= 0) then 
    print *, str, status
  end if
end subroutine check



integer function calc_dstar(btd8, btd11, dstar) result(status)
  implicit none 
  
  real, dimension(:,:), intent(in)    ::  btd8 
  real, dimension(:,:), intent(in)    ::  btd11
  real, dimension(:,:), intent(inout) ::  dstar
  integer, dimension(2)               ::  dims 
  integer                             ::  i, j 

! -- intialize thermal constants 
  real(kind=8) , parameter   ::  A = -0.05
  real(kind=8) , parameter   ::  B = 10.0  
  real(kind=8)               ::  btd8btd11 
  
  dims  = shape(btd8)
! --calculate D* 
!   where (btd11 > -900 .AND. btd8 > -900.0)
!     dstar = exp(((btd11 - A)/ (btd8 - B))) 
!   end where

! --calculate D*  : modified to avoid potential error when (btd11(i,j) - A)/ (btd8(i,j) - B) is very low
  do i = 1, dims(1)
  do j = 1, dims(2)
    if (btd11(i,j) < -900 .or. btd8(i,j) < -900.0) cycle
    btd8btd11  = ((btd11(i,j) - A)/ (btd8(i,j) - B))
    if (btd8btd11 <= -80) dstar(i,j)=0.
    if (btd8btd11 >= 50) dstar(i,j)=exp(50.) 
    if (btd8btd11 > -80 .and. btd8btd11 < 50)dstar(i,j) = exp(((btd11(i,j) - A)/ (btd8(i,j) - B))) 
  end do
  end do

  status = 0

  return
  
end function calc_dstar


integer function rescale_array(orig_arr, res_arr, factor) result(status)
  implicit none
  
  real, dimension(:,:), intent(in)      ::  orig_arr
  real, dimension(:,:), intent(inout)   ::  res_arr
  integer, intent(in)                   ::  factor
  real, dimension(:,:), allocatable     ::  dummy_arr
  integer, dimension(2)                 ::  dims2
  real                                  ::  sum_arr
  
  integer                               ::  i, j, fi, fj
  integer                               ::  ii, jj, narr
  
  dims2 = (/size(orig_arr,1), size(orig_arr,2)/)
  if (mod(dims2(1), factor) /= 0 .OR. mod(dims2(2), factor) /= 0) then
    print *, "ERROR: Array dimensions must be evenly divisible by factor: ", factor
    status = -1
    return
  end if
  
  allocate(dummy_arr(factor,factor))
  narr  = 0

  dims2 = (/size(res_arr,1), size(res_arr,2)/) 
  do j = 1, dims2(2)
    do i = 1, dims2(1)
      
      ii = ((i-1) * factor) + 1
      jj = ((j-1) * factor) + 1
      dummy_arr = orig_arr(ii:ii+(factor-1),jj:jj+(factor-1))
    
      narr    = 0
      sum_arr = 0
        do fi =1, factor
        do fj =1, factor   
          if (dummy_arr(fi, fj) > -900) then 
            narr    = narr+1
            sum_arr = sum_arr+dummy_arr(fi, fj)
          end if
        end do 
        end do
      
!       res_arr(i,j) = sum(orig_arr(ii:ii+(factor-1),jj:jj+(factor-1))) / (factor*factor)
      if (narr > 0) res_arr(i,j)  = sum_arr/narr
    end do
  end do
  
  status = 0
  return
end function rescale_array  
  
end module

