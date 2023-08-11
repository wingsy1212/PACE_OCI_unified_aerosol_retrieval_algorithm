module viirs_db_utils
use modis_surface, only:  check_status  
private

! -- types
public  ::  viirs_gas_correction
public  ::  viirs_gas_correction_fullwv
public  ::  viirs_gas_transmittance
public  ::  viirs_calib_correction
public  ::  load_dbdt_region_table
public  ::  create_dbdt_aot550,create_dbdt_pace

! -- functions
public  ::  create_viirs_l2
public  ::  calc_minmax_ratio
public  ::  calc_smoke_mask
public  ::  calc_pyrocb_mask					
public  ::  calc_high_alt_smoke_mask	
public  ::  calc_smoke_ae_mask
public  ::  calc_gas_correction
public  ::  calc_gas_correction_fullwv
public  ::  calc_gas_transmittance
public  ::  calc_calibration_correction
public  ::  ocean_cloud_filter
public  ::  extrapolate_uv

public  ::  latitude, longitude, ocld_mask

integer, parameter  ::  i8  = selected_int_kind(3)
integer, parameter  ::  i64 = selected_int_kind(18)
integer, parameter	:: 	r64 = selected_real_kind(p=10)

type  ::  viirs_gas_correction
  
  real, dimension(:,:), allocatable   ::  m01
  real, dimension(:,:), allocatable   ::  m02
  real, dimension(:,:), allocatable   ::  m03
  real, dimension(:,:), allocatable   ::  m04
  real, dimension(:,:), allocatable   ::  m05
  real, dimension(:,:), allocatable   ::  m06
  real, dimension(:,:), allocatable   ::  m07
  real, dimension(:,:), allocatable   ::  m08
  real, dimension(:,:), allocatable   ::  m10
  real, dimension(:,:), allocatable   ::  m11
  
end type viirs_gas_correction

type  ::  viirs_gas_correction_fullwv

  real, dimension(:,:), allocatable   ::  m11

end type viirs_gas_correction_fullwv

type  ::  viirs_gas_transmittance

  real, dimension(:,:), allocatable   ::  m01
  real, dimension(:,:), allocatable   ::  m02
  real, dimension(:,:), allocatable   ::  m03
  real, dimension(:,:), allocatable   ::  m04
  real, dimension(:,:), allocatable   ::  m05
  real, dimension(:,:), allocatable   ::  m06
  real, dimension(:,:), allocatable   ::  m07
  real, dimension(:,:), allocatable   ::  m08
  real, dimension(:,:), allocatable   ::  m09
  real, dimension(:,:), allocatable   ::  m10
  real, dimension(:,:), allocatable   ::  m11

end type viirs_gas_transmittance

type  ::  viirs_calib_correction

  real   ::  m01
  real   ::  m02
  real   ::  m03
  real   ::  m04
  real   ::  m05
  real   ::  m06
  real   ::  m07
  real   ::  m08
  real   ::  m10
  real   ::  m11
  
end type viirs_calib_correction

real, dimension(:,:), allocatable     ::  latitude, longitude
integer, dimension(:,:), allocatable  ::  ocld_mask
integer, dimension(:,:), allocatable      ::  dbdt_region

contains


integer function extrapolate_uv(dbdt_aot,uvaod) result(status)
 implicit none

 real, dimension(:,:,:), intent(inout)   ::  dbdt_aot
 real, dimension(:,:,:), intent(inout)::  uvaod
 real                                 ::  wv480, wv550, wv670
 real                                 ::  dd,coefa, coefb, coefc, coefd, a0, a1, a2,f, angexp,aodrat
 integer                              ::  i, j!, ilat, ilon
 wv480  = 466.  !480.
 wv550  = 554. !550.
 wv670  = 645. ! 670.
 dd      = alog(wv480/wv670)
 
! print *, 'start UV AOD extrapolation'

 status = -1
 do i = 1, size(dbdt_aot,1)
  do j = 1, size(dbdt_aot,2)
   if (dbdt_aot(i,j,1) <= 0.0  .or. dbdt_aot(i,j,2) <= 0.0 .or. dbdt_aot(i,j,3) <= 0.0 ) cycle 

!    dbdt_aot(i,j,1)  = 2.691315521
!    dbdt_aot(i,j,2)  = 2.203356908
!    dbdt_aot(i,j,3)  = 1.800963554
   
   coefa = (log(dbdt_aot(i,j,2)) - log(dbdt_aot(i,j,1)))/(log(wv550) - log(wv480))										
   coefb = (log(dbdt_aot(i,j,3)) - log(dbdt_aot(i,j,2)))/(log(wv670) - log(wv550))										
   coefc = ((log(wv550))**2.-(log(wv480))**2.)/(log(wv550) - log(wv480))																	
   coefd = ((log(wv670))**2.-(log(wv550))**2.)/(log(wv670) - log(wv550))	

   aodrat = dbdt_aot(i,j,1) /dbdt_aot(i,j,3)
   angexp = alog(aodrat)
   angexp = -1.*angexp/dd     
   f = 0		
   if (angexp > 1) f = log(1.05+(0.10*log(dbdt_aot(i,j,2))))

   a2 = (coefa - coefb)/(coefc - coefd)										
   a1= coefb-(coefd*a2)										
   a0 = log(dbdt_aot(i,j,3)) - (a1*log(wv670)) - (a2*(log(wv670))**2.)+f

   uvaod(i,j,1) = exp(a0 + (a1*log(388.)) + (a2*(log(388.))**2))
   uvaod(i,j,2) = exp(a0 + (a1*log(354.)) + (a2*(log(354.))**2))   
  end do
 end do
   
 status = 0
 return

end function extrapolate_uv



integer function create_dbdt_pace(lat, lon, db488, db550, db670, db_lqa,db_type, dt_laot, dt_lqa, &
&                                 dbdt_aot, dbdt_qa, dbdt_alg,dtspec,lreflc_mean,dtfmf,dbdtfmf,&
&                                 dt_cldmsk,lcld_mask,lat_ori, lon_ori, dbdt_refl, sza ) result(status)
 implicit none

 real(kind=4), dimension(:,:), intent(in)::  lat, lon,lat_ori, lon_ori 
 real, dimension(:,:), intent(in)        ::  db488, db550, db670, sza
 real, dimension(:,:,:), intent(in)      ::  dt_laot,lreflc_mean
 integer, dimension(:,:), intent(in)     ::  db_lqa,db_type
 real, dimension(:,:), intent(in)        ::  dt_lqa
 real, dimension(:,:,:), intent(inout)   ::  dbdt_aot,dtspec,dbdt_refl
 integer, dimension(:,:), intent(inout)  ::  dbdt_qa,dt_cldmsk,lcld_mask
 integer, dimension(:,:), intent(inout)  ::  dbdt_alg
 real, dimension(:,:), intent(inout)     ::  dtfmf,dbdtfmf
 real, dimension(:,:), allocatable       ::  to_iof
 real, parameter                 ::  d2r = 3.14159/180.0   ! convert degrees to radians


 integer       ::  i, j, ilat, ilon
 integer       ::  region
 real          ::  meanfmf
 integer, dimension(2)      :: dims
 
! print *, 'start DB/DT AOD merging'
 status = -1
 dims = shape(sza)
 allocate(to_iof(dims(1),dims(2)), stat=status)
 if (status /= 0) then
   print *, "ERROR: Unable to allocate I/F conversion array: ", status
   stop
 end if
 
 !dbdt merged fmf temporal
 dims = shape(dbdtfmf)
 to_iof  = (cos(sza*d2r)/3.14159)
 
 dbdt_refl(:,:,1)  = lreflc_mean(9, :,:)/to_iof
 dbdt_refl(:,:,2)  = lreflc_mean(10, :,:)/to_iof 
 dbdt_refl(:,:,3)  = lreflc_mean(2, :,:)/to_iof
 dbdt_refl(:,:,4)  = lreflc_mean(3, :,:)/to_iof
 dbdt_refl(:,:,5)  = lreflc_mean(4, :,:)/to_iof
 dbdt_refl(:,:,6)  = lreflc_mean(5, :,:)/to_iof
 dbdt_refl(:,:,7)  = lreflc_mean(6, :,:)/to_iof
 dbdt_refl(:,:,8)  = lreflc_mean(7, :,:)/to_iof
 dbdt_refl(:,:,9)  = lreflc_mean(8, :,:)/to_iof
 
!      dbdt_refl(:,:,1)  = dtspec(:,:,1)
!      dbdt_refl(:,:,2)  = dtspec(:,:,2) 
!      dbdt_refl(:,:,3)  = dtspec(:,:,3)
!      dbdt_refl(:,:,4)  = dtspec(:,:,4)
!      dbdt_refl(:,:,5)  = dtspec(:,:,5)
!      dbdt_refl(:,:,6)  = dtspec(:,:,6)
!      dbdt_refl(:,:,7)  = dtspec(:,:,7)
!      dbdt_refl(:,:,8)  = dtspec(:,:,8)
!      dbdt_refl(:,:,9)  = dtspec(:,:,9) 

   
 !merge AOD
 do i = 1, size(db550,1)
  do j = 1, size(db550,2)

   dbdt_aot(i,j,:) = -9999
   dbdt_qa(i,j) = 0
   dbdt_alg(i,j) = -999   ! 1=db 2=dt
   if (lat(i,j) < -900.0 .OR. lon(i,j) < -900.0) cycle

   if (dt_lqa(i,j)>=3) dbdt_alg(i,j)= 2 !DT
   if (db_lqa(i,j)>=2) dbdt_alg(i,j)= 1 !DB
   if (dt_lqa(i,j)>=3 .and. db_lqa(i,j)>=2) then 
    dbdt_alg(i,j)= 1
    if (lat(i,j) > -30.0 .and. lat(i,j) < 45.0 .and. lon(i,j) > -96.0 .and. lon(i,j) < -58.0) dbdt_alg(i,j)= 2
    if (lat(i,j) > 45.0 .and. lon(i,j) > -20.0 .and. lon(i,j) < 5.0) dbdt_alg(i,j)= 2
    if (lat(i,j) > 35.0 .and. lon(i,j) > 5.0 .and. lon(i,j) < 36.0) dbdt_alg(i,j)= 2
    if (lat(i,j) < 3.0 .and. lon(i,j) > 5.0 .and. lon(i,j) < 52.0) dbdt_alg(i,j)= 2
    if (lat(i,j) > -10.0 .and. lat(i,j) < 30.0 .and. lon(i,j) > 82.0 .and. lon(i,j) < 112.0) dbdt_alg(i,j)= 2
    if (lat(i,j) > -10.0 .and. lat(i,j) < 19.0 .and. lon(i,j) > 112.0 .and. lon(i,j) < 155.0) dbdt_alg(i,j)= 2
    if (lat(i,j) > -30.0 .and. lat(i,j) < 15.0 .and. lon(i,j) > 145.0 .and. lon(i,j) < 160.0) dbdt_alg(i,j)= 2
   endif
   
   if (dbdt_alg(i,j) < -900.0) cycle
   if (dt_lqa(i,j)>=3 .or. db_lqa(i,j)>=2) dbdt_qa(i,j) = 3
   if (dbdt_alg(i,j) == 2) dbdt_aot(i,j,:) = dt_laot(i,j,1:3) 
   if (dbdt_alg(i,j) == 1) then 
    dbdt_aot(i,j,1) = db488(i,j)
    dbdt_aot(i,j,2) = db550(i,j) 
    dbdt_aot(i,j,3) = db670(i,j) 
   endif
   
   if (dbdt_alg(i,j) == 1) dbdtfmf(i,j)=db_type(i,j)+100.

   if (dbdt_alg(i,j) == 2) then 
     dbdtfmf(i,j)=dtfmf(i,j)
     dbdt_refl(i,j,1)  = dtspec(i,j,1)
     dbdt_refl(i,j,2)  = dtspec(i,j,2) 
     dbdt_refl(i,j,3)  = dtspec(i,j,3)
     dbdt_refl(i,j,4)  = dtspec(i,j,4)
     dbdt_refl(i,j,5)  = dtspec(i,j,5)
     dbdt_refl(i,j,6)  = dtspec(i,j,6)
     dbdt_refl(i,j,7)  = dtspec(i,j,7)
     dbdt_refl(i,j,8)  = dtspec(i,j,8)
     dbdt_refl(i,j,9)  = dtspec(i,j,9)     
   endif
   
  end do
 end do

!  do i = 1, size(db550,1)
!   do j = 1, size(db550,2)
!    if (dbdt_alg(i,j) < -900.0) dtspec(i,j,:)=-999.0
!    if (dbdt_alg(i,j) == 1) then 
! !      dtspec(i,j,1)=-999.0  !354nm
! !      dtspec(i,j,2)=-999.0 !388nm
!      dtspec(i,j,3)=lreflc_mean(2,i,j)  !480nm
!      dtspec(i,j,4)=lreflc_mean(3,i,j)  !550nm
!      dtspec(i,j,5)=lreflc_mean(4,i,j)  !670nm
!      dtspec(i,j,6)=lreflc_mean(5,i,j)  !870nm
!      dtspec(i,j,7)=lreflc_mean(6,i,j)  !1240nm
!      dtspec(i,j,8)=lreflc_mean(7,i,j)  !1.64micron
!      dtspec(i,j,9)=lreflc_mean(8,i,j)  !2.25micron
!    endif
!   end do
!  end do 
 
 !merge cloud mask
 do i = 1, size(lon_ori,1)
  do j = 1, size(lon_ori,2)

   if (lat_ori(i,j) < -900.0 .OR. lon_ori(i,j) < -900.0) cycle

   if (lat_ori(i,j) > -30.0 .and. lat_ori(i,j) < 45.0 .and. lon_ori(i,j) > -96.0 .and. lon_ori(i,j) < -58.0) &
   &  lcld_mask(i,j)= dt_cldmsk(i,j)
   if (lat_ori(i,j) > -10.0 .and. lat_ori(i,j) < 30.0 .and. lon_ori(i,j) > 82.0 .and. lon_ori(i,j) < 112.0) &
   &  lcld_mask(i,j)= dt_cldmsk(i,j)
   if (lat_ori(i,j) > -10.0 .and. lat_ori(i,j) < 19.0 .and. lon_ori(i,j) > 112.0 .and. lon_ori(i,j) < 155.0) &
   &  lcld_mask(i,j)= dt_cldmsk(i,j)
   if (lat_ori(i,j) > -30.0 .and. lat_ori(i,j) < 15.0 .and. lon_ori(i,j) > 145.0 .and. lon_ori(i,j) < 160.0) &
   &  lcld_mask(i,j)= dt_cldmsk(i,j)  
   
   if (lat_ori(i,j) > 45.0 .and. lon_ori(i,j) > -20.0 .and. lon_ori(i,j) < 5.0) lcld_mask(i,j)= dt_cldmsk(i,j)
   if (lat_ori(i,j) > 35.0 .and. lon_ori(i,j) > 5.0 .and. lon_ori(i,j) < 36.0) lcld_mask(i,j)= dt_cldmsk(i,j)
   if (lat_ori(i,j) < 3.0 .and. lon_ori(i,j) > 5.0 .and. lon_ori(i,j) < 52.0) lcld_mask(i,j)= dt_cldmsk(i,j)


  end do
 end do 
 !substitute dt_cldmsk to lcld_mask, final output will be dt_cldmsk which is DBDT merged
!  dt_cldmsk  = lcld_mask
   
 status = 0
 return

end function create_dbdt_pace

! -- Performs the actual mixing of the two data sets based on the dbdt_region array and
! --  confidence flags.
! -- Returns -1 on failure, otherwise 0.
!
! -- Over ocean, dark target data is used.  Over land, dark target or deep blue values
! --  are used based on the dbdt region map. In areas designated as mixed areas, the
! --  following table explains which value will be used as a function of QA.
!                 DB QA
!                 1   2   3
!                 ------------
!               1 FV  DB  DB
!       DT QA   2 FV  DB  DB
!               3 DT  DT  AVG
!
! --  where DB=Deep Blue, DT=Dark Target, AVG=mean of DB and DT values, FV=fill value.
!-----------------------------------------------------------------------------------------
  integer function create_dbdt_aot550(lat, lon, db_laot550, db_lqa, dt_laot550, dt_lqa, dt_oaot550, &
  &                                   dt_oqa, dbdt_aot550, dbdt_qa, dbdt_alg) result(status)
    implicit none

    real(kind=4), dimension(:,:), intent(in)        ::  lat, lon
    real, dimension(:,:), intent(in)     ::  db_laot550
    real, dimension(:,:), intent(in)     ::  dt_laot550, dt_oaot550
    integer, dimension(:,:), intent(in)     ::  db_lqa
    integer, dimension(:,:), intent(in)     ::  dt_lqa, dt_oqa
    real, dimension(:,:), intent(inout)  ::  dbdt_aot550
    integer, dimension(:,:), intent(inout)  ::  dbdt_qa
    integer, dimension(:,:), intent(inout)  ::  dbdt_alg

    integer       ::  i, j, ilat, ilon
    integer       ::  region
    
!    print *, 'start DB/DT AOD merging'
    status = -1
    do i = 1, size(db_laot550,1)
      do j = 1, size(db_laot550,2)

        dbdt_aot550(i,j) = -9999
        dbdt_qa(i,j) = 0
        dbdt_alg(i,j) = -999
        if (lat(i,j) < -900.0 .OR. lon(i,j) < -900.0) cycle

!       -- use Dark Target values over oceans
        if (dt_oaot550(i,j) /= -9999 .AND. dt_oqa(i,j) > 0) then
          dbdt_aot550(i,j)  = dt_oaot550(i,j)
          dbdt_qa(i,j)      = dt_oqa(i,j)
          dbdt_alg(i,j)     = 0       ! Dark Target used
        else
          region = get_dbdt_region(lat(i,j), lon(i,j), status)
          if (status /= 0) then
            print *, "ERROR: Failed to get DBDT region: ", status
            status = -1
            return
          end if

          select case (region)
!           -- Dark Target region only
            case (0)
              if (dt_laot550(i,j) /= -9999) then
                dbdt_aot550(i,j) = dt_laot550(i,j)
                dbdt_qa(i,j)     = dt_lqa(i,j)
                dbdt_alg(i,j)    = 0

              end if
!           -- Deep Blue region only
            case (1)
              if (db_laot550(i,j) /= -9999) then
                dbdt_aot550(i,j) = db_laot550(i,j)
                dbdt_qa(i,j)     = db_lqa(i,j)
                dbdt_alg(i,j)    = 1
              end if
!           -- Transition, best QA or mix
            case (2)
              if ((dt_lqa(i,j) == 3 .AND. db_lqa(i,j) == 3) .AND.   &
              &     (dt_laot550(i,j) /= -9999 .AND. db_laot550(i,j) /= -9999)) then
                dbdt_aot550(i,j)  = (db_laot550(i,j) + dt_laot550(i,j)) / 2
                dbdt_qa(i,j)      = dt_lqa(i,j)
                dbdt_alg(i,j)     = 2
              else if ((db_lqa(i,j) == 3 .AND. dt_lqa(i,j) < 3) .AND. db_laot550(i,j) /= -9999) then
                dbdt_aot550(i,j)  = db_laot550(i,j)
                dbdt_qa(i,j)      = db_lqa(i,j)
                dbdt_alg(i,j)     = 1
              else if ((dt_lqa(i,j) == 3 .AND. db_lqa(i,j) < 3) .AND. dt_laot550(i,j) /= -9999) then
                dbdt_aot550(i,j)  = dt_laot550(i,j)
                dbdt_qa(i,j)      = dt_lqa(i,j)
                dbdt_alg(i,j)     = 0
              else if ((db_lqa(i,j) == 2 .AND. dt_lqa(i,j) <= 2) .AND. db_laot550(i,j) /= -9999) then
                dbdt_aot550(i,j)  = db_laot550(i,j)
                dbdt_qa(i,j)      = db_lqa(i,j)
                dbdt_alg(i,j)     = 1
              end if

!           -- No region information in table.
            case (-999)
              cycle

!           -- Uh oh.
            case default
              print *, "ERROR: Invalid DBDT region index specified: ", dbdt_region(i,j)
              status = -1
              return

            end select
          end if
!          print *, 'final: ', i, j, dbdt_aot550(i,j), dbdt_qa(i,j), dbdt_alg(i,j)
        end do
      end do
      
      status = 0
      return

  end function create_dbdt_aot550
 
 ! -- Load input table describing mixing regions, dbdt_file, into dbdt_region array 
! --  for month month.
! -- Returns -1 on failure, otherwise 0.
!-----------------------------------------------------------------------------------------
  integer function load_dbdt_region_table(dbdt_file, month) result(status)
    implicit none

    include 'hdf.inc'
    include 'dffunc.inc'
        
    character(len=*), intent(in)      ::  dbdt_file
    integer, intent(in)               ::  month

    character(len=255)                ::  sds_name
    integer                           ::  sd_id, sds_index, sds_id
    integer, dimension(2)             ::  start2, stride2, edges2

    character(len=2)                  ::  str_month
    logical                           ::  file_exists

    character(len=255), parameter     ::  func_name = "load_dbdt_region_table"

    status = -1

!   -- check input parameters make sense and files exist.
    if (month > 12 .OR. month < 1) then
      print *, trim(func_name) // ": Invalid month specified: ", month
      status = -1
      return
    end if

    inquire(file=trim(dbdt_file), exist=file_exists, iostat=status)
    if (status /= 0) then
      print *, "ERROR: Unable to inquire about dbdt file: ", status
      status = -1
      return
    end if
    if (file_exists .eqv. .false.) then
      print *, "ERROR: DBDT file does not exist. Please check the config file."
      status = -1
      return
    end if

!   -- allocate array and open, read file.
                allocate(dbdt_region(1440, 720), stat=status)
                if (status /= 0) then
      print *, "ERROR: Unable to allocate dbdt array: ", status
      return
    end if 
  
                 sd_id = sfstart(dbdt_file, DFACC_READ)
    if (sd_id == FAIL ) then
        print *,"ERROR: failed to start SDS interface on dbdt file: ", sd_id
      status = -1
      return
    end if

    write(str_month, '(I0.2)') month
    sds_name = str_month // '_dbdt_regions'
    sds_index = sfn2index(sd_id, sds_name)
    if (sds_index == FAIL) then
      print *, "ERROR: Unable to find index of "//trim(sds_name)//" SDS: ", sds_index
      status = -1
      return
    end if
    sds_id = sfselect(sd_id, sds_index)
    if (sds_id == FAIL) then
      print *, "ERROR: Unable to select "//trim(sds_name)//" SDS: ", sds_id
      stop
    end if

    start2  = (/0,0/)
    stride2 = (/1,1/)
    edges2  = (/1440,720/)
    status = sfrdata(sds_id, start2, stride2, edges2, dbdt_region)
    if (status == FAIL) then
      print *, "ERROR: Unable to read "//trim(sds_name)//" SDS: ", status
      status = -1
      return
    end if

    status = sfendacc(sds_id)
    if (status == FAIL) then
      print *, "ERROR: Unable to end access to "//trim(sds_name)//" SDS: ", status
      return
    end if

    status = sfend(sd_id)
    if (status /= 0) then
      print *, "ERROR: Unable to close dbdt file: ", status
      return
    end if

    status = 0
    return
  end function load_dbdt_region_table
   
! -- Helper function for converting a lat, lon into the appropriate index into
! --  dbdt_regions and returns value from dbdt_regions array.
! -- Returns -1 on failure, otherwise 0.
!-----------------------------------------------------------------------------------------
  integer function get_dbdt_region(lat, lon, status) result(region)
    implicit none

    real, intent(in)        ::  lat, lon
    integer, intent(inout)  ::  status

    integer                 ::  ilat, ilon

    status = -1
    region = -999
    ilat = lat2index(lat, status)
    if (status /= 0) then
      print *, "ERROR: Unable to convert latitude to index into DBDT region table: ", lat
      status = -1
      return
    end if
    ilon = lon2index(lon, status)
    if (status /= 0) then
      print *, "ERROR: Unable to convert longitude to index into DBDT region table: ", lon
      status = -1
      return
    end if

    region = dbdt_region(ilon, ilat)
    status = 0
    return

  end function get_dbdt_region

! -- Helper function for get_dbdt_region. Converts longitude value into an index into
! --  dbdt_regions array.
! -- Returns -1 on failure, otherwise 0.
!-----------------------------------------------------------------------------------------
  integer function lon2index(lon, status) result(index)
    implicit none

    real, intent(in)        ::  lon
    integer, intent(inout)  ::  status

    status = -1
    index = -999

    if (lon < -180.0 .OR. lon > 180.0) then
      print *, "ERROR: Longitude is out of range: ", lon
      status = -1
      return
    end if

    index = floor((lon + 180.0) * 4.0) + 1
    status = 0
    return

  end function lon2index

! -- Helper function for get_dbdt_region. Converts latitude value into an index into
! --  dbdt_regions array.
! -- Returns -1 on failure, otherwise 0.
!-----------------------------------------------------------------------------------------
  integer function lat2index(lat, status) result(index)
    implicit none

    real, intent(in)        ::  lat
    integer, intent(inout)  ::  status

    status = -1
    index = -999

    if (lat < -90.0 .OR. lat > 90.0) then
      print *, "ERROR: Latitude is out of range: ", lat
      status = -1
      return
    end if

    index = floor((lat + 90.0) * 4.0) + 1
    status = 0
    return

  end function lat2index









                     
integer function calc_btd8(m14_rad, m15_rad, btd8) result(status)
  implicit none
  
  real, dimension(:,:), intent(in)    ::  m14_rad
  real, dimension(:,:), intent(in)    ::  m15_rad
  real, dimension(:,:), intent(inout) ::  btd8

  where(m14_rad > -900.0 .AND. m15_rad > -900.0) 
    btd8=0.0144*(1/(8.6e-6*log(25.354e+2/m14_rad+1)))- &
    &          0.0144*(1/(11.02e-6*log(7.3388e+2/m15_rad+1)))
  end where
  
  status = 0
  return
          
end function calc_btd8

integer function calc_btd11(m15_rad, m16_rad, btd11) result(status)
  implicit none
  
  real, dimension(:,:), intent(in)    ::  m15_rad
  real, dimension(:,:), intent(in)    ::  m16_rad
  real, dimension(:,:), intent(inout) ::  btd11
  
  where(m15_rad > -900.0 .AND. m16_rad > -900.0) 
    btd11=0.0144*(1/(11.02e-6*log(7.3388e+2/m15_rad+1)))- &
    &       0.0144*(1/(12.02e-6*log(4.7534e+2/m16_rad+1)))      
  end where
  
  status = 0
  return
  
end function calc_btd11


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Creates L2 output file and fills !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer function put_attribute_l2(name, dset_name, val, id) result(status)
  use h5lt
  use hdf5
  character(len=255)                      ::  name, val, dset_name
  integer(hid_t)                          ::  id

  call h5ltset_attribute_string_f(id, dset_name, trim(name), trim(val), status)
  if (status < 0) then
    print *, "ERROR: Unable to create "//trim(name)//" attribute on "//trim(dset_name)//": ", status
    stop
  end if  
  
end function put_attribute_l2

integer function put_attribute_l2_float(name, dset_name, val, id, size) result(status)
  use h5lt
  use hdf5
  character(len=255)                      ::  name, dset_name
  integer(hid_t)                          ::  id
  real, dimension(1)                      ::  val
  integer(size_t)                         ::  size

  call h5ltset_attribute_float_f(id, dset_name, name, val, size, status)
  if (status /= 0) then
    print *, "ERROR: Failed to create fill value attribute on "//trim(dset_name)//": ", status
  return
  end if 
end function put_attribute_l2_float

integer function put_attribute_l2_float_arr(name, dset_name, val, id, size) result(status)
  use h5lt
  use hdf5
  character(len=255)                      ::  name, dset_name
  integer(hid_t)                          ::  id
  real, dimension(2)                      ::  val
  integer(size_t)                         ::  size

  call h5ltset_attribute_float_f(id, dset_name, name, val, size, status)
  if (status /= 0) then
    print *, "ERROR: Failed to create fill value attribute on "//trim(dset_name)//": ", status
  return
  end if 
end function put_attribute_l2_float_arr

integer function put_attribute_l2_int_arr(name, dset_name, val, id, size) result(status)
  use h5lt
  use hdf5
  character(len=255)                      ::  name, dset_name
  integer(hid_t)                          ::  id
  integer, dimension(2)                   ::  val
  integer(size_t)                         ::  size

  call h5ltset_attribute_int_f(id, dset_name, name, val, size, status)
  if (status /= 0) then
    print *, "ERROR: Failed to create fill value attribute on "//trim(dset_name)//": ", status
  return
  end if 
end function put_attribute_l2_int_arr

subroutine make_dimension(dimtype, dset_name, dset,dset2, hid, attr_val,dimn_name,&
    &   lname,unit_name, unit_val,did )
    
  use h5lt
  use hdf5
  use h5ds
        implicit none

  character(len=255) , intent(in)         ::  dset_name,dimtype
  real, dimension(:), intent(in)          ::  dset
  integer, dimension(:), intent(in)       ::  dset2
  integer(hid_t)                          ::  hid, dspace_id
  integer(hid_t), intent(out)             ::  did
  character(len=255)                      ::  attr_val, dimn_name,lname,unit_name, unit_val
  integer(hsize_t), dimension(1)          ::  dims1
  integer                                 ::  status
  
  if (dimtype .eq. "real") then
     dims1 = (/size(dset,1)/)  
  endif
  if (dimtype .eq. "integer") then
     dims1 = (/size(dset2,1)/)  
  endif

  call h5screate_simple_f(1, dims1,  dspace_id, status)
  if (status /= 0) then
    print *, "ERROR: Unable to create dataspace for "//trim(dset_name)//": ", status
    stop
  end if
 
  call h5dcreate_f(hid, dset_name, H5T_NATIVE_REAL,  dspace_id, did, status)
  if (status /= 0) then
    print *, "ERROR: Unable to create dataset for "//trim(dset_name)//" : ", status
    stop
  end if

  if (dimtype .eq. "real") then
    call h5dwrite_f(did, H5T_NATIVE_REAL, dset, dims1, status)
    if (status /= 0) then
      print *, "ERROR: Unable to write dataset "//trim(dset_name)//" : ", status
      stop
    end if
  endif
 
  if (dimtype .eq. "integer") then
    call h5dwrite_f(did, H5T_NATIVE_integer, dset2, dims1, status)
    if (status /= 0) then
      print *, "ERROR: Unable to write dataset "//trim(dset_name)//" : ", status
      stop
    end if
  endif
  
  call h5sclose_f(dspace_id, status)
  if (status /= 0) then
    print *, "ERROR: Unable to close dataspace for "//trim(dset_name)//" : ", status
    stop
  end if
    
  status = put_attribute_l2(lname ,dset_name,attr_val , hid)
  status = put_attribute_l2(unit_name ,dset_name,unit_val , hid)

  ! Make it into a dimension scale
  call H5DSset_scale_f(did, status, dimn_name)
  if (status < 0) then
    print *, "ERROR: Unable to set "//trim(dset_name)//" as a dimension scale: ", status
    stop
  end if
  
  return
end subroutine make_dimension

subroutine make_dset_i(dimtype, dset_name, dset,dset2,dset3, hid, attr_val,&
    &   lname,unit_name, unit_val,did,dimid1,  dimid2,dimid3)
    
  use h5lt
  use hdf5
  use h5ds
  implicit none

  character(len=255) , intent(in)         ::  dset_name,dimtype
  integer, dimension(:), intent(in)       ::  dset
  integer, dimension(:,:), intent(in)     ::  dset2
  integer, dimension(:,:,:), intent(in)   ::  dset3
  integer(hid_t)                          ::  hid, dspace_id,dimid1, dimid2,dimid3
  integer(hid_t), intent(out)             ::  did
  character(len=255)                      ::  attr_val, lname,unit_name, unit_val,geo_name, geo_val
  integer(hsize_t), dimension(1)          ::  dims1
  integer(hsize_t), dimension(2)          ::  dims2
  integer(hsize_t), dimension(3)          ::  dims3
  integer                                 ::  status
  
  if (dimtype .eq. "1") then
     dims1 = (/size(dset,1)/)  
  endif
  if (dimtype .eq. "2") then
     dims2 = (/size(dset2,1), size(dset2,2)/)
  endif
  if (dimtype .eq. "3") then
     dims3 = (/size(dset3,1), size(dset3,2), size(dset3,3)/)
  endif

  if (dimtype .eq. "1") then
   call h5screate_simple_f(1, dims1,  dspace_id, status)
  endif
  if (dimtype .eq. "2") then
   call h5screate_simple_f(2, dims2,  dspace_id, status)
  endif
  if (dimtype .eq. "3") then
   call h5screate_simple_f(3, dims3,  dspace_id, status)
  endif  
     
  if (status /= 0) then
    print *, "ERROR: Unable to create dataspace for "//trim(dset_name)//": ", status
    stop
  end if
 
  call h5dcreate_f(hid, dset_name, H5T_NATIVE_integer,  dspace_id, did, status)
  if (status /= 0) then
    print *, "ERROR: Unable to create dataset for "//trim(dset_name)//" : ", status
    stop
  end if

  if (dimtype .eq. "1") then
    call h5dwrite_f(did, H5T_NATIVE_integer, dset, dims1, status)
  endif
  if (dimtype .eq. "2") then
    call h5dwrite_f(did, H5T_NATIVE_integer, dset2, dims2, status)
  endif
  if (dimtype .eq. "3") then
    call h5dwrite_f(did, H5T_NATIVE_integer, dset3, dims3, status)
  endif  
        
  if (status /= 0) then
    print *, "ERROR: Unable to write dataset "//trim(dset_name)//" : ", status
    stop
  end if
  
  call h5sclose_f(dspace_id, status)
  if (status /= 0) then
    print *, "ERROR: Unable to close dataspace for "//trim(dset_name)//" : ", status
    stop
  end if
    
  status = put_attribute_l2(lname ,dset_name,attr_val , hid)
  status = put_attribute_l2(unit_name ,dset_name,unit_val , hid)

  ! Attach dimension scales
  call H5DSattach_scale_f(did, dimid1, 1, status)
  if (status < 0) then
    print *, "ERROR: Unable to set dimension scale for "//trim(dset_name)//": ", status
    stop
  end if
  
  if (dimtype .eq. "2" .or. dimtype .eq. "3" ) then
    call H5DSattach_scale_f(did, dimid2, 2, status)
    if (status < 0) then
      print *, "ERROR: Unable to set dimension scale for "//trim(dset_name)//": ", status
      stop
    end if !dsid_lat , dsid_at
  endif

  if (dimtype .eq. "3" ) then
    call H5DSattach_scale_f(did, dimid3, 3, status)
    if (status < 0) then
      print *, "ERROR: Unable to set dimension scale for "//trim(dset_name)//": ", status
      stop
    end if !dsid_lat , dsid_lat
  endif
    
  geo_name = "coordinates"
  geo_val  = "Longitude Latitude"
  call h5ltset_attribute_string_f(hid, trim(dset_name), trim(geo_name), trim(geo_val), status)
  if (status < 0) then
    print *, "ERROR: Unable to create "//trim(geo_name)//" attribute on "//trim(dset_name)//": ", status
    stop
  end if
     
  return
end subroutine make_dset_i

subroutine make_dset_f(dimtype, dset_name, dset,dset2,dset3, hid, attr_val,&
    &   lname,unit_name, unit_val,did,dimid1,  dimid2,dimid3)
    
  use h5lt
  use hdf5
  use h5ds
  implicit none

  character(len=255) , intent(in)         ::  dset_name,dimtype
  real, dimension(:), intent(in)          ::  dset
  real, dimension(:,:), intent(in)        ::  dset2
  real, dimension(:,:,:), intent(in)      ::  dset3
  integer(hid_t)                          ::  hid, dspace_id,dimid1, dimid2,dimid3
  integer(hid_t), intent(out)             ::  did
  character(len=255)                      ::  attr_val, lname,unit_name, unit_val,geo_name, geo_val
  integer(hsize_t), dimension(1)          ::  dims1
  integer(hsize_t), dimension(2)          ::  dims2
  integer(hsize_t), dimension(3)          ::  dims3
  integer                                 ::  status,skip_cor
  
  skip_cor  = 0
  if (dset_name .eq. "/Longitude") skip_cor = 1
  if (dset_name .eq. "/Latitude") skip_cor = 1
  
  if (dimtype .eq. "1") then
     dims1 = (/size(dset,1)/)  
  endif
  if (dimtype .eq. "2") then
     dims2 = (/size(dset2,1), size(dset2,2)/)
  endif
  if (dimtype .eq. "3") then
     dims3 = (/size(dset3,1), size(dset3,2), size(dset3,3)/)
  endif

  if (dimtype .eq. "1") then
   call h5screate_simple_f(1, dims1,  dspace_id, status)
  endif
  if (dimtype .eq. "2") then
   call h5screate_simple_f(2, dims2,  dspace_id, status)
  endif
  if (dimtype .eq. "3") then
   call h5screate_simple_f(3, dims3,  dspace_id, status)
  endif  
     
  if (status /= 0) then
    print *, "ERROR: Unable to create dataspace for "//trim(dset_name)//": ", status
    stop
  end if
 
  call h5dcreate_f(hid, dset_name, H5T_NATIVE_REAL,  dspace_id, did, status)
  if (status /= 0) then
    print *, "ERROR: Unable to create dataset for "//trim(dset_name)//" : ", status
    stop
  end if

  if (dimtype .eq. "1") then
    call h5dwrite_f(did, H5T_NATIVE_REAL, dset, dims1, status)
  endif
  if (dimtype .eq. "2") then
    call h5dwrite_f(did, H5T_NATIVE_REAL, dset2, dims2, status)
  endif
  if (dimtype .eq. "3") then
    call h5dwrite_f(did, H5T_NATIVE_REAL, dset3, dims3, status)
  endif  
        
  if (status /= 0) then
    print *, "ERROR: Unable to write dataset "//trim(dset_name)//" : ", status
    stop
  end if
  
  call h5sclose_f(dspace_id, status)
  if (status /= 0) then
    print *, "ERROR: Unable to close dataspace for "//trim(dset_name)//" : ", status
    stop
  end if
    
  status = put_attribute_l2(lname ,dset_name,attr_val , hid)
  status = put_attribute_l2(unit_name ,dset_name,unit_val , hid)

  ! Attach dimension scales
  call H5DSattach_scale_f(did, dimid1, 1, status)
  if (status < 0) then
    print *, "ERROR: Unable to set dimension scale for "//trim(dset_name)//": ", status
    stop
  end if
  
  if (dimtype .eq. "2" .or. dimtype .eq. "3" ) then
    call H5DSattach_scale_f(did, dimid2, 2, status)
    if (status < 0) then
      print *, "ERROR: Unable to set dimension scale for "//trim(dset_name)//": ", status
      stop
    end if !dsid_lat , dsid_at
  endif

  if (dimtype .eq. "3" ) then
    call H5DSattach_scale_f(did, dimid3, 3, status)
    if (status < 0) then
      print *, "ERROR: Unable to set dimension scale for "//trim(dset_name)//": ", status
      stop
    end if !dsid_lat , dsid_lat
  endif
  
  if (skip_cor .eq. 0) then   
    geo_name = "coordinates"
    geo_val  = "Longitude Latitude"
    call h5ltset_attribute_string_f(hid, trim(dset_name), trim(geo_name), trim(geo_val), status)
    if (status < 0) then
      print *, "ERROR: Unable to create "//trim(geo_name)//" attribute on "//trim(dset_name)//": ", status
      stop
    end if
  end if
     
  return
end subroutine make_dset_f


integer function create_viirs_l2(output_file, time, lat, lon, sza, vza, raa, &
&         sca, aot550, aot550_best, naot550, aot412, aot488, aot670, aot550_sd, &
&         ae, ae_best, ssa, qa_flag, alg, oalg, oaot,oaot550, oaot550_best, &
&         oae, oae_best, ofmf, ofmf_best, osd550, onaot550, &
&         oss, oqa_flag, omodel_flag, caot550, caot550_best, cae, cae_best, n_valid_pixels, &
&         ws, oz, wv, wd, ndvi_avg, rcndvi_avg, ndvi, sr_avg, dstar_avg, btd11_avg, &
&         turbid_res, elev_avg,  ltype, ctype,elev_avg_land,&
&         elev_avg_ocean, smoke_count,oreflc_mean, lreflc_mean,hires,alg_old,alg_old2,platform,&
&         opt_err,DBDTaot,DBDTqa, DBDTflag,dbdtspec,dbdtfmf) result(status)

  use h5lt
  use hdf5
  use h5ds
    
!   use viirs_ocean_aero, only: viirs_ocean_cell_output
  
  implicit none
  
  character(len=*), intent(in)            ::  output_file, platform
  real(kind=8), dimension(:), intent(in) ::  time
  real, dimension(:,:), intent(in)        ::  lat
  real, dimension(:,:), intent(in)        ::  lon
  real, dimension(:,:), intent(in)        ::  sza
  real, dimension(:,:), intent(in)        ::  vza
  real, dimension(:,:), intent(in)        ::  raa
  real, dimension(:,:), intent(in)        ::  sca
  real, dimension(:,:), intent(in)        ::  aot550
  real, dimension(:,:), intent(in)        ::  aot550_best
  integer, dimension(:,:), intent(in)     ::  naot550
  real, dimension(:,:), intent(in)        ::  aot412
  real, dimension(:,:), intent(in)        ::  aot488
  real, dimension(:,:), intent(in)        ::  aot670
  real, dimension(:,:), intent(in)        ::  opt_err
  
  real, dimension(:,:), intent(in)        ::  aot550_sd
  real, dimension(:,:), intent(in)        ::  ae
  real, dimension(:,:), intent(in)        ::  ae_best
  real, dimension(:,:,:), intent(in)      ::  ssa
  integer, dimension(:,:), intent(in)     ::  qa_flag
  integer, dimension(:,:), intent(in)     ::  alg, alg_old,alg_old2
  integer, dimension(:,:), intent(in)     ::  oalg
  real, dimension(:,:,:), intent(in)      ::  oaot
  real, dimension(:,:), intent(in)        ::  oaot550_best,oaot550
  real, dimension(:,:),  intent(in)       ::  oae, oae_best
  real, dimension(:,:),  intent(in)       ::  ofmf, ofmf_best
  real, dimension(:,:), intent(in)        ::  osd550
  integer, dimension(:,:), intent(in)     ::  onaot550
  real, dimension(:,:), intent(in)        ::  oss
  integer, dimension(:,:), intent(in)     ::  oqa_flag
  integer, dimension(:,:), intent(in)     ::  omodel_flag, ltype, ctype
  real, dimension(:,:), intent(in)        ::  caot550           ! combined land/ocean AOT550
  real, dimension(:,:), intent(in)        ::  caot550_best      ! combined land/ocean AOT550, best estimate
  real, dimension(:,:), intent(in)        ::  cae               ! combined land/ocean AE
  real, dimension(:,:), intent(in)        ::  cae_best          ! combined land/ocean AE, best estimate
  integer, dimension(:,:), intent(in)     ::  n_valid_pixels
  real, dimension(:,:), intent(in)        ::  ws              ! wind speed
  real, dimension(:,:), intent(in)        ::  oz              ! ozone
  real, dimension(:,:), intent(in)        ::  wv              ! precipitable water vapor
  real, dimension(:,:), intent(in)        ::  wd              ! wind direction
  real, dimension(:,:), intent(in)        ::  ndvi_avg
  real, dimension(:,:), intent(in)        ::  rcndvi_avg
  real, dimension(:,:), intent(in)        ::  elev_avg,elev_avg_land,elev_avg_ocean  ! Cell average elevation
!   real, dimension(:,:), intent(in)        ::  chl_avg         ! Cell average climatological Chl
  real, dimension(:,:,:), intent(in)      ::  sr_avg          ! surface reflectivity
  real, dimension(:,:), intent(in)        ::  dstar_avg
  real, dimension(:,:), intent(in)        ::  btd11_avg
  integer, dimension(:,:), intent(in)     ::  smoke_count
  real, dimension(:,:), intent(in)        ::  turbid_res 
  real, dimension(:,:), intent(in)        ::  ndvi 
  real, dimension(:,:,:), intent(in)      ::  oreflc_mean, lreflc_mean
  real, dimension(:,:,:), intent(in)      ::  DBDTaot,dbdtspec
  real, dimension(:,:), intent(in)        :: dbdtfmf
  integer, dimension(:,:), intent(in)     ::  DBDTqa, DBDTflag
  

  real, dimension(1)                      ::  rdum1
  real, dimension(2,2)                    ::  rdum2
  real, dimension(2,2,2)                  ::  rdum3
  integer, dimension(1)                   ::  idum1
  integer, dimension(2,2)                 ::  idum2  
  integer, dimension(2,2,2)               ::  idum3

  character(len=255)                      ::  dset_name
  character(len=255)                      ::  attr_name, unit_name,units,lname,fname, valname
  character(len=255)                      ::  grp_name

  integer(hid_t)                          ::  h5f_id
  integer(hid_t)                          ::  h5g_id
  integer(hid_t)                          ::  dset_id
  integer(hid_t)                          ::  dspace_id
  integer(hsize_t), dimension(1)          ::  dims1
  integer(hsize_t), dimension(2)          ::  dims2
  integer(hsize_t), dimension(3)          ::  dims3  
  integer(hid_t)                          ::  dsid_lbands, dsid_obands, dsid_lon 
  integer(hid_t)                          ::  dsid_lat , dsid, dumid,dsid_oref,dsid_dbdtref
                                              ! Dimension scale IDs
  character(len=255)                      ::  dimn_lbands, dimn_obands, dimn_lon, dimn_lat,dimn_name,dimtype  
                                              ! Dimension scale names                                              
  
  integer(size_t)                         :: attr_size,iattr_size,val_size
  integer, dimension(:), allocatable      :: nxtrack_arr, natrack_arr,ndbdtwl_arr
  integer                                 :: nxt,nat, nx, na
  integer                                 :: dum_idx, ndim
  real(kind=8), dimension(:,:), allocatable  ::  time_arr

!   TYPE(C_PTR) :: f_ptr

  real, dimension(1)                      ::  fv
  real, dimension(2)                      ::  val_range
  real, dimension(9)                      ::  dbdt_wl
  integer, dimension(2)                   ::  ival_range 
  integer, dimension(1)                   ::  ifv,zerofv
  real*8, dimension(1)                    ::  dfv
  character(len=140)                      ::  long_name
  character(len=255)                      ::  attr_val, unit_val,attr_val2
  real(kind=r64), dimension(:), allocatable, target :: tmp_target
  
  real, dimension(1)                      ::  nxtrack     !number of cells across track
  real, dimension(1)                      ::  natrack     !number of cells along track
  real, dimension(3)                      ::  land_bands  !wavelengths of spectral AOD in land retrieval
  real, dimension(7)                      ::  ocean_bands !wavelengths of spectral AOD in ocean retrieval
  real, dimension(8)                      ::  ref_bands   !wavelengths of spectral reflectance

  real, dimension(:,:,:), allocatable     ::  aot

  logical                                 :: extended_output, hires, normal_output
  extended_output = .false.
  
  normal_output = .true.
  if (hires) normal_output = .false.


  natrack = size(lat,2) 
  nxtrack = size(lat,1)
  na = size(lat,2) 
  nx = size(lat,1)
  land_bands  = [412., 488., 670.]
  ocean_bands = [488.0,	550.0,	670.0,	865.0,	1240.0,	1640.0,	2250.0]
  ref_bands =   [412.0, 488.0,	550.0,	670.0,	865.0,	1240.0,	1640.0,	2250.0]

  call h5open_f (status)
  if (status < 0) then 
    print *, "ERROR: Unable to start HDF interface: ", status
    return
  end if
  
  call h5fcreate_f(trim(output_file), H5F_ACC_TRUNC_F, h5f_id, status)
  if (status < 0) then 
    print *, "ERROR: Unable to create VIIRS DB output file: ", status
    print *, "File: ", trim(output_file)
    return
  end if
  
  unit_name = "units"
  unit_val = "1" 
  lname = "long_name" 
  fname = "_FillValue" 
  valname = "valid_range" 
  fv = (/-999.0/)
  attr_size = size(fv, 1)
  ifv = (/-999/)
  iattr_size = size(ifv, 1)
  dfv = (/-999.0/)
  zerofv  = 0
  val_range =(/-999.0, -999.0/)
  ival_range =(/-999, -999/)
  val_size = size(val_range, 1)
  
!Global attribute
  dset_name = "/"

  attr_name = "title"
  attr_val = "VIIRS Deep Blue aerosol (VNPAERDB)" 
  status = put_attribute_l2( attr_name ,dset_name,attr_val , h5f_id)
 
  attr_name = "processing_level"
  attr_val = "L2"
  status = put_attribute_l2( attr_name ,dset_name,attr_val , h5f_id) 

  attr_name = "cdm_data_type"
  attr_val = "swath"
  status = put_attribute_l2( attr_name ,dset_name,attr_val , h5f_id) 

  attr_name = "keywords_vocabulary"
  attr_val = "NASA Global Change Master Directory (GCMD) Science Keywords"
  status = put_attribute_l2( attr_name ,dset_name,attr_val , h5f_id) 

  attr_name = "Keywords"    !very weird error. 'keywords' doesn't work on SIPS machine,
  ! but 'Keywords' works. Couldn't find the reason (W.Kim Jan 2018)
!   attr_val = "aerosol optical depth thickness angstrom exponent land ocean deep blue viirs"
  attr_val = "EARTH SCIENCE > ATMOSPHERE > AEROSOLS > AEROSOL OPTICAL DEPTH/THICKNESS > ANGSTROM EXPONENT, &
  &EARTH SCIENCE > ATMOSPHERE > AEROSOLS > AEROSOL OPTICAL DEPTH/THICKNESS"
  
  status = put_attribute_l2( attr_name ,dset_name,attr_val , h5f_id) 
  
  attr_name = "license"
  attr_val = "http://science.nasa.gov/earth-science/earth-science-data/data-information-policy/"
  status = put_attribute_l2( attr_name ,dset_name,attr_val , h5f_id) 

  attr_name = "stdname_vocabulary"
  attr_val = "NetCDF Climate and Forecast (CF) Metadata Convention"
  status = put_attribute_l2( attr_name ,dset_name,attr_val , h5f_id) 

  attr_name = "naming_authority"
  attr_val = "gov.nasa.gsfc.sci.atmos"
  status = put_attribute_l2( attr_name ,dset_name,attr_val , h5f_id)   

  !Set dimension scales
  !land_bands
    dset_name   = "Land_Bands"  
    attr_val    = "Wavelengths of spectral AOD in land retrieval"
    dimn_name   = "dim_land_bands"
    dimtype     = "real"
    attr_val2   = "nm"
    call make_dimension(dimtype,dset_name, land_bands,nxtrack_arr, h5f_id, attr_val,&
                & dimn_name,lname,unit_name, attr_val2,dsid_lbands)

!   if (normal_output) then 
  !ocean_bands
    dset_name   = "Ocean_Bands"  
    attr_val    = "Wavelengths of spectral AOD in ocean retrieval"
    dimn_name   = "dim_ocean_bands"
    dimtype     = "real"
    attr_val2   = "nm"
    call make_dimension(dimtype,dset_name, ocean_bands,nxtrack_arr, h5f_id, attr_val,&
                & dimn_name,lname,unit_name, attr_val2, dsid_obands)
!   endif !normal_output

  !reflectance_bands
    dset_name   = "Reflectance_Bands"  
    attr_val    = "Wavelengths of spectral reflectance"
    dimn_name   = "dim_reflectance_bands"
    dimtype     = "real"
    attr_val2   = "nm"
    call make_dimension(dimtype,dset_name, ref_bands,nxtrack_arr, h5f_id, attr_val,&
                & dimn_name,lname,unit_name, attr_val2, dsid_oref)

  !reflectance_bands
    dset_name   = "DBDT Reflectance_Bands"  
    attr_val    = "Wavelengths of spectral reflectance"
    dimn_name   = "dim_reflectance_bands"
    dimtype     = "real"
    attr_val2   = "nm"
    dbdt_wl=(/354, 388, 480, 550, 670, 870, 1240, 1640, 2250/)
    call make_dimension(dimtype,dset_name, dbdt_wl,ndbdtwl_arr,& 
          & h5f_id, attr_val, dimn_name,lname,unit_name, attr_val2, dsid_dbdtref)                
  !nxtrack  
    allocate(nxtrack_arr(nx) , stat=status)    
    if (status /= 0) then
      print *, "ERROR: Failed to allocate array for spectral nxtrack_arr: ", status
      return
    end if 
    nxtrack_arr(:) = 0    
    do nxt  = 1, nx
       nxtrack_arr(nxt)  = nxt
    end do 
    dset_name   = "Idx_Xtrack"  
    attr_val    = "Index of cells across track"
    dimn_name   = "dim_lon"
    dimtype     = "integer"    
    call make_dimension(dimtype,dset_name,ocean_bands, nxtrack_arr, h5f_id, attr_val,&
                & dimn_name,lname,unit_name, unit_val,dsid_lon)

  !natrack
    allocate(natrack_arr(na) , stat=status)    
    if (status /= 0) then
      print *, "ERROR: Failed to allocate array for spectral natrack_arr: ", status
      return
    end if 
    natrack_arr(:) = 0  
    do nat  = 1, na
       natrack_arr(nat)  = nat
    end do 
    dset_name   = "Idx_Atrack"  
    attr_val    = "Index of cells along track"
    dimn_name   = "dim_lat"
    dimtype     = "integer"    
    call make_dimension(dimtype,dset_name,ocean_bands, natrack_arr, h5f_id, attr_val,&
                & dimn_name,lname,unit_name, unit_val,dsid_lat)

! -- since this dataset uses a larger integer, we have to drop down to the full
! -- HDF5 library API rather than the Lite API. 
    if (normal_output) then 
    dims1 = (/size(time,1)/)
    allocate(tmp_target(dims1(1)), stat=status)
    tmp_target(:) = time

    dset_name = "/Scan_Start_Time"
    dims2 = (/size(lon,1), size(lon,2)/)
    allocate(time_arr(dims2(1),dims2(2)), stat=status)
    time_arr(:,:) = -999.0    
    do nxt  = 1, nx
    do nat  = 1, na
       time_arr(nxt, nat)  = tmp_target(nat)
!        print * , nxt, nat
    end do
    end do     

    if (status /= 0) then
      print *, "ERROR: Failed to allocate tmp target array for "//trim(dset_name)//": ", status
      return
    end if
 
    call h5screate_simple_f(2, dims2 , dspace_id, status)
    if (status /= 0) then
      print *, "ERROR: Failed to create dataspace for "//trim(dset_name)//": ", status
      return
    end if
  
    call H5Dcreate_f(h5f_id, dset_name, h5kind_to_type(r64,H5_real_KIND),  dspace_id, dset_id, status)
    if (status /= 0) then
      print *, "ERROR: Failed to create dataset for "//trim(dset_name)//": ", status
      return
    end if
    
    call h5dwrite_f(dset_id, h5kind_to_type(r64,H5_real_KIND), time_arr, dims2, status)
    if (status /= 0) then
      print *, "ERROR: Failed to write dataset for "//trim(dset_name)//": ", status
      return
    end if  

    call H5DSattach_scale_f(dset_id, dsid_lat, 1, status)
    if (status < 0) then
      print *, "ERROR: Unable to set dimension scale for "//trim(dset_name)//": ", status
      stop
    end if
    
    call H5DSattach_scale_f(dset_id, dsid_lon, 2, status)
    if (status < 0) then
      print *, "ERROR: Unable to set dimension scale for "//trim(dset_name)//": ", status
      stop
    end if 
    
    attr_val  = "Scan start time (TAI93)";
    status = put_attribute_l2(lname ,dset_name,attr_val ,h5f_id)
    attr_val  = "seconds"
    status = put_attribute_l2(unit_name ,dset_name,attr_val , h5f_id)
!     status = put_attribute_l2_float(fname, dset_name, dfv, h5f_id, attr_size)  
    call h5ltset_attribute_double_f(h5f_id, dset_name, fname, dfv, attr_size, status)
    if (status /= 0) then
      print *, "ERROR: Failed to create fill value attribute on "//trim(dset_name)//": ", status
    return
    end if
       
    attr_val  = "coordinates"
    attr_val2 = "Longitude Latitude"
    call h5ltset_attribute_string_f(h5f_id, trim(dset_name), trim(attr_val), trim(attr_val2), status)
    if (status < 0) then
      print *, "ERROR: Unable to create "//trim(attr_val)//" attribute on "//trim(dset_name)//": ", status
      stop
    end if
    endif !normal_output 
    

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Geolocation and geometry !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Latitude
    dset_name = "/Latitude"
    attr_val = "Latitude"    
    attr_val2 = "degrees"
    dimtype  = "2"
    call make_dset_f(dimtype,dset_name, rdum1,lat, rdum3, h5f_id, attr_val,&
         & lname,unit_name, attr_val2, dsid, dsid_lat, dsid_lon, dumid)
    val_range =(/-90.0, 90.0/) 
    status = put_attribute_l2_float_arr(valname, dset_name, val_range, h5f_id, val_size) 
    
  ! Longitude
    dset_name = "/Longitude"
    attr_val = "Longitude"    
    attr_val2 = "degrees"
    dimtype  = "2"     
    call make_dset_f(dimtype,dset_name, rdum1,lon, rdum3, h5f_id, attr_val,&
         & lname,unit_name, attr_val2, dsid, dsid_lat, dsid_lon, dumid)
    val_range =(/-180.0, 180.0/) 
    status = put_attribute_l2_float_arr(valname, dset_name, val_range, h5f_id, val_size)          

    
  if (normal_output) then 
  ! Solar zenith angle
    dset_name = "/Solar_Zenith_Angle"
    attr_val = "Solar zenith angle"    
    attr_val2 = "degrees"
    dimtype  = "2"
    call make_dset_f(dimtype,dset_name, rdum1,sza, rdum3, h5f_id, attr_val,&
         & lname,unit_name, attr_val2, dsid, dsid_lat, dsid_lon, dumid)
    status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size)
    val_range =(/0.0, 90.0/) 
    status = put_attribute_l2_float_arr(valname, dset_name, val_range, h5f_id, val_size)    
    
  ! Viewing (sensor) zenith angle
    dset_name = "/Viewing_Zenith_Angle"
    attr_val = "Viewing zenith angle"    
    attr_val2 = "degrees"
    dimtype  = "2"    
    call make_dset_f(dimtype,dset_name, rdum1,vza, rdum3, h5f_id, attr_val,&
         & lname,unit_name, attr_val2, dsid, dsid_lat, dsid_lon, dumid)
    status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size)
    val_range =(/0.0, 90.0/) 
    status = put_attribute_l2_float_arr(valname, dset_name, val_range, h5f_id, val_size) 
        
  ! Relative azimuth angle
    dset_name = "/Relative_Azimuth_Angle"
    attr_val = "Relative azimuth angle (Gordon convention)"
    attr_val2 = "degrees"
    dimtype  = "2"    
    call make_dset_f(dimtype,dset_name, rdum1,raa, rdum3, h5f_id, attr_val,&
         & lname,unit_name, attr_val2, dsid, dsid_lat, dsid_lon, dumid)
    status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size)
    val_range =(/0.0, 180.0/) 
    status = put_attribute_l2_float_arr(valname, dset_name, val_range, h5f_id, val_size)  
           
  ! Scattering angle
    dset_name = "/Scattering_Angle"
    attr_val = "Scattering angle"
    attr_val2 = "degrees"
    dimtype  = "2"    
    call make_dset_f(dimtype,dset_name, rdum1,sca, rdum3, h5f_id, attr_val,&
         & lname,unit_name, attr_val2, dsid, dsid_lat, dsid_lon, dumid)
    status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size)
    val_range =(/0.0, 180.0/) 
    status = put_attribute_l2_float_arr(valname, dset_name, val_range, h5f_id, val_size)  
  endif !normal_output
 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! DBDT merged product !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! DBDT spectral
    dset_name = "/DBDT_spectral_AOD"
    attr_val = "DBDT merged spectral AOD,0.48, 0.55, 0.67"
    dimtype  = "3" 
    call make_dset_f(dimtype,dset_name, rdum1,rdum2, DBDTaot, h5f_id, attr_val,&
         & lname,unit_name, unit_val, dsid, dsid_lbands, dsid_lat , dsid_lon)
    status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size) 
    val_range =(/0.0, 5.0/) 
    status = put_attribute_l2_float_arr(valname, dset_name, val_range, h5f_id, val_size) 

  ! DBDT spectral
    dset_name = "/DBDT_spectral_reflectance"
    attr_val = "DBDT merged spectral reflectance, 354, 388, 480, 550, 670, 870, 1240, 1640, 2250 nm"
    dimtype  = "3" 
    call make_dset_f(dimtype,dset_name, rdum1,rdum2, dbdtspec, h5f_id, attr_val,&
         & lname,unit_name, unit_val, dsid, dsid_dbdtref, dsid_lat , dsid_lon)
    status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size) 
    val_range =(/0.0, 5.0/) 
    status = put_attribute_l2_float_arr(valname, dset_name, val_range, h5f_id, val_size) 

  ! DBDT QA
    dset_name = "DBDT_QA"
    attr_val = "DBDT merged QA"
    dimtype  = "2"   
    call make_dset_i(dimtype,dset_name, idum1, DBDTqa, idum3, h5f_id, attr_val,&
         & lname,unit_name, unit_val, dsid, dsid_lat , dsid_lon, dumid)
    call h5ltset_attribute_int_f(h5f_id, dset_name, fname, zerofv, attr_size, status)
    if (status /= 0) then
      print *, "ERROR: Failed to create fill value attribute on "//trim(dset_name)//": ", status
    return
    end if
    ival_range =(/0, 3/) 
    status = put_attribute_l2_int_arr(valname, dset_name, ival_range, h5f_id, val_size) 
    
  ! DBDT flag
    dset_name = "DBDT_flag"
    attr_val = "DBDT merged flag"
    dimtype  = "2"   
    call make_dset_i(dimtype,dset_name, idum1, DBDTflag, idum3, h5f_id, attr_val,&
         & lname,unit_name, unit_val, dsid, dsid_lat , dsid_lon, dumid)
    call h5ltset_attribute_int_f(h5f_id, dset_name, fname, zerofv, attr_size, status)
    if (status /= 0) then
      print *, "ERROR: Failed to create fill value attribute on "//trim(dset_name)//": ", status
    return
    end if
    ival_range =(/0, 3/) 
    status = put_attribute_l2_int_arr(valname, dset_name, ival_range, h5f_id, val_size)         

  ! DBDT fmf
    dset_name = "/DBDT_fmf"
    attr_val = "DBDT merged fmf"
    dimtype  = "2"    
    call make_dset_f(dimtype,dset_name, rdum1,dbdtfmf, rdum3, h5f_id, attr_val,&
         & lname,unit_name, unit_val, dsid, dsid_lat, dsid_lon, dumid) 
    status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size) 
    val_range =(/0.0, 120.0/) 
    status = put_attribute_l2_float_arr(valname, dset_name, val_range, h5f_id, val_size) 
          
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Land and ocean combined retrievals !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! AOD at 550 nm for land and ocean
    dset_name = "/Aerosol_Optical_Thickness_550_Land_Ocean"
    attr_val = "Deep Blue/SOAR aerosol optical thickness at 550 nm over land &
  &and ocean"
    dimtype  = "2"    
    call make_dset_f(dimtype,dset_name, rdum1,caot550, rdum3, h5f_id, attr_val,&
         & lname,unit_name, unit_val, dsid, dsid_lat, dsid_lon, dumid) 
    status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size) 
    val_range =(/0.0, 5.0/) 
    status = put_attribute_l2_float_arr(valname, dset_name, val_range, h5f_id, val_size) 

  if (normal_output) then       
    ! AOD at 550 nm for land and ocean, QA-filtered
    dset_name = "/Aerosol_Optical_Thickness_550_Land_Ocean_Best_Estimate"
    attr_val = "Deep Blue/SOAR aerosol optical thickness at 550 nm over land &
   &and ocean, QA-filtered"
    dimtype  = "2"    
    call make_dset_f(dimtype,dset_name, rdum1, caot550_best, rdum3, h5f_id, attr_val,&
         & lname,unit_name, unit_val, dsid, dsid_lat, dsid_lon, dumid)         
    status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size)  
    val_range =(/0.0, 5.0/) 
    status = put_attribute_l2_float_arr(valname, dset_name, val_range, h5f_id, val_size) 
    
  ! AE for land and ocean
    dset_name = "/Angstrom_Exponent_Land_Ocean"
    attr_val = "Deep Blue/SOAR Angstrom exponent over land and ocean"
    dimtype  = "2"    
    call make_dset_f(dimtype,dset_name, land_bands, cae ,aot, h5f_id, attr_val,&
         & lname,unit_name, unit_val, dsid,dsid_lat , dsid_lon,dsid_lbands) 
    status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size)  
    val_range =(/-0.5, 2.5/) 
    status = put_attribute_l2_float_arr(valname, dset_name, val_range, h5f_id, val_size) 
    
  ! AE for land and ocean, QA-filtered
    dset_name = "/Angstrom_Exponent_Land_Ocean_Best_Estimate"
    attr_val = "Deep Blue/SOAR Angstrom exponent over land and ocean, QA-filtered"
    dimtype  = "2"    
    call make_dset_f(dimtype,dset_name, rdum1, cae_best, rdum3, h5f_id, attr_val,&
         & lname,unit_name, unit_val, dsid, dsid_lat, dsid_lon, dumid)
    status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size) 
    val_range =(/-0.5, 2.5/) 
    status = put_attribute_l2_float_arr(valname, dset_name, val_range, h5f_id, val_size) 
  endif !normal_output
    
  !!!!!!!!!!!!!!!!!!!
  ! Land retrievals !
  !!!!!!!!!!!!!!!!!!!
  ! AOD at 550 nm over land
    dset_name = "/Aerosol_Optical_Thickness_550_Land"
    attr_val = "Deep Blue aerosol optical thickness at 550 nm over land"
    dimtype  = "2"    
    call make_dset_f(dimtype,dset_name, rdum1, aot550, rdum3, h5f_id, attr_val,&
         & lname,unit_name, unit_val, dsid, dsid_lat, dsid_lon, dumid)
    status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size) 
    val_range =(/0.0, 5.0/) 
    status = put_attribute_l2_float_arr(valname, dset_name, val_range, h5f_id, val_size) 
  
  if (normal_output) then 
  ! AOD at 550 nm over land, QA-filtered
    dset_name = "/Aerosol_Optical_Thickness_550_Land_Best_Estimate"
    attr_val = "Deep Blue aerosol optical thickness at 550 nm over land, QA-filtered"
    dimtype  = "2"    
    call make_dset_f(dimtype,dset_name, rdum1, aot550_best, rdum3,  h5f_id, attr_val,&
         & lname,unit_name, unit_val, dsid,dsid_lat , dsid_lon,dumid) 
    status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size) 
    val_range =(/0.0, 5.0/) 
    status = put_attribute_l2_float_arr(valname, dset_name, val_range, h5f_id, val_size) 
      
  ! Spectral AOD over land, QA-filtered
    dset_name = "/Spectral_Aerosol_Optical_Thickness_Land"
    allocate(aot(size(aot412,1), size(aot412,2),3), stat=status)
    if (status /= 0) then
      print *, "ERROR: Failed to allocate array for spectral AOT: ", status
      return
    end if
    aot(:,:,1) = aot412
    aot(:,:,2) = aot488
    aot(:,:,3) = aot670
    attr_val = "Deep Blue spectral aerosol optical thickness at 412, 488, and 670 nm over land"
    dimtype  = "3" 
    call make_dset_f(dimtype,dset_name, rdum1,rdum2, aot, h5f_id, attr_val,&
         & lname,unit_name, unit_val, dsid, dsid_lbands, dsid_lat , dsid_lon)
    status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size) 
    val_range =(/0.0, 5.0/) 
    status = put_attribute_l2_float_arr(valname, dset_name, val_range, h5f_id, val_size) 
    
  ! AE over land
    dset_name = "/Angstrom_Exponent_Land"
    attr_val = "Deep Blue Angstrom exponent over land; 412/488 nm when Algorithm_Flag_Land=0, &
  & 488/670 nm when Algorithm_Flag_Land>0"
    dimtype  = "2"   
    call make_dset_f(dimtype,dset_name, rdum1, ae, rdum3, h5f_id, attr_val,&
         & lname,unit_name, unit_val, dsid, dsid_lat , dsid_lon, dumid)
    status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size) 
    val_range =(/-0.5, 2.5/) 
    status = put_attribute_l2_float_arr(valname, dset_name, val_range, h5f_id, val_size) 
    
  ! AE over land, QA-filtered
    dset_name = "/Angstrom_Exponent_Land_Best_Estimate"
    attr_val = "Deep Blue Angstrom exponent over land; 412/488 nm when Algorithm_Flag_Land=0, &
  & 488/670 nm when Algorithm_Flag_Land>0, QA-filtered"
    dimtype  = "2"   
    call make_dset_f(dimtype,dset_name, rdum1, ae_best, rdum3, h5f_id, attr_val,&
         & lname,unit_name, unit_val, dsid, dsid_lat , dsid_lon, dumid)
    status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size) 
    val_range =(/-0.5, 2.5/) 
    status = put_attribute_l2_float_arr(valname, dset_name, val_range, h5f_id, val_size) 

  ! Spectral SSA over land
    dset_name = "/Spectral_Single_Scattering_Albedo_Land"
    attr_val = "Deep Blue single scattering albedo over land; 412/488/670 nm"
    dimtype  = "3"   
    call make_dset_f(dimtype,dset_name, rdum1, rdum2, ssa, h5f_id, attr_val,&
         & lname,unit_name, unit_val, dsid, dsid_lbands,dsid_lat , dsid_lon)
    status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size) 
    val_range =(/0.0, 1.0/) 
    status = put_attribute_l2_float_arr(valname, dset_name, val_range, h5f_id, val_size) 

  ! Spectral ocean reflectance
    dset_name = "/Spectral_TOA_Reflectance_Ocean"
    attr_val = "Cell-averaged cloud-screened spectral TOA reflectance over ocean"
    dimtype  = "3"  
    attr_val2 = "I/F" 
    call make_dset_f(dimtype,dset_name, rdum1, rdum2, oreflc_mean, h5f_id, attr_val,&
         & lname,unit_name, attr_val2, dsid, dsid_lat , dsid_lon, dsid_oref)
    status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size) 
    val_range =(/0.0, 1.0/) 
    status = put_attribute_l2_float_arr(valname, dset_name, val_range, h5f_id, val_size) 

  ! Spectral land reflectance
    dset_name = "/Spectral_TOA_Reflectance_Land"
    attr_val = "Cell-averaged cloud-screened spectral TOA reflectance over land"
    dimtype  = "3"   
    attr_val2 = "I/F"
    call make_dset_f(dimtype,dset_name, rdum1, rdum2, lreflc_mean, h5f_id, attr_val,&
         & lname,unit_name,  attr_val2 , dsid, dsid_lat , dsid_lon, dsid_oref)
    status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size) 
    val_range =(/0.0, 1.0/) 
    status = put_attribute_l2_float_arr(valname, dset_name, val_range, h5f_id, val_size) 
        
  ! Land QA flag
    dset_name = "/Aerosol_Optical_Thickness_QA_Flag_Land"
    attr_val = "Deep Blue quality assurance flag over land. 0=no retrieval, 1=poor, &
  & 2=moderate, 3=good"
    dimtype  = "2"   
    call make_dset_i(dimtype,dset_name, idum1, qa_flag, idum3, h5f_id, attr_val,&
         & lname,unit_name, unit_val, dsid, dsid_lat , dsid_lon, dumid)
    call h5ltset_attribute_int_f(h5f_id, dset_name, fname, zerofv, attr_size, status)
    if (status /= 0) then
      print *, "ERROR: Failed to create fill value attribute on "//trim(dset_name)//": ", status
    return
    end if
    ival_range =(/0, 3/) 
    status = put_attribute_l2_int_arr(valname, dset_name, ival_range, h5f_id, val_size) 
         
  ! Land algorithm flag
    dset_name = "/Algorithm_Flag_Land"
    attr_val = "Deep Blue algorithm flag over land. 0=arid DB, 1=vegetated, 2=mixed"
    dimtype  = "2"   
    call make_dset_i(dimtype,dset_name, idum1, alg_old, idum3, h5f_id, attr_val,&
         & lname,unit_name, unit_val, dsid, dsid_lat , dsid_lon, dumid)
!     status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size) 
    call h5ltset_attribute_int_f(h5f_id, dset_name, fname, ifv, attr_size, status)
    if (status /= 0) then
      print *, "ERROR: Failed to create fill value attribute on "//trim(dset_name)//": ", status
    return
    end if
    ival_range =(/0, 2/) 
    status = put_attribute_l2_int_arr(valname, dset_name, ival_range, h5f_id, val_size) 

  ! Land algorithm flag(db/veg/2.2/mix)
!     dset_name = "/Algorithm_Flag_Land_for_test2"
!     attr_val = "Deep Blue algorithm flag over land. 0=arid DB, 1=vegetated, 2=2.2micron, 3=mixed"
!     dimtype  = "2"   
!     call make_dset_i(dimtype,dset_name, idum1, alg_old2, idum3, h5f_id, attr_val,&
!          & lname,unit_name, unit_val, dsid, dsid_lat , dsid_lon, dumid)
! !     status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size) 
!     call h5ltset_attribute_int_f(h5f_id, dset_name, fname, ifv, attr_size, status)
!     if (status /= 0) then
!       print *, "ERROR: Failed to create fill value attribute on "//trim(dset_name)//": ", status
!     return
!     end if
!     ival_range =(/0, 3/) 
!     status = put_attribute_l2_int_arr(valname, dset_name, ival_range, h5f_id, val_size) 
    
  ! Land algorithm flag
!     dset_name = "/Algorithm_Flag_Land_for_test"
!     attr_val = "Deep Blue algorithm flag only for algorithm development. 1=DB/AERONET BRDF,"//&
!        &"2=DB/BRDF coefficient, 3=DB/minimum table, 11=VEG, 12=2.2micron, 20=mix"
!     dimtype  = "2"   
!     call make_dset_i(dimtype,dset_name, idum1, alg, idum3, h5f_id, attr_val,&
!          & lname,unit_name, unit_val, dsid, dsid_lat , dsid_lon, dumid)
! !     status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size) 
!     call h5ltset_attribute_int_f(h5f_id, dset_name, fname, ifv, attr_size, status)
!     if (status /= 0) then
!       print *, "ERROR: Failed to create fill value attribute on "//trim(dset_name)//": ", status
!     return
!     end if
!     ival_range =(/0, 20/) 
!     status = put_attribute_l2_int_arr(valname, dset_name, ival_range, h5f_id, val_size)     
    
  ! retrieval condition 
!     dset_name = "/Retrieval Condition "
!     attr_val = "TBD"
!     dimtype  = "2"   
!     call make_dset_i(dimtype,dset_name, idum1, rcond, idum3, h5f_id, attr_val,&
!          & lname,unit_name, unit_val, dsid, dsid_lat , dsid_lon, dumid)
!     status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size) 
    
    endif !normal_output    
 
  !!!!!!!!!!!!!!!!!!!!
  ! Ocean retrievals !
  !!!!!!!!!!!!!!!!!!!!
  ! AOD at 550 nm, ocean
    dset_name = "/Aerosol_Optical_Thickness_550_Ocean"
    attr_val = "SOAR aerosol optical thickness at 550 nm over water"
    dimtype  = "2"    
    call make_dset_f(dimtype,dset_name, rdum1, oaot550, rdum3,  h5f_id, attr_val,&
     & lname,unit_name, unit_val, dsid,dsid_lat , dsid_lon,dumid) 
    ! we used to use oaot(:,:,2) for 550nm AOD but oaot(:,:,2)a and oaot550 is slightly different
    status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size) 
    val_range =(/0.0, 5.0/) 
    status = put_attribute_l2_float_arr(valname, dset_name, val_range, h5f_id, val_size) 
    
  if (normal_output) then 
  ! AOD at 550 nm, ocean, QA-filtered
    dset_name = "/Aerosol_Optical_Thickness_550_Ocean_Best_Estimate"
    attr_val = "SOAR aerosol optical thickness at 550 nm over water, QA-filtered"
    dimtype  = "2"    
    call make_dset_f(dimtype,dset_name, rdum1, oaot550_best, rdum3,  h5f_id, attr_val,&
         & lname,unit_name, unit_val, dsid,dsid_lat , dsid_lon,dumid) 
    status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size) 
    val_range =(/0.0, 5.0/) 
    status = put_attribute_l2_float_arr(valname, dset_name, val_range, h5f_id, val_size) 
    
  ! AE, ocean
    dset_name = "/Angstrom_Exponent_Ocean"
    attr_val = "SOAR Angstrom exponent (550/865 nm) over water"
    dimtype  = "2"    
    call make_dset_f(dimtype,dset_name, rdum1, oae, rdum3,  h5f_id, attr_val,&
         & lname,unit_name, unit_val, dsid,dsid_lat , dsid_lon,dumid) 
    status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size) 
    val_range =(/-0.5, 2.5/) 
    status = put_attribute_l2_float_arr(valname, dset_name, val_range, h5f_id, val_size) 
    
  ! AE, ocean, QA-filtered
    dset_name = "/Angstrom_Exponent_Ocean_Best_Estimate"
    attr_val = "SOAR Angstrom exponent (550/865 nm) over water, QA-filtered"
    dimtype  = "2"    
    call make_dset_f(dimtype,dset_name, rdum1, oae_best, rdum3,  h5f_id, attr_val,&
         & lname,unit_name, unit_val, dsid,dsid_lat , dsid_lon,dumid) 
    status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size) 
    val_range =(/-0.5, 2.5/) 
    status = put_attribute_l2_float_arr(valname, dset_name, val_range, h5f_id, val_size) 
    
  ! Spectral AOD, ocean
    dset_name = "/Spectral_Aerosol_Optical_Thickness_Ocean"
    attr_val = "SOAR spectral aerosol optical thickness at &
       & 488, 550, 670, 865, 1240, 1610, 2250 nm over water"
    dimtype  = "3"    
    call make_dset_f(dimtype,dset_name, rdum1, rdum2, oaot,  h5f_id, attr_val,&
         & lname,unit_name, unit_val, dsid,dsid_obands,dsid_lat , dsid_lon) 
    status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size) 
    val_range =(/0.0, 10.0/) 
    status = put_attribute_l2_float_arr(valname, dset_name, val_range, h5f_id, val_size) 
    
  ! QA flag, ocean
    dset_name = "/Aerosol_Optical_Thickness_QA_Flag_Ocean"
    attr_val = "SOAR quality assurance flag over water. 0=no retrieval, 1=poor, &
   & 3=good, 2=only in NRT product"
    dimtype  = "2"   
    call make_dset_i(dimtype,dset_name, idum1, oqa_flag, idum3, h5f_id, attr_val,&
         & lname,unit_name, unit_val, dsid, dsid_lat , dsid_lon, dumid)
    call h5ltset_attribute_int_f(h5f_id, dset_name, fname, zerofv, attr_size, status)
    if (status /= 0) then
      print *, "ERROR: Failed to create fill value attribute on "//trim(dset_name)//": ", status
    return
    end if
    ival_range =(/0, 3/) 
    status = put_attribute_l2_int_arr(valname, dset_name, ival_range, h5f_id, val_size) 
        
  ! Chosen aerosol model flag, ocean
    dset_name = "/Aerosol_Type_Ocean"
    attr_val = "SOAR retrieved aerosol optical model over water. -999=no retrieval, 1=dust, &
    &2=fine dominated, 3=maritime, 4=mixed"
    dimtype  = "2"   
    call make_dset_i(dimtype,dset_name, idum1, omodel_flag, idum3, h5f_id, attr_val,&
         & lname,unit_name, unit_val, dsid, dsid_lat , dsid_lon, dumid)
    call h5ltset_attribute_int_f(h5f_id, dset_name, fname, ifv, attr_size, status)
    if (status /= 0) then
      print *, "ERROR: Failed to create fill value attribute on "//trim(dset_name)//": ", status
    return
    end if
    ival_range =(/1, 4/) 
    status = put_attribute_l2_int_arr(valname, dset_name, ival_range, h5f_id, val_size)  
            
  ! Chosen aerosol model flag, land
    dset_name = "/Aerosol_Type_Land"
    attr_val = "Aerosol optical model over land. -999=no retrieval, 0=dust, &
    &1=smoke, 2=high altitude smoke, 3=pyrocumulonimbus clouds, 4=non-smoke fine mode, &
    &5=mixed, 6=background"
    dimtype  = "2"   
    call make_dset_i(dimtype,dset_name, idum1, ltype, idum3, h5f_id, attr_val,&
         & lname,unit_name, unit_val, dsid, dsid_lat , dsid_lon, dumid)
    call h5ltset_attribute_int_f(h5f_id, dset_name, fname, ifv, attr_size, status)
    if (status /= 0) then
      print *, "ERROR: Failed to create fill value attribute on "//trim(dset_name)//": ", status
    return
    end if
    ival_range =(/0, 6/) 
    status = put_attribute_l2_int_arr(valname, dset_name, ival_range, h5f_id, val_size) 
       
  ! Chosen aerosol model flag, land_ocean
    dset_name = "/Aerosol_Type_Land_Ocean"
    attr_val = "Aerosol optical model over land and water. -999=no retrieval, 0=dust(land+ocean), &
    &1=smoke, 2=high altitude smoke, 3=pyrocumulonimbus clouds, 4=non-smoke fine mode, &
    &5=mixed(land+ocean), 6=background(land+ocean maritime), 7=fine dominated"
    dimtype  = "2"   
    call make_dset_i(dimtype,dset_name, idum1, ctype, idum3, h5f_id, attr_val,&
         & lname,unit_name, unit_val, dsid, dsid_lat , dsid_lon, dumid)
    call h5ltset_attribute_int_f(h5f_id, dset_name, fname, ifv, attr_size, status)
    if (status /= 0) then
      print *, "ERROR: Failed to create fill value attribute on "//trim(dset_name)//": ", status
    return
    end if
    ival_range =(/0, 7/) 
    status = put_attribute_l2_int_arr(valname, dset_name, ival_range, h5f_id, val_size) 
        
  ! Retrieval algorithm flag, ocean
    dset_name = "/Algorithm_Flag_Ocean"
    attr_val = "SOAR algorithm flag over water. -999=no retrieval, 0=full &
        &retrieval, 1=turbid/shallow, 2=mixed"
    dimtype  = "2"   
    call make_dset_i(dimtype,dset_name, idum1, oalg, idum3, h5f_id, attr_val,&
         & lname,unit_name, unit_val, dsid, dsid_lat , dsid_lon, dumid)
!     status = put_attribute_l2_float(fname, dset_name, ifv, h5f_id, attr_size) 
    call h5ltset_attribute_int_f(h5f_id, dset_name, fname, ifv, attr_size, status)
    if (status /= 0) then
      print *, "ERROR: Failed to create fill value attribute on "//trim(dset_name)//": ", status
    return
    end if
    ival_range =(/0, 2/) 
    status = put_attribute_l2_int_arr(valname, dset_name, ival_range, h5f_id, val_size) 
    
  ! Fine mode fraction of AOD at 550 nm, ocean
    dset_name = "Fine_Mode_Fraction_550_Ocean"
    attr_val = "SOAR fine mode fraction of aerosol optical thickness at 550 nm over water"
    dimtype  = "2"    
    call make_dset_f(dimtype,dset_name, rdum1, ofmf, rdum3,  h5f_id, attr_val,&
         & lname,unit_name, unit_val, dsid,dsid_lat , dsid_lon,dumid) 
    status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size) 
    val_range =(/0.0, 1.0/) 
    status = put_attribute_l2_float_arr(valname, dset_name, val_range, h5f_id, val_size) 
    
  ! Fine mode fraction of AOD at 550 nm, ocean, QA-filtered
    dset_name = "Fine_Mode_Fraction_550_Ocean_Best_Estimate"
    attr_val = "SOAR fine mode fraction of aerosol optical thickness at 550 nm over water, QA-filtered"
    dimtype  = "2"    
    call make_dset_f(dimtype,dset_name, rdum1, ofmf_best, rdum3,  h5f_id, attr_val,&
         & lname,unit_name, unit_val, dsid,dsid_lat , dsid_lon,dumid) 
    status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size) 
    val_range =(/0.0, 1.0/) 
    status = put_attribute_l2_float_arr(valname, dset_name, val_range, h5f_id, val_size) 

    endif !normal_output
        
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Debug info (which most users don't want so  !
  ! can be switched off in final L2 production) !
  ! Still need to decide if any of this should  !
  ! make it into final L2 files.                !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (normal_output) then 

  ! Number of non-fill pixels in cell
    dset_name = "/Number_Valid_Pixels"
    attr_val = "Number of non-fill L1b pixels in cell"
    dimtype  = "2"    
    call make_dset_i(dimtype,dset_name, idum1, n_valid_pixels, idum3,  h5f_id, attr_val,&
         & lname,unit_name, unit_val, dsid,dsid_lat , dsid_lon,dumid) 
!     status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size)
    call h5ltset_attribute_int_f(h5f_id, dset_name, fname, ifv, attr_size, status)
    if (status /= 0) then
      print *, "ERROR: Failed to create fill value attribute on "//trim(dset_name)//": ", status
    return
    end if
    ival_range =(/0, 100/) 
    status = put_attribute_l2_int_arr(valname, dset_name, ival_range, h5f_id, val_size) 
        
  ! Number of pixels used for retrieval in cell, land
    dset_name = "/Number_Of_Pixels_Used_Land"
    attr_val = "Deep Blue number of pixels used within cell for retrieval over land"
    dimtype  = "2"    
    call make_dset_i(dimtype,dset_name, idum1, naot550, idum3,  h5f_id, attr_val,&
         & lname,unit_name, unit_val, dsid,dsid_lat , dsid_lon,dumid)    
!     status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size) 
    call h5ltset_attribute_int_f(h5f_id, dset_name, fname, ifv, attr_size, status)
    if (status /= 0) then
      print *, "ERROR: Failed to create fill value attribute on "//trim(dset_name)//": ", status
    return
    end if
    ival_range =(/0, 100/) 
    status = put_attribute_l2_int_arr(valname, dset_name, ival_range, h5f_id, val_size) 
     
  ! Number of pixels used for retrieval in cell, ocean
    dset_name = "/Number_Of_Pixels_Used_Ocean"
    attr_val = "SOAR number of pixels used within cell for retrieval over water"
    dimtype  = "2"    
    call make_dset_i(dimtype,dset_name, idum1, onaot550, idum3,  h5f_id, attr_val,&
         & lname,unit_name, unit_val, dsid,dsid_lat , dsid_lon,dumid)    
!     status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size) 
    call h5ltset_attribute_int_f(h5f_id, dset_name, fname, ifv, attr_size, status)
    if (status /= 0) then
      print *, "ERROR: Failed to create fill value attribute on "//trim(dset_name)//": ", status
    return
    end if
    ival_range =(/0, 100/) 
    status = put_attribute_l2_int_arr(valname, dset_name, ival_range, h5f_id, val_size) 
    
  ! Standard deviation of AOD in cell, land
    dset_name = "/Aerosol_Optical_Thickness_550_STDV_Land"
    attr_val = "Deep Blue standard deviation of aerosol optical thickness at 550 nm &
   & within cell over land"
    dimtype  = "2"    
    call make_dset_f(dimtype,dset_name, rdum1, aot550_sd, rdum3,  h5f_id, attr_val,&
         & lname,unit_name, unit_val, dsid,dsid_lat , dsid_lon,dumid)  
    status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size) 
    val_range =(/0.0, 5.0/) 
    status = put_attribute_l2_float_arr(valname, dset_name, val_range, h5f_id, val_size) 
    
  ! Standard deviation of AOD in cell, ocean
    dset_name = "/Aerosol_Optical_Thickness_550_STDV_Ocean"
    attr_val = "SOAR standard deviation of aerosol optical thickness at 550 nm &
  & within cell over water"
    dimtype  = "2"    
    call make_dset_f(dimtype,dset_name, rdum1, osd550, rdum3,  h5f_id, attr_val,&
         & lname,unit_name, unit_val, dsid,dsid_lat , dsid_lon,dumid)  
    status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size) 
    val_range =(/0.0, 5.0/) 
    status = put_attribute_l2_float_arr(valname, dset_name, val_range, h5f_id, val_size) 

  ! Ocean retrieval erro, ocean
!     dset_name = "/Retrieval_Error_Ocean"
!     attr_val = "SOAR retrieval error of aerosol optical thickness at 550 nm &
!   & within cell over water"
!     dimtype  = "2"    
!     call make_dset_f(dimtype,dset_name, rdum1, opt_err, rdum3,  h5f_id, attr_val,&
!          & lname,unit_name, unit_val, dsid,dsid_lat , dsid_lon,dumid)  
!     status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size) 
!     val_range =(/-100.0, 100.0/) 
!     status = put_attribute_l2_float_arr(valname, dset_name, val_range, h5f_id, val_size) 
           
  ! Spectral surface reflectance averaged for cell, land
    dset_name = "/Spectral_Surface_Reflectance"
    attr_val = "Deep Blue spectral (412, 488, 670 nm) surface reflectance over land"
    dimtype  = "3"    
    call make_dset_f(dimtype,dset_name, rdum1, rdum2, sr_avg,  h5f_id, attr_val,&
         & lname,unit_name, unit_val, dsid, dsid_lbands, dsid_lat , dsid_lon)  
    status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size) 
    val_range =(/0.0, 1.0/) 
    status = put_attribute_l2_float_arr(valname, dset_name, val_range, h5f_id, val_size)   
    
  ! TOA NDVI
    dset_name = "/TOA_NDVI"
    attr_val = "Average normalized difference vegetation index"
    dimtype  = "2"    
    call make_dset_f(dimtype,dset_name, rdum1, ndvi_avg, rdum3,  h5f_id, attr_val,&
         & lname,unit_name, unit_val, dsid,dsid_lat , dsid_lon,dumid)  
    status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size) 
    val_range =(/-1.0, 1.0/) 
    status = put_attribute_l2_float_arr(valname, dset_name, val_range, h5f_id, val_size) 
      
  ! Mean over-ocean sum of square residuals
    dset_name = "/Ocean_Sum_Squares"
    attr_val = "Average sum of square residuals for SOAR ocean retrieval"
    dimtype  = "2"    
    call make_dset_f(dimtype,dset_name, rdum1, oss, rdum3,  h5f_id, attr_val,&
         & lname,unit_name, unit_val, dsid,dsid_lat , dsid_lon,dumid)  
    status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size) 
    val_range =(/0.0, 10000000.0/) 
    status = put_attribute_l2_float_arr(valname, dset_name, val_range, h5f_id, val_size) 
    
  ! Ingested wind speed
    dset_name = "/Wind_Speed"
    attr_val = "Ancillary wind speed"
    attr_val2 = "meters per second"
    dimtype  = "2"    
    call make_dset_f(dimtype,dset_name, rdum1, ws, rdum3,  h5f_id, attr_val,&
         & lname,unit_name, attr_val2, dsid,dsid_lat , dsid_lon,dumid)  
    status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size)
    val_range =(/0.0, 100.0/) 
    status = put_attribute_l2_float_arr(valname, dset_name, val_range, h5f_id, val_size) 
         
  ! Ingested wind direction
    dset_name = "/Wind_Direction"
    attr_val = "Ancillary wind direction"
    attr_val2 = "degrees"
    dimtype  = "2"    
    call make_dset_f(dimtype,dset_name, rdum1, wd, rdum3,  h5f_id, attr_val,&
         & lname,unit_name, attr_val2, dsid,dsid_lat , dsid_lon,dumid) 
    status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size)
    val_range =(/0.0, 360.0/) 
    status = put_attribute_l2_float_arr(valname, dset_name, val_range, h5f_id, val_size) 
        
  ! Ingested total column ozone
    dset_name = "/Total_Column_Ozone"
    attr_val = "Ancillary total column ozone amount"
    attr_val2 = "atm/cm"
    dimtype  = "2"    
    call make_dset_f(dimtype,dset_name, rdum1, oz, rdum3,  h5f_id, attr_val,&
         & lname,unit_name, attr_val2, dsid,dsid_lat , dsid_lon,dumid) 
    status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size)
    val_range =(/0.0, 1.0/) 
    status = put_attribute_l2_float_arr(valname, dset_name, val_range, h5f_id, val_size)  
       
  ! Ingested column water vapour
    dset_name = "/Precipitable_Water"
    attr_val = "Ancillary total column water vapor amount"
    attr_val2 = "cm"
    dimtype  = "2"    
    call make_dset_f(dimtype,dset_name, rdum1, wv, rdum3,  h5f_id, attr_val,&
         & lname,unit_name, attr_val2, dsid,dsid_lat , dsid_lon,dumid) 
    status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size)
    val_range =(/0.0, 100.0/) 
    status = put_attribute_l2_float_arr(valname, dset_name, val_range, h5f_id, val_size)    

!  ! Average pixel elevation in cell
!     dset_name = "/Cell_Average_Elevation"
!     attr_val = "Cell-averaged elevation above sea level"
!     attr_val2 = "m"
!     dimtype  = "2"    
!     call make_dset_f(dimtype,dset_name, rdum1, elev_avg, rdum3,  h5f_id, attr_val,&
!          & lname,unit_name, attr_val2, dsid,dsid_lat , dsid_lon,dumid) 
!     status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size) 

  ! Average pixel elevation in cell
    dset_name = "/Cell_Average_Elevation_Land"
    attr_val = "Cell-averaged elevation above sea level, land pixels"
    attr_val2 = "m"
    dimtype  = "2"
    call make_dset_f(dimtype,dset_name, rdum1, elev_avg_land, rdum3,  h5f_id, attr_val,&
         & lname,unit_name, attr_val2, dsid,dsid_lat , dsid_lon,dumid)
    status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size)
    val_range =(/-500.0, 10000.0/) 
    status = put_attribute_l2_float_arr(valname, dset_name, val_range, h5f_id, val_size) 
    
  ! Average pixel elevation in cell
    dset_name = "/Cell_Average_Elevation_Ocean"
    attr_val = "Cell-averaged elevation above sea level, ocean pixels"
    attr_val2 = "m"
    dimtype  = "2"
    call make_dset_f(dimtype,dset_name, rdum1, elev_avg_ocean, rdum3,  h5f_id, attr_val,&
         & lname,unit_name, attr_val2, dsid,dsid_lat , dsid_lon,dumid)
    status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size)
    val_range =(/-500.0, 10.0/) 
    status = put_attribute_l2_float_arr(valname, dset_name, val_range, h5f_id, val_size) 
        
  ! Average Chl concentration from climatology
!     dset_name = "/Cell_Average_Chl"
!     attr_val = "Cell-averaged climatological Chlorophyll concentration from &
!         &ancillary data"
!     attr_val2 = "mg/m3"
!     dimtype  = "2"    
!     call make_dset_f(dimtype,dset_name, rdum1, chl_avg, rdum3,  h5f_id, attr_val,&
!          & lname,unit_name, attr_val2, dsid,dsid_lat , dsid_lon,dumid) 
!     status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size) 
!     val_range =(/0.0, 100.0/) 
!     status = put_attribute_l2_float_arr(valname, dset_name, val_range, h5f_id, val_size) 
    endif !normal_output            
    !===============!
    !extended output!
    !===============!

    if (extended_output) then 

    ! Rayleigh-corrected NDVI
    dset_name = "/RC_NDVI"
    attr_val = "Average normalized difference vegetation index, corrected &
   & for Rayleigh scattering"
    dimtype  = "2"
    call make_dset_f(dimtype,dset_name, rdum1, rcndvi_avg, rdum3,  h5f_id, attr_val,&
         & lname,unit_name, unit_val, dsid,dsid_lat , dsid_lon,dumid)
    status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size)

    ! Dstar
    dset_name = "/Dstar"
    attr_val = "Average Dstar value"
    dimtype  = "2"
    call make_dset_f(dimtype, dset_name, rdum1, dstar_avg, rdum3, h5f_id, attr_val,&
         & lname,unit_name, unit_val, dsid, dsid_lat, dsid_lon, dumid)
    status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size)

    ! BTD11
    dset_name = "/BTD11"
    attr_val = "Average brightness temperature difference between 11 and 12 microns"
    dimtype  = "2"
    call make_dset_f(dimtype, dset_name, rdum1, btd11_avg, rdum3, h5f_id, attr_val,&
         & lname,unit_name, unit_val, dsid, dsid_lat, dsid_lon, dumid)
    status = put_attribute_l2_float(fname, dset_name, fv, h5f_id, attr_size)

    ! smoke count
    dset_name = "/Number_Smoke_Pixels"
    attr_val = "Number of smoke pixels within cell"
    dimtype  = "2"
    call make_dset_i(dimtype, dset_name, idum1, smoke_count, idum3, h5f_id, attr_val,&
         & lname,unit_name, unit_val, dsid, dsid_lat, dsid_lon, dumid)

    end if ! end extended_output
 
  call h5fclose_f(h5f_id, status)
  if (status /= 0) then
    print *, "ERROR: Failed to close VIIRS DB output file: ", status
    return
  end if
  
  call h5close_f(status)
  if (status /= 0) then
    print *, "ERROR: Failed to close HDF interface: ", status
    return
  end if
  
  deallocate(nxtrack_arr)
  if (status /= 0) then
    print *, "ERROR: Failed to allocate array for spectral nxtrack_arr: ", status
    return
  end if
  
  deallocate(natrack_arr)
  if (status /= 0) then
    print *, "ERROR: Failed to allocate array for spectral nxtrack_arr: ", status
    return
  end if
  
  if (normal_output) then 
   deallocate(tmp_target)
   if (status /= 0) then
     print *, "ERROR: Failed to allocate array for spectral tmp_target: ", status
     return
   end if
  
   deallocate(time_arr)
   if (status /= 0) then
     print *, "ERROR: Failed to allocate array for spectral time_arr: ", status
     return
   end if

!   deallocate(aot)
!   if (status /= 0) then
!     print *, "ERROR: Failed to allocate array for spectral aot: ", status
!     return
!   endif !normal_output
  end if  

end function create_viirs_l2

integer function calc_minmax_ratio(refl412, mm_ratio) result(status)
  implicit none
  
  real, dimension(:,:), intent(in)    ::  refl412
  real, dimension(:,:), intent(inout) ::  mm_ratio

  real        ::  xmin, xmax
  integer     ::  i, j, i1, i2, j1, j2
  
  status = -1
  
  do j = 1, size(refl412,2)
    do i = 1, size(refl412,1)
      if (refl412(i,j) < -900.0) continue
      
      i1 = max(i-1,1)
      i2 = min(i+1,size(refl412,1))
      j1 = max(j-1,1)
      j2 = min(j+1,size(refl412,2))
      xmin = -999.0 
      xmax = -999.0
      if (count(refl412(i1:i2,j1:j2) > -900.0) > 0) then
        xmin = minval(refl412(i1:i2,j1:j2), refl412(i1:i2,j1:j2) > -900.0)
        xmax = maxval(refl412(i1:i2,j1:j2), refl412(i1:i2,j1:j2) > -900.0)
        mm_ratio(i,j) = xmax/xmin
        if (xmin == 0) mm_ratio(i,j) = xmin/xmin  !0/0=nan
        if (xmin /= 0) mm_ratio(i,j) = xmax/xmin           
      end if    
      
    end do
  end do
  
  status = 0
  return

end function calc_minmax_ratio

! -- calculate gas correction factors for each band for VIIRS.
! -- SZA and VZA are in degrees. ozone is in atm/cm, water vapor should be in 
! -- cm, and surface pressure should be in atms.
! -- good status = 0, bad status = -1.
! -- Coefficients from Dr. Sayer's IDL script, viirs_new_gas_correctoin.pro. 
! -- Derived from gas correction expressions from paper by Dr. Robert Levy:
! -- "The Collection 6 MODIS aerosol products over land and ocean"
! -- http://www.atmos-meas-tech.net/6/2989/2013/amt-6-2989-2013.html
!
! -- Water vapor reduced by half below based on Tanre paper:
! -- D. Tanre, B. N. Holben and Y. J. Kaufman, "Atmospheric correction against 
! -- algorithm for NOAA-AVHRR products: theory and application," in IEEE 
! -- Transactions on Geoscience and Remote Sensing, vol. 30, no. 2, pp. 231-248, Mar 1992.
! -- doi: 10.1109/36.134074

type(viirs_gas_correction) function calc_gas_correction(sza, vza, ozone, wvapor, &
&                                      platform, sfcprs, status, mask) result(vgt)
  implicit none
  
  real, dimension(:,:), intent(in)      ::  sza             ! in degrees
  real, dimension(:,:), intent(in)      ::  vza             ! in degrees
  real, dimension(:,:), intent(in)      ::  ozone           ! in atm/cm
  real, dimension(:,:), intent(in)      ::  wvapor          ! in cm
  real, dimension(:,:), intent(in)      ::  sfcprs          ! in atms
  character(len=*), intent(in)          ::  platform
  integer, intent(inout)                ::  status   
  integer, dimension(:,:), intent(in), optional ::  mask
  
  real, dimension(4)        ::  amf_coeffs_oz, amf_coeffs_w,amf_coeffs_cs
  real                      ::  amf_oz          ! air mass factors for ozone
  real                      ::  amf_wv          ! water vapor, and constant
  real                      ::  amf_cs          ! species.
  
  real                      ::  trans_oz1, trans_oz2        ! ozone transmittance
  real                      ::  trans_wv1, trans_wv2        ! water vapor transmittance
  real                      ::  trans_cs1, trans_cs2        ! constant species transmittance
  
  integer, parameter        ::  nbands = 10
  real, dimension(2,nbands) ::  oz_coeffs
  real, dimension(3,nbands) ::  wv_coeffs                                  
  real, dimension(nbands)   ::  cs_coeffs
  
  real, dimension(:,:), allocatable       ::  ozone_du
  real, dimension(:,:,:), allocatable    ::  tmp_gt
  integer, dimension(:,:), allocatable   ::  run_mask
  
  real                      ::  xsza
  real                      ::  xvza
  real                      ::  xoz
  real                      ::  xwv
  real                      ::  xsp
  
  real                      ::  mu
  integer, dimension(2)     ::  dims
  
  real      ::  amf
  real      ::  x
  integer   ::  i, j, k
  
  real, parameter ::  d2r = 3.14159265/180.0
  logical   ::  do_print
  
!  oz_coeffs =   (/  0.000285210,    1     412
!                    0.00287980,     2     445
!                    0.0180350,      3     488
!                    0.0838500,      4     555
!                    0.0433130,      5     670
!                    0.0106730,      6     745
!                    7.67350e-5,     7     865
!                    1.52580e-8,     8     1.2
!                    0.000000,       9     1.38
!                    0.000000,       10    1.6
!                    0.000000/)      11    2.2
!  
!  wv_coeffs =       4.04370e-5      -0.000986480    -7.37470e-6  
!                    -7.23950e-7     -0.000124690    7.14210e-8
!                    6.77590e-6      -0.000372640    -1.22700e-6
!                    -0.000122860    -0.000247090    2.07450e-5
!                    -0.000517040    -3.06490e-5     7.73180e-5
!                    -0.00533640     0.00186690      0.000872150
!                    -0.00251020     0.000712850     0.000381480
!                    -0.00377030     0.00238370      0.000591240
!                    0.000000        0.000000        0.000000
!                    -0.00115360     0.000863490     0.000137830
!                    -0.00162120     0.00101020      0.000265270
!                    
!  cs_coeffs =       -0.000280560    0.00116490      0.000281710     -0.00111620     7.43100e-5    -0.000304890
!                    -2.83280e-5     0.000103750     2.90410e-5      -0.000102150    7.52440e-6    -2.70540e-5
!                    -0.000117540    0.000366230     0.000120750     -0.000375200    3.12710e-5    -9.67470e-5
!                    -9.96060e-5     0.000311280     0.000102420     -0.000322650    2.64560e-5    -8.17780e-5
!                    -0.00198180     0.00846380      0.00177870      -0.00954910     0.000519320   -0.00231570
!                    -0.00183480     0.00397870      0.00209930      -0.00513400     0.000496360   -0.00107000
!                    -2.75520e-5     0.00112460      8.43890e-6      0.000202290     2.69090e-6    -9.68680e-6
!                    -0.000904070    0.00737160      1.24250e-5      -0.000592510    0.000146410   -0.00118650
!                    0.000000        0.000000        0.000000        0.000000        0.000000      0.000000
!                    -0.0209480      0.00393730      0.00301690      0.0403560       0.00425260    0.00454670
!                    -0.0470690      0.0398200       -0.0126610      -0.0422850      0.00771930    -0.0136530

!  Band Index                 1         2         3         4         5         6         7         8         9        10
!  Bands:                   M01       M02       M03       M04       M05       M06       M07       M08       M10       M11
  oz_coeffs(1,:)  = (/-3.09E-07,-5.72E-05,-1.25E-04,-4.75E-05,-4.79E-05,-7.72e-05, 4.18E-07, 1.19E-07, 1.19E-07,-2.61E-08/)
  oz_coeffs(2,:)  = (/ 4.71E-07, 2.99E-06, 1.98E-05, 9.08E-05, 4.37E-05, 1.09e-05, 2.24E-06, 5.17E-26, 1.03E-25, 3.28E-09/)
  
  wv_coeffs(1,:)  = (/-9.61E+00,-8.52E+00,-9.65E+00,-7.50E+00,-7.69E+00,-5.74E+00,-6.05E+00,-5.16E+00,-6.43E+00,-5.85E+00/)
  wv_coeffs(2,:)  = (/ 9.16E-01, 9.90E-01, 9.87E-01, 9.84E-01, 9.95E-01, 1.00E+00, 9.65E-01, 9.59E-01, 1.02E+00, 1.28E+00/)
  wv_coeffs(3,:)  = (/-2.01E-02,-1.49E-03, 1.80E-04,-3.87E-03,-1.10E-02,-1.88E-02,-1.53E-02,-2.67E-02,-3.60E-03,-5.04E-03/)
  
  cs_coeffs(:)    = (/ 2.47E-04, 3.77E-04, 1.84E-03, 8.34E-04, 1.44E-03, 5.21e-04, 2.45E-05, 1.19E-02, 2.13E-02, 5.32E-02/)
  
  status = -1 
  
  dims = shape(sza)
  allocate(tmp_gt(dims(1),dims(2),nbands), run_mask(dims(1),dims(2)), stat=status)
  if (status /= 0) then
    print *, "ERROR: Unable to allocate tmp gas transmittance array: ", status
    return
  end if
  
! -- initialize run_mask to all 1 unless mask input parameter was used. Then use that.
  run_mask(:,:) = 1
  if (present(mask)) then
    run_mask = mask
  end if
  
  allocate(vgt%m01(dims(1),dims(2)), vgt%m02(dims(1),dims(2)), vgt%m03(dims(1),dims(2)), &
  &       vgt%m04(dims(1),dims(2)), vgt%m05(dims(1),dims(2)), vgt%m06(dims(1),dims(2)),  &
  &       vgt%m07(dims(1),dims(2)), vgt%m08(dims(1),dims(2)),  &
  &       vgt%m10(dims(1),dims(2)), vgt%m11(dims(1),dims(2)), stat=status)
  if (status /= 0) then 
    print *, "ERROR: Failed to allocate gas transmittance object arrays: ", status
    return
  end if 
  
! -- convert ozone values from atm/cm to dobson units (DU)
  allocate(ozone_du(dims(1),dims(2)), stat=status)
  if (status /= 0) then
    print *, "ERROR: Failed to allocate ozone array: ", status
    return
  end if
  ozone_du = ozone * 1000.0

  do k = 1, nbands
    if (k == 9) do_print = .true.
    amf_coeffs_oz  = (/268.45,0.5,115.42,-3.2922/)
    amf_coeffs_w   = (/0.0311,0.1,92.471,-1.3814/)
    amf_coeffs_cs  = (/0.4567,0.07,96.484,-1.697/)   
    
    if (platform .eq. 'GOES' .or. platform .eq. 'AHI') then
      if (k .eq. 1 .or. k .eq. 2 .or. k .eq. 6 .or. k .eq. 8) cycle
    end if
    if (platform .eq. 'GOES' .and. k .eq. 4 ) cycle
    
    do j = 1, dims(2)
      do i = 1, dims(1)
       select case (run_mask(i,j))
       case (1)        
        xvza = vza(i,j)
        xsza = sza(i,j)
        xwv  = wvapor(i,j)   ! cut water vapor by half, see Tanre (1992), moved to deep_blue.f90
        xoz  = ozone_du(i,j)
        xsp  = sfcprs(i,j)
        
!       -- calculate tranmittance correction factor for VZA
!       ---------------------------------------------------
        trans_oz1 = -999.0
        trans_wv1 = -999.0
        trans_cs1 = -999.0

        mu = cos(d2r*xvza)

!       -- calculate ozone air mass factor and transmittance correction factor    
        amf_oz      = 1.0 / (mu + amf_coeffs_oz(1) * ((xvza**amf_coeffs_oz(2)) * &
              &  ((amf_coeffs_oz(3)-xvza)**amf_coeffs_oz(4))))       
        trans_oz1    = exp(oz_coeffs(1,k)+oz_coeffs(2,k)*amf_oz*xoz)
        
!       -- calculate water vapor air mass factor
        amf_wv      = 1.0 / (mu + amf_coeffs_w(1) * ((xvza**amf_coeffs_w(2)) * &
              & ((amf_coeffs_w(3)-xvza)**amf_coeffs_w(4))))
        trans_wv1    = exp(exp(wv_coeffs(1,k) + wv_coeffs(2,k)*log(amf_wv*xwv) + &
        &               wv_coeffs(3,k)*((log(amf_wv*xwv))**2)))
        
!       -- calculate constant species air mass factor
        amf_cs      = 1.0 / (mu + amf_coeffs_cs(1) * ((xvza**amf_coeffs_cs(2)) * &
              & ((amf_coeffs_cs(3)-xvza)**amf_coeffs_cs(4))))
        trans_cs1    = exp(cs_coeffs(k) * amf_cs * xsp)
        
!       -- calculate tranmittance correction factor for SZA
!       ---------------------------------------------------
        trans_oz2 = -999.0
        trans_wv2 = -999.0
        trans_cs2 = -999.0

        mu = cos(d2r*xsza)

!       -- calculate ozone air mass factor and transmittance correction factor
        amf_oz      = 1.0 / (mu + amf_coeffs_oz(1) * ((xsza**amf_coeffs_oz(2)) * &
              & ((amf_coeffs_oz(3)-xsza)**amf_coeffs_oz(4))))       
        trans_oz2    = exp(oz_coeffs(1,k)+oz_coeffs(2,k)*amf_oz*xoz)
        
!       -- calculate water vapor air mass factor
        amf_wv      = 1.0 / (mu + amf_coeffs_w(1) * ((xsza**amf_coeffs_w(2)) * &
              & ((amf_coeffs_w(3)-xsza)**amf_coeffs_w(4))))        
        trans_wv2    = exp(exp(wv_coeffs(1,k) + wv_coeffs(2,k)*log(amf_wv*xwv) + &
        &               wv_coeffs(3,k)*((log(amf_wv*xwv))**2)))
        
!       -- calculate constant species air mass factor
        amf_cs      = 1.0 / (mu + amf_coeffs_cs(1) * ((xsza**amf_coeffs_cs(2)) * &
              & ((amf_coeffs_cs(3)-xsza)**amf_coeffs_cs(4))))     
        trans_cs2    = exp(cs_coeffs(k) * amf_cs * xsp)
        tmp_gt(i,j,k) = trans_oz1 * trans_wv1 * trans_cs1 * trans_oz2 * trans_wv2 * trans_cs2
       case default
        cycle
       end select 
      end do
    end do
  end do

! -- rewrite transmittances to object
  vgt%m01 = tmp_gt(:,:,1)
  vgt%m02 = tmp_gt(:,:,2)
  vgt%m03 = tmp_gt(:,:,3)
  vgt%m04 = tmp_gt(:,:,4)
  vgt%m05 = tmp_gt(:,:,5)
  vgt%m06 = tmp_gt(:,:,6)
  vgt%m07 = tmp_gt(:,:,7)
  vgt%m08 = tmp_gt(:,:,8)
  vgt%m10 = tmp_gt(:,:,9)
  vgt%m11 = tmp_gt(:,:,10)
  
  deallocate(tmp_gt, stat=status)
  if (status /= 0) then
    print *, "ERROR: Unable to deallocate tmp gas transmittance array: ", status
    return
  end if
    
  status = 0

  return
  
end function calc_gas_correction

type(viirs_gas_correction_fullwv) function calc_gas_correction_fullwv(sza, vza, ozone, wvapor, &
&                                       sfcprs, status, mask) result(vgt)
  implicit none

  real, dimension(:,:), intent(in)      ::  sza             ! in degrees
  real, dimension(:,:), intent(in)      ::  vza             ! in degrees
  real, dimension(:,:), intent(in)      ::  ozone           ! in atm/cm
  real, dimension(:,:), intent(in)      ::  wvapor          ! in cm
  real, dimension(:,:), intent(in)      ::  sfcprs          ! in atms
  integer, intent(inout)                ::  status
  integer, dimension(:,:), intent(in), optional ::  mask

  real, dimension(4)        ::  amf_coeffs
  real                      ::  amf_oz          ! air mass factors for ozone
  real                      ::  amf_wv          ! water vapor, and constant
  real                      ::  amf_cs          ! species.

  real                      ::  trans_oz        ! ozone transmittance
  real                      ::  trans_wv        ! water vapor transmittance
  real                      ::  trans_cs        ! constant species transmittance

  integer, parameter        ::  nbands = 1
  real, dimension(2,nbands) ::  oz_coeffs
  real, dimension(3,nbands) ::  wv_coeffs
  real, dimension(nbands)   ::  cs_coeffs

  real, dimension(:,:), allocatable       ::  ozone_du
  real, dimension(:,:,:), allocatable    ::  tmp_gt
  integer, dimension(:,:), allocatable   ::  run_mask

  real                      ::  xsza
  real                      ::  xvza
  real                      ::  xoz
  real                      ::  xwv
  real                      ::  xsp

  real                      ::  mu
  integer, dimension(2)     ::  dims

  real      ::  amf
  real      ::  x
  integer   ::  i, j, k

  real, parameter ::  d2r = 3.14159265/180.0
  logical   ::  do_print

!  Band Index                 1         2         3         4         5        6         7         8         9        10
!  Bands:                   M01       M02       M03       M04       M05        M06       M07       M08       M10       M11
  oz_coeffs(1,:)  = (/-2.61E-08/)
  oz_coeffs(2,:)  = (/ 3.28E-09/)

  wv_coeffs(1,:)  = (/-5.85E+00/)
  wv_coeffs(2,:)  = (/ 1.28E+00/)
  wv_coeffs(3,:)  = (/-5.04E-03/)

  cs_coeffs(:)    = (/ 5.32E-02/)

  status = -1

  dims = shape(sza)
  allocate(tmp_gt(dims(1),dims(2),nbands), run_mask(dims(1),dims(2)), stat=status)
  if (status /= 0) then
    print *, "ERROR: Unable to allocate tmp gas transmittance array: ", status
    return
  end if

! -- initialize run_mask to all 1 unless mask input parameter was used. Then use
! that.
  run_mask(:,:) = 1
  if (present(mask)) then
    run_mask = mask
  end if

  allocate(vgt%m11(dims(1),dims(2)), stat=status)
  if (status /= 0) then
    print *, "ERROR: Failed to allocate gas transmittance object arrays: ", status
    return
  end if

! -- convert ozone values from atm/cm to dobson units (DU)
  allocate(ozone_du(dims(1),dims(2)), stat=status)
  if (status /= 0) then
    print *, "ERROR: Failed to allocate ozone array: ", status
    return
  end if
  ozone_du = ozone * 1000.0

  do k = 1, nbands
    if (k == 9) do_print = .true.

    do j = 1, dims(2)
      do i = 1, dims(1)

        if (run_mask(i,j) == 0) cycle   ! -- skip pixel if run_mask is 0.

        xvza = vza(i,j)
        xsza = sza(i,j)
        xwv  = wvapor(i,j)   ! cut water vapor by half, see Tanre (1992), moved to deep_blue.f90
        xoz  = ozone_du(i,j)
        xsp  = sfcprs(i,j)

!       -- calculate tranmittance correction factor for VZA
!       ---------------------------------------------------
        trans_oz = -999.0
        trans_wv = -999.0
        trans_cs = -999.0

        mu = cos(d2r*xvza)

!       -- calculate ozone air mass factor and transmittance correction factor
        amf_coeffs  = (/268.45,0.5,115.42,-3.2922/)
        amf_oz      = mu + amf_coeffs(1) * ((xvza**amf_coeffs(2)) * ((amf_coeffs(3)-xvza)**amf_coeffs(4)))
        amf_oz      = 1.0 / amf_oz

        trans_oz    = exp(oz_coeffs(1,k)+oz_coeffs(2,k)*amf_oz*xoz)

!       -- calculate water vapor air mass factor
        amf_coeffs  = (/0.0311,0.1,92.471,-1.3814/)
        amf_wv      = mu + amf_coeffs(1) * ((xvza**amf_coeffs(2)) * ((amf_coeffs(3)-xvza)**amf_coeffs(4)))
        amf_wv      = 1.0 / amf_wv

        trans_wv    = exp(exp(wv_coeffs(1,k) + wv_coeffs(2,k)*log(amf_wv*xwv) + &
        &               wv_coeffs(3,k)*((log(amf_wv*xwv))**2)))

!       -- calculate constant species air mass factor
        amf_coeffs  = (/0.4567,0.07,96.484,-1.697/)
        amf_cs      = mu + amf_coeffs(1) * ((xvza**amf_coeffs(2)) * ((amf_coeffs(3)-xvza)**amf_coeffs(4)))
        amf_cs      = 1.0 / amf_cs

        trans_cs    = exp(cs_coeffs(k) * amf_cs * xsp)

        tmp_gt(i,j,k) = trans_oz * trans_wv * trans_cs

!       -- calculate tranmittance correction factor for SZA
!       ---------------------------------------------------
        trans_oz = -999.0
        trans_wv = -999.0
        trans_cs = -999.0

        mu = cos(d2r*xsza)

!       -- calculate ozone air mass factor and transmittance correction factor
        amf_coeffs  = (/268.45,0.5,115.42,-3.2922/)
        amf_oz      = mu + amf_coeffs(1) * ((xsza**amf_coeffs(2)) * ((amf_coeffs(3)-xsza)**amf_coeffs(4)))
        amf_oz      = 1.0 / amf_oz

        trans_oz    = exp(oz_coeffs(1,k)+oz_coeffs(2,k)*amf_oz*xoz)

!       -- calculate water vapor air mass factor
        amf_coeffs  = (/0.0311,0.1,92.471,-1.3814/)
        amf_wv      = mu + amf_coeffs(1) * ((xsza**amf_coeffs(2)) * ((amf_coeffs(3)-xsza)**amf_coeffs(4)))
        amf_wv      = 1.0 / amf_wv

        trans_wv    = exp(exp(wv_coeffs(1,k) + wv_coeffs(2,k)*log(amf_wv*xwv) + &
        &               wv_coeffs(3,k)*((log(amf_wv*xwv))**2)))

!       -- calculate constant species air mass factor
        amf_coeffs  = (/0.4567,0.07,96.484,-1.697/)
        amf_cs      = mu + amf_coeffs(1) * ((xsza**amf_coeffs(2)) * ((amf_coeffs(3)-xsza)**amf_coeffs(4)))
        amf_cs      = 1.0 / amf_cs

        trans_cs    = exp(cs_coeffs(k) * amf_cs * xsp)

        tmp_gt(i,j,k) = tmp_gt(i,j,k) * trans_oz * trans_wv * trans_cs

      end do
    end do
  end do

! -- rewrite transmittances to object
  vgt%m11 = tmp_gt(:,:,1)

  deallocate(tmp_gt, stat=status)
  if (status /= 0) then
    print *, "ERROR: Unable to deallocate tmp gas transmittance array: ", status
    return
  end if

  status = 0

  return

end function calc_gas_correction_fullwv

type(viirs_gas_transmittance) function calc_gas_transmittance(sza, vza, ozone, wvapor, &
&                                       sfcprs, status, mask) result(vgt)
  implicit none

  real, dimension(:,:), intent(in)      ::  sza
  real, dimension(:,:), intent(in)      ::  vza
  real, dimension(:,:), intent(in)      ::  ozone
  real, dimension(:,:), intent(in)      ::  wvapor
  real, dimension(:,:), intent(in)      ::  sfcprs
  integer, intent(inout)                ::  status
  integer, dimension(:,:), intent(in), optional ::  mask

  real                      ::  trans_oz        ! ozone transmittance
  real                      ::  trans_wv        ! water vapor transmittance
  real                      ::  trans_cs        ! constant species transmittance

  real, dimension(11)       ::  oz_coeffs
  real, dimension(3,11)     ::  wv_coeffs
  real, dimension(6,11)     ::  cs_coeffs

  real, dimension(:,:,:), allocatable    ::  tmp_gt
  integer, dimension(:,:), allocatable   ::  run_mask

  integer, dimension(2)     ::  dims

  real      ::  amf
  real      ::  x
  integer   ::  i, j, k

  real, parameter ::  d2r = 3.14159265/180.0

  oz_coeffs =   (/ 0.000285210,0.00287980,0.0180350,0.0838500,0.0433130,0.0106730,  &
  &                7.67350e-5,1.52580e-8,0.000000,0.000000,0.000000/)
  wv_coeffs = transpose(reshape((/4.04370e-5,-7.23950e-7,6.77590e-6,-0.000122860,-0.000517040,&
  &   -0.00533640,-0.00251020,-0.00377030,0.000000,-0.00115360,-0.00162120, &
  &   -0.000986480,-0.000124690,-0.000372640,-0.000247090,-3.06490e-5,      &
  &   0.00186690,0.000712850,0.00238370,0.000000,0.000863490,0.00101020,    &
  &   -7.37470e-6,7.14210e-8,-1.22700e-6,2.07450e-5,7.73180e-5,0.000872150, &
  &   0.000381480,0.000591240,0.000000,0.000137830,0.000265270/), (/11,3/)))
  cs_coeffs = transpose(reshape((/-0.000280560,-2.83280e-5,-0.000117540,-9.96060e-5,-0.00198180, &
  &   -0.00183480,-2.75520e-5,-0.000904070,0.000000,-0.0209480,-0.0470690,  &
  &   0.00116490,0.000103750,0.000366230,0.000311280,0.00846380,0.00397870, &
  &   0.00112460,0.00737160,0.000000,0.00393730,0.0398200,0.000281710,      &
  &   2.90410e-5,0.000120750,0.000102420,0.00177870,0.00209930,8.43890e-6,  &
  &   1.24250e-5,0.000000,0.00301690,-0.0126610,-0.00111620,-0.000102150,   &
  &   -0.000375200,-0.000322650,-0.00954910,-0.00513400,0.000202290,        &
  &   -0.000592510,0.000000,0.0403560,-0.0422850,7.43100e-5,7.52440e-6,     &
  &   3.12710e-5,2.64560e-5,0.000519320,0.000496360,2.69090e-6,0.000146410, &
  &   0.000000,0.00425260,0.00771930,-0.000304890,-2.70540e-5,-9.67470e-5,  &
  &   -8.17780e-5,-0.00231570,-0.00107000,-9.68680e-6,-0.00118650,0.000000, &
  &   0.00454670,-0.0136530/), (/11,6/)))

  status = -1

  dims = shape(sza)
  allocate(tmp_gt(dims(1),dims(2),11), run_mask(dims(1),dims(2)), stat=status)
  if (status /= 0) then
    print *, "ERROR: Unable to allocate tmp gas transmittance array: ", status
    return
  end if

! -- initialize run_mask to all 1 unless mask input parameter was used. Then use that.
  run_mask(:,:) = 1
  if (present(mask)) then
    run_mask = mask
  end if

  allocate(vgt%m01(dims(1),dims(2)), vgt%m02(dims(1),dims(2)), vgt%m03(dims(1),dims(2)), &
  &       vgt%m04(dims(1),dims(2)), vgt%m05(dims(1),dims(2)), vgt%m06(dims(1),dims(2)),  &
  &       vgt%m07(dims(1),dims(2)), vgt%m08(dims(1),dims(2)), vgt%m09(dims(1),dims(2)),  &
  &       vgt%m10(dims(1),dims(2)), vgt%m11(dims(1),dims(2)), stat=status)
  if (status /= 0) then
    print *, "ERROR: Failed to allocate gas transmittance object arrays: ", status
    return
  end if

  do k = 1, 11
    do j = 1, dims(2)
      do i = 1, dims(1)

        if (run_mask(i,j) == 0) cycle   ! -- skip pixel if run_mask is 0.

        amf = (1.0/cos(d2r*sza(i,j))) + (1.0/cos(d2r*vza(i,j)))
        x = -1*amf * ozone(i,j)
        trans_oz = exp(x*oz_coeffs(k))

        x = amf * wvapor(i,j)
        trans_wv = exp(x*wv_coeffs(1,k) + log(x)*wv_coeffs(2,k) + x*log(x)*wv_coeffs(3,k))

        x = sfcprs(i,j)
        trans_cs = exp((x*cs_coeffs(1,k) + cs_coeffs(2,k)*log(x))*amf) *        &
                    exp((x*cs_coeffs(3,k) + cs_coeffs(4,k)*log(x))*log(amf)) *  &
                    exp((x*cs_coeffs(5,k) + cs_coeffs(6,k)*log(x))*amf*log(amf))

        tmp_gt(i,j,k) = trans_oz * trans_wv * trans_cs

      end do
    end do
  end do

! -- rewrite transmittances to object
  vgt%m01 = tmp_gt(:,:,1)
  vgt%m02 = tmp_gt(:,:,2)
  vgt%m03 = tmp_gt(:,:,3)
  vgt%m04 = tmp_gt(:,:,4)
  vgt%m05 = tmp_gt(:,:,5)
  vgt%m06 = tmp_gt(:,:,6)
  vgt%m07 = tmp_gt(:,:,7)
  vgt%m08 = tmp_gt(:,:,8)
  vgt%m09 = tmp_gt(:,:,9)
  vgt%m10 = tmp_gt(:,:,10)
  vgt%m11 = tmp_gt(:,:,11)

  deallocate(tmp_gt, stat=status)
  if (status /= 0) then
    print *, "ERROR: Unable to deallocate tmp gas transmittance array: ", status
    return
  end if

  status = 0

  return

end function calc_gas_transmittance

! -- Detect and flag smoke over the whole globe. 
! -- Created using smoke near Fresno on September 5-9, 2013.
integer function calc_smoke_mask(ler412, ler488, ler670, toa1380, toa2250, &
&                                bt11, gzone, lc, sr650, mask) result(status)
  implicit none
  
  real, dimension(:,:), intent(in)    ::  ler412
  real, dimension(:,:), intent(in)    ::  ler488
  real, dimension(:,:), intent(in)    ::  ler670
  real, dimension(:,:), intent(in)    ::  toa1380
  real, dimension(:,:), intent(in)    ::  toa2250
  real, dimension(:,:), intent(in)    ::  bt11
  integer, dimension(:,:), intent(in) ::  gzone
  integer, dimension(:,:), intent(in) ::  lc
  real, dimension(:,:), intent(in)    ::  sr650
  
  integer, dimension(:,:), intent(inout) ::  mask
  
  real                                ::  rat488670
  real                                ::  rat412488
  integer     ::  i, j, i1, i2, j1, j2
  
  status = -1

  mask(:,:) = 0  
  do j = 1, size(mask,2)
    do i = 1, size(mask,1)
      if (toa1380(i,j) < -900.0) continue
      
      i1 = max(i-1,1)
      i2 = min(i+1,size(mask,1))
      j1 = max(j-1,1)
      j2 = min(j+1,size(mask,2))

      rat488670 = ler488(i,j) / ler670(i,j)
      rat412488 = ler412(i,j) / ler488(i,j)
    
!       if (ler412(i,j) > 12.0 .AND. count(bt11(i1:i2,j1:j2) > 286.0) == 9) then
      if (ler412(i,j) > 12.0 ) then     
        if ((sr650(i,j) < 0.08 .AND. count(toa2250(i1:i2,j1:j2) < 0.035) == 9) .OR. &
        & (sr650(i,j) >= 0.08 .AND. count(toa2250(i1:i2,j1:j2) < 0.065) == 9)) then 
          if (toa2250(i,j) < 0.025) then
            mask(i,j) = 1
          else 
            if (ler412(i,j) > 20.0 .AND. toa2250(i,j) < 0.04) then
              mask(i,j) = 1
            else
              if (ler488(i,j) > 0.0 .AND. ler670(i,j) > 0.0) then
                rat488670 = ler488(i,j) / ler670(i,j)
                if (rat488670 > 0.88) then
                  mask(i,j) = 1
                else 
                  if (rat488670 > 0.7 .AND. rat488670 < 0.88 .AND. toa2250(i,j) < 0.04  &
                  &   .AND. toa1380(i,j) > 0.0015) then
                    mask(i,j) = 1
                  end if
                end if
              end if
            end if    
          end if
        end if
      end if

     if (sr650(i,j) >= 0.12 .OR. ler670(i,j) > 50.0) mask(i,j) = 0
     if (lc(i,j) < 1 .OR. lc(i,j) > 4) mask(i,j) = 0
     if (rat412488/rat488670 > 0.94) mask(i,j) = 0
    end do
  end do
  
  status = 0
  return

end function calc_smoke_mask

integer function calc_pyrocb_mask(ler412, ler470, ler650, toa1380, toa2100,bt8, mask) result(status)	
  implicit none
  
  real, dimension(:,:), intent(in)    ::  ler412
  real, dimension(:,:), intent(in)    ::  ler470
  real, dimension(:,:), intent(in)    ::  ler650
  real, dimension(:,:), intent(in)    ::  toa1380
  real, dimension(:,:), intent(in)    ::  toa2100
  real, dimension(:,:), intent(in)    ::  bt8
  integer, dimension(:,:), intent(inout) ::  mask
  
  real                                ::  rat412470
  integer     ::  i, j
  
  status = -1
  
  mask(:,:) = 0  
  do j = 1, size(mask,2)
    do i = 1, size(mask,1)
     if (toa1380(i,j) < -900.0) continue
    
    rat412470 = ler412(i,j) / ler470(i,j)
    if (toa1380(i,j) > 0.06 .AND. toa2100(i,j) > 0.2 .AND. bt8(i,j) < 268.0 .AND. &
    &   ler650(i,j) > 24.0 .AND. rat412470 < 0.7) then
      mask(i,j) = 1
    end if  
    
    end do
  end do
  
  status = 0
  return

end function calc_pyrocb_mask

integer function calc_high_alt_smoke_mask(ler412, ler470, ler650, toa1380, toa2100, bt11, gzone, sr650, lc, mask) result(status) 
  implicit none
  
  real, dimension(:,:), intent(in)    ::  ler412
  real, dimension(:,:), intent(in)    ::  ler470
  real, dimension(:,:), intent(in)    ::  ler650
  real, dimension(:,:), intent(in)    ::  toa1380
  real, dimension(:,:), intent(in)    ::  toa2100
  real, dimension(:,:), intent(in)    ::  bt11
  real, dimension(:,:), intent(in)    ::  sr650
  
  integer, dimension(:,:), intent(in) ::  gzone
  integer, dimension(:,:), intent(inout) ::  lc
  integer, dimension(:,:), intent(inout) ::  mask
  
  real                                ::  rat470650
  real                                ::  rat412470
  
  integer     ::  i, j, i1, i2, j1, j2
  
  status = -1
  
  mask(:,:) = 0  
  do j = 1, size(mask,2)
    do i = 1, size(mask,1)
      if (toa1380(i,j) < -900.0) continue
      
      i1 = max(i-1,1)
      i2 = min(i+1,size(mask,1))
      j1 = max(j-1,1)
      j2 = min(j+1,size(mask,2))
      
      rat470650 = ler470(i,j) / ler650(i,j)
      rat412470 = ler412(i,j) / ler470(i,j)
      
     if (ler412(i,j) > 12.0 .AND. count(bt11(i1:i2,j1:j2) > 280.0) == 9) then
       if ((sr650(i,j) < 0.08 .AND. count(toa2100(i1:i2,j1:j2) < 0.25) == 9) .OR. &
       & (sr650(i,j) >= 0.08 .AND. count(toa2100(i1:i2,j1:j2) < 0.25) == 9)) then 
      
          if (ler470(i,j) > 18.0 .AND. (rat470650 > 0.7 .AND. rat470650 < 0.88) .AND. &
          &   toa2100(i,j) < 0.22 .AND. toa1380(i,j) > 0.0032) then
            mask(i,j) = 1
          end if

          if (ler470(i,j) > 18.0 .AND. (rat470650 > 0.7 .AND. rat470650 < 0.88) .AND. &
          &   toa2100(i,j) < 0.22 .AND. rat412470 < 0.8) then
            mask(i,j) = 1
          end if

          if (ler470(i,j) > 18.0 .AND. toa1380(i,j) > 0.0015 .AND. &
          &   bt11(i,j) > 268.0 .AND. rat412470 < 0.8) then
            mask(i,j) = 1
          end if

!     -- smoke detected only over Yuma (zone 31)
          if (gzone(i,j) == 31 ) then
          if (ler470(i,j) > 18.0 .AND. (rat470650 > 0.9 .AND. rat470650 < 1.06) .AND. &
          &   (rat412470 > 0.7 .AND. rat412470 < 0.88) .AND. &
          &   toa2100(i,j) < 0.22 .AND. toa1380(i,j) > 0.0032) then
            mask(i,j) = 1
          end if

          endif

          
        end if
      end if

     if (sr650(i,j) >= 0.12 .OR. toa2100(i,j) > 0.06 .OR. ler650(i,j)> 50.0 .OR. &
     &   toa1380(i,j) < 0.0005 .OR. toa1380(i,j) > 0.013) mask(i,j) = 0
     if (rat412470/rat470650 > 0.94) mask(i,j) = 0
     if (lc(i,j) < 1 .OR. lc(i,j) > 4) mask(i,j) = 0
    end do
  end do
  
  status = 0
  return
  
end function calc_high_alt_smoke_mask

! -- Detect and flag smoke over Australia and ConUS.
! -- Should be run after the cloud masking -- these tests don't discriminate
! -- between cloud and smoke and assume the clouds have already been removed. 
! -- Only meant to be used to reset the AE to 1.8
integer function calc_smoke_ae_mask(ler412, ler488, ler670, toa1380,  &
&                                dstar, gzone, mask) result(status)
  implicit none
  
  real, dimension(:,:), intent(in)    ::  ler412
  real, dimension(:,:), intent(in)    ::  ler488
  real, dimension(:,:), intent(in)    ::  ler670
  real, dimension(:,:), intent(in)    ::  toa1380
  real, dimension(:,:), intent(in)    ::  dstar
  integer, dimension(:,:), intent(in) ::  gzone
  
  integer, dimension(:,:), intent(inout) ::  mask
  
  real                                ::  rat488670
  integer     ::  i, j
  
  status = -1
  
  mask(:,:) = 0  
  do j = 1, size(mask,2)
    do i = 1, size(mask,1)
      if (toa1380(i,j) < -900.0) continue
      
!     -- smoke detected only over ConUS (zone 13 and 18)
      if (gzone(i,j) == 13 .OR. gzone(i,j) == 18) then 

!       -- LER 488/670 ratio test over ConUS
        if (ler488(i,j) > 0.0 .AND. ler670(i,j) > 0.0 .AND. dstar(i,j) > -900.0) then
          rat488670 = ler488(i,j) / ler670(i,j)
          if (rat488670 > 0.82 .AND. dstar(i,j) < 1.06) then
            mask(i,j) = 1
          end if
        end if
      end if
      
!     -- smoke test over Australia
      if (gzone(i,j) == 12) then  
        if (ler488(i,j) > 0.0 .AND. ler670(i,j) > 0.0 .AND. dstar(i,j) > -900.0) then
          rat488670 = ler488(i,j) / ler670(i,j)
          if (rat488670 > 0.80 .AND. dstar(i,j) < 1.06) then
            mask(i,j) = 1
          end if
        end if
      end if  
    end do
  end do
  
  status = 0
  return

end function calc_smoke_ae_mask


! -- calculate calibration correction factors for each band for VIIRS.
! -- good status = 0, bad status = -1.
type(viirs_calib_correction) function calc_calibration_correction(yy, mo, dd, status) result(vcc)
  use calendars
  
  implicit none
  
  integer, intent(in)   ::  yy    ! year
  integer, intent(in)   ::  mo    ! month
  integer, intent(in)   ::  dd    ! day of month  
  integer, intent(inout)::  status
  
  type(gdatetime)       ::  epoc_date 
  type(gdatetime)       ::  curr_date 
  type(datetime)        ::  diff
  type(datetime)        ::  dt1, dt2
  
  real                  ::  nyears
  
  status = -1
  
  epoc_date = gdatetime(2010, 1, 1, 0, 0, 0, 0, 0)
  dt1 = fixed_from_gregorian(epoc_date)
  
  curr_date = gdatetime(yy, mo, dd, 0, 0, 0, 0, 0)
  dt2 = fixed_from_gregorian(curr_date)
  
  diff = dt2 - dt1
  nyears = diff%day / 365.25
  print *, 'nyears: ', nyears

  ! Set all to 1 to turn calibration corr off
!  vcc%m01 = 1.000
!  vcc%m02 = 1.000
!  vcc%m03 = 1.000
!  vcc%m04 = 1.000
!  vcc%m05 = 1.000
!  vcc%m06 = 1.000
!  vcc%m07 = 1.000
!  vcc%m08 = 1.000
!  vcc%m10 = 1.000
!  vcc%m11 = 1.000

  ! Uncomment to turn calibration correction on. See Sayer et al AMT 2017 for source of numbers.
  ! Table 3 http://www.atmos-meas-tech.net/10/1425/2017/amt-10-1425-2017.pdf 
  vcc%m01 = 0.995
  vcc%m02 = 1.000
  vcc%m03 = 0.992
  vcc%m04 = 0.956
  vcc%m05 = 0.941
  vcc%m06 = 0.966
  vcc%m07 = 0.9544 + 0.0018*nyears
  vcc%m08 = 1.003  + 0.0019*nyears
  vcc%m10 = 0.9646 + 0.0035*nyears
  vcc%m11 = 0.931

  status = 0
  return
  
end function  

integer function rescale_array(orig_arr, res_arr, factor) result(status)
  implicit none
  
  real, dimension(:,:), intent(in)      ::  orig_arr
  real, dimension(:,:), intent(inout)   ::  res_arr
  integer, intent(in)                   ::  factor
  
  integer, dimension(2)                 ::  dims2
  
  integer                               ::  i, j
  integer                               ::  ii, jj
  
  dims2 = (/size(orig_arr,1), size(orig_arr,2)/)
  if (mod(dims2(1), factor) /= 0 .OR. mod(dims2(2), factor) /= 0) then
    print *, "ERROR: Array dimensions must be evenly divisible by factor: ", factor
    status = -1
    return
  end if

  dims2 = (/size(res_arr,1), size(res_arr,2)/) 
  do j = 1, dims2(2)
    do i = 1, dims2(1)
      
      ii = ((i-1) * factor) + 1
      jj = ((j-1) * factor) + 1
      
      res_arr(i,j) = sum(orig_arr(ii:ii+(factor-1),jj:jj+(factor-1))) / (factor*factor)
    end do
  end do
  
  status = 0
  return
end function rescale_array  
  
end module

