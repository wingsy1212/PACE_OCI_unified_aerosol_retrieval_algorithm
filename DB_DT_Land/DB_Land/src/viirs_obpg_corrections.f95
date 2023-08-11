! --  module to contain the various functions related to the Ocean Biology Processing
! --  group corrections to the MODIS TOA reflectance data. Currently only
! --  applied to the 412nm and 470nm bands. 

module viirs_obpg_corrections
  integer, parameter  :: R_DBL 	= kind(1.0d0)
  integer, parameter  :: I_SNGL 	= selected_int_kind(4)
  
  real, parameter, private     :: d2r = 0.017453292519943295
  real, parameter, private     :: pi=3.141592653589793238462643383279502884197
  
  type  ::  xcal_table
    integer                                       ::  ncoeffs
    integer                                       ::  ndets
    integer                                       ::  nmirror
    integer                                       ::  ndates
    real(R_DBL), dimension(:), allocatable        ::  secs
    real(R_DBL), dimension(:,:,:,:), allocatable  ::  m11, m12, m13
  end type xcal_table

  type  ::  gain_table
    integer                                       ::  ndates
    real(R_DBL), dimension(:), allocatable        ::  secs
    real, dimension(:), allocatable               ::  gain
  end type gain_table
  
  type  ::  rayl_table
    integer                               ::  nstokes
    integer                               ::  nsza
    integer                               ::  nvza
    integer                               ::  nraa
    real, dimension(:), allocatable       ::  sza_nodes 
    real, dimension(:), allocatable       ::  vza_nodes
    real, dimension(:), allocatable       ::  raa_nodes
    real, dimension(:,:,:,:), allocatable ::  raylrad
  end type rayl_table
  
  type(xcal_table)                          ::  xcal412, xcal470
  type(gain_table)                          ::  gain412, gain470
  type(rayl_table)                          ::  rayl412, rayl470, rayl650
  
  contains
  
! -- load_obpg_xcal()
! -- read OBPG cross-calibration tables and store in xcal_table type.
! -- 0 = success, -1 = failure
  integer function load_obpg_xcal_tables(xcal412_tbl, xcal470_tbl) result(status)
    implicit none
    	
    character(len=*), intent(in)  ::  xcal412_tbl, xcal470_tbl
    
    logical                       ::  file_exists
    integer                       ::  i, j
    
    status = -1
    
    inquire(file=trim(xcal412_tbl), exist=file_exists, iostat=status)
    if (status /= 0) then
      print *, "ERROR: Unable to inquire about OBPG xcal table: ", status
      print *, "Missing table: ", trim(xcal412_tbl)
      return
    end if
    if (file_exists .eqv. .false.) then
      print *, "ERROR: OBPG xcal table does not exist. Please check the config file."
      print *, "Missing table: ", trim(xcal412_tbl)
      status = -1
      return
    end if
  
    inquire(file=trim(xcal470_tbl), exist=file_exists, iostat=status)
    if (status /= 0) then
      print *, "ERROR: Unable to inquire about OBPG xcal table: ", status
      print *, "Missing table: ", trim(xcal470_tbl)
      return
    end if
    if (file_exists .eqv. .false.) then
      print *, "ERROR: OBPG xcal table does not exist. Please check the config file."
      print *, "Missing table: ", trim(xcal470_tbl)
      status = -1
      return
    end if
    
    xcal412 = read_xcal_table(xcal412_tbl, status)
    if (status /= 0) then
      print *, "ERROR: Failed to read 412 xcal table: ", status
      return
    end if
    
    xcal470 = read_xcal_table(xcal470_tbl, status)  
    if (status /= 0) then
      print *, "ERROR: Failed to read 470 xcal table: ", status
      return
    end if

    status = 0
    return
    
  end function load_obpg_xcal_tables
  
! -- load_obpg_gain()
! -- read OBPG-based vicarious gain tables and store in gain_table type.
! -- 0 = success, -1 = failure
! -- NOT CURRENTLY USED, NO POLARIZATION CORRECTION YET FROM OBPG.
  integer function load_obpg_gain_tables(gain412_tbl, gain470_tbl) result(status)
    implicit none
    	
    character(len=*), intent(in)  ::  gain412_tbl, gain470_tbl
    
    logical                       ::  file_exists
    integer                       ::  i, j
    
    status = -1
    
    inquire(file=trim(gain412_tbl), exist=file_exists, iostat=status)
    if (status /= 0) then
      print *, "ERROR: Unable to inquire about OBPG gain table: ", status
      print *, "Missing table: ", trim(gain412_tbl)
      return
    end if
    if (file_exists .eqv. .false.) then
      print *, "ERROR: OBPG gain table does not exist. Please check the config file."
      print *, "Missing table: ", trim(gain412_tbl)
      status = -1
      return
    end if
  
    inquire(file=trim(gain470_tbl), exist=file_exists, iostat=status)
    if (status /= 0) then
      print *, "ERROR: Unable to inquire about OBPG gain table: ", status
      print *, "Missing table: ", trim(gain470_tbl)
      return
    end if
    if (file_exists .eqv. .false.) then
      print *, "ERROR: OBPG gain table does not exist. Please check the config file."
      print *, "Missing table: ", trim(gain470_tbl)
      status = -1
      return
    end if
    
    gain412 = read_gain_table(gain412_tbl, status)
    if (status /= 0) then
      print *, "ERROR: Failed to read 412 gain table: ", status
      return
    end if
    
    gain470 = read_gain_table(gain470_tbl, status)  
    if (status /= 0) then
      print *, "ERROR: Failed to read 470 gain table: ", status
      return
    end if
    
    status = 0
    return
    
  end function load_obpg_gain_tables

! --  load_obpg_rayl_tables()
! --  read in the Rayleigh radiance tables create by MJ for the OBPG polarization
! --  correction. Data is stored in the rayl_table objects.
! -- NOT CURRENTLY USED, NO POLARIZATION CORRECTION YET FROM OBPG.
  integer function load_obpg_rayl_tables(rayl412_tbl, rayl470_tbl, rayl650_tbl) &
    result(status)
    
    implicit none
    
    character(len=*), intent(in)            ::  rayl412_tbl
    character(len=*), intent(in)            ::  rayl470_tbl
    character(len=*), intent(in)            ::  rayl650_tbl
    
    status = -1
    
    status = check_file_exists(trim(rayl412_tbl))
    if (status /= 0) then
      print *, "ERROR: Input OBPG 412 nm Rayleigh table does not exist: ", status
      return
    end if
    
    status = check_file_exists(trim(rayl470_tbl))
    if (status /= 0) then
      print *, "ERROR: Input OBPG 470 nm Rayleigh table does not exist: ", status
      return
    end if
    
    status = check_file_exists(trim(rayl650_tbl))
    if (status /= 0) then
      print *, "ERROR: Input OBPG 650 nm Rayleigh table does not exist: ", status
      return
    end if
    
    rayl412 = read_rayl_table(rayl412_tbl, status)
    if (status /= 0) then
      print *, "ERROR: Failed to read in 412nm Rayleigh radiance table: ", status
      return
    end if
    
    rayl470 = read_rayl_table(rayl470_tbl, status)
    if (status /= 0) then
      print *, "ERROR: Failed to read in 470nm Rayleigh radiance table: ", status
      return
    end if
    
    rayl650 = read_rayl_table(rayl650_tbl, status)
    if (status /= 0) then
      print *, "ERROR: Failed to read in 650nm Rayleigh radiance table: ", status
      return
    end if
    
    status = 0
    return
  end function load_obpg_rayl_tables

! --  calc_rvs412_corr()
! --  Uses the OBPG 412nm cross-calibration tables to calculate the RVS correction
! --  factor. Returned via rvs412.
!--   0 = success, -1 = failure.
  integer function calc_rvs412_corr(mirror_side, ev_cntr_time, rvs412, p_offset) result(status)
    implicit none
		
		integer, dimension(:), intent(in) 	::	mirror_side
		real(R_DBL), dimension(:), intent(in)			  ::	ev_cntr_time
		real(R_DBL), dimension(:,:), intent(inout)	::	rvs412 
  	integer, dimension(:), intent(in), optional ::  p_offset

		status = calc_rvs_corr(xcal412, mirror_side, ev_cntr_time, rvs412, p_offset)
		if (status /= 0) then
		  print *, "ERROR: Failed to calculate 412 nm RVS correction: ", status
		  return
		end if
		
		return
		
  end function calc_rvs412_corr
  
! --  calc_rvs488_corr()
! --  Uses the OBPG 470nm cross-calibration tables to calculate the RVS correction
! --  factor. Returned via rvs470.
!--   0 = success, -1 = failure.
  integer function calc_rvs488_corr(mirror_side, ev_cntr_time, rvs488, p_offset) result(status)
    implicit none
		
		integer, dimension(:), intent(in) 	::	mirror_side
		real(R_DBL), dimension(:), intent(in)			  ::	ev_cntr_time
		real(R_DBL), dimension(:,:), intent(inout)	::	rvs488 
  	integer, dimension(:), intent(in), optional ::  p_offset


		status = calc_rvs_corr(xcal470, mirror_side, ev_cntr_time, rvs488, p_offset)
		if (status /= 0) then
		  print *, "ERROR: Failed to calculate 488 nm RVS correction: ", status
		  return
		end if
		
		return
		
  end function calc_rvs488_corr
  
! --  calc_pol412_corr()
! --  Uses the OBPG 412nm cross-calibration tables to calculate the polarization 
! --  correction factor. Returned via pol412.
!--   0 = success, -1 = failure.
! -- NOT CURRENTLY USED, NO POLARIZATION CORRECTION YET FROM OBPG.
  integer function calc_pol412_corr(lat, lon, mirror_side, ev_cntr_time, sza, saa, vza, &
  &    vaa, raa, t_inst2ecr, toa412, polz, p_offset) result(status)
    implicit none
    
    real, dimension(:,:), intent(in)            ::  lat, lon
    integer, dimension(:), intent(in)  	::	mirror_side
		real(R_DBL), dimension(:), intent(in)	  		::	ev_cntr_time
    real, dimension(:,:), intent(in)            ::  sza, saa
    real, dimension(:,:), intent(in)            ::  vza, vaa
    real, dimension(:,:), intent(in)            ::  raa
    real(R_DBL), dimension(:,:), intent(in)     ::  t_inst2ecr
    real, dimension(:,:), intent(in)            ::  toa412
		real, dimension(:,:), intent(inout)	        ::	polz
		integer, dimension(:), intent(in), optional ::  p_offset
		
    status = -1
    
    status = calc_pol_corr(xcal412, rayl412, lat, lon, mirror_side, ev_cntr_time, sza, saa, vza, &
    &    vaa, raa, t_inst2ecr, toa412, polz, p_offset)
  	if (status /= 0) then
  	  print *, "ERROR: Failed to calculate polarization correction factor for 412nm: ", status
  	  return
  	end if
  	
  	return
		
  end function calc_pol412_corr

! --  calc_pol88_corr()
! --  Uses the OBPG 488nm cross-calibration tables to calculate the polarization 
! --  correction factor. Returned via pol470.
!--   0 = success, -1 = failure.
! -- NOT CURRENTLY USED, NO POLARIZATION CORRECTION YET FROM OBPG.
  integer function calc_pol488_corr(lat, lon, mirror_side, ev_cntr_time, sza, saa, vza, &
  &    vaa, raa, t_inst2ecr, toa470, polz, p_offset) result(status)
    implicit none
    
    real, dimension(:,:), intent(in)            ::  lat, lon
    integer, dimension(:), intent(in)  	::	mirror_side
		real(R_DBL), dimension(:), intent(in)	  		::	ev_cntr_time
    real, dimension(:,:), intent(in)            ::  sza, saa
    real, dimension(:,:), intent(in)            ::  vza, vaa
    real, dimension(:,:), intent(in)            ::  raa
    real(R_DBL), dimension(:,:), intent(in)     ::  t_inst2ecr
    real, dimension(:,:), intent(in)            ::  toa470
		real, dimension(:,:), intent(inout)	        ::	polz
		integer, dimension(:), intent(in), optional ::  p_offset
		
    status = -1
    
    status = calc_pol_corr(xcal470, rayl470, lat, lon, mirror_side, ev_cntr_time, sza, saa, vza, &
    &    vaa, raa, t_inst2ecr, toa470, polz, p_offset)
  	if (status /= 0) then
  	  print *, "ERROR: Failed to calculate polarization correction factor for 470nm: ", status
  	  return
  	end if
  	
  	return
		
  end function calc_pol488_corr
  
  real function get_obpg_gain412(ev_cntr_time, status) result(gain)
    implicit none
    
    real(R_DBL), dimension(:), intent(in)	  	::	ev_cntr_time
    integer, intent(inout)                    ::  status
    
    status = -1
    
    gain = get_obpg_gain(gain412, ev_cntr_time, status)
    if (status /= 0) then
      print *, "ERROR: Failed to OBGP gain coefficient for 412nm: ", status
      return
    end if
    
    status = 0
    return
  end function get_obpg_gain412

  real function get_obpg_gain470(ev_cntr_time, status) result(gain)
    implicit none
    
    real(R_DBL), dimension(:), intent(in)	  	::	ev_cntr_time
    integer, intent(inout)                    ::  status
    
    status = -1
    
    gain = get_obpg_gain(gain470, ev_cntr_time, status)
    if (status /= 0) then
      print *, "ERROR: Failed to OBGP gain coefficient for 470nm: ", status
      return
    end if
        
    status = 0
    return
  end function get_obpg_gain470
  
  real function get_obpg_gain(gain_tbl, ev_cntr_time, status) result(gain)
    implicit none
      
    type(gain_table), intent(in)            ::  gain_tbl
    real(R_DBL), dimension(:), intent(in)	  ::	ev_cntr_time
    integer, intent(inout)                  ::  status
    
 		integer				                          ::	t_index1, t_index2
 		real                                    ::  wt

    status = -1
    gain = -999.0
    
!   -- use the first time of the granule to determine gain
!   -- if xcal_locate status is non-zero:
!   -- -1, date is before any available dates in the file. Set gain to 1.0 and reset
!   -- status.
!   -- 1, date is after any dates in the table. Set gain to last value available and return.
!   -- else, bomb out.
    wt = -999.0
    t_index1 = -999
    t_index2 = -999
    t_index1 = xcal_locate(gain_tbl%secs, ev_cntr_time(1), status)
    t_index2 = t_index1 + 1
    if (status == 0) then
      wt = (ev_cntr_time(1)-gain_tbl%secs(t_index1))/(gain_tbl%secs(t_index2)-gain_tbl%secs(t_index1))
		  gain = gain_tbl%gain(t_index1) + wt*(gain_tbl%gain(t_index2) - gain_tbl%gain(t_index1))		
      status = 0    
    else
      if (status == -1) then    
        gain = 1.0
        status = 0
      else if (status == 1) then
        gain = gain_tbl%gain(size(gain_tbl%gain))
        status = 0
      else
        print *, "ERROR: Failed to find time interpolation index: ", status
        status = -1
      end if
    end if
        
    return
    
  end function get_obpg_gain
  
! --  read_xcal_table()
! --  read in the cross calibration tables from OBPG into the xcal_table 
! --  object xcal and return.
! --  0 = success, -1 = failure.
  type(xcal_table) function read_xcal_table(tbl_file, status) result(xcal)
    implicit none
    
    include 'hdf.inc'
    include 'dffunc.inc'
    
    character(len=*), intent(in)          ::  tbl_file
    integer, intent(inout)                ::  status
    
    integer(I_SNGL), dimension(:), allocatable    ::  tmp_years, tmp_days
    
    character(len=255)      ::  sds_name
    integer                 ::  sd_id, sds_index, sds_id
    integer, dimension(1)   ::  start1, stride1, edges1
    integer, dimension(4)   ::  start4, stride4, edges4

    integer                 ::  rank, nattrs, ntype
    integer, dimension(255) ::  dimsizes
    
    integer                 ::  i
    real                    ::  xhour
    
		sd_id = sfstart(tbl_file, DFACC_READ)
    if (sd_id == FAIL ) then
    	print *, "ERROR: failed to start SDS interface on xcal file: ", sd_id
      status = -1
      return
    end if
        
    sds_name = 'M11'
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
    
    status = sfginfo(sds_id, sds_name, rank, dimsizes, ntype, nattrs) 
    if (rank /= 4) then
      print *, "ERROR: Unexpected rank for xcal tables: ", rank
      status = -1
      return
    end if
    
    xcal%ncoeffs = dimsizes(4)
    xcal%ndets   = dimsizes(3)
    xcal%nmirror = dimsizes(2)
    xcal%ndates  = dimsizes(1)
    allocate( xcal%m11(xcal%ndates,xcal%nmirror,xcal%ndets,xcal%ncoeffs),     &
    &         xcal%m12(xcal%ndates,xcal%nmirror,xcal%ndets,xcal%ncoeffs),     &
    &         xcal%m13(xcal%ndates,xcal%nmirror,xcal%ndets,xcal%ncoeffs),     &
    &         xcal%secs(xcal%ndates), stat=status)
    if (status /= 0) then
      print *, "ERROR: Failed to allocate xcal table objects: ", status
      return
    end if
    
!   -- read in data from table
    start4  = (/0,0,0,0/)
    stride4 = (/1,1,1,1/)
    edges4  = dimsizes(1:rank)
    status = sfrdata(sds_id, start4, stride4, edges4, xcal%m11)
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
            
    sds_name = 'm12'
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
    
    start4  = (/0,0,0,0/)
    stride4 = (/1,1,1,1/)
    edges4  = dimsizes(1:rank)
    status = sfrdata(sds_id, start4, stride4, edges4, xcal%m12)
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
    
    sds_name = 'm13'
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
    
    start4  = (/0,0,0,0/)
    stride4 = (/1,1,1,1/)
    edges4  = dimsizes(1:rank)
    status = sfrdata(sds_id, start4, stride4, edges4, xcal%m13)
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
    
    allocate(tmp_years(xcal%ndates), tmp_days(xcal%ndates), stat=status)
    if (status /= 0) then
      print *, "ERROR: Unable to allocate years and days arrays: ", status
      return
    end if
    
    sds_name = 'year'
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
    
    start1  = (/0/)
    stride1 = (/1/)
    edges1  = xcal%ndates
    status = sfrdata(sds_id, start1, stride1, edges1, tmp_years)
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
    
    sds_name = 'day'
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
    
    start1  = (/0/)
    stride1 = (/1/)
    edges1  = xcal%ndates
    status = sfrdata(sds_id, start1, stride1, edges1, tmp_days)
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

!   -- convert years, days to seconds since Jan 1, 1993.
    xhour=12.0  ! assume m12, m13, rvs were derived at 12UTC
    do i = 1, xcal%ndates
      
      xcal%secs(i)  = yds_to_sec(tmp_years(i), tmp_days(i), xhour, status)
      if (status /= 0) then
        print *, "ERROR: Failed to convert year, day to seconds: ", &
        &   tmp_years(i), tmp_days(i), status
        return
      end if
    
    end do
    
    deallocate(tmp_years, tmp_days, stat=status)
    if (status /= 0) then
      print *, "ERROR: Failed to deallocate years and days arrays: ", status
      return
    end if
    
    status = sfend(sd_id)
    if (status /= 0) then
      print *, "ERROR: Unable to close land aerosol model file: ", status
      return
    end if
    
    return
  end function read_xcal_table
  
! --  read_gain_table()
! --  read in the cross calibration tables from OBPG into the xcal_table 
! --  object gain and return.
! --  0 = success, -1 = failure.
! -- NOT CURRENTLY USED, NO POLARIZATION CORRECTION YET FROM OBPG.
  type(gain_table) function read_gain_table(tbl_file, status) result(gain)
    implicit none
    
    include 'hdf.inc'
    include 'dffunc.inc'
    
    character(len=*), intent(in)          ::  tbl_file
    integer, intent(inout)                ::  status
    
    integer(I_SNGL), dimension(:), allocatable    ::  tmp_years, tmp_days
    
    character(len=255)      ::  sds_name
    integer                 ::  sd_id, sds_index, sds_id
    integer, dimension(1)   ::  start1, stride1, edges1

    integer                 ::  rank, nattrs, ntype
    integer, dimension(255) ::  dimsizes
    
    integer                 ::  i
    real                    ::  xhour
    
		sd_id = sfstart(tbl_file, DFACC_READ)
    if (sd_id == FAIL ) then
    	print *, "ERROR: failed to start SDS interface on gain file: ", sd_id
      status = -1
      return
    end if
    
    sds_name = 'year'
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
    
    status = sfginfo(sds_id, sds_name, rank, dimsizes, ntype, nattrs) 
    if (rank /= 1) then
      print *, "ERROR: Unexpected rank for gain tables: ", rank
      status = -1
      return
    end if

    gain%ndates  = dimsizes(1)
    
    allocate(tmp_years(gain%ndates), tmp_days(gain%ndates), stat=status)
    if (status /= 0) then
      print *, "ERROR: Unable to allocate years and days arrays: ", status
      return
    end if
    
    allocate( gain%secs(gain%ndates), gain%gain(gain%ndates), stat=status)
    if (status /= 0) then
      print *, "ERROR: Failed to allocate gain table objects: ", status
      return
    end if
    
    start1  = (/0/)
    stride1 = (/1/)
    edges1  = gain%ndates
    status = sfrdata(sds_id, start1, stride1, edges1, tmp_years)
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
    
    sds_name = 'day'
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
    
    start1  = (/0/)
    stride1 = (/1/)
    edges1  = gain%ndates
    status = sfrdata(sds_id, start1, stride1, edges1, tmp_days)
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

!   -- convert years, days to seconds since Jan 1, 1993.
    xhour=12.0  ! assume m12, m13, rvs were derived at 12UTC
    do i = 1, gain%ndates
      
      gain%secs(i)  = yds_to_sec(tmp_years(i), tmp_days(i), xhour, status)
      if (status /= 0) then
        print *, "ERROR: Failed to convert year, day to seconds: ", &
        &   tmp_years(i), tmp_days(i), status
        return
      end if
    
    end do
    
    deallocate(tmp_years, tmp_days, stat=status)
    if (status /= 0) then
      print *, "ERROR: Failed to deallocate years and days arrays: ", status
      return
    end if
    
    sds_name = 'gain'
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
    
    start1  = (/0/)
    stride1 = (/1/)
    edges1  = gain%ndates
    status = sfrdata(sds_id, start1, stride1, edges1, gain%gain)
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
      print *, "ERROR: Unable to close land aerosol model file: ", status
      return
    end if
    
    return
  end function read_gain_table
! -----------------------------------------------------------
!   Purpose: To convert year and DoY to seconds since
!           1993/01/01 00:00:00
!   Inputs: 
!     iyear: Year from xcal data file (xcal_modist.hdf) 
!     ijday: Day of year from xcal data file 
!     xhr  : hour (min, sec in decimal)
! 
!   Output:
!     xcalsec: seconds since 1993/01/01 
! 
!   ver. 1.0  Written by MJ (Aug 2008) 
!   ver. 2.0  Updated by Corey (Aug 2013)      
  real function yds_to_sec(yr, dy, hr, status) result(sec)
    implicit none
    
    integer(I_SNGL), intent(in)         ::  yr
    integer(I_SNGL), intent(in)         ::  dy
    real, intent(in)                    ::  hr
    integer, intent(inout)              ::  status
    
    
    integer, parameter          ::  start_year = 1993
      
    integer                     ::  is_leap, leap
    integer                     ::  ndays
    integer                     ::  iyr
        
!   -- check that year is later than 1993. None of this works for years < 1993.
    if (yr < 1993) then
      print *, "ERROR: Invalid year. Must be >= 1993: ", yr
      status = -1
      return
    end if
    
!   -- loop through years-1 adding up days. Exclude last year since it isn't
!   --  complete.
    sec   = 0.0
    ndays = 0      
    do iyr = start_year, yr-1
      is_leap = mod(iyr,4)
      if (is_leap == 0) then 
         leap=1
      else
         leap=0
      end if 

      ndays = ndays + 365 + leap
    end do

!   -- add in the days from the current year and convert to seconds.
    ndays = ndays + dy - 1
    
    sec = 86400.0D0*ndays + 3600.0D0*hr
    
    return

  end function yds_to_sec 
  
! --  read_rayl_table()
! --  Generic read routine for MJ's binary Rayleigh radiance tables for the OBPG
! --  polarization corrections. Data is returned in a rayl_table object, rayl.
! -- NOT CURRENTLY USED, NO POLARIZATION CORRECTION YET FROM OBPG.
  type(rayl_table) function read_rayl_table(tbl_file, status) result(rayl)
    implicit none
    
    character(len=*), intent(in)      ::  tbl_file
    integer, intent(inout)            ::  status

    integer                           ::  i
    
    status = -1
    
!   --  initialize our rayl_table object with node values. Hard-coded 
!   --  since these aren't documented directly in the files.
    rayl%nstokes  = 3
    rayl%nsza     = 10
    rayl%nvza     = 46
    rayl%nraa     = 31
    allocate(rayl%sza_nodes(rayl%nsza), rayl%vza_nodes(rayl%nvza),  &
    & rayl%raa_nodes(rayl%nraa), stat=status)
    if (status /= 0) then
      print *, "ERROR: Failed to allocate rayl_table members: ", status
      return
    end if
    
    allocate(rayl%raylrad(rayl%nstokes, rayl%nsza, rayl%nvza, rayl%nraa), stat=status)
    if (status /= 0) then
      print *, "ERROR: Failed to allocate Rayleigh radiance array: ", status
      return
    end if

    rayl%sza_nodes = (/0.0,8.0,16.0,24.0,32.0,40.0,48.0,56.0,64.0,72.0/)
    
    do i = 1, rayl%nvza
      rayl%vza_nodes(i) = 2.0*(i-1)
    end do

    do i = 1, rayl%nraa
      rayl%raa_nodes(i) = 6.0*(i-1)
    end do
  
!   -- read in the actual table.
    open(1,file=tbl_file, status='old') 
    read(1,*) rayl%raylrad
    close(1)
    
    status = 0
    return
  end function read_rayl_table
  
! --  check_file_exists()
! --  Helper function that checks if the file exists.
! --  0 = success, -1 = failure.
  integer function check_file_exists(in_file) result(status)
    implicit none
    
    character(len=*), intent(in)      ::  in_file
    
    logical                           ::  file_exists

    status = -1
    
    inquire(file=trim(in_file), exist=file_exists, iostat=status)
    if (status /= 0) then
      print *, "ERROR: Unable to inquire about file: ", status
      return
    end if
    
    if (file_exists .eqv. .false.) then
      status = -1
      print *, "ERROR: File does not exist: ", status
      return
    else
      status = 0
    end if
    
    return
  end function check_file_exists

! --  calc_pol_corr()
! --  Generic function to calculate polarization correction factors using the tables 
! --  in the specified xcal xcal_table object.
! --  0 = success, -1 = failure  
! -- NOT CURRENTLY USED, NO POLARIZATION CORRECTION YET FROM OBPG.
  integer function calc_pol_corr(xcal, rayl, lat, lon, mirror_side, ev_cntr_time, sza, saa, vza, &
  &    vaa, raa, t_inst2ecr, toa, polz, p_offset) result(status)
    implicit none
    
    type(xcal_table), intent(in)                ::  xcal
    type(rayl_table), intent(in)                ::  rayl
    real, dimension(:,:), intent(in)            ::  lat, lon
    integer, dimension(:), intent(in)  	        ::	mirror_side
		real(R_DBL), dimension(:), intent(in)	  		::	ev_cntr_time
    real, dimension(:,:), intent(in)            ::  sza, saa
    real, dimension(:,:), intent(in)            ::  vza, vaa
    real, dimension(:,:), intent(in)            ::  raa
    real(R_DBL), dimension(:,:), intent(in)     ::  t_inst2ecr
    real, dimension(:,:), intent(in)            ::  toa
		real, dimension(:,:), intent(inout)	        ::	polz
		integer, dimension(:), intent(in), optional ::  p_offset
		
		real(R_DBL)                                 ::  pn1, pn2, pn3, pn4, pn5
		real(R_DBL)                                 ::  x1, x2
		real(R_DBL)                                 ::  wt
    real(R_DBL), dimension(xcal%ncoeffs)				::	m12_1, m12_2
    real(R_DBL), dimension(xcal%ncoeffs)				::	m13_1, m13_2
    real(R_DBL)                                 ::  m12, m13
		real                                        ::  lq, lu, lx
		
		real(R_DBL)                                 ::  alpha
		
		real                                        ::  toa_nosza
		real                                        ::  cossza
		integer				                              ::	iscan, icell, iframe
		integer                                     ::  idet, ims, icoef
		integer				                              ::	t_index1, t_index2
    
    status = -1
		
		do iscan = 1, size(polz,2)
		  icell = (iscan-1)/16 + 1
			if (mirror_side(icell) < 0 .OR. mirror_side(icell) > 1) then 
     		print *, "WARNING: Invalid mirror side value: ", mirror_side(icell), iscan
     		cycle
      end if
      
!     --  find indexes in the xcal_table%sec array to interpolate between.
!     --  if ev_cntr_time > 1993 but before any correction nodes, 
!     --  (xcal_locate returns -1), set correction to 1.0 for all pixels in
!     --  this scan line. If xcal_locate returns 1, time is after all correction
!     --  nodes, use last available corrections!
      wt = -999.0
      t_index1 = -999
      t_index2 = -999
      t_index1 = xcal_locate(xcal%secs, ev_cntr_time(icell), status)
      t_index2 = t_index1 + 1
      if (status == 0) then 
   			wt = (ev_cntr_time(icell)-xcal%secs(t_index1))/(xcal%secs(t_index2)-xcal%secs(t_index1))
      else
        if (status == -1) then    
          polz(:,iscan) = 1.0
          cycle
        else if (status == 1) then
          t_index1 = size(xcal%m12,1) - 1
          t_index2 = size(xcal%m12,1)
          wt = 1.0
        else
          print *, "ERROR: Failed to find time interpolation index: ", status
          return
        end if
      end if
			
			idet = mod(iscan-1, 16) + 1    
			ims = mirror_side(icell) + 1		! convert mirror side value to array index
			m12_1(:) = -999.0
			m12_2(:) = -999.0
			m13_1(:) = -999.0
			m13_2(:) = -999.0
			
			do icoef = 1, xcal%ncoeffs

				m12_1(icoef)=xcal%m12(t_index1,ims,idet,icoef)
				m12_2(icoef)=xcal%m12(t_index2,ims,idet,icoef)

				m13_1(icoef)=xcal%m13(t_index1,ims,idet,icoef)
				m13_2(icoef)=xcal%m13(t_index2,ims,idet,icoef)
	              				
			end do
      
			do iframe = 1, size(polz,1)
				pn1 = float(iframe)-1.0D0 ! according to Bryan Franz 
        if (present(p_offset)) then   ! have to recalculate original pixel number for subsets
				  pn1 = float(iframe+p_offset(1)-1)-1.0D0
  			end if
				pn2 = pn1 * pn1
				pn3 = pn2 * pn1
				pn4 = pn3 * pn1
				pn5 = pn4 * pn1			
							               ! 				
				x1 = m12_1(1)+pn1*m12_1(2)+pn2*m12_1(3)+pn3*m12_1(4)+pn4*m12_1(5)+pn5*m12_1(6)
				x2 = m12_2(1)+pn1*m12_2(2)+pn2*m12_2(3)+pn3*m12_2(4)+pn4*m12_2(5)+pn5*m12_2(6)
				m12 = x1 + wt*(x2 - x1)	
				
				x1 = m13_1(1)+pn1*m13_1(2)+pn2*m13_1(3)+pn3*m13_1(4)+pn4*m13_1(5)+pn5*m13_1(6)
				x2 = m13_2(1)+pn1*m13_2(2)+pn2*m13_2(3)+pn3*m13_2(4)+pn4*m13_2(5)+pn5*m13_2(6)
				m13 = x1 + wt*(x2 - x1)	
        
!       -- calculate rotation angle between instrument reference plane and Stokes vector
!       -- See Gordon et al, "Atmospheric correction of ocean color sensors: 
!       -- analysis of the effects of residual instrument polarization sensitivity", 
!       -- http://dx.doi.org/10.1364/AO.36.006938 				
				alpha = calc_rotation_angle(lon(iframe,iscan),lat(iframe,iscan),vza(iframe,iscan), &
				&                           vaa(iframe,iscan),t_inst2ecr(:,icell),status)
				if (status /= 0) then
				  print *, "ERROR: Failed to calculate rotation angle: ", status
				  return
				end if

!       -- retrieve the Rayleigh radiance and Q and U components (lq, lu) of the 
!       --  polarization from the LUT's. Must multiply by PI to match the MODIS data. 
				status = get_rayl_radiance(rayl, sza(iframe,iscan), vza(iframe,iscan),  &
				& raa(iframe,iscan), lq, lu, lx)
				if (status /= 0) then
				  print *, "ERROR: Failed to get Rayleigh-only radiances from LUTs: ", status
				  return
				end if
								
				lq = lq * pi
				lu = lu * pi
          
!       -- Remove the cos(sza) factor from the MODIS TOA reflectances. Then plug
!       --  everything into the formula (see Meister et al, "Cross calibration of 
!       -- ocean-color bands from Moderate Resolution Imaging Spectradiometer on Terra
!       -- platorm.", Applied Optics, Volume 47, No. 36, December 20, 2008.
				cossza=cos(d2r*sza(iframe,iscan))
				toa_nosza = toa(iframe,iscan) * cossza
				polz(iframe,iscan) = calc_polz_factor(saa(iframe,iscan), vaa(iframe,iscan), m12,&
				&   m13, alpha, lq, lu, toa_nosza, status)
				if (status /= 0) then
				  print *, "ERROR: Failed to calculate polarization factor: ", status
				  return
				end if
				
			end do
		end do

		status = 0
		return
		
  end function calc_pol_corr

! --  calc_rvs_corr()
! --  Generic function to calculate RVS correction factors using the tables 
! --  in the specified xcal xcal_table object. Straight-forward process of fetching
! --  the appropriate coefficients from the table, calculating the RVS correction 
! --  factor, and interpolate that factor in time.
! --  0 = success, -1 = failure  
  integer function calc_rvs_corr(xcal, mirror_side, ev_cntr_time, rvs, p_offset)  &
  &   result(status)
    implicit none
    
		type(xcal_table), intent(in)                ::  xcal     
    integer, dimension(:), intent(in)  	::	mirror_side
		real(R_DBL), dimension(:), intent(in)	  		::	ev_cntr_time
		real(R_DBL), dimension(:,:), intent(inout)	::	rvs
		integer, dimension(:), intent(in), optional ::  p_offset        ! pixel offset from 
		                                                                ! extracts.
		
		real(R_DBL)                                 ::  pn1, pn2, pn3, pn4, pn5
		real(R_DBL)                                 ::  x1, x2
		real(R_DBL)                                 ::  wt
    real(R_DBL), dimension(xcal%ncoeffs)				::	rvs_1, rvs_2
		
		integer				::	iscan, icell, iframe, idet, ims, icoef
		integer				::	t_index1, t_index2
    
    status = -1
		
		do iscan = 1, size(rvs,2)
		  icell = (iscan-1)/16 + 1
!		  print *, 'iscan, icell, mirror side: ', iscan, icell, mirror_side(icell)
			if (mirror_side(icell) < 0 .OR. mirror_side(icell) > 1) then 
     		print *, "WARNING: Invalid mirror side value: ", mirror_side(icell), iscan
     		cycle
      end if
      
!     --  find indexes in the xcal_table%sec array to interpolate between.
!     --  if ev_cntr_time > 1993 but before any correction nodes, 
!     --  (xcal_locate returns -1), set correction to 1.0 for all pixels in
!     --  this scan line. If xcal_locate returns 1, time is after all correction
!     --  nodes, abort!
      wt = -999.0
      t_index1 = -999
      t_index2 = -999
      t_index1 = xcal_locate(xcal%secs, ev_cntr_time(icell), status)
      t_index2 = t_index1 + 1
      if (status == 0) then
        wt = (ev_cntr_time(icell)-xcal%secs(t_index1))/(xcal%secs(t_index2)-xcal%secs(t_index1))
      else
        if (status == -1) then    
          rvs(:,iscan) = 1.0
          cycle
        else if (status == 1) then
          t_index1 = size(xcal%m11,1) - 1
          t_index2 = size(xcal%m11,1)
          wt = 1.0
        else
          print *, "ERROR: Failed to find time interpolation index: ", status
          return
        end if
      end if
      
			idet = mod(iscan-1, 16) + 1    
			ims = mirror_side(icell) + 1		! convert mirror side value to array index
			rvs_1(:) = -999.0
			rvs_2(:) = -999.0
			do icoef = 1, xcal%ncoeffs
				rvs_1(icoef)=xcal%m11(t_index1,ims,idet,icoef)
				rvs_2(icoef)=xcal%m11(t_index2,ims,idet,icoef)
			end do
			
			do iframe = 1, size(rvs,1)
				pn1 = float(iframe)-1.0D0 ! according to Bryan Franz 
  			if (present(p_offset)) then   ! have to recalculate original pixel number for subsets
				  pn1 = float(iframe+p_offset(1)-1)-1.0D0
  			end if
				pn2 = pn1 * pn1
				pn3 = pn2 * pn1
				pn4 = pn3 * pn1
				pn5 = pn4 * pn1			
							          
				x1 = rvs_1(1)+pn1*rvs_1(2)+pn2*rvs_1(3)+pn3*rvs_1(4)+pn4*rvs_1(5)+pn5*rvs_1(6)
				x2 = rvs_2(1)+pn1*rvs_2(2)+pn2*rvs_2(3)+pn3*rvs_2(4)+pn4*rvs_2(5)+pn5*rvs_2(6)
				rvs(iframe,iscan) = x1 + wt*(x2 - x1)				  
			end do
		end do

		status = 0
		return
		
  end function calc_rvs_corr
  
! -- perform a binary search of array xx for j such that x lies between xx(j) and xx(j+1).
! -- see Numerical Recipes in Fortran, Second Edition, p.110
! -- returns  values:
! --     0 if xx(1) > x and size(xx) if xx(size(xx)) < x
! --    -1, 1 if x < xx(1) or x > xx(size(xx)) respectively.
! --    -2 if xx is not monotonic.
  integer function xcal_locate(xx,x,status) result(j)
    implicit none

    real(R_DBL), dimension(:), intent(in)    ::  xx
    real(R_DBL), intent(in)                  ::  x
    integer, intent(inout)                   ::  status

    integer                           ::  jl, jm, ju
    integer                           ::  i, n

    n = size(xx)
    status = 0

!		-- check that x is within xx array
		if (x < xx(1)) then
		  j = 1
			status = -1
			return
		end if
		
		if (x > xx(n)) then
		  j = n
		  status = 1
		  return
		end if
		
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
      end if
    end do

    j = jl

    return

  end function xcal_locate

! -----------------------------------------------------------------------------
! calc_rotation_angle() - computes the frame rotation angle for polarization correction 
!                                                                              
! adapted from Seadas src in 'l1_modis_hdf.c'  (<--Seadas5.1)                                 
!                                                                              
! arguments:	lon: longitude
! 	lat: latitude
!		senz:sensor zenith angle
!		sena:sensor azimuth angle
!		mnorm:taken from MODIS geolocation T_inst2ECR
! output: alpha
! modified by JCW
! version 1: keep same way as original file
!
! modified by MJ
! The sign of alpha was changed to conform Gordon/Meister papers 
! and to conform "l1_modis_hdf.c" in seadas5.2).
! -----------------------------------------------------------------------------
! -- NOT CURRENTLY USED, NO POLARIZATION CORRECTION YET FROM OBPG.
  real(R_DBL) function calc_rotation_angle (lon, lat, senz, sena, mnorm, status) result(alpha)
      implicit none

      real,  intent(in):: lon, lat, senz, sena
      real(R_DBL), dimension(3), intent(in):: mnorm
      integer, intent(inout)  ::  status
      
      real(R_DBL), parameter :: d2r_d = 0.017453292519943295
      
! local variables 
      integer		:: i
      real	:: slon, clon
      real	:: slat, clat
      real	:: szen, czen
      real	:: sazi, cazi

      real, dimension(3) :: e, n, v, r, s
      real    	 :: sdvcr, vdr, sdr, sdv 
      
      status = -1
      
!    /* invert mirror normal */
      do i = 1, 3
	 		r(i) = -mnorm(i)
      enddo

      slon = sin(d2r*lon) 
      clon = cos(d2r*lon) 
      slat = sin(d2r*lat) 
      clat = cos(d2r*lat) 
      szen = sin(d2r*senz) 
      czen = cos(d2r*senz) 
      sazi = sin(d2r*sena) 
      cazi = cos(d2r*sena) 

!     /* pixel coordinate system (north, east, vertical) in ECR */

      e(1) = -slon 
      e(2) =  clon 
      e(3) =  0.0 

      n(1) = -slat * clon 
      n(2) = -slat * slon 
      n(3) =  clat 

      v(1) =  clat * clon 
      v(2) =  clat * slon 
      v(3) =  slat 

!     /* sensor view vector in ECR */

      do  i = 1, 3
        s(i) = e(i) * szen * sazi         &
               + n(i) * szen * cazi        &
               + v(i) * czen 
      end do

!     /* compute rotation angle (alpha) from pixel normal (v) to mirror */
!     /* normal (r) about sensor view vector (s)  (Wertz, p. 757)       */

      sdvcr = s(1) * (v(2)*r(3)-v(3)*r(2)) &
             + s(2) * (v(3)*r(1)-v(1)*r(3)) &
             + s(3) * (v(1)*r(2)-v(2)*r(1)) 

      vdr   = v(1)*r(1) + v(2)*r(2) + v(3)*r(3) 
      sdr   = s(1)*r(1) + s(2)*r(2) + s(3)*r(3) 
      sdv   = v(1)*s(1) + v(2)*s(2) + v(3)*s(3) 

!      /* negated to be consistent with Gordon et al. */
!...   alpha = (1/d2r_d) * atan2(sdvcr,(vdr-sdr*sdv)) - 90.0  ! <--Seadas5.1

      alpha = -1.0D0 * ((1.0/d2r_d) * atan2(sdvcr,(vdr-sdr*sdv)) - 90.0) 
      
      status = 0
      return
      
    end function calc_rotation_angle

! --  calc_polz_factor()
! --  Straight-forward calculation of the polarization correction factor. Implements
! --  equation #4 in M-J Jeong et al, "Impacts of Cross-Platform Vicarious Calibration
! --  on the Deep Blue Aerosol Retrievals for Moderate Resolution Imaging Spectradiometer 
! --  Aboard Terra," 10.1109/TGRS.2011.2153205
! -- NOT CURRENTLY USED, NO POLARIZATION CORRECTION YET FROM OBPG.
  real function calc_polz_factor(saa, vaa, m12, m13, alpha, lq, lu, lx, status) result(polz_factor)
    implicit none
    
    real, intent(in)          ::  saa       ! solar azimuth angle
    real, intent(in)          ::  vaa       ! viewing azimuth angle
    real(R_DBL), intent(in)   ::  alpha	    ! polarization frame rotation angle
    real, intent(in)          ::  lq        ! Rayleigh radiance, q-component 
    real, intent(in)          ::  lu        ! Rayleigh radiance, u-component 
    real, intent(in)          ::  lx        ! normalization radiance
    real(R_DBL), intent(in)  :: 	m12, m13  ! coeff m12, m13
    integer, intent(inout)    ::  status

    real    ::  azimuth, tmp_lu    
    real    ::  lqu, luq

!   -- check azimuth angle and change sign on Lu if needed.
    azimuth = saa-vaa
		if (azimuth < 0.0) azimuth = 360.0 + azimuth
		
		if (azimuth < 180.0) then
		  tmp_lu = -1.0*lu
		else
		  tmp_lu = lu
		end if

!	Compute the polarization correction for each band.
    if (lx > 0.0) then
      lqu = lq*cos(2.0*alpha*d2r) + tmp_lu*sin(2.0*alpha*d2r)  
      luq = tmp_lu*cos(2.0*alpha*d2r) - lq*sin(2.0*alpha*d2r)  
      polz_factor = 1.0 / (1.0 - (m12*lqu)/lx - (m13*luq)/lx)  
	  else
      polz_factor = 1.0
	  end if
    
    status = 0
    return
    
	end function calc_polz_factor  

! --  get_rayl_radiance()
! --  Perform a table lookup using the tables from the rayl rayl_table input object
! --  and return the Q and U components (lq, lu) of the Rayleigh-only TOA radiances.  
! --  lx is the actual Rayleigh radiance.
! -- NOT CURRENTLY USED, NO POLARIZATION CORRECTION YET FROM OBPG.
  integer function get_rayl_radiance(rayl, sza, vza, raa, lq, lu, lx) result(status)
    implicit none
    
    type(rayl_table)          ::  rayl
    real, intent(in)          ::  sza, vza, raa
    real, intent(inout)       ::  lq, lu, lx
    
    real, dimension(3)        ::  yy
    real                      ::  y, dy
    
    integer                   ::  ia
    real                      ::  raa_fixed

!   --  interpolation code doesn't like relative azimuth angles == 180.0.
    raa_fixed = raa
    if(raa >= 180.0) raa_fixed = 179.9999
    do ia = 1, 3 
      call intep(rayl%sza_nodes, rayl%vza_nodes, rayl%raa_nodes, &
      &     rayl%raylrad, rayl%nsza, rayl%nvza, rayl%nraa, ia,sza,vza,raa_fixed,y,dy)
      yy(ia) = y/pi
    end do
    lx = yy(1)      
    lq = yy(2)
    lu = yy(3)
      
    status = 0
    return
    
    contains

!   --  intep()
!   --  helper function performs a 3D polynomial interpolation of input ya array.
!   --  For reference, check out Numerical Recipes.
    subroutine intep(x1a,x2a,x3a,ya,m,n,l,ia,x1,x2,x3,y,dy)
      real, dimension(:)    ::  x1a
      real, dimension(:)    ::  x2a
      real, dimension(:)    ::  x3a
      real, dimension(:,:,:,:)  ::  ya
      integer               ::  m, n, l, ia
      
      real                  ::  x1, x2, x3, y, dy
      
      integer, parameter    ::  nmax = 46
      integer, parameter    ::  mmax = 10
      integer, parameter    ::  lmax = 31
      
      real, dimension(4)    ::  xx2a, xx1a
      real, dimension(4)    ::  yntmp, ymtmp
      real, dimension(lmax) ::  yltmp
      
      integer               ::  i, j, k, o, ii
      real                  ::  dif, dift, frac
      integer               ::  mbeg, nbeg
      integer               ::  nsm, nsn
      
      call search3(x3,x3a,l,ii)
      frac = (x3-x3a(ii))/(x3a(ii+1)-x3a(ii))

      nsm=1
      dif=x1-x1a(1)
      do i=1,m
        dift=x1-x1a(i)
        if (dift.gt.0. .and. dift.lt.dif) then
          nsm=i
          dif=dift
        endif
      enddo
      mbeg = nsm - 2
      if (mbeg.le.0) then
        mbeg = 0
      else if (mbeg.gt.m-4) then
        mbeg = m-4
      endif

      nsn=1
      dif=x2-x2a(1)
      do i=1,n
        dift=x2-x2a(i)
        if (dift.gt.0. .and. dift.lt.dif) then
          nsn=i
          dif=dift
        endif
      enddo
      nbeg = nsn - 2
      if (nbeg.le.0) then
        nbeg = 0
      else if (nbeg.gt.n-4) then
        nbeg = n-4
      endif

      do 12 j=1,4
        do 11 k=1,4
          do 10 o=1,l
           yltmp(o)=ya(ia,j+mbeg,k+nbeg,o)
10        continue
          yntmp(k) = yltmp(ii)*(1.-frac) + yltmp(ii+1)*frac
          xx2a(k) = x2a(k+nbeg)
11      continue
        call polint2(xx2a,yntmp,4,x2,ymtmp(j),dy)
        xx1a(j) = x1a(j+mbeg)
12    continue
      call polint2(xx1a,ymtmp,4,x1,y,dy)
      return
      end subroutine intep

!     -- polint2()
!     -- polynomial interpolation routine, see Numerical Recipes.
      subroutine polint2(xa,ya,n,x,y,dy)
        real, dimension(:)    ::  xa, ya
        integer               ::  n
        real                  ::  x, y, dy
        
        integer, parameter    ::  nmax = 50
        real, dimension(nmax) ::  c, d
        
        real      ::  den, dif, dift, ho, hp, w
        integer   ::  i, m, ns
        
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
      end subroutine polint2

      subroutine  search3(xbar,x,n,i)
!
!     purpose
!       locate position in table of point at which interpolation is
!       required
!
!     usage
!       call  search3 (xbar, x, n, i)
!
!     description of parameters
!       xbar   - point at which interpolation is required
!       x      - array containing independent variable
!       n      - number of points in x array
!       i      - index specifying segment containing xbar
!
      real                ::  xbar
      real, dimension(:)  ::  x
      integer             ::  n, i
      
      real, parameter     ::  b = 0.69314718
      integer             ::  k, m
      
      if (n.lt.2) go to 15
      if(x(1).gt.x(2)) go to 17
      m = int((log(float(n)))/b)
      i=2**m
      k=i
   10 k=k/2
      if (xbar.ge.x(i).and.xbar.lt.x(i+1))return
      if (xbar.gt.x(i)) go to 12
      i = i-k
      go to 10
   12 i = i+k
      if (i.le.n) go to 10
      i=n
      go to 10
   15 print *, "Search n is less than 2."
      return
   17 print *, "Search table is not in increasing order."
      return
  100 format(5x,6a8)
      end subroutine search3
      
  end function get_rayl_radiance 

! -- NOT CURRENTLY USED, NO POLARIZATION CORRECTION YET FROM OBPG.
  subroutine calc_sensor_orientation(pos, vel, att, rm, coef)
    
    implicit none
        
    real, dimension(3), intent(in)        ::  pos
    real, dimension(3), intent(in)        ::  vel
    real, dimension(3), intent(in)        ::  att
    real, dimension(3,3), intent(inout)   ::  rm
    real, dimension(10), intent(inout)    ::  coef
    
    real(kind=R_DBL)                        ::  pm
    real(kind=R_DBL)                        ::  omf2p
    real(kind=R_DBL)                        ::  pxy
    real(kind=R_DBL)                        ::  rd
    real(kind=R_DBL)                        ::  yprim
    real(kind=R_DBL)                        ::  temp
    
    real, dimension(3)          ::  vc
    real, dimension(3)          ::  xpri
    real, dimension(3)          ::  ypri
    real, dimension(3)          ::  zpri
    real, dimension(3,3)        ::  yrp
    real, dimension(3,3)        ::  rn
    
    real(kind=R_DBL), parameter             ::  omf2 = (1.0d0 - (1.0/298.257d0))**2
    real(kind=R_DBL), parameter             ::  omegae = 7.29211585494e-05
    real(kind=R_DBL), parameter             ::  rem = 6371.0
    real(kind=R_DBL), parameter             ::  re = 6378.137
    
    vc(1) = vel(1) - omegae*pos(2)
    vc(2) = vel(2) - omegae*pos(1)
    vc(3) = vel(3)
!    print *, 'vc: ', vc
    
    pm = sqrt(dble(pos(1)*pos(1) + pos(2)*pos(2) + pos(3)*pos(3)))
    omf2p = (omf2*rem + pm - rem) / pm
    pxy = pos(1)*pos(1) + pos(2)*pos(2)
    temp = sqrt(dble(pos(3)*pos(3) + omf2p*omf2p*pxy))
    zpri(1) = -omf2p*pos(1) / temp
    zpri(2) = -omf2p*pos(2) / temp
    zpri(3) = -pos(3) / temp
    call crossp(vc, zpri, ypri)
    
    yprim = sqrt(dble(ypri(1)*ypri(1) + ypri(2)*ypri(2) + ypri(3)*ypri(3)))
    ypri(1) = -ypri(1) / yprim
    ypri(2) = -ypri(2) / yprim
    ypri(3) = -ypri(3) / yprim
    call crossp(ypri, zpri, xpri)
    
    rn(1,:) = xpri(:)
    rn(2,:) = ypri(:)
    rn(3,:) = zpri(:)
!    print *, 'xpri: ', xpri
!    print *, 'ypri: ', ypri
!    print *, 'zpri: ', zpri
    
    call oceuler(att, yrp)
!    print *, 'yrp: ', yrp
    
    call matmpy(yrp, rn, rm)
!    print *, 'rm: ', rm
    
    rd = 1.d0 / omf2
    coef(1) = 1.d0 + (rd-1.d0)*rm(1,3)*rm(1,3)
    coef(2) = 1.d0 + (rd-1.d0)*rm(2,3)*rm(2,3)
    coef(3) = 1.d0 + (rd-1.d0)*rm(3,3)*rm(3,3)
    coef(4) = (rd-1.d0)*rm(1,3)*rm(2,3)*2.d0
    coef(5) = (rd-1.d0)*rm(1,3)*rm(3,3)*2.d0
    coef(6) = (rd-1.d0)*rm(2,3)*rm(3,3)*2.d0
    coef(7) = (rm(1,1)*pos(1) + rm(1,2)*pos(2) + rm(1,3)*pos(3)*rd)*2.d0
    coef(8) = (rm(2,1)*pos(1) + rm(2,2)*pos(2) + rm(2,3)*pos(3)*rd)*2.d0
    coef(9) = (rm(3,1)*pos(1) + rm(3,2)*pos(2) + rm(3,3)*pos(3)*rd)*2.d0
    coef(10)= pos(1)*pos(1) + pos(2)*pos(2) + pos(3)*pos(3)*rd - re*re
!    print *, 'coef: ', coef
    
    return
    
  end subroutine calc_sensor_orientation  
  
  subroutine crossp(v1, v2, v3)
  
    implicit none
    
    real, dimension(3), intent(in)    ::  v1
    real, dimension(3), intent(in)    ::  v2
    real, dimension(3), intent(inout) ::  v3
    
    v3(1) = v1(2)*v2(3) - v1(3)*v2(2)
    v3(2) = v1(3)*v2(1) - v1(1)*v2(3)
    v3(3) = v1(1)*v2(2) - v1(2)*v2(1)
    
    return
  
  end subroutine crossp
  
  subroutine oceuler(a, xm)
    
    implicit none
    
    real, dimension(3,3), intent(inout)     ::  xm
    real, dimension(3), intent(in)          ::  a
    
    real, dimension(3,3)    :: xm1
    real, dimension(3,3)    :: xm2
    real, dimension(3,3)    :: xm3
    real, dimension(3,3)    :: xmm
      
    real                    ::  c1
    real                    ::  c2
    real                    ::  c3
    real                    ::  s1
    real                    ::  s2
    real                    ::  s3
    
    real(kind=R_DBL), parameter   ::  radeg = 180.0 / 3.14159265359
    
    xm1(:,:) = 0.0
    xm2(:,:) = 0.0
    xm3(:,:) = 0.0
    
    c1 = cos(a(1) / radeg)
    s1 = sin(a(1) / radeg)
    c2 = cos(a(2) / radeg)
    s2 = sin(a(2) / radeg)
    c3 = cos(a(3) / radeg)
    s3 = sin(a(3) / radeg)
    
    xm1(1,1) = 1.d0
    xm1(2,2) = c1
    xm1(3,3) = c1
    xm1(2,3) = s1
    xm1(3,2) = -s1
    xm2(2,2) = 1.d0
    xm2(1,1) = c2
    xm2(3,3) = c2
    xm2(3,1) = s2
    xm2(1,3) = -s2
    xm3(3,3) = 1.d0
    xm3(2,2) = c3
    xm3(1,1) = c3
    xm3(1,2) = s3
    xm3(2,1) = -s3
    
    call matmpy(xm2, xm3, xmm)
    call matmpy(xm1,xmm, xm)
    
    return
  
  end subroutine oceuler

! -- NOT CURRENTLY USED, NO POLARIZATION CORRECTION YET FROM OBPG.
  subroutine matmpy(xm1, xm2, xm3)
    
    implicit none
    
    real, dimension(3,3), intent(in)    :: xm1
    real, dimension(3,3), intent(in)    :: xm2
    real, dimension(3,3), intent(inout) :: xm3
    
    integer   ::  i, j, k
    
    xm3(:,:) = 0.0
    do i = 1, 3
      do j = 1, 3
        do k = 1, 3
          xm3(i,j) = xm3(i,j) + xm1(i,k)*xm2(k,j)
        end do
      end do
    end do
    
    return
    
    end subroutine matmpy
end module
