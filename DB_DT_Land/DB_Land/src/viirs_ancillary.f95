!-----------------------------------------------------------------
!  Module with a subroutine to read in wind speed, total ozone,
!  and precipitable water from an ncep GDAS grib file and
!  functions to get values from these fields at input  lat and lon 
!-----------------------------------------------------------------
module viirs_ancillary
  implicit none
  
  private 
  
  public  ::  load_geos5_data
  public  ::  get_winds
  public  ::  get_pwat
  public  ::  get_tozne
  public  ::  get_wind_dir
  public  ::  get_ps
 
  integer, parameter      ::  rdbl = selected_real_kind(p=13)
  
  type  ::  ncep_gdas_grib
    integer                             ::  ni
    integer                             ::  nj
    real(kind=rdbl), dimension(:,:), allocatable   ::  tozne
    real(kind=rdbl), dimension(:,:), allocatable   ::  pwat
    real(kind=rdbl), dimension(:,:), allocatable   ::  wndspd
    real(kind=rdbl), dimension(:,:), allocatable   ::  u_wndspd
    real(kind=rdbl), dimension(:,:), allocatable   ::  v_wndspd
    real(kind=rdbl), dimension(:,:), allocatable   ::  ps    
  end type ncep_gdas_grib
  
  type(ncep_gdas_grib)  ::  gdas
    
  contains
  
  
  integer function load_geos5_data(filename1, filename2, hr, min) result(status)
    implicit none
    
    character(len=255), intent(in)      ::  filename1
    character(len=255), intent(in)      ::  filename2
    integer, intent(in)                 ::  hr
    integer, intent(in)                 ::  min
    
    type(ncep_gdas_grib)                ::  gdas1
    type(ncep_gdas_grib)                ::  gdas2
                 
    real                                ::  hrs
    real                                ::  t1,time_interval
    
    integer                             ::  tmp
    
    status = -1
    time_interval =6.0
!   -- read in both pre-granule and post-granule NCEP GDAS GFS data.
!   -- read_ncep_gdas converts GDAS 1D GRIB data into 2D lat/lon array.
    gdas1 = read_geos5(filename1, status)    
    if (status /= 0) then
      print *, "ERROR: Failed to read GEOS5 ancillary data: ", status
      print *, "File: ", trim(filename1)
      return
    end if    
    
    gdas2 = read_geos5(filename2, status)
    if (status /= 0) then
      print *, "ERROR: Failed to read GEOS5 ancillary data: ", status
      print *, "File: ", trim(filename2)
      return
    end if
            
!   -- interpolate between the pre- and post-granule data based on hr and min values.
    hrs = hr + (min / 60.0)
    select case (hr)
    case (0:5)
      t1 = 0
    case (6:11)
      t1 = 6
    case (12:17)
      t1 = 12
    case (18:23)
      t1 = 18
    case default
      print *, "ERROR: Invalid hour detected: ", hr
      status = -1
      return
    end select
    
!   -- for NRT retrieval
    if (filename1(19:22) == 'fp.f') then
      print *, 'start NRT Ancillary process'
      time_interval =1.0      
      if (min >= 30) then
         t1 = hr + 0.5
      else
         t1 = hr - 0.5
      end if     
    endif
    
    allocate(gdas%u_wndspd(gdas1%ni, gdas1%nj), gdas%v_wndspd(gdas1%ni, gdas1%nj),  &
    & gdas%wndspd(gdas1%ni, gdas1%nj),gdas%ps(gdas1%ni, gdas1%nj), stat=status)
    if (status /= 0) then
      print *, "ERROR: Failed to allocate interpolated windspeed arrays: ", status
      return
    end if

    allocate(gdas%tozne(gdas1%ni, gdas1%nj), stat=status)
    if (status /= 0) then
      print *, "ERROR: Failed to allocate interpolated ozone array: ", status
      return
    end if
    
    allocate(gdas%pwat(gdas1%ni, gdas1%nj), stat=status)
    if (status /= 0) then
      print *, "ERROR: Failed to allocate interpolated preciptable water array: ", status
      return
    end if
    
    gdas%ni       = gdas1%ni
    gdas%nj       = gdas1%nj
    gdas%tozne    = gdas1%tozne + (gdas2%tozne-gdas1%tozne) * (hrs-t1)/time_interval
    gdas%pwat     = gdas1%pwat + (gdas2%pwat-gdas1%pwat) * (hrs-t1)/time_interval
    gdas%u_wndspd = gdas1%u_wndspd + (gdas2%u_wndspd-gdas1%u_wndspd) * (hrs-t1)/time_interval
    gdas%v_wndspd = gdas1%v_wndspd + (gdas2%v_wndspd-gdas1%v_wndspd) * (hrs-t1)/time_interval
    gdas%wndspd   = gdas1%wndspd + (gdas2%wndspd-gdas1%wndspd) * (hrs-t1)/time_interval
    gdas%ps       = gdas1%ps + (gdas2%ps-gdas1%ps) * (hrs-t1)/time_interval
    return
    
  end function load_geos5_data
!
! Getter method to fetch surface pressure from GRIB data for given latitude and longitude, lat and lon.
! Data is interpolated linearly from the surrounding 4 points from the GRIB sample data.
!
! Returns the windspeed, status: -1 = fail, 0 = success.
!---------------------------------------------------------------------------------------------------      
  real function get_ps(rlat, rlon, status,nrt) result(ps)
    implicit none
    
    real, intent(in)        :: rlat
    real, intent(in)        :: rlon
    integer, intent(out)    :: status
    integer, intent(in)     :: nrt
    
    real, dimension(4)      ::  f
    real                    ::  x1, x2
    real                    ::  y1, y2
    
    integer :: i1, i2
    integer :: j1, j2

    status = -1
    if(rlat > 90.0 .or. rlat < -90.0) then
      print*, 'ERROR: lat out of range must be 90 to -90:', rlat
      return
    end if
    if (rlon < -180.0 .OR. rlon > 180.0) then
      print*, 'ERROR: lon out of range must be [-180,180):', rlon
      return
    end if

    status = get_interp_indexes(rlat, rlon, i1, i2, j1, j2,nrt)
    if (status /= 0) then
      print *, "ERROR: Failed to get indexes of surrounding data for interpolation: ", status
      return
    end if
    
    if (i1 > 576 .OR. i2 > 576) then
      print *, 'i1, i2, lon: ', i1, i2, rlon
    end if    
    
!   -- perform 2D bilinear interpolation (http://en.wikipedia.org/wiki/Bilinear_interpolation)
    f = (/gdas%ps(i1,j1), gdas%ps(i2,j1), gdas%ps(i2,j2), gdas%ps(i1,j2)/)
    x1 = index2lon(i1, status,nrt)
    x2 = index2lon(i2, status,nrt)
    y1 = index2lat(j1, status,nrt)
    y2 = index2lat(j2, status,nrt)

    ps = f(1)*(x2-rlon)*(y2-rlat) + f(2)*(rlon-x1)*(y2-rlat) +  &
    &     f(3)*(rlon-x1)*(rlat-y1) + f(4)*(x2-rlon)*(rlat-y1)
    ps = 1.0 / ((x2-x1)*(y2-y1)) * ps
    
    status = 0
    return
    
  end function get_ps    
!
! Getter method to fetch windspeed from GRIB data for given latitude and longitude, lat and lon.
! Data is interpolated linearly from the surrounding 4 points from the GRIB sample data.
!
! Returns the windspeed, status: -1 = fail, 0 = success.
!---------------------------------------------------------------------------------------------------      
  real function get_winds(rlat, rlon, status,nrt) result(ws)
    implicit none
    
    real, intent(in)        :: rlat
    real, intent(in)        :: rlon
    integer, intent(out)    :: status
    integer, intent(in)     :: nrt
    
    real, dimension(4)      ::  f
    real                    ::  x1, x2
    real                    ::  y1, y2
    
    integer :: i1, i2
    integer :: j1, j2

    status = -1
    if(rlat > 90.0 .or. rlat < -90.0) then
      print*, 'ERROR: lat out of range must be 90 to -90:', rlat
      return
    end if
    if (rlon < -180.0 .OR. rlon > 180.0) then
      print*, 'ERROR: lon out of range must be [-180,180):', rlon
      return
    end if

    status = get_interp_indexes(rlat, rlon, i1, i2, j1, j2,nrt)
    if (status /= 0) then
      print *, "ERROR: Failed to get indexes of surrounding data for interpolation: ", status
      return
    end if
    
    if (i1 > 576 .OR. i2 > 576) then
      print *, 'i1, i2, lon: ', i1, i2, rlon
    end if    
    
!   -- perform 2D bilinear interpolation (http://en.wikipedia.org/wiki/Bilinear_interpolation)
    f = (/gdas%wndspd(i1,j1), gdas%wndspd(i2,j1), gdas%wndspd(i2,j2), gdas%wndspd(i1,j2)/)
    x1 = index2lon(i1, status,nrt)
    x2 = index2lon(i2, status,nrt)
    y1 = index2lat(j1, status,nrt)
    y2 = index2lat(j2, status,nrt)

    ws = f(1)*(x2-rlon)*(y2-rlat) + f(2)*(rlon-x1)*(y2-rlat) +  &
    &     f(3)*(rlon-x1)*(rlat-y1) + f(4)*(x2-rlon)*(rlat-y1)
    ws = 1.0 / ((x2-x1)*(y2-y1)) * ws
    
    status = 0
    return
    
  end function get_winds
  
!
! Getter method to calculate windspeed directions from windspeed u/v data for given 
! latitude and longitude, lat and lon.
!
! Data is interpolated linearly from the surrounding 4 points from the GRIB sample data.
!
! Returns the windspeed, status: -1 = fail, 0 = success.
!---------------------------------------------------------------------------------------------------      
  real function get_wind_dir(rlat, rlon, status,nrt) result(wd)
    implicit none
    
    real, intent(in)        :: rlat
    real, intent(in)        :: rlon
    integer, intent(out)    :: status
    integer, intent(in)     :: nrt
    
    real, dimension(4)      ::  f
    real(kind=8)                    ::  x1, x2
    real(kind=8)                    ::  y1, y2
    
    real(kind=8)    :: u, v
    integer :: i1, i2
    integer :: j1, j2
    
    status = -1
    if(rlat >= 90.0 .or. rlat < -90.0) then
      print*, 'ERROR: lat out of range must be 90 to -90:', rlat
      return
    end if
    if (rlon < -180.0 .OR. rlon >= 180.0) then
      print*, 'ERROR: lon out of range must be [-180,180):', rlon
      return
    end if

    status = get_interp_indexes(rlat, rlon, i1, i2, j1, j2,nrt)
    if (status /= 0) then
      print *, "ERROR: Failed to get indexes of surrounding data for interpolation: ", status
      return
    end if

!   -- perform 2D bilinear interpolation (http://en.wikipedia.org/wiki/Bilinear_interpolation)
!   -- u-direction
    f = (/gdas%u_wndspd(i1,j1), gdas%u_wndspd(i2,j1), gdas%u_wndspd(i2,j2), gdas%u_wndspd(i1,j2)/)
    x1 = index2lon(i1, status,nrt)
    x2 = index2lon(i2, status,nrt)
    y1 = index2lat(j1, status,nrt)
    y2 = index2lat(j2, status,nrt)
    
    u = f(1)*(x2-rlon)*(y2-rlat) + f(2)*(rlon-x1)*(y2-rlat) +  &
    &     f(3)*(rlon-x1)*(rlat-y1) + f(4)*(x2-rlon)*(rlat-y1)
    u = 1.0 / ((x2-x1)*(y2-y1)) * u

!   -- v-direction    
    f = (/gdas%v_wndspd(i1,j1), gdas%v_wndspd(i2,j1), gdas%v_wndspd(i2,j2), gdas%v_wndspd(i1,j2)/)
    x1 = index2lon(i1, status,nrt)
    x2 = index2lon(i2, status,nrt)
    y1 = index2lat(j1, status,nrt)
    y2 = index2lat(j2, status,nrt)
    
    v = f(1)*(x2-rlon)*(y2-rlat) + f(2)*(rlon-x1)*(y2-rlat) +  &
    &     f(3)*(rlon-x1)*(rlat-y1) + f(4)*(x2-rlon)*(rlat-y1)
    v = (1 / (x2-x1)*(y2-y1)) * v
    
    wd = (atan2(-1.0*u,-1.0*v) * 180.0/3.1415926535)
    if (wd < 0) wd = wd + 360.0

    status = 0
    return
    
  end function get_wind_dir
  
! Getter method to fetch preciptable water from GRIB data for given latitude and longitude, lat and lon.
! Data is interpolated linearly from the surrounding 4 points from the GRIB sample data.
!
! Returns the precipitable water value, status: -1 = fail, 0 = success.
!---------------------------------------------------------------------------------------------------      
  real(8) function get_pwat(rlat, rlon, status,nrt) result(pw)
    implicit none
    
    real, intent(in)        :: rlat
    real, intent(in)        :: rlon
    integer, intent(out)    :: status
    integer, intent(in)     :: nrt
    
    real, dimension(4)      ::  f
    real                    ::  x1, x2
    real                    ::  y1, y2
    
    integer :: i1, i2
    integer :: j1, j2

    status = -1
    if(rlat >= 90.0 .or. rlat < -90.0) then
      print*, 'ERROR: lat out of range must be 90 to -90:', rlat
      return
    end if
    if (rlon < -180.0 .OR. rlon >= 180.0) then
      print*, 'ERROR: lon out of range must be [-180,180):', rlon
      return
    end if

!   -- perform 2D bilinear interpolation (http://en.wikipedia.org/wiki/Bilinear_interpolation)
    status = get_interp_indexes(rlat, rlon, i1, i2, j1, j2,nrt)
    if (status /= 0) then
      print *, "ERROR: Failed to get indexes of surrounding data for interpolation: ", status
      return
    end if
    
    f = (/gdas%pwat(i1,j1), gdas%pwat(i2,j1), gdas%pwat(i2,j2), gdas%pwat(i1,j2)/)
    x1 = index2lon(i1, status,nrt)
    x2 = index2lon(i2, status,nrt)
    y1 = index2lat(j1, status,nrt)
    y2 = index2lat(j2, status,nrt)
    
    pw = f(1)*(x2-rlon)*(y2-rlat) + f(2)*(rlon-x1)*(y2-rlat) +  &
    &     f(3)*(rlon-x1)*(rlat-y1) + f(4)*(x2-rlon)*(rlat-y1)
    pw = 1.0 / ((x2-x1)*(y2-y1)) * pw
    
    status = 0
    return
    
  end function get_pwat
 
! Getter method to fetch total column ozone from GRIB data for given latitude and longitude, lat and lon.
! Data is interpolated linearly from the surrounding 4 points from the GRIB sample data.
!
! Returns the total colum ozone, status: -1 = fail, 0 = success.
!---------------------------------------------------------------------------------------------------      
  real(8) function get_tozne(rlat, rlon, status, nrt) result(oz)
    implicit none
    
    real, intent(in)        :: rlat
    real, intent(in)        :: rlon
    integer, intent(out)    :: status
    integer, intent(in)     :: nrt
    
    real, dimension(4)      ::  f
    real                    ::  x1, x2
    real                    ::  y1, y2
    
    integer :: i1, i2
    integer :: j1, j2

    status = -1
    if(rlat >= 90.0 .or. rlat < -90.0) then
      print*, 'ERROR: lat out of range must be 90 to -90:', rlat
      return
    end if
    if (rlon < -180.0 .OR. rlon >= 180.0) then
      print*, 'ERROR: lon out of range must be [-180,180):', rlon
      return
    end if

    status = get_interp_indexes(rlat,rlon, i1, i2, j1, j2,nrt)
    if (status /= 0) then
      print *, "ERROR: Failed to get indexes of surrounding data for interpolation: ", status
      return
    end if

!   -- perform 2D bilinear interpolation (http://en.wikipedia.org/wiki/Bilinear_interpolation)
    f = (/gdas%tozne(i1,j1), gdas%tozne(i2,j1), gdas%tozne(i2,j2), gdas%tozne(i1,j2)/)
    x1 = index2lon(i1, status,nrt)
    x2 = index2lon(i2, status,nrt)
    y1 = index2lat(j1, status,nrt)
    y2 = index2lat(j2, status,nrt)
    
    oz = f(1)*(x2-rlon)*(y2-rlat) + f(2)*(rlon-x1)*(y2-rlat) +  &
    &     f(3)*(rlon-x1)*(rlat-y1) + f(4)*(x2-rlon)*(rlat-y1)
    oz = 1.0 / ((x2-x1)*(y2-y1)) * oz
    
    status = 0
    return
    
  end function get_tozne

  
!
! Returns a ncep_gdas_grib object.
!---------------------------------------------------------------------------------------------------    
  type(ncep_gdas_grib) function read_geos5(filename1, status) result(gdas0)
    ! load the wind speed, precip water, and total ozone data from the
    ! file into arrays
    use netcdf
  
    implicit none
  
    character(len=*),intent(in) :: filename1
    integer, intent(inout)        ::  status
    
    character(len=255)    :: dset_name
    
    integer               ::  nc_id
    integer               ::  dim_id
    integer               ::  dset_id
    
    integer               ::  nlat
    integer               ::  nlon
    integer               ::  ntime
    
    real,dimension(:,:,:), allocatable ::  tmp_anc
    
    status = -1
        
    status = nf90_open(filename1, nf90_nowrite, nc_id)
    if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to open geolocation file: ", status
      return
    end if
    
    dset_name = 'lat'
    status = nf90_inq_dimid(nc_id, dset_name, dim_id)
    if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get size of dimension "//trim(dset_name)//": ", status
      return
    end if
    
    status = nf90_inquire_dimension(nc_id, dim_id, len=nlat)
    if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to size of dimension "//trim(dset_name)//": ", status
      return
    end if
    gdas0%nj = nlat
    
    dset_name = 'lon'
    status = nf90_inq_dimid(nc_id, dset_name, dim_id)
    if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get size of dimension "//trim(dset_name)//": ", status
      return
    end if
    
    status = nf90_inquire_dimension(nc_id, dim_id, len=nlon)
    if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to size of dimension "//trim(dset_name)//": ", status
      return
    end if
    gdas0%ni = nlon
    
    dset_name = 'time'
    status = nf90_inq_dimid(nc_id, dset_name, dim_id)
    if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get size of dimension "//trim(dset_name)//": ", status
      return
    end if
    
    status = nf90_inquire_dimension(nc_id, dim_id, len=ntime)
    if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to size of dimension "//trim(dset_name)//": ", status
      return
    end if
    
    allocate(gdas0%u_wndspd(nlon,nlat), gdas0%v_wndspd(nlon,nlat),  &
    &       gdas0%pwat(nlon,nlat), gdas0%tozne(nlon,nlat),gdas0%ps(nlon,nlat), stat=status)
    if (status /= 0) then 
      print *, "ERROR: Failed to allocate GEOS5 data arrays: ", status
      return
    end if
   
    allocate(tmp_anc(nlon,nlat,ntime), stat=status)
    if (status /= 0) then 
      print *, "ERROR: Failed to allocate tmp data arrays: ", status
      return
    end if
   
    dset_name = 'U10M'
    status = nf90_inq_varid(nc_id, dset_name, dset_id)
    if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
    end if
 
    status = nf90_get_var(nc_id, dset_id, tmp_anc)
    if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
    end if            
 
    gdas0%u_wndspd = tmp_anc(:,:,1)

    dset_name = 'V10M'
    status = nf90_inq_varid(nc_id, dset_name, dset_id)
    if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
    end if
 
    status = nf90_get_var(nc_id, dset_id, tmp_anc)
    if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
    end if             

    gdas0%v_wndspd = tmp_anc(:,:,1)
    
    dset_name = 'TQV'
    status = nf90_inq_varid(nc_id, dset_name, dset_id)
    if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
    end if
 
    status = nf90_get_var(nc_id, dset_id, tmp_anc)
    if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
    end if            
 
    gdas0%pwat = tmp_anc(:,:,1)
   
    dset_name = 'TO3'
    status = nf90_inq_varid(nc_id, dset_name, dset_id)
    if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
    end if
 
    status = nf90_get_var(nc_id, dset_id, tmp_anc)
    if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
    end if            
 
    gdas0%tozne = tmp_anc(:,:,1)

    dset_name = 'PS'
    status = nf90_inq_varid(nc_id, dset_name, dset_id)
    if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
      return
    end if
 
    status = nf90_get_var(nc_id, dset_id, tmp_anc)
    if (status /= NF90_NOERR) then
      print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
      return
    end if            
 
    gdas0%ps = tmp_anc(:,:,1)   
    
    deallocate(tmp_anc, stat=status)
    if (status /= 0) then
      print *, "WARNING: Failed to deallocate tmp data array: ", status
    end if
    
    status = nf90_close(nc_id)
    if (status /= NF90_NOERR) then
      print *, 'ERROR: Failed to close GEOS-5 ancillary data file: ', status
      return
    end if
    
!   -- check that all we read all of data that we need.    
    if (.NOT. allocated(gdas0%u_wndspd)) then  
      print *, 'ERROR: No u wind read from file.'
      status = -1
      return
    end if
    if (.NOT. allocated(gdas0%v_wndspd)) then  
      print *, 'ERROR: No v wind read from file.'
      status = -1
      return
    end if
    if (.NOT. allocated(gdas0%pwat)) then  
      print *, 'ERROR: No preciptable water read from file.'
      status = -1
      return
    end if
    if (.NOT. allocated(gdas0%tozne)) then   
      print *, 'ERROR: No ozone read from file.'
      status = -1
      return
    end if

!   -- all data exists, so let's calculate actual total windspeed.
    allocate(gdas0%wndspd(nlon, nlat), stat=status)
    if (status /= 0) then
      print *, 'ERROR: Failed to allocate array for GDAS windspeeds: ', status
      return
    end if
  
    gdas0%wndspd = sqrt(gdas0%v_wndspd**2 + gdas0%u_wndspd**2) 
    
    status = 0
    return
    
  end function read_geos5

!
! Helper function converts latitudes into indexes into the 2D lat/lon data arrays.
! Assumes a 1x1 degree resolution. Index returned is for the lower-left corner of the grid box
! containing the location specified.
!
! Returns the index, status: -1 = fail, 0 = success.
!---------------------------------------------------------------------------------------------------     
  integer function lat2index(lat, status,nrt) result(index)
    implicit none
    
    real, intent(in)      ::  lat
    integer, intent(in)  ::  nrt    
    integer, intent(out)  ::  status
    
    if (lat < -90.0 .OR. lat > 90.0) then
      print *, "ERROR: Invalid latitude: ", lat
      status = -1
      return
    end if
    
    if (nrt == 0) then   
      index = floor((lat + 90.0) / 0.5) + 1
    end if
    if (nrt == 1) then     
      index = floor((lat + 90.0) / 0.25) + 1
    end if
    
    status = 0
    return
    
  end function lat2index
  
!
! Helper function converts indexes into longitudes. Assumes a 1x1 degree resolution.
!
! Returns the latitude, status: -1 = fail, 0 = success.
!---------------------------------------------------------------------------------------------------      
  real function index2lat(index, status, nrt) result(lat)
    implicit none
    
    integer, intent(in)   ::  index, nrt
    integer, intent(out)  ::  status
    
    if (index < 1 .OR. index > gdas%nj) then
      print *, 'ERROR: Index is out of bounds of windspeed array: ', index, gdas%nj
      status = -1
      return
    end if
    
    if (nrt == 0) then      
      lat = -90.0 + 0.5*(index-1) 
    end if
    if (nrt == 1) then 
      lat = -90.0 + 0.25*(index-1) 
    endif

    status = 0
    return
    
  end function index2lat

!
! Helper function converts longitude into indexes into the 2D lat/lon data arrays.
! Assumes a 1x1 degree resolution. Index returned is for the lower-left corner of the grid box
! containing the location specified.
!
! Returns the index, status: -1 = fail, 0 = success.
!---------------------------------------------------------------------------------------------------      
  integer function lon2index(lon, status,nrt) result(index)
    implicit none
    
    real, intent(in)      ::  lon
    integer, intent(in)   ::  nrt    
    integer, intent(out)  ::  status
    
    if (lon < -180.0 .OR. lon >= 180.0) then
      print *, "ERROR: Invalid longitude: ", lon
      status = -1
      return
    end if
    
    if (nrt == 0) then    
      index = floor((lon + 180.0)/0.625) + 1
      if (index > gdas%ni) index = index - gdas%ni    ! loop around if we're close to the dateline
    end if
    if (nrt == 1) then    
      index = floor((lon + 180.0)/0.3125) + 1
      if (index > gdas%ni) index = index - gdas%ni    ! loop around if we're close to the dateline
    end if
    status = 0
    return
    
  end function lon2index

!
! Helper function converts index into longitudes. Assumes a 1x1 degree resolution.
!
! Returns the index, status: -1 = fail, 0 = success.
!---------------------------------------------------------------------------------------------------      
  real function index2lon(index, status,nrt) result(lon)
    implicit none
    
    integer, intent(in)   ::  index,nrt
    integer, intent(out)  ::  status
    
    if (index < 1 .OR. index > gdas%ni) then
      print *, 'ERROR: Index is out of bounds of windspeed array: ', index, gdas%ni
      status = -1
      return
    end if
    
    if (nrt == 0) then
      lon = -180.0 + 0.625*(index-1)
    end if
    if (nrt == 1) then
      lon = -180.0 + 0.3125*(index-1)
    end if
    
    status = 0
    return
    
  end function index2lon

!
! Based on values of lat and lon, returns indexes of surrounding points for 2D linear interpolation.
! i1 and i2 are the indexes to the west and east of the given longitude, lon.
! j1 and j2 are the indexes to the south and north of the given latitude, lat.
! Corner cases involving dateline and north pole (exactly 90.0N) should be handled correctly. :D
!
! Returns status: -1 = fail, 0 = success.
!---------------------------------------------------------------------------------------------------      
  integer function get_interp_indexes(lat, lon, i1, i2, j1, j2,nrt) result(status)
    implicit none
    
    real, intent(in)        ::  lat
    real, intent(in)        ::  lon
    integer, intent(in)   ::  nrt
    integer, intent(out)  ::  i1
    integer, intent(out)  ::  i2
    integer, intent(out)  ::  j1
    integer, intent(out)  ::  j2
    
    status = -1
    
!   -- specifically exclude locations where lat==90.0 exactly. This would break
!   -- the interpolation below.
    if(lat >= 90.0 .or. lat < -90.0) then
      print*, 'ERROR: lat out of range must be [-90,90):', lat
      return
    end if
    if (lon < -180.0 .OR. lon >= 180.0) then
      print*, 'ERROR: lon out of range must be [-180,180):', lon
      return
    end if
    
  !   -- get indexes into the windspeed table for lower left point.
    i1 = lon2index(lon, status,nrt)
    if (status /= 0) then
      print *, "ERROR: Failed to convert longitude to index into ancillary data array: ", status
      return
    end if
    
    j1 = lat2index(lat, status,nrt)
    if (status /= 0) then
      print *, "ERROR: Failed to convert latitude to index into ancillary data array: ", status
      return
    end if
    
!   -- simply calculate the next indexes. Watch for instances where we're 
!   -- crossing the dateline (i2 == gdas%ni+1) and circle back to index 1.
    i2 = i1 + 1
    if (i2 == gdas%ni+1) i2 = 1
    j2 = j1 + 1
  
    status = 0
    return
    
  end function get_interp_indexes
  
 end module viirs_ancillary
