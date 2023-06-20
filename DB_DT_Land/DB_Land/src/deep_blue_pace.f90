!-----------------------------------------------------------------------
! !F90
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:  none
!
! !OUTPUT PARAMETERS:  none
!
! !REVISION HISTORY:
!
!
! !TEAM-UNIQUE HEADER:
!
! !REFERENCES AND CREDITS
!
! !DESIGN NOTES:
!
!
!   Externals:
!
!
! !END
!-----------------------------------------------------------------------
program deep_blue
  use landcover
  use calendars, only: season_from_doy

  use viirs_config, only:                               &
                            viirs_config_type,          &
                            load_viirs_config
                            
  use viirs_db_utils, only:                               &
!                             viirs_db_svm,                 &
                            viirs_gas_correction,         &
                            viirs_gas_correction_fullwv,  &
                            viirs_calib_correction,       &
!                             load_viirs_db_data,           &
!                             load_viirs_db_data_nasa,      &
!                             load_h8_data,                 &
                            calc_minmax_ratio,            &
                            calc_gas_correction,          &
                            calc_gas_correction_fullwv,   &
                            calc_calibration_correction,  &
                            calc_smoke_mask,              &
                            calc_smoke_ae_mask,           &
        									  calc_pyrocb_mask,             &		
								            calc_high_alt_smoke_mask,     &                           
                            ocld_mask
                            
  use viirs_ocean_aero, only:                           &
                            load_viirs_ocean_aerosol_luts,  &
                            unload_viirs_ocean_aerosol_luts, &
                            ocean_cloud_filter, &
                            ocean_cloud_filter_ahi, &                            
                            load_bathymetry_lut, &
                            load_chl_lut

  use modis_surface, only:                              &
                            load_terrainflg_tables,     &
                            load_seasonal_desert,       &
                            load_brdf,                  &
                            unload_brdf,                &
                            set_limits,                 &
                            set_limits6,                &
                            load_hdfLER,                &
                            get_LER412,                 &
                            get_LER470,                 &
                            get_LER650,                 &
                            latlon_to_index_ler,        &
                            get_geographic_zone,        &
                            get_sfc_elev_std,           &
                            load_swir_coeffs
  
  use viirs_aerosol_luts, only:                         &
                            load_viirs_aerosol_luts,    &
                            unload_viirs_aerosol_luts
  
  use viirs_ler_luts
                            
  use db_debug
  use read_l1b
  use viirs_ancillary, only:                            &
                            load_geos5_data,            &
                            get_pwat,                   &
                            get_tozne,                  &
                            get_winds,                  &
                            get_wind_dir,               &
                            get_ps
 
  use seawifs_surface_pressure, only :                  & 
                            load_surfterr_table,        & 
                            get_elevation,              &
                            get_surface_pressure
 
  use screening
                            
  implicit none

  include 'newaottbl90.inc'
  include 'sfc21tbl90.inc'
  
  integer, parameter  :: R_DBL 	= kind(1.0d0)


  type (viirs_db_svm)             ::  viirs_data
  type (viirs_config_type)        ::  config

  character(len=255)              ::  config_file
  
  integer                         ::  args
  character(len=255)              ::  usage
  character(len=5)                ::  platform
  type(viirs_calib_correction)          ::  vcc
 
  integer, dimension(:,:), allocatable  ::  lcld_mask, tmp_mask!,ocld_mask
  integer, dimension(:,:), allocatable  ::  lskip_mask
  integer, dimension(:,:), allocatable  ::  oskip_mask
  integer, dimension(:,:), allocatable  ::  smoke_mask
  integer, dimension(:,:), allocatable  ::  smoke_ae_mask
	integer, dimension(:,:), allocatable  ::  pyrocb_mask		
  integer, dimension(:,:), allocatable  ::  high_alt_smoke_mask	
  integer, dimension(:,:), allocatable  ::  snow_mask, snow_mask2, month
  real, dimension(:,:), allocatable     ::  sr650
  real, dimension(:,:), allocatable     ::  minmaxbt
  real, dimension(:,:), allocatable     ::  sfcstd
  integer, dimension(:,:), allocatable  ::  gzflg
  real, dimension(:,:), allocatable     ::  wv
  real, dimension(:,:), allocatable     ::  oz
  real, dimension(:,:), allocatable     ::  windsp
  real, dimension(:,:), allocatable     ::  winddir
  integer, dimension(:,:), allocatable  ::  lc                  ! land cover
  integer, dimension(:,:), allocatable  ::  bathy ! bathymetry
  real, dimension(:,:), allocatable     ::  chl ! log-10 Chl climatology

  real, dimension(:,:), allocatable     ::  ler412
  real, dimension(:,:), allocatable     ::  ler488
  real, dimension(:,:), allocatable     ::  ler670
  real, dimension(:,:), allocatable     ::  qdf412
  real, dimension(:,:), allocatable     ::  qdf488
  real, dimension(:,:), allocatable     ::  qdf670, test_arr
  real                                ::  rat488670
  real                                ::  rat412488
 
  type(viirs_gas_correction)         ::  gas_corr
  type(viirs_gas_correction_fullwv)         ::  gas_corr_full
    
  real, dimension(:,:), allocatable     ::  m01_nogc,m03_nogc,m05_nogc,m11_oldgc,m09_noiof,m07_noiof

  integer                               ::  ilat, ilon, nrt

  integer, dimension(:,:), allocatable  ::  n_total_pixels
  integer                               ::  cell_resolution,min_flag
  integer, dimension(2)                 ::  cell_dims
  
  real                                  ::  ndsi
  real                                  ::  ndvi_lower, ndvi_upper
  real                                  ::  dd
  integer                         ::  i,j 
  integer                         ::  i1, i2, j1, j2
  integer                         ::  status
  integer                         ::  season
  integer, dimension(2)           ::  dims2
  real, dimension(:,:), allocatable    ::  to_iof
  real, parameter                 ::  d2r = 3.14159/180.0   ! convert degrees to radians

  type(output_file)   ::  debug_file
  type(dataset), dimension(:), allocatable  :: all_output_datasets
  logical ::  hires
  
  hires  = .false.      ! set hires .true. to use different aggregation method
  status = -1
  
! -- when using gfortran        
  usage = "./deep_blue <CONFIG_FILE>"
  args = command_argument_count()
  if ( args < 1 ) then
   write(*,*) usage
   stop
  end if

  call get_command_argument(1, value=config_file, status=status)
  if (status /= 0) then
    print *, "ERROR: Failed to get config file from command arguments: ", status
    stop
  end if

  config = load_viirs_config(config_file, status)
  if (status /= 0) then
    print *, "ERROR: Failed to read in VIIRS configuration file: ", status
    stop
  end if
  platform  =  config%platform    !'VIIRS' or 'AHI' or 'GOES' 
  print *, 'Start Deep_blue.f90 :', platform

!-----------------------------------------------------------------------------------------
! -- debug output, part 1
!-----------------------------------------------------------------------------------------
!-- setup our storage for the debug data output arrays.
 allocate(all_output_datasets(5), stat=status)
 if (status /= 0) then 
  	print *, "ERROR: Failed to allocate output dataset array: ", status
  	stop
 end if
 
 do i = 1, size(all_output_datasets,1)  
   allocate(all_output_datasets(i)%values(viirs_data%xscan, viirs_data%scan), stat=status)
   if (status /= 0) then
     print *, 'ERROR: Failed to allocate all_output_datasets(i)%values: ', i, status
     stop
   end if   
 end do
  
!-----------------------------------------------------------------------------------------
! -- read all input data (reflectances, geolocation, geometry, land/sea, etc...
!-----------------------------------------------------------------------------------------

! For pace
   viirs_data = load_pace_data(config%geo_m, config%l1b_m,config%year, config%month, config%day, status)
   if (status /= 0) then
     print *, "ERROR: Failed to load NASA L1B VIIRS data: ", status
     stop
   endif
! For Himawari
  if (platform .eq. 'AHI') then 
   viirs_data = load_h8_data(config%geo_m, status)
   if (status /= 0) then
     print *, "ERROR: Failed to load Himawari8 data: ", status
     stop
   endif
  endif
! For GOES
  if (platform .eq. 'GOES') then 
!    viirs_data = load_goes_data(config%geo_m, config%segment, status)
   viirs_data = load_goes_data(config%geo_m,config%goes_1,config%goes_2,config%goes_3,config%goes_4,&
      &   config%goes_5,config%goes_6,config%goes_7,config%goes_11,config%goes_13,config%goes_14,&
      &   config%goes_15,config%segment, status)         !load_goes_data(config%segment, status) 
   if (status /= 0) then
     print *, "ERROR: Failed to load Himawari8 data: ", status
     stop
   endif
  endif

! -- apply calibration correction factors to TOA reflectance

  vcc = calc_calibration_correction(viirs_data%yr, viirs_data%mo, viirs_data%dy, status)
  if (status /= 0) then
    print *, 'ERROR: Failed to calculate VIIRS calibration correction factors: ', status
    stop
  end if

! -- load auxiliary input files
! -- aerosol inversion LUTs
  print *, 'Loading aerosol LUTs...'
  status = load_viirs_aerosol_luts(config%aerosol_land_file, config%aerosol_dust_file)
  if (status /= 0) then
    print *, "ERROR: Failed to load VIIRS aerosol LUTS: ", status
    stop
  end if
  print *, 'done.'

! -- load ocean aerosol inversion LUTs
  print *, 'Loading ocean LUTs...'
  print *, trim(config%aerosol_ocean_dust_file)
  print *, trim(config%aerosol_ocean_fine_file)
  print *, trim(config%aerosol_ocean_mari_file)
  print *, trim(config%aerosol_ocean_mix_file)
  status = load_viirs_ocean_aerosol_luts(config%aerosol_ocean_dust_file,  &
  & config%aerosol_ocean_fine_file, config%aerosol_ocean_mari_file, &
  & config%aerosol_ocean_mix_file,platform)
  if (status /= 0) then
    print *, "ERROR: Failed to load VIIRS ocean aerosol LUT: ", status
    stop
  end if
  
! -- load surface elevation and pressure table
  if (platform .eq. 'AHI' .or. platform .eq. 'GOES') then   
  status = load_surfterr_table(config%surfpressure_file)
  if (status /= 0) then
    print *, "ERROR: Failed to load surface pressure and elevation input file: ", status
    stop
  end if
  end if
  
! -- load bathymetry data base.
! A Sayer 01 Aug 2017
  print *, 'Loading bathymetry...'
  allocate(bathy(viirs_data%xscan,viirs_data%scan), stat=status)
  status = load_bathymetry_lut(config%bathymetry_lut_file, &
  & viirs_data%lat,viirs_data%lon,viirs_data%xscan,viirs_data%scan,bathy)
  if (status /= 0) then
    print *, "ERROR: Failed to load and use bathymetry LUT: ", status
    stop
  end if

! -- load Chl climatology base.
! A Sayer 25 Sep 2017
  print *, 'Loading climatological Chl...'
  allocate(chl(viirs_data%xscan,viirs_data%scan), stat=status)
  status = load_chl_lut(config%chl_lut_file, &
  & viirs_data%lat,viirs_data%lon,viirs_data%xscan, &
  & viirs_data%scan,config%month,chl)
  if (status /= 0) then
    print *, "ERROR: Failed to load and use Chl LUT: ", status
    stop
  end if

! -- LER luts
  call load_viirs_ler_luts(config%ler_lut_file, status)
  if (status /= 0) then
    print *, "ERROR: Failed to load VIIRS LER LUTS: ", status
    stop
  end if

! -- land cover
  status = load_landcover(trim(config%landcover_file))
  if (status /= 0) then
    print *, "ERROR: Unable to load land cover input file: ", status
    stop
  end if
        
! -- geo zones
  season = season_from_doy(2005,viirs_data%doy)
  print *, "season: ", season
  status = load_terrainflg_tables(config%geozone_file,season)
  if (status /= 0) then
    print *, "ERROR: Unable to load geozone table: ", status
    stop
  end if
     
! -- seasonal desert data
  status = load_seasonal_desert(config%seasonal_deserts_file)
  if (status /= 0) then
    print *, "ERROR: Unable to load seasonal deserts file: ", status
  end if

! -- 670nm BRDF data
  status = load_brdf(config%brdfbase_file)
  if (status /= 0) then
    print *, "ERROR: Unable to load BRDF base input file: ", status
    stop
  end if

! -- surface database and BRDF coefficients
  dims2 = (/viirs_data%xscan, viirs_data%scan/)
  status = set_limits(dims2, viirs_data%lat, viirs_data%lon)
  if (status /= 0) then
    print *, "ERROR: Failure to set limits on surface coefficients table: ", status
    stop
  end if

  status = load_hdfLER(config%modis_surfdb_file, config%viirs_surfdb_file, config%surfcoeffs_file)
  if (status /= 0) then
     print *, "ERROR: Unable to load surface BRDF coefficients: ", status
     stop
  end if

! -- swir vs. vis surface coeffs
  dims2 = (/viirs_data%xscan, viirs_data%scan/)
  status = set_limits6(dims2, viirs_data%lat, viirs_data%lon)
  if (status /= 0) then
    print *, "ERROR: Failure to set limits on 2.2 um surface coefficients table: ", status
    stop
  end if

  status = load_swir_coeffs(config%swir_vis_surfcoeffs_file)
  if (status /= 0) then
    print *, "ERROR: Unable to load swir vs. vis surface coeffs file: ", status
  end if

! -- vegetated retrieval landcover 
  call get_lut_igbp_land_cover(config%veg_landcover_file, status)
  if (status /= 0) then
    print *, "ERROR: Failed to read in landcover input for vegetated retrieval: ", status
    stop
  end if
  
  call get_lut_211sfc(config%veg_sfc21_file, status)
  if (status /= 0) then
    print *, "ERROR: Failed to read in landcover input for vegetated 2.1um sfc table: ", status
    stop
  end if
  
  allocate(wv(viirs_data%xscan,viirs_data%scan), oz(viirs_data%xscan,viirs_data%scan),&
  &        windsp(viirs_data%xscan,viirs_data%scan),&
  &        winddir(viirs_data%xscan,viirs_data%scan), lc(viirs_data%xscan,viirs_data%scan), stat=status)
  if (status /= 0) then
    print *, "ERROR: Failed to allocate water vapor, ozone arrays, and wind speed arrays: ", status
    stop
  end if
  
  lc(:,:) =-999

  print *, 'ancillary files: '
  print *, 'pre: ', trim(config%gdas_file1)
  print *, 'post: ', trim(config%gdas_file2)
  print *, 'hr: ', viirs_data%hr
  print *, 'min: ', viirs_data%min
  
  status = load_geos5_data(config%gdas_file1, config%gdas_file2, viirs_data%hr, viirs_data%min)
  
  if (status /= 0) then
    print *, "ERROR: Failed to load GEOS5 data: ", status
    stop
  end if
  print *, 'done loading ancillary data.'

! -- check NRT retrieval
    nrt = 0
    if (config%gdas_file1(20:23) == 'fp.f') then
      print *, 'start NRT process'
      nrt = 1
    endif
    if (config%gdas_file1(20:23) /= 'fp.f' .and. config%gdas_file1(20:23) /= 'fpit') then
      print *, 'ERROR: check GEOS ancillary data'
      stop
    endif
    
! -- interpolating ancillary data
   if (platform .eq. 'AHI'  .or. platform .eq. 'GOES') then 
     viirs_data%land_mask(:,:) =-999
   end if
   do j = 1, viirs_data%scan
    do i = 1, viirs_data%xscan
	
      if (viirs_data%lat(i,j) < -900.0 .OR. viirs_data%lon(i,j) < -900.0) cycle
!       if (viirs_data%m03_refl(i,j) < -900.0) cycle

!     -- 10 converts to cm.
      wv(i,j) = get_pwat(viirs_data%lat(i,j),(viirs_data%lon(i,j)),status,nrt) / 10.0
      if (status /= 0) then
        print *, "ERROR: Failed to fill water vapor array: ", status
	      print *, "lat", viirs_data%lat(i,j),"lon", viirs_data%lon(i,j)
	  	  print *, "i:", i, "j:", j
	  	  stop
      end if

!     -- 0.001 converts DU to atm/cm.
	    oz(i,j) = get_tozne(viirs_data%lat(i,j),(viirs_data%lon(i,j)),status,nrt) * 0.001
	    if (status /= 0) then
        print *, "ERROR: Failed to fill ozone array: ", status
	  	  print *, "i:", i, "j:", j
	  	  stop
      end if
     
      windsp(i,j) = get_winds(viirs_data%lat(i,j),(viirs_data%lon(i,j)),status,nrt)
      if (status /= 0) then
        print *, "ERROR: Failed to fill wind array: ", status
        print *, "i:", i, "j:", j
        stop
      end if

      winddir(i,j) = get_wind_dir(viirs_data%lat(i,j),(viirs_data%lon(i,j)),status,nrt)
      if (status /= 0) then
        print *, "ERROR: Failed to fill wind direction array: ", status
        print *, "i:", i, "j:", j
        stop
      end if
      
      lc(i,j) = get_landcover(viirs_data%lat(i,j),viirs_data%lon(i,j), status)
      if (status /= 0) then
        print *, "ERROR: Failed to get land cover value: ", status
        print *, "i:", i, "j:", j
        stop
      end if      

      ! For Himawari & GOES16 (VIIRS has land type information in L1b file
      
      if (platform .eq. 'AHI'  .or. platform .eq. 'GOES') then 
        if (lc(i,j)==0) then 
          viirs_data%land_mask(i,j)  = 3  !ocean
        else
          viirs_data%land_mask(i,j)  = 1  !land
          if (lc(i,j)==6) viirs_data%land_mask(i,j)  = 0  !dessert
        end if
      
        viirs_data%elev(i,j) = get_elevation(viirs_data%lat(i,j),viirs_data%lon(i,j), status)
        if (status /= 0) then
          print *, "ERROR: Failed to get surface elevation: ", status
          print *, "i: ", i, "j: ", j
          stop
        end if  
        
        viirs_data%ps(i,j) = get_ps(viirs_data%lat(i,j),(viirs_data%lon(i,j)),status,nrt)/100.
        if (status /= 0) then
          print *, "ERROR: Failed to fill GEOS5 PS array: ", status
          print *, "i:", i, "j:", j
          stop
        end if        
      end if      

    end do
  end do
  print *, 'done interpolating ancillary data.'

!-----------------------------------------------------------------------------------------
! -- calculate and store valid pixels count, assume 8x8 pixel aggregation
!-----------------------------------------------------------------------------------------
  cell_resolution = 8      ! number of pixels to aggregate
  if (platform .eq. 'AHI'  .or. platform .eq. 'GOES') cell_resolution = 4 ! For Himawari
  if (hires) then
   cell_resolution = 2 
  endif
     
  cell_dims = (/ceiling(viirs_data%xscan/float(cell_resolution)), ceiling(viirs_data%scan/float(cell_resolution))/)
  allocate(n_total_pixels(cell_dims(1), cell_dims(2)), stat=status)
  if (status /= 0) then
    print *, "ERROR: Failed to allocate array for valid pixel counts: ", status
    stop
  end if
  n_total_pixels(:,:) = 0
  
  print *, 'start n_total_pixels calc...'
  
  do j = 1, cell_dims(2)
    do i = 1, cell_dims(1)
      i1 = (cell_resolution * (i-1)) + 1
      i2 = min((cell_resolution * i), viirs_data%xscan)
      j1 = (cell_resolution * (j-1)) + 1
      j2 = min((cell_resolution * j), viirs_data%scan)
      n_total_pixels(i,j) = count(viirs_data%m01_refl(i1:i2, j1:j2) > -900.0 &
      & .AND. viirs_data%land_mask(i1:i2, j1:j2) > -900)
      if (platform .eq. 'AHI'  .or. platform .eq. 'GOES') n_total_pixels(i,j) = count(viirs_data%m03_refl(i1:i2, j1:j2) > -900.0 &
      & .AND. viirs_data%land_mask(i1:i2, j1:j2) > -900)
    end do
  end do
  print *, 'done.'
  
!-----------------------------------------------------------------------------------------
! -- cloud filtering
!-----------------------------------------------------------------------------------------  
   print *, 'start cloud filtering...'
    status= cloud_screen(viirs_data,config, wv,platform,lcld_mask, ocld_mask,&
    &  month,gzflg, snow_mask, snow_mask2,sfcstd, sr650, minmaxbt)

   if (status /= 0) then 
     print *, "ERROR: cloud screening: ", status
     stop
   end if  
  
!-----------------------------------------------------------------------------------------
! -- convert reflectances to I/F
!-----------------------------------------------------------------------------------------
  print *, 'start I/F conversion...'
  allocate(to_iof(viirs_data%xscan,viirs_data%scan),m09_noiof(viirs_data%xscan,viirs_data%scan),&
          &       m07_noiof(viirs_data%xscan,viirs_data%scan), stat=status)
  if (status /= 0) then
    print *, "ERROR: Unable to allocate I/F conversion array: ", status
    stop
  end if
  
  if (platform .ne. 'GOES') to_iof  = (cos(viirs_data%sza*d2r)/3.14159)
  !reflectance unit of GOES16 is I*pi/F0, other satellites are I*pi/(F*cos(sza))
  if (platform .eq. 'GOES') to_iof  = (1./3.14159)      
  
!   if (platform .eq. 'VIIRS') then  
    where (viirs_data%m01_refl > -900.0) viirs_data%m01_refl = viirs_data%m01_refl*to_iof
!     where (viirs_data%m02_refl > -900.0) viirs_data%m02_refl = viirs_data%m02_refl*to_iof
    where (viirs_data%m08_refl > -900.0) viirs_data%m08_refl = viirs_data%m08_refl*to_iof
!   end if
!   if (platform .ne. 'GOES') then
    where (viirs_data%m04_refl > -900.0) viirs_data%m04_refl = viirs_data%m04_refl*to_iof
!   end if  
!   if (platform .ne. 'AHI') then
    m09_noiof(:,:) = viirs_data%m09_refl
    where (viirs_data%m09_refl > -900.0) viirs_data%m09_refl = viirs_data%m09_refl*to_iof 
!   end if 
  where (viirs_data%m03_refl > -900.0) viirs_data%m03_refl = viirs_data%m03_refl*to_iof
  where (viirs_data%m05_refl > -900.0) viirs_data%m05_refl = viirs_data%m05_refl*to_iof
  m07_noiof(:,:) = viirs_data%m07_refl
  where (viirs_data%m07_refl > -900.0) viirs_data%m07_refl = viirs_data%m07_refl*to_iof
  where (viirs_data%m10_refl > -900.0) viirs_data%m10_refl = viirs_data%m10_refl*to_iof
  where (viirs_data%m11_refl > -900.0) viirs_data%m11_refl = viirs_data%m11_refl*to_iof
  
  deallocate(to_iof, stat=status)
  if (status /= 0) then
    print *, 'ERROR: Unable to deallocate I/F conversion array: ', status
    stop
  end if
  
!-----------------------------------------------------------------------------------------
! -- trace gas correction
!-----------------------------------------------------------------------------------------
  print *, 'start trace gas correction...'
  allocate(tmp_mask(viirs_data%xscan,viirs_data%scan), stat=status)
  if (status /= 0) then
    print *, "ERROR: Failed to allocate valid pixel mask for trace gas correction calculation: ", status
    stop
  end if
  tmp_mask(:,:) = 1
  if (platform .eq. 'VIIRS') where (viirs_data%m01_refl < -900) tmp_mask = 0
  if (platform .eq. 'AHI'  .or. platform .eq. 'GOES') where (viirs_data%m03_refl < -900) tmp_mask = 0

!-----------------------------------------------------------------------------------------
! -- no gas corrections
! -----------------------------------------------------------------------------------------
  allocate(m01_nogc(viirs_data%xscan,viirs_data%scan),  &
           m03_nogc(viirs_data%xscan,viirs_data%scan), &
           m05_nogc(viirs_data%xscan,viirs_data%scan), stat=status)

  if (status /= 0) then
    print *, "ERROR: Failed to allocate arrays for new gas correction: ", status
    stop
  end if
  m01_nogc(:,:) = viirs_data%m01_refl
  m03_nogc(:,:) = viirs_data%m03_refl
  m05_nogc(:,:) = viirs_data%m05_refl

!-----------------------------------------------------------------------------------------
! -- new gas corrections (full wv, called '_oldgc')
! -----------------------------------------------------------------------------------------
  print *, 'start full wv correction...'
  allocate(m11_oldgc(viirs_data%xscan,viirs_data%scan), stat=status)

  if (status /= 0) then
    print *, "ERROR: Failed to allocate arrays for old gas correction: ", status
    stop
  end if

  m11_oldgc(:,:) = -999.0

  ! only for m11, function should be modified to include other bands
  gas_corr_full = calc_gas_correction_fullwv(viirs_data%sza, viirs_data%vza, oz, &
  &                              wv, viirs_data%ps/1013.25, status, mask=tmp_mask)
  if (status /= 0) then
    print *, "ERROR: Failed to calculate trace gas corrections: ", status
    stop
  end if
  print *, 'corrections calculated, applying...'
  
!   if (platform .eq. 'VIIRS') then
    where (viirs_data%m11_refl > -900.0) m11_oldgc = viirs_data%m11_refl*gas_corr_full%m11
!   else
!     m11_oldgc = viirs_data%m11_refl
!   end if

  deallocate(gas_corr_full%m11, stat=status)
  if (status /= 0) then
    print *, "ERROR: Unable to deallocate gaseous transmittance arrays: ", status
    stop
  end if
  print *, 'finish full wv correction...' 
  
!-----------------------------------------------------------------------------------------
! -- new gas corrections (half wv)
! -----------------------------------------------------------------------------------------

  print *, 'start half wv correction...'
  gas_corr = calc_gas_correction(viirs_data%sza, viirs_data%vza, oz, &
  &                   wv/2.0, platform, viirs_data%ps/1013.25, status, mask=tmp_mask)
  if (status /= 0) then 
    print *, "ERROR: Failed to calculate trace gas corrections: ", status
    stop
  end if
  print *, 'corrections calculated, applying...'
  
  if (platform .eq. 'VIIRS') then
   where (viirs_data%m01_refl > -900.0) viirs_data%m01_refl = viirs_data%m01_refl*gas_corr%m01
!    where (viirs_data%m02_refl > -900.0) viirs_data%m02_refl = viirs_data%m02_refl*gas_corr%m02
   where (viirs_data%m08_refl > -900.0) viirs_data%m08_refl = viirs_data%m08_refl*gas_corr%m08
   where (viirs_data%m04_refl > -900.0) viirs_data%m04_refl = viirs_data%m04_refl*gas_corr%m04
   where (viirs_data%m03_refl > -900.0) viirs_data%m03_refl = viirs_data%m03_refl*gas_corr%m03
   where (viirs_data%m05_refl > -900.0) viirs_data%m05_refl = viirs_data%m05_refl*gas_corr%m05
   where (viirs_data%m07_refl > -900.0) viirs_data%m07_refl = viirs_data%m07_refl*gas_corr%m07
   where (viirs_data%m10_refl > -900.0) viirs_data%m10_refl = viirs_data%m10_refl*gas_corr%m10
   where (viirs_data%m11_refl > -900.0) viirs_data%m11_refl = viirs_data%m11_refl*gas_corr%m11   
  else
    print*, 'non VIIRS retrieaval, no gas correction'
  end if
 
  deallocate(gas_corr%m01, gas_corr%m02, gas_corr%m03, gas_corr%m04, gas_corr%m05, &
  &         gas_corr%m06, gas_corr%m07, gas_corr%m08, gas_corr%m10,   &
  &         gas_corr%m11, stat=status)
  if (status /= 0) then
    print *, "ERROR: Unable to deallocate gaseous transmittance arrays: ", status
    stop
  end if
 
  deallocate(tmp_mask, stat=status)
  if (status /= 0) then
    print *, "ERROR: Unable to deallocate valid pixel mask for gas transmittance calculation: ", status
    stop
  end if
  print *, 'finish half wv correction...'
  print *, 'finish trace gas correction...'

  if (platform .eq. 'VIIRS') then
    where (viirs_data%m01_refl > -900.0)
       viirs_data%ndvi = (viirs_data%m07_refl - viirs_data%m05_refl) /  &
      &         (viirs_data%m07_refl + viirs_data%m05_refl)
    end where 
  end if 
  
  if (platform .eq. 'AHI'  .or. platform .eq. 'GOES') then
    where (viirs_data%m03_refl > -900.0)
       viirs_data%ndvi = (viirs_data%m07_refl - viirs_data%m05_refl) /  &
      &         (viirs_data%m07_refl + viirs_data%m05_refl)
    end where 
  end if 
!-----------------------------------------------------------------------------------------
! -- calculate LER
!-----------------------------------------------------------------------------------------
  print *, 'start LER calculation...'
  allocate(ler412(viirs_data%xscan,viirs_data%scan), ler488(viirs_data%xscan,viirs_data%scan),  &
  &         ler670(viirs_data%xscan,viirs_data%scan), stat=status)
  if (status /= 0) then
    print *, "ERROR: Failed to allocate LER arrays: ", status
    stop
  end if
  ler412(:,:) = -999.0
  ler488(:,:) = -999.0
  ler670(:,:) = -999.0
  
  allocate(qdf412(viirs_data%xscan,viirs_data%scan), qdf488(viirs_data%xscan,viirs_data%scan),  &
  &         qdf670(viirs_data%xscan,viirs_data%scan), stat=status)
  if (status /= 0) then
    print *, "ERROR: Failed to allocate LER arrays: ", status
    stop
  end if
  qdf412(:,:) = -999.0
  qdf488(:,:) = -999.0
  qdf670(:,:) = -999.0
  
  do j = 1, viirs_data%scan
    do i = 1, viirs_data%xscan  
    
      xlat = viirs_data%lat(i,j)
      xlong = viirs_data%lon(i,j)
      sza = viirs_data%sza(i,j)
      xthet = viirs_data%vza(i,j)
      xphi = viirs_data%raa(i,j)
      
      !     -- move read_data operations here.
      if (platform .eq. 'VIIRS') then
        if (m05_nogc(i,j) <= 0.0 .OR. m03_nogc(i,j) <= 0.0 .OR. m01_nogc(i,j) <= 0.0) cycle
      end if
      if (platform .eq. 'AHI'  .or. platform .eq. 'GOES') then 
        if (m05_nogc(i,j) < 0.0 .OR. m03_nogc(i,j) < 0.0) cycle      
      end if
      if (xlat < -900 .OR. xlong < -900) cycle
      
      xnvalm(1) = -1.0*alog10(m05_nogc(i,j))
      xnvalm(2) = -1.0*alog10(m05_nogc(i,j))
      xnvalm(3) = -1.0*alog10(m03_nogc(i,j))
      xnvalm(4) = -1.0*alog10(m03_nogc(i,j))
      xnvalm(6) = -1.0*alog10(m05_nogc(i,j))
      if (platform .eq. 'VIIRS') xnvalm(5) = -1.0*alog10(m01_nogc(i,j))
      
      isnow = 0
      pcloud = 0.7
      
      ilat = floor(xlat*10.0) + 900 + 1
      if (ilat > 1800) ilat = 1800
      if (ilat < 1)    ilat = 1
      
      ilon = floor(xlong*10.0) + 1800 + 1
      if (ilon > 3600) ilon = 3600
      if (ilon < 1)    ilon = 1
      
      sfref412 = get_LER412(ilat, ilon, viirs_data%ndvi(i,j), viirs_data%sca(i,j), &
                            viirs_data%raa(i,j),min_flag)/100.0
      sfref470 = get_LER470(ilat, ilon, viirs_data%ndvi(i,j), viirs_data%sca(i,j), &
                            viirs_data%raa(i,j),min_flag)/100.0
      sfref650 = get_LER650(ilat, ilon, viirs_data%ndvi(i,j), viirs_data%sca(i,j), &
                            viirs_data%raa(i,j),min_flag)/100.0
      pteran   = viirs_data%ps(i,j)/1013.25 

      call total
            
      ler412(i,j) = realbuf(22)
      ler488(i,j) = realbuf(23)
      ler670(i,j) = realbuf(6)
      
      qdf412(i,j) = realbuf(18)
      qdf488(i,j) = realbuf(19)
      qdf670(i,j) = realbuf(20)
      
    end do
  end do
  if (platform .eq. 'AHI'  .or. platform .eq. 'GOES') ler412(:,:) = -999.0
  print *, 'end LER calculation'
    
!-----------------------------------------------------------------------------------------
! -- cloud filtering2
!-----------------------------------------------------------------------------------------  
    print *, 'start cloud filtering2...'    
    
    status= cloud_screen2(viirs_data,platform,ler412,ler488,ler670, sr650, sfcstd, snow_mask, &
    & snow_mask2,ocld_mask, gzflg,lcld_mask,month,lc,lskip_mask, oskip_mask,smoke_mask, &
    & smoke_ae_mask,pyrocb_mask	,high_alt_smoke_mask,m09_noiof,m07_noiof)  
   print *, 'cloud filtering2 complete.'		

   if (status /= 0) then 
     print *, "ERROR: cloud screening: ", status
     stop
   end if 
   
   allocate(test_arr(viirs_data%xscan,viirs_data%scan), stat=status)
   test_arr= viirs_data%land_mask

! modified by CH 2/18/2020

  do j = 1, viirs_data%scan
  do i = 1, viirs_data%xscan

   if (gzflg(i,j) == 12 .and. ler412(i,j) > 12.0 .and. ler412(i,j)< 50.0) then
   dd = (ler412(i,j)*ler670(i,j)) / (ler488(i,j)*ler488(i,j))
    if (dd < 0.9 .and. minmaxbt(i,j) < 1.04) then
        lcld_mask(i,j) = 0
        lskip_mask(i,j) = 0
        smoke_mask(i,j) = 1
        smoke_ae_mask(i,j) = 1
    endif
   endif

  end do
  end do

!-----------------------------------------------------------------------------------------
! -- debug output, part 2start half wv correction
!-----------------------------------------------------------------------------------------
 !  all_output_datasets(1)  = dataset('lat', viirs_data%lat)
!   all_output_datasets(2)  = dataset('lon', viirs_data%lon)
!   all_output_datasets(3) = dataset('ler488', ler488)
!   all_output_datasets(4) = dataset('ler670', ler670) 
!   all_output_datasets(5)  = dataset('dstar', viirs_data%dstar) 
! !   all_output_datasets(1)  = dataset('ocld_mask', ocld_mask) 
! !   all_output_datasets(3)  = dataset("m03",viirs_data%m03_refl)
! !   all_output_datasets(4)  = dataset("sza",0.11*cos(viirs_data%sza*(3.1415926535 / 180.0))) 
! !     all_output_datasets(3)  = dataset('land_mask', test_arr)
! !     all_output_datasets(4)  = dataset('ocld_mask', ocld_mask)
! !   all_output_datasets(3)  = dataset('land_mask', test_arr)
! ! !   all_output_datasets(3)  = dataset("ndvi", viirs_data%ndvi)
! !   all_output_datasets(3)  = dataset("solar_zenith_angle",viirs_data%sza)    
! !   all_output_datasets(4) = dataset("viewing_zenith_angle",viirs_data%vza)  
! !   all_output_datasets(5)  = dataset("raa",viirs_data%raa)
! !   all_output_datasets(6)  = dataset("vaa",viirs_data%vaa)
! !   all_output_datasets(7)  = dataset("saa",viirs_data%saa)
! ! ! !  
! !   all_output_datasets(8)  = dataset("m03",viirs_data%m03_refl)
! !   all_output_datasets(9)  = dataset("m05",viirs_data%m05_refl)
! !   all_output_datasets(10)  = dataset("m07",viirs_data%m07_refl)
! !   all_output_datasets(11)  = dataset("m10",viirs_data%m10_refl)
! !   all_output_datasets(12)  = dataset("m11",viirs_data%m11_refl)
! !   all_output_datasets(13)  = dataset('land_mask', test_arr)
! !     all_output_datasets(14)  = dataset('ocld_mask', ocld_mask)
! 
! ! 
! ! !   all_output_datasets(3) = dataset('ref2.2', viirs_data%m11_refl)
! ! !   all_output_datasets(4) = dataset('ler670', ler670)
! ! !   all_output_datasets(5) = dataset('lat', viirs_data%lat)
! ! !   all_output_datasets(6) = dataset('lon', viirs_data%lon)
! ! !   all_output_datasets(1)  = dataset("lskip_mask", lskip_mask)
! ! !   all_output_datasets(2)  = dataset("oskip_mask",oskip_mask)
! ! !   all_output_datasets(3)  = dataset("lcld_mask",lcld_mask)
! ! !   all_output_datasets(4)  = dataset("ndvi", viirs_data%ndvi)
! ! !   all_output_datasets(3)  = dataset("ocld_mask",ocld_mask)
! ! !   all_output_datasets(6)  = dataset('ler412', ler412)
! ! !   all_output_datasets(7) = dataset('ler488', ler488)
! ! !   all_output_datasets(8) = dataset('ler670', ler670)  
! ! !   all_output_datasets(9)  = dataset("solar_zenith_angle",viirs_data%sza)    
! ! !   all_output_datasets(10)  = dataset("viewing_zenith_angle",viirs_data%vza)    
! ! 
! ! !   all_output_datasets(5)  = dataset("ps", viirs_data%ps)
! ! !   all_output_datasets(6)  = dataset("geos ps", ps)
! ! !   all_output_datasets(7)  = dataset("ndvi", viirs_data%ndvi)
! ! !   all_output_datasets(8)  = dataset("scattering_angle",viirs_data%sca)
! ! !   all_output_datasets(9)  = dataset('ler412', ler412)
! ! !   all_output_datasets(10) = dataset('ler488', ler488)
! ! !   all_output_datasets(11) = dataset('ler670', ler670)
! ! 
!   print *, 'about to call debug...'
!   debug_file = output_file(trim(config%output_l2)//".debug.hdf", all_output_datasets)
!   status = db_debug_output(debug_file)
!   if (status /= 0) then
!     print *, "ERROR: Failed to dump debug output: ", status
!     stop
!   end if
! 
!   do i = 1, size(all_output_datasets,1)
!     deallocate(all_output_datasets(i)%values, stat=status)
!     if (status /= 0) then
!       print *, 'ERROR: Failed to deallocate all_output_datasets(i)%values: ', i, status
!       stop
!     end if
!   end do
! 
!   deallocate(all_output_datasets)
!   print *, 'After debug...'
! -----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! -- clean up unnecessary variables
!-----------------------------------------------------------------------------------------

  deallocate(sr650, gzflg,sfcstd,snow_mask, snow_mask2,month, stat=status)
  if (status /= 0) then
    print *, "ERROR: Unable to deallocate sr650, gzflg arrays: ", status
    stop
  end if  
  
  deallocate(ocld_mask, lcld_mask,lc, stat=status)
  if (status /= 0) then 
    print *, "ERROR: Unable to deallocate cloud mask arrays: ", status
    stop
  end if  
  
!-----------------------------------------------------------------------------------------
! -- process data
!-----------------------------------------------------------------------------------------
  print *, 'dimension :',viirs_data%xscan, viirs_data%scan, cell_dims(1), cell_dims(2)
  print *, 'about to call modis...'
  call modis(viirs_data%doy, viirs_data%xscan, viirs_data%scan, cell_dims(1), cell_dims(2), &
  &   viirs_data%scan_time, viirs_data%lat, viirs_data%lon, viirs_data%sza, viirs_data%vza, &
  &   viirs_data%raa, viirs_data%sca, viirs_data%m01_refl, viirs_data%m03_refl,             &
  &   viirs_data%m04_refl, viirs_data%m05_refl, viirs_data%m07_refl, viirs_data%m08_refl,   &
  &   viirs_data%m10_refl, viirs_data%m11_refl, m01_nogc, m03_nogc, m05_nogc, m11_oldgc,    &
  &   viirs_data%land_mask, viirs_data%dstar, viirs_data%ndvi, viirs_data%elev, viirs_data%ps, &
  &   lskip_mask, oskip_mask, smoke_mask,  &
  &   smoke_ae_mask,pyrocb_mask, high_alt_smoke_mask,  n_total_pixels, config%output_l2,    &
  &   windsp, winddir, oz, wv, ler412, ler488, ler670, qdf412, qdf488, qdf670, bathy, chl, vcc, &
  &   viirs_data%m15_bt, viirs_data%btd11, hires, cell_resolution,nrt,platform)
  print *, 'back from modis...'
  
!-----------------------------------------------------------------------------------------
! -- clean up
!-----------------------------------------------------------------------------------------

  deallocate(windsp, winddir, oz, wv, stat=status)
  if (status /= 0) then
    print *, "ERROR: Failed to deallocate tmp ozone and water vapor arrays1: ", status
    stop
  end if
  
  deallocate(lskip_mask, oskip_mask, smoke_mask, smoke_ae_mask, stat=status)
  if (status /= 0) then 
    print *, "ERROR: Unable to deallocate skip and smoke mask arrays2: ", status
    stop
  end if
  
  if (platform .eq. 'VIIRS') then  
!     deallocate(viirs_data%m01_refl, viirs_data%m02_refl, viirs_data%m03_refl, &
!     & viirs_data%m04_refl, viirs_data%m05_refl, viirs_data%m07_refl, viirs_data%m08_refl, &
!     & viirs_data%m09_refl, viirs_data%m10_refl, viirs_data%m11_refl, stat=status)
!     if (status /= 0) then
!       print *, "ERROR: Unable to deallocate VIIRS DB reflectance data array3: ", status
!       stop
!     end if
  end if 
  if (platform .eq. 'AHI' ) then  
    deallocate(viirs_data%m03_refl, viirs_data%m04_refl, viirs_data%m05_refl, viirs_data%m07_refl, &
    & viirs_data%m10_refl, viirs_data%m11_refl, stat=status)
    if (status /= 0) then
      print *, "ERROR: Unable to deallocate VIIRS DB reflectance data array4: ", status
      stop
    end if
  end if   
  if ( platform .eq. 'GOES') then  
    deallocate(viirs_data%m03_refl, viirs_data%m05_refl, viirs_data%m07_refl, &
    & viirs_data%m10_refl, viirs_data%m11_refl, stat=status)
    if (status /= 0) then
      print *, "ERROR: Unable to deallocate VIIRS DB reflectance data array5: ", status
      stop
    end if
  end if   
  
  deallocate(viirs_data%lat, viirs_data%lon, viirs_data%sza, viirs_data%vza, &
  & viirs_data%raa, viirs_data%sca, stat=status)
  if (status /= 0) then
    print *, "ERROR: Unable to dellocate VIIRS DB geolocation data arrays6: ", status
    stop
  end if
  
!   deallocate(viirs_data%btd8, viirs_data%btd11, viirs_data%dstar, viirs_data%ndvi, &
!   & viirs_data%land_mask, viirs_data%elev, viirs_data%ps, stat=status)
!   if (status /= 0) then
!     print *, "ERROR: Unable to deallocate VIIRS DB ancillary data arrays: ", status
!     stop
!   end if

! J. Lee  3 October 2017  Removed some of the TOA reflectance variables
  deallocate(m01_nogc, m03_nogc, m05_nogc, stat=status)
  if (status /= 0) then
    print *, "ERROR: Unable to deallocate VIIRS new gas-corrected DB reflectance data array: ", status
    stop
  end if

  deallocate(m11_oldgc, stat=status)
  if (status /= 0) then
    print *, "ERROR: Unable to deallocate VIIRS old gas-corrected DB reflectance data array: ", status
    stop
  end if
! end J. Lee removed some of the TOA reflectance variables
  
  deallocate(n_total_pixels, stat=status)
  if (status /= 0) then 
    print *, "ERROR: Unable to deallocate pixel count array: ", status
    stop
  end if

! -- unload various input data
  
  call unload_landcover(status)
  if (status /= 0) then
    print *, "ERROR: Unable to unload landcover data. Continuing: ", status
  end if
  
  call unload_brdf(status)
  if (status /= 0) then
    print *, "ERROR: Unable to unload BRDF input file. Continuing: ", status
  end if
  
  call unload_viirs_aerosol_luts(status)
  if (status /= 0) then
    print *, "ERROR: Unable to unload VIIRS aerosol LUTS. Continuing: ", status
  end if
  
!-----------------------------------------------------------------------------------------
! -- done
!-----------------------------------------------------------------------------------------
  print *, 'Done : deep_blue.f90'
end program deep_blue

