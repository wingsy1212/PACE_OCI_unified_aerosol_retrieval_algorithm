module modis_surface
   !-----------------------------------------------------------------------------------------
   ! Module to organize functions related the surface reflectivity calculations and tables
   !   for Aqua MODIS.
   !
   ! Corey Bettenhausen
   ! Science Systems and Applications, Inc
   ! NASA Goddard Space Flight Center
   !
   ! Jeremy Warner
   ! Science Systems and Applications, Inc
   ! NASA Goddard Space Flight Center
   !-----------------------------------------------------------------------------------------
   implicit none

   private
  
   public  ::  load_terrainflg_tables, load_seasonal_desert
   public  ::  load_swir_coeffs
   public  ::  get_swir_coeffs412, get_swir_coeffs470
   public  ::  get_swir_stderr412, get_swir_stderr470
   public  ::  get_swir_range
   public  ::  get_brdfcorr_sr
   public  ::  load_brdf, unload_brdf
   public  ::  get_aot500
   public  ::  terrain_flag_new, terrain_flag, sfc_elev_std
   public  ::  set_limits, load_hdfLER, set_limits6
   public  ::  get_LER412, get_LER470, get_LER650, get_LER865
   public  ::  get_modis_LER412, get_modis_LER470, get_modis_LER650, get_modis_LER865
   public  ::  get_viirs_LER412, get_viirs_LER488, get_viirs_LER670
   public  ::  get_viirs_modisbrdf_LER412, get_viirs_modisbrdf_LER488, get_viirs_modisbrdf_LER670
   public  ::  latlon_to_index_ler,check_status
   public  ::  get_geographic_zone, get_sfc_elev_std, get_background_aod! 9 January 2018 JLee
   private ::  readLER5, readLER2
  
   !   -- brdf650    = summer surface reflectivity at 650nm, used to transfer BRDF from
   !                     AERONET site to pixel location.
   !   -- aero_sites = list of AERONET site names as listed on the AERONET website
   !   -- aero_zones = geographical zone index of site based upon seawifs_terrainflag_*.hdf input
   !   -- aero_types = land cover type of AERONET site based upon landcover_*.hdf input
   !   -- aero_elev  = elevation of AERONET site in meters from surface_pressure_*.hdf input
   !   -- aero_sr*   = surface reflectance of AERONET site based on surface database at 135 degrees in
   !                     winter(1), spring(2), summer(3), and fall(4).
   real, dimension(:,:), allocatable               ::  brdf650
   character(len=255), dimension(:), allocatable   ::  aero_sites
   integer, dimension(:), allocatable              ::  aero_zones
   integer, dimension(:), allocatable              ::  aero_types
   integer, dimension(:), allocatable              ::  aero_elev
   real, dimension(:,:), allocatable               ::  aero_sr412, aero_sr470, aero_sr650, aero_bgaod

   integer,dimension(3600,1800)  ::  terrain_flag_new
   real, dimension(3600,1800)    ::  terrain_flag, sfc_elev_std
   real, dimension(360,180)    ::  bg_aod! 9 January 2018 JLee

   integer                                  :: LERstart(2), LERedge(2), dateline
   integer                                  :: LERstart6(2), LERedge6(2), dateline6
   real, dimension(:,:,:,:), allocatable  :: coefs650_fwd, coefs470_fwd, coefs412_fwd
   real, dimension(:,:,:,:), allocatable  :: coefs650_all, coefs470_all, coefs412_all

   real, dimension(:,:), allocatable        :: gref412_all, gref412_fwd
   real, dimension(:,:), allocatable        :: gref470_all, gref470_fwd
   real, dimension(:,:), allocatable        :: gref650_all, gref650_fwd
   real, dimension(:,:), allocatable        :: gref865_all
  
   ! -- VIIRS, all-angle surface database
   real, dimension(:,:), allocatable        :: vgref412_all
   real, dimension(:,:), allocatable        :: vgref488_all
   real, dimension(:,:), allocatable        :: vgref670_all

   ! -- 2.2 um surface database
   real, dimension(:,:,:), allocatable      :: swir_coeffs412, swir_coeffs470
   real, dimension(:,:), allocatable        :: swir_stderr412, swir_stderr470
   real, dimension(:,:), allocatable        :: swir_min, swir_max
  
   real, parameter   ::  NDVI1_CUTOFF = 0.18
   real, parameter   ::  NDVI2_CUTOFF = 0.35

contains

   ! @TODO: should rename this to load_geozone_table or something even more descriptive.
   integer function load_terrainflg_tables(tflg_file, season) result(status)
      use netcdf
      USE OCIUAAER_Config_Module
      implicit none

      character(len=255), intent(in)  ::  tflg_file
      integer, intent(in)    ::  season
    
      real, dimension(:,:), allocatable   ::  tmptfn
      real, dimension(:,:,:), allocatable ::  tmpaod
      integer                             ::  i, j

      ! HDF vars
      character(len=255)    ::  sds_name
      character(len=255)    ::  dset_name
      character(len=255)    ::  attr_name
      character(len=255)    ::  group_name

      integer               ::  nc_id
      integer               ::  dim_id
      integer               ::  dset_id
      integer               ::  grp_id
      integer               ::  sd_id, sds_index, sds_id
      integer, dimension(2) ::  start2, stride2, edges2
      integer, dimension(3) ::  start3, stride3, edges3
 
      status = -1
    
      !   -- allocate our tmp array
      allocate(tmptfn(3600,1800), stat=status)
      if (status /= 0) then
         print *, "ERROR: Unable to allocate tmp array for geo zone data: ", status
         return
      end if

      !   -- allocate our tmp array
      allocate(tmpaod(360,180,4), stat=status)
      if (status /= 0) then
         print *, "ERROR: Unable to allocate tmp array for background aod data: ", status
         return
      end if

      status = nf90_open(cfg%db_nc4, nf90_nowrite, nc_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to open deepblue lut_nc4 file: ", status
         return
      end if

      group_name = 'geozone'
      status = nf90_inq_ncid(nc_id, group_name, grp_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of group "//trim(group_name)//": ", status
         return
      end if

      start2  = (/1,1/)
      stride2 = (/1,1/)
      edges2  = (/3600,1800/)
      dset_name = 'geographical_zone_flag'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      status = nf90_get_var(grp_id, dset_id, tmptfn, start=start2, &
         stride=stride2, count=edges2)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if

      dset_name = 'surface_elevation_stddev'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      status = nf90_get_var(grp_id, dset_id, sfc_elev_std, start=start2, &
         stride=stride2, count=edges2)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if

      start3  = (/1,1,1/)
      stride3 = (/1,1,1/)
      edges3  = (/360,180,4/)
      dset_name = 'background_aod'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      status = nf90_get_var(grp_id, dset_id, tmpaod, start=start3, &
         stride=stride3, count=edges3)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if

      status = nf90_close(nc_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to close lut_nc4 file: ", status
         return
      end if

      terrain_flag_new = int(tmptfn)
      bg_aod(1:360,1:180) = tmpaod(1:360,1:180,season)
      !print *, "aod test", bg_aod(1,68), tmpaod(1,68,season)
 
      !   -- clean up tmptfn
      deallocate(tmptfn, stat=status)
      if (status /= 0) then
         print *, "ERROR: Unable to deallocate tmp array for geo zone data: ", status
         return
      end if

      !   -- clean up tmpaod
      deallocate(tmpaod, stat=status)
      if (status /= 0) then
         print *, "ERROR: Unable to deallocate tmp array for geo zone data: ", status
         return
      end if

      status = 0
      return
   end function load_terrainflg_tables
  
   integer function load_seasonal_desert(file) result(status)
    
      !   include 'hdf.f90'
      !   include 'dffunc.f90'
      use netcdf
      USE OCIUAAER_Config_Module

      implicit none

      character(len=255), intent(in)  ::  file
    
      real, dimension(:,:), allocatable   ::  tmptfn
      integer                             ::  i, j

      ! HDF vars
      character(len=255)    ::  sds_name
      character(len=255)    ::  dset_name
      character(len=255)    ::  attr_name
      character(len=255)    ::  group_name

      integer               ::  nc_id
      integer               ::  dim_id
      integer               ::  dset_id
      integer               ::  grp_id
      integer               ::  sd_id, sds_index, sds_id
      integer, dimension(2) ::  start2, stride2, edges2
    
      status = -1
    
      status = nf90_open(cfg%db_nc4, nf90_nowrite, nc_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to open deepblue lut_nc4 file: ", status
         return
      end if

      group_name = 'seasonal_deserts'
      status = nf90_inq_ncid(nc_id, group_name, grp_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of group "//trim(group_name)//": ", status
         return
      end if

      start2  = (/1,1/)
      stride2 = (/1,1/)
      edges2  = (/3600,1800/)
      dset_name = 'seasonal_desert_flag'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
      status = nf90_get_var(grp_id, dset_id, terrain_flag, start=start2, &
         stride=stride2, count=edges2)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if

      status = nf90_close(nc_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to close lut_nc4 file: ", status
         return
      end if

      status = 0
      return
   end function load_seasonal_desert

   ! -- initialize the AERONET site and surface variables and load the BRDF base
   ! --    reflectivity file into brdf650 array.
   !-----------------------------------------------------------------------------------------
   integer function load_brdf(brdffile) result(status)
    
      !   include 'hdf.f90'
      !   include 'dffunc.f90'
      use netcdf
      USE OCIUAAER_Config_Module

      implicit none

      character(len=255), intent(in)    ::  brdffile
    
      integer, parameter              ::  nsites = 30
    
      character(len=255)    ::  sds_name
      character(len=255)    ::  dset_name
      character(len=255)    ::  attr_name
      character(len=255)    ::  group_name

      integer               ::  nc_id
      integer               ::  dim_id
      integer               ::  dset_id
      integer               ::  grp_id
      integer               ::  sd_id, sds_index, sds_id
      integer, dimension(2) ::  start2, stride2, edges2, dims2
    
      !   -- assume a successful return.
      status = 0
    
      !   -- allocate and fill AERONET-related arrays.
      allocate(aero_sites(nsites), aero_zones(nsites), aero_types(nsites), aero_elev(nsites), &
         & stat=status)
      if (status /= 0) then
         print *, "ERROR: Unable to allocate AERONET site info arrays: ", status
         return
      end if
    
      allocate(aero_sr412(nsites,4), aero_sr470(nsites,4), aero_sr650(nsites,4), aero_bgaod(nsites,4), stat=status)
      if (status /= 0) then
         print *, "ERROR: Unable to allocate AERONET SR arrays: ", status
         return
      end if

      !   -- aero_types = land cover type over AERONET site
      !   --    0. ocean
      !   --    1. forest
      !   --    2. grasslands
      !   --    3. croplands
      !   --    4. urban
      !   --    5. snow/ice
      !   --    6. barren/desert

      !   -- Banizoumbou
      aero_sites(1)   = 'Banizoumbou'
      aero_types(1)   = 2
      aero_zones(1)   = 5
      aero_elev(1)    = 250
      aero_sr412(1,:) = (/7.08923, 7.71880, 8.48224, 6.62584/)
      aero_sr470(1,:) = (/10.5942, 11.6695, 12.4470, 9.9028/)
      aero_sr650(1,:) = (/28.7862, 31.9045, 31.5499, 25.0119/)
      aero_bgaod(1,:) = (/0.15000, 0.22800, 0.26300, 0.18500/) !MODIS
      !    aero_bgaod(1,:) = (/0.18200, 0.41348, 0.37300, 0.19649/) !MISR
    
      !   -- Tinga Tingana
      aero_sites(2)   = 'Tinga_Tingana'
      aero_types(2)   = 6
      aero_zones(2)   = 12
      aero_elev(2)    = 38
      aero_sr412(2,:) = (/8.3397, 9.3348, 8.3018, 10.4549/)
      aero_sr470(2,:) = (/11.5649, 12.8902, 11.5255, 13.5749/)
      aero_sr650(2,:) = (/29.0277, 31.3340, 27.9479, 30.2249/)
      aero_bgaod(2,:) = (/0.02300, 0.01900, 0.01800, 0.02100/)
      !    aero_bgaod(2,:) = (/0.08670, 0.06223, 0.05391, 0.10703/)
 
      !   -- Zinder_Airport
      aero_sites(3)   = 'Zinder_Airport'
      aero_types(3)   = 2
      aero_zones(3)   = 1
      aero_elev(3)    = 456
      aero_sr412(3,:) = (/7.59892, 8.63954, 8.38775, 6.58359/)
      aero_sr470(3,:) = (/11.3537, 12.8399, 12.3052, 9.87239/)
      aero_sr650(3,:) = (/26.7511, 30.1220, 27.5969, 21.1145/)
      aero_bgaod(3,:) = (/0.15400, 0.29600, 0.30200, 0.14300/)
      !    aero_bgaod(3,:) = (/0.15523, 0.37014, 0.35720, 0.17437/)
    
      !   -- Moldova
      aero_sites(4)   = 'Moldova'
      aero_types(4)   = 4
      aero_zones(4)   = 17
      aero_elev(4)    = 205
      aero_sr412(4,:) = (/-999.000, 6.10670, 5.23876, 5.76111/)
      aero_sr470(4,:) = (/-999.000, 7.08618, 6.13133, 6.63834/)
      aero_sr650(4,:) = (/-999.000, 8.11562, 7.51764, 7.96273/)
      aero_bgaod(4,:) = (/0.01900, 0.05900, 0.07700, 0.04300/)
      !    aero_bgaod(4,:) = (/0.07450, 0.12211, 0.12758, 0.07907/)

      !   -- Beijing
      aero_sites(5)   = 'Beijing'
      aero_types(5)   = 4
      aero_zones(5)   = 16
      aero_elev(5)    = 92
      aero_sr412(5,:) = (/5.04678, 6.91650, 6.04189, 5.75009/)
      aero_sr470(5,:) = (/7.59624, 8.45864, 7.31190, 7.05383/)
      aero_sr650(5,:) = (/11.7573, 11.3573, 9.2663, 8.54484/)
      aero_bgaod(5,:) = (/0.13800, 0.18100, 0.14400, 0.11400/)
      !    aero_bgaod(5,:) = (/0.12670, 0.19608, 0.20208, 0.11582/)
    
      !   -- Kanpur, India except for urban areas and Thar Desert
      aero_sites(6)    = 'Kanpur'
      aero_types(6)    = 3
      aero_zones(6)    = 15
      aero_elev(6)     = 123
      !   use fall BRDF for summer, 26 January 2018 JLee, TEST
      aero_sr412(6,:)  = (/8.76996,6.24308,8.49987,7.49987/)
      aero_sr470(6,:)  = (/6.26996,5.74308,7.99987,5.99987/)
      aero_sr650(6,:)  = (/10.25542,10.1785,11.75790,10.75790/)
      aero_bgaod(6,:) = (/0.27400, 0.24800, 0.27600, 0.27900/)
      !    aero_bgaod(6,:) = (/0.31233, 0.28328, 0.44091, 0.36488/)
    
      !   -- Modena
      aero_sites(7)   = 'Modena'
      aero_types(7)   = 4
      aero_zones(7)   = 17
      aero_elev(7)    = 56
      aero_sr412(7,:) = (/3.7846,5.30996,5.72852,5.69932/)
      aero_sr470(7,:) = (/5.3279,6.31288,7.11941,6.54509/)
      aero_sr650(7,:) = (/5.5748,9.44464,9.73280,10.1099/)
      aero_bgaod(7,:) = (/0.02100, 0.09200, 0.09900, 0.03900/)
      !    aero_bgaod(7,:) = (/0.06867, 0.13446, 0.15744, 0.09253/)

      !   -- Palencia
      aero_sites(8)   = 'Palencia'
      aero_types(8)   = 3
      aero_zones(8)   = 17
      aero_elev(8)    = 750
      aero_sr412(8,:) = (/-999.000,5.03951,4.76740,-999.000/)
      aero_sr470(8,:) = (/-999.000,6.17346,6.62762,-999.000/)
      aero_sr650(8,:) = (/-999.000,8.70520,10.7275,-999.000/)
      aero_bgaod(8,:) = (/0.02500, 0.04000, 0.02700, 0.03300/)
      !    aero_bgaod(8,:) = (/0.05151, 0.09216, 0.10017, 0.06591/)
    
      !   -- Lecce_University
      aero_sites(9)   = 'Lecce_University'
      aero_types(9)   = 2
      aero_zones(9)   = 17
      aero_elev(9)    = 30
      aero_sr412(9,:) = (/4.68698,4.01772,6.14852,5.86577/)
      aero_sr470(9,:) = (/5.16447,5.68885,8.15023,6.82434/)
      aero_sr650(9,:) = (/10.1211,10.6153,11.0301,11.3503/)
      aero_bgaod(9,:) = (/0.05700, 0.08800, 0.06200, 0.06800/)
      !    aero_bgaod(9,:) = (/0.06377, 0.12146, 0.13148, 0.08287/)

      !   -- Fresno_2
      !   This controls N. America urban areas. Fresno AERONET validation is affected by
      !   both Fresno_2 and Fresno_GZ18. Currently optimized for general urban areas,
      !   as thinking of making separate geozone for Fresno if validation is not
      !   acceptable at Fresno. Only baseline change would be needed for the new zone.
      aero_sites(10)   = 'Fresno_2'
      aero_types(10)   = 4
      aero_zones(10)   = 13
      aero_elev(10)    = 0.0
      aero_sr412(10,:) = (/5.68700,4.64569,4.42003,4.78884/)
      aero_sr470(10,:) = (/6.91660,6.96638,6.80781,6.73356/)
      aero_sr650(10,:) = (/11.5361,12.2151,12.2892,12.6795/)
      aero_bgaod(10,:) = (/0.05300, 0.11500, 0.08500, 0.08100/)
      !    aero_bgaod(10,:) = (/0.08402, 0.13809, 0.14728, 0.10320/)
    
      !   -- Fresno (Central Valley)
      aero_sites(11)   = 'Fresno_GZ18'
      aero_types(11)   = 2
      aero_zones(11)   = 18
      aero_elev(11)    = 0.0
      aero_sr412(11,:) = (/6.18700,5.14569,4.92003,5.28884/)
      aero_sr470(11,:) = (/7.41660,7.46638,7.30781,7.23356/)
      aero_sr650(11,:) = (/11.5361,12.2151,12.2892,12.6795/)
      aero_bgaod(11,:) = (/0.05300, 0.11500, 0.08500, 0.08100/)
      !    aero_bgaod(11,:) = (/0.08402, 0.13809, 0.14728, 0.10320/)

      !   -- IER_Cinzana
      aero_sites(12)   = 'IER_Cinzana'
      aero_types(12)   = 2
      aero_zones(12)   = 5
      aero_elev(12)    = 285
      aero_sr412(12,:) = (/5.33969,6.89590,7.78313,5.45146/)
      aero_sr470(12,:) = (/8.17876,10.2201,11.0532,7.67885/)
      aero_sr650(12,:) = (/18.6043,21.9242,19.8147,13.6748/)
      aero_bgaod(12,:) = (/0.16800, 0.24200, 0.12900, 0.17400/)
      !    aero_bgaod(12,:) = (/0.14072, 0.32845, 0.29905, 0.18086/)
    
      !   -- Agoufou
      aero_sites(13)   = 'Agoufou'
      aero_types(13)   = 2
      aero_zones(13)   = -1 !5
      aero_elev(13)    = 305
      aero_sr412(13,:) = (/6.33764,7.20075,7.12166,5.88014/)
      aero_sr470(13,:) = (/10.3036,11.2734,10.7413,9.34117/)
      aero_sr650(13,:) = (/26.6428,30.4116,27.0584,21.6639/)
      aero_bgaod(13,:) = (/0.11800, 0.20500, 0.19900, 0.13200/)
      !    aero_bgaod(13,:) = (/0.12562, 0.29455, 0.39775, 0.17390/)

      !   -- Saada -- leave disabled, decided not to use it. troublesome site.
      aero_sites(14)   = 'Saada'
      aero_types(14)   = 3
      aero_zones(14)   = -1 !2
      aero_elev(14)    = 420
      aero_sr412(14,:) = (/7.30339, 5.90723, 6.37791, 6.20939/)
      aero_sr470(14,:) = (/8.68933, 7.76850, 8.46196, 8.15088/)
      aero_sr650(14,:) = (/14.1430, 14.5881, 16.7061, 15.5649/)
      aero_bgaod(14,:) = (/0.08300, 0.06400, 0.08800, 0.08700/)
      !    aero_bgaod(14,:) = (/0.03898, 0.06964, 0.08859, 0.07119/)

      !   -- Trelew (S. America)
      aero_sites(15)   = 'Trelew'
      aero_types(15)   = 6
      aero_zones(15)   = 14
      aero_elev(15)    = 15
      aero_sr412(15,:) = (/5.29937, 5.30638, 6.01197, 5.75946/)
      aero_sr470(15,:) = (/8.20220, 7.37385, 7.43250, 7.71553/)
      aero_sr650(15,:) = (/14.0610, 11.7312, 11.2763, 12.9785/)
      aero_bgaod(15,:) = (/0.02200, 0.01900, 0.01700, 0.01900/)
      !    aero_bgaod(15,:) = (/0.06490, 0.03365, 0.03397, 0.05836/)

      !   -- Carpentras
      aero_sites(16)   = 'Carpentras'
      aero_types(16)   = 3
      aero_zones(16)   = 17
      aero_elev(16)    = 100
      aero_sr412(16,:) = (/-999.000,4.27180,3.84850,3.60839/)
      aero_sr470(16,:) = (/-999.000,5.77824,5.63915,5.02537/)
      aero_sr650(16,:) = (/-999.000,9.71739,9.57229,8.67115/)
      aero_bgaod(16,:) = (/0.01900, 0.03200, 0.03500, 0.02100/)
      !    aero_bgaod(16,:) = (/0.03448, 0.06474, 0.05331, 0.03388/)

      !   -- 25km BRDF
      !    aero_sr412(16,:) = (/-999.000,4.72180,4.64850,4.20839/)
      !    aero_sr470(16,:) = (/-999.000,6.22824,6.43915,5.62537/)
      !    aero_sr650(16,:) = (/-999.000,9.71739,9.57229,8.67115/)

      !        -- Pune, India urban areas
      aero_sites(17) = 'Pune'
      aero_types(17)   = 4
      aero_zones(17)   = 19
      aero_elev(17)    = 559
      aero_sr412(17,:) = (/4.49376,6.22264,4.81305,7.31305/)
      aero_sr470(17,:) = (/5.42197,8.08891,5.49410,7.99410/)
      aero_sr650(17,:) = (/8.40501,11.6605,6.73313,9.23313/)
      aero_bgaod(17,:) = (/0.20400, 0.16400, 0.07500, 0.17100/)
      !    aero_bgaod(17,:) = (/0.14888, 0.21124, 0.27165, 0.19066/)
    
      !        -- Evora, Spain
      aero_sites(18)     = 'Evora'
      aero_types(18)   = 3
      aero_zones(18)   = 22
      aero_elev(18)    = 293
      aero_sr412(18,:) = (/4.95347,4.48004,4.75238,5.64016/)
      aero_sr470(18,:) = (/5.60902,5.80674,7.54495,7.83002/)
      aero_sr650(18,:) = (/6.80235,6.94325,13.3975,11.9871/)
      aero_bgaod(18,:) = (/0.01900, 0.03400, 0.02600, 0.02500/)
      !    aero_bgaod(18,:) = (/0.03544, 0.06548, 0.09886, 0.05251/)

      !        -- Blida, N. Africa
      aero_sites(19)     = 'Blida'
      aero_types(19)   = 3
      aero_zones(19)   = -1 !2
      aero_elev(19)    = 230
      aero_sr412(19,:) = (/-999.000,5.20722,5.84409,-999.000/)
      aero_sr470(19,:) = (/-999.000,7.35584,7.89343,-999.000/)
      aero_sr650(19,:) = (/-999.000,11.1594,13.5330,-999.000/)
      aero_bgaod(19,:) = (/0.04400, 0.04900, 0.07700, 0.05800/)
      !    aero_bgaod(19,:) = (/0.07116, 0.10362, 0.15463, 0.09386/)

      !        -- Blida, N. Africa
      aero_sites(20)     = 'Blida_High'
      aero_types(20)   = 3
      aero_zones(20)   = -1 !2
      aero_elev(20)    = 600
      aero_sr412(20,:) = (/-999.000,5.20722,5.84409,-999.000/)
      aero_sr470(20,:) = (/-999.000,7.35584,7.89343,-999.000/)
      aero_sr650(20,:) = (/-999.000,11.1594,13.5330,-999.000/)
      aero_bgaod(20,:) = (/0.04400, 0.04900, 0.07700, 0.05800/)
      !    aero_bgaod(20,:) = (/0.07116, 0.10362, 0.15463, 0.09386/)

      !   -- GZ24_Only, covers Taklimakan Desert, only AOT models are used here.
      aero_sites(21)   = 'GZ24_Only'
      aero_types(21)   = -1
      aero_zones(21)   = 24
      aero_elev(21)    = -1
      aero_sr412(21,:) = (/-999.0,-999.0,-999.0,-999.0/)      ! new from surf. coeffs.
      aero_sr470(21,:) = (/-999.0,-999.0,-999.0,-999.0/)
      aero_sr650(21,:) = (/-999.0,-999.0,-999.0,-999.0/)
      aero_bgaod(21,:) = (/ -999.0,  -999.0,  -999.0,  -999.0/)
      !    aero_bgaod(21,:) = (/ -999.0,  -999.0,  -999.0,  -999.0/)
   
      !   -- Ilorin
      aero_sites(22)   = 'Ilorin'
      aero_types(22)   = 2
      aero_zones(22)   = 26
      aero_elev(22)    = 350
      aero_sr412(22,:) = (/4.79848, 4.13429, -999.000, -999.000/)
      aero_sr470(22,:) = (/5.73108, 5.07124, -999.000, -999.000/)
      aero_sr650(22,:) = (/10.0571, 9.28994, -999.000, -999.000/)
      aero_bgaod(22,:) = (/0.34500, 0.29700, 0.17300, 0.16400/)
      !    aero_bgaod(22,:) = (/0.32718, 0.32547, -999.00000, 0.22936/)

      !   -- CCNY
      aero_sites(23)   = 'CCNY'
      aero_types(23)   = 4
      aero_zones(23)   = 25
      aero_elev(23)    = 0.0
      aero_sr412(23,:) = (/5.7380,6.3655,8.7437,5.3349/)
      aero_sr470(23,:) = (/7.0723,7.5391,8.8168,6.8278/)
      aero_sr650(23,:) = (/10.1025,10.7149,10.1311,10.5906/)
      aero_bgaod(23,:) = (/0.04800, 0.06600, 0.13500, 0.06100/)
      !    aero_bgaod(23,:) = (/0.06593, 0.09566, 0.09318, 0.05144/)

      !   -- Ilorin
      aero_sites(24)   = 'Ilorin_Transition'
      aero_types(24)   = 2
      !    aero_zones(24)   = -1 !27
      aero_zones(24)   = 27
      aero_elev(24)    = 350
      aero_sr412(24,:) = (/4.79848, 4.13429, -999.000, -999.000/)
      aero_sr470(24,:) = (/5.73108, 5.07124, -999.000, -999.000/)
      aero_sr650(24,:) = (/10.0571, 9.28994, -999.000, -999.000/)
      aero_bgaod(24,:) = (/0.34500, 0.29700, 0.17300, 0.16400/)
      !    aero_bgaod(24,:) = (/0.32718, 0.32547, -999.00000, 0.22936/)

      !   -- SACOL
      aero_sites(25)   = 'SACOL'
      aero_types(25)   = 2
      aero_zones(25)   = 28
      aero_elev(25)    = 1965
      aero_sr412(25,:) = (/6.57751, 5.85782, 4.26251, 5.79214/)
      aero_sr470(25,:) = (/8.5020, 8.2185, 5.56137, 6.24013/)
      aero_sr650(25,:) = (/16.6909, 16.8518, 11.5214, 12.5133/)
      aero_bgaod(25,:) = (/0.03700, 0.05400, 0.05700, 0.03400/)
      !    aero_bgaod(25,:) = (/0.12611, 0.20761, 0.17829, 0.11342/)
    
      !   -- Mexico_City
      aero_sites(26)   = 'Mexico_City'
      aero_types(26)   = 4
      aero_zones(26)   = 29
      aero_elev(26)    = 2268.0
      aero_sr412(26,:) = (/6.73461, 6.20030, -999.000, 8.10955 /)
      aero_sr470(26,:) = (/7.50571, 7.88785, -999.000, 9.46562/)
      aero_sr650(26,:) = (/7.7320, 10.2994, -999.000, 11.9709/)
      aero_bgaod(26,:) = (/0.01900, 0.02100, 0.03900, 0.02600/)
      !    aero_bgaod(26,:) = (/0.07039, 0.10752, 0.11487, 0.09446/)
    
      !   -- Solar Village
      aero_sites(27)   = 'Solar_Village'
      aero_types(27)   = 6
      aero_zones(27)   = 10
      aero_elev(27)    = 764.0
      aero_sr412(27,:) = (/10.4297, 10.8623, 10.7472, 11.9705/)
      aero_sr470(27,:) = (/15.0892, 16.1351, 16.0690, 17.0390/)
      aero_sr650(27,:) = (/32.0747, 34.5677, 35.3692, 34.6681/)
      aero_bgaod(27,:) = (/0.10100, 0.09800, 0.16700, 0.10900/)
      !    aero_bgaod(27,:) = (/0.14651, 0.27687, 0.34036, 0.23912/)
    
      !   -- Jaipur, Thar Desert
      aero_sites(28)   = 'Jaipur'
      aero_types(28)   = 4
      aero_zones(28)   = 20
      aero_elev(28)    = 450.0
      aero_sr412(28,:) = (/6.46991, 7.40196, 7.28651, 5.22799/)
      aero_sr470(28,:) = (/8.49850, 9.42026, 9.49201, 7.03474/)
      aero_sr650(28,:) = (/11.3653, 12.0653, 15.2039, 10.3618/)
      aero_bgaod(28,:) = (/0.07100, 0.10500, 0.09500, 0.07100/)
      !    aero_bgaod(28,:) = (/0.16877, 0.23449, 0.43655, 0.22099/)

      !   -- NW_India_Desert
      aero_sites(29)   = 'NW_India_Desert'
      aero_types(29)   = 4
      aero_zones(29)   = 30
      aero_elev(29)    = 450.0
      aero_sr412(29,:) = (/7.09280, 5.90470, 6.97091, 4.24017/)
      aero_sr470(29,:) = (/8.31369, 7.69160, 8.73495, 5.41526/)
      aero_sr650(29,:) = (/11.8653, 13.1653, 15.0039, 10.9618/)
      aero_bgaod(29,:) = (/0.07100, 0.10500, 0.09500, 0.07100/)
      !    aero_bgaod(29,:) = (/0.16877, 0.23449, 0.43655, 0.22099/)
    
      !   --Yuma
      aero_sites(30)   = 'Yuma'
      aero_types(30)   = 6
      aero_zones(30)   = 31
      aero_elev(30)    = 63
      aero_sr412(30,:) = (/6.7668, 6.9406, 8.3705, 7.4484/)
      aero_sr470(30,:) = (/9.8905, 9.9898, 10.8432, 10.9740/)
      aero_sr650(30,:) = (/24.7466, 24.4755, 25.5649, 25.4377/)
      aero_bgaod(30,:) = (/0.07100, 0.12300, 0.11500, 0.08600/)
      !    aero_bgaod(30,:) = (/0.05629, 0.12101, 0.13732, 0.07340/)

      !   -- read in base BRDF reflectivitiy @ 650nm from infile.
      allocate(brdf650(3600,1800), stat=status)
      if (status /= 0) then
         print *, "ERROR: Unable to allocate array for BRDF base data: ", status
         return
      end if

      status = nf90_open(cfg%db_nc4, nf90_nowrite, nc_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to open deepblue lut_nc4 file: ", status
         return
      end if

      group_name = 'brdfbase'
      status = nf90_inq_ncid(nc_id, group_name, grp_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of group "//trim(group_name)//": ", status
         return
      end if

      dset_name = 'brdf_base_650'
      status = nf90_inq_varid(grp_id, dset_name, dset_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
         return
      end if
    
      start2  = (/1,1/)
      stride2 = (/1,1/)
      edges2  = (/3600,1800/)
      status = nf90_get_var(grp_id, dset_id, brdf650, start=start2, &
         stride=stride2, count=edges2)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
         return
      end if

      status = nf90_close(nc_id)
      if (status /= NF90_NOERR) then
         print *, "ERROR: Failed to close lut_nc4 file: ", status
         return
      end if
    
      return
   end function load_brdf

   ! -- deallocate brdf650 array and AERONET site and surface reflectivity variables.
   subroutine unload_brdf(status)
      implicit none

      integer, intent(inout)    ::  status
    
      deallocate(brdf650, aero_sites, aero_types, aero_zones, stat=status)
      if (status /= 0) then
         print *, "ERROR: Unable to deallocate BRDF data: ", status
         return
      end if
    
      deallocate(aero_sr412, aero_sr470, aero_sr650, aero_bgaod, stat=status)
      if (status /= 0) then
         print *, "ERROR: Unable to deallocate AERONET SR arrays: ", status
         return
      end if
    
      return
    
   end subroutine unload_brdf
  
! -- calculate and return the BRDF-corrected surface reflectivity for pixel at lat,lon
! --  based on matching AERONET BRDF's.
! -- return codes:
! --  0 = success
! -- -1 = failure due to no AERONET BRDF
! -- -2 = failure due to no baseline surface reflectivity values.
  integer function get_brdfcorr_sr(lat, lon, ra, sa, vza, amf, elev, month, ndvi, stdv, gzone, lc_type, bgaod, &
                                 & sr412, sr470, sr650, use_alternate_brdf, debug) result(status)
    
    implicit none
    
    character(len=50), parameter    ::  func_name = "get_brdfcorr_sr"
    
    real, intent(in)          ::  lat
    real, intent(in)          ::  lon
    real, intent(in)          ::  ra
    real, intent(in)          ::  sa              ! -- scattering angle
    real, intent(in)          ::  vza             ! -- viewing zenith angle
    real, intent(in)          ::  amf             ! -- air mass factor
    real, intent(in)          ::  elev            ! -- surface elevation in meters
    integer, intent(in)       ::  month
    real, intent(in)          ::  ndvi
    real, intent(in)          ::  stdv            ! std. dev. of TOA412, 0.2deg radius
    integer, intent(in)       ::  gzone
    integer, intent(in)       ::  lc_type
    real, intent(in)          ::  bgaod
    real, intent(inout)       ::  sr412
    real, intent(inout)       ::  sr470
    real, intent(inout)       ::  sr650
    logical, intent(in), optional :: use_alternate_brdf        ! use alternate BRDF fits
    logical, intent(in), optional :: debug
    
    real                      ::  refsr650
    integer                   ::  ilat, ilon
    integer                   ::  m
    integer                   ::  fillcnt,min_flag
    
    character(len=255)                  ::  asite
    real, dimension(:), allocatable     ::  maero412, maero470, maero650
    real, dimension(:), allocatable     ::  mbgaod
    integer, dimension(:), allocatable  ::  msiteindx
    integer, dimension(:), allocatable  ::  sorted
    real                                ::  ab412, ab470, ab650     ! AERONET BRDF SR's
    real                                ::  ac412, ac470, ac650     ! AERONET BRDF constants
    real                                ::  xnorm_412_f1, xnorm_412_f2
    real                                ::  xnorm_470_f1, xnorm_470_f2
    real                                ::  frac, normfrac
    real                                ::  aod_corr_factor
    integer                             ::  i, ii, jj, cnt
    integer                             ::  season
     
    logical                             ::  dflag
        
    real                  ::  mb_sr412, mb_sr470, mb_sr650
    real                  ::  m_sr412, m_sr470, m_sr650
    real                  ::  v_sr412, v_sr488, v_sr670
      
    dflag = .false.
    min_flag  = 0
    if (present(debug)) then
      dflag = debug
    end if
    
    if (dflag) then
      print *, trim(func_name)//', lat, lon, raa, scat, elev, month, ndvi, gzone, lc, sr412, sr470, sr650: ' &
      & , lat, lon, ra, sa, elev, month, ndvi, gzone, lc_type, sr412, sr470, sr650
    end if
    
!   -- convert geolocation into array indices.
    ilat = floor(lat*10.0) + 900 + 1
    ilon = floor(lon*10.0) + 1800 + 1
    
    if (ilat > 1800) ilat = 1800
    if (ilon > 3600) ilon = 3600
    if (dflag) print *, trim(func_name)//', lat, lon, ilat, ilon: ', lat, lon, ilat, ilon
    
!   -- convert month to season
    select case (month)
      case (12,1,2)
        season = 1
      case(3:5)
        season = 2
      case (6:8)
        season = 3
      case (9:11)
        season = 4
      case default
        print *, "ERROR: Invalid month specified: ", month
        status = -1
        return
    end select
    
!   -- set up our reference surface reflectance
    refsr650 = brdf650(ilon,ilat)
    
!   -- do we have an AERONET site in the same zone with the same land cover type?
    m = 0
    m = count(aero_zones == gzone .AND. aero_types == lc_type .AND. (elev < 500 .EQV. aero_elev < 500))
    
!   -- create exception for China, europe, morocco, spain -- only match zone
    if (gzone == 16 .OR. (gzone == 17 .OR. gzone == 2) .OR. gzone == 22) then
      m = count(aero_zones == gzone .AND. (elev < 500 .EQV. aero_elev < 500))
    end if
    
!     -- create exception for Fresno Valley, Australia, tropical Sahel, Mexico_City - match by region only.
      if (gzone == 18 .OR. gzone == 12 .OR. (gzone == 26 .OR. gzone == 27) .OR. gzone == 29) then
         m = count(aero_zones == gzone)
      end if

!     -- create exception for high elevation Tibet/China zone - match by region only.
      if (gzone == 28) then
           m = count(aero_zones == gzone)
      end if

!     -- create exception for Jaipur zone - match by region only.
      if (gzone == 20) then
           m = count(aero_zones == gzone)
      end if

      !     -- create exception for NW_India_Desert zone - match by region only.
      if (gzone == 30) then
           m = count(aero_zones == gzone)
      end if

!     -- create exception for Pune - match by region only.
      if (gzone == 19) then
         m = count(aero_zones == gzone)
      end if

!   -- create exception for Kanpur - match by region only, 9 January 2018 JLee
    if (gzone == 15) then
      m = count(aero_zones == gzone)
    end if

!   -- create exception for Sahel, geozone =5, landcover=2 to ignore elevation.
!   JLee added gzone 13 (N. America, urban)
    if ((gzone == 5 .AND. lc_type == 2) .OR. gzone == 1 .or. gzone == 13) then
      m = count(aero_zones == gzone .AND. aero_types == lc_type)
    end if

!   -- create exception for Barren North America, geozone =31 to ignore elevation above 750m.
    if (gzone == 31 .AND. elev < 750) then
      m = count(aero_zones == gzone)
    end if

!   -- if barren, use surface tables except in zone 2 (morocco), 12 (australia), and 14 (S. America)
    if (gzone /= 2 .AND. gzone /= 12 .AND. gzone /= 14 .AND. gzone /= 28 .AND. &
    &   gzone /= 10 .AND. gzone /= 20 .AND. gzone /= 30 .AND. gzone /= 31) then
      if (lc_type == 6) m = 0   ! reset over barren surfaces to force use of surface tables.
    end if
    
    if (m > 0) then
!     -- allocate our arrays to store the matching AERONET data.
!     -- no explicit deallocate() as these should automatically be
!     -- deallocated when the function ends.
      if (allocated(maero412)) deallocate(maero412, stat=status)
      if (allocated(maero470)) deallocate(maero470, stat=status)
      if (allocated(maero650)) deallocate(maero650, stat=status)
      if (allocated(mbgaod)) deallocate(mbgaod, stat=status)
      if (allocated(msiteindx)) deallocate(msiteindx, stat=status)
      if (allocated(sorted)) deallocate(sorted, stat=status)
      allocate(maero412(m), maero470(m), maero650(m), mbgaod(m), msiteindx(m), sorted(m), stat=status)
      if (status /= 0) then
        print *, "ERROR: Failed to allocate AERONET 650 SR match arrays: ", status
        return
      end if
      
      cnt = 0
    
!     -- get and store base table SR values at each matching AERONET site at 412, 470, and 650.
      do i = 1, size(aero_sites)        ! i = AERONET site index
        select case (gzone)
          case (2, 16, 17, 22)   ! -- china, europe, spain only match by zone and elevation -- no land cover.
            if (aero_zones(i) == gzone .AND. (elev < 500 .EQV. aero_elev(i) < 500)) then
              cnt = cnt + 1
              maero412(cnt) = aero_sr412(i,season)
              maero470(cnt) = aero_sr470(i,season)
              maero650(cnt) = aero_sr650(i,3)        ! < -- always use summer for 650nm to match refsr650
              mbgaod(cnt)   = aero_bgaod(i,season)
              msiteindx(cnt) = i
             
              if (dflag) then
                print '(A,A,A,I4,I4)', trim(func_name), ', matching site: ', trim(aero_sites(i)), aero_zones(i), aero_types(i)
                print '(A,A,3(F11.6,1X))', trim(func_name), ', AERONET Baseline SR: ', maero412(cnt), maero470(cnt), maero650(cnt)
              end if
            end if
          
          case (18, 12, 20, 26, 27, 28, 29, 30, 31)                  ! Fresno Valley, only match by region
            if (aero_zones(i) == gzone) then
               cnt = cnt + 1
               maero412(cnt) = aero_sr412(i,season)
              maero470(cnt) = aero_sr470(i,season)
              maero650(cnt) = aero_sr650(i,3)        ! < -- always use summer for 650nm to match refsr650
              mbgaod(cnt)   = aero_bgaod(i,season)
              msiteindx(cnt) = i
              
              if (dflag) then
                print '(A,A,A,I4,I4)', trim(func_name), ', matching site: ', trim(aero_sites(i)), aero_zones(i), aero_types(i)
                print '(A,A,3(F11.6,1X))', trim(func_name), ', AERONET Baseline SR: ', maero412(cnt), maero470(cnt), maero650(cnt)
              end if
            end if
            
          case (15, 19)                ! Pune, only match by region, added Kanpur and India high elevation 31 January 2018 JLee
            if (aero_zones(i) == gzone) then
               cnt = cnt + 1
               maero412(cnt) = aero_sr412(i,season)
              maero470(cnt) = aero_sr470(i,season)
              maero650(cnt) = aero_sr650(i,3)        ! < -- always use summer for 650nm to match refsr650
              mbgaod(cnt)   = aero_bgaod(i,season)
              msiteindx(cnt) = i
              
              if (dflag) then
                print '(A,A,A,I4,I4)', trim(func_name), ', matching site: ', trim(aero_sites(i)), aero_zones(i), aero_types(i)
                print '(A,A,3(F11.6,1X))', trim(func_name), ', AERONET Baseline SR: ', maero412(cnt), maero470(cnt), maero650(cnt)
              end if
            end if
            
          case (1, 5, 10, 13)                ! N. Africa, Solar Villge (Saudi Arabia)
            if (aero_zones(i) == gzone .AND. aero_types(i) == lc_type) then
              cnt = cnt + 1
              maero412(cnt) = aero_sr412(i,season)
              maero470(cnt) = aero_sr470(i,season)
              maero650(cnt) = aero_sr650(i,3)        ! < -- always use summer for 650nm to match refsr650
              mbgaod(cnt)   = aero_bgaod(i,season)
              msiteindx(cnt) = i
              
              if (dflag) then
                print '(A,A,A,I4,I4)', trim(func_name), ', matching site: ', trim(aero_sites(i)), aero_zones(i), aero_types(i)
                print '(A,A,3(F11.6,1X))', trim(func_name), ', AERONET Baseline SR: ', maero412(cnt), maero470(cnt), maero650(cnt)
              end if
            end if
          
          case default ! -- everywhere else match by land cover type and geozone and elevation.
            if (aero_zones(i) == gzone .AND. aero_types(i) == lc_type .AND. (elev < 500 .EQV. aero_elev(i) < 500)) then
              cnt = cnt + 1
              maero412(cnt) = aero_sr412(i,season)
              maero470(cnt) = aero_sr470(i,season)
              maero650(cnt) = aero_sr650(i,3)        ! < -- always use summer for 650nm to match refsr650
              mbgaod(cnt)   = aero_bgaod(i,season)
              msiteindx(cnt) = i
              
              if (dflag) then
                print '(A,A,A,I4,I4)', trim(func_name), ', matching site: ', trim(aero_sites(i)), aero_zones(i), aero_types(i)
                print '(A,A,3(F11.6,1X))', trim(func_name), ', AERONET Baseline SR: ', maero412(cnt), maero470(cnt), maero650(cnt)
              end if
            end if
          end select
      end do
      
!     -- can we interpolate between zone's AERONET sites?
      call sortrx(m, maero650, sorted)
      if (refsr650 >= minval(maero650) .AND. refsr650 < maxval(maero650)) then

!       -- find where refsr650 fits in maero650() and interpolate between the two sites
        do i = 1, m-1
          if (refsr650 >= maero650(sorted(i)) .AND. refsr650 < maero650(sorted(i+1))) then
            ii = sorted(i)
            asite = aero_sites(msiteindx(ii))
            status = get_aeronet_brdf_sr(trim(asite), month, ra, sa, vza, ndvi, stdv, ab412,  &
            &                            ab470, ab650, ac412, ac470, ac650, use_alternate_brdf, debug=dflag)
            if (status /= 0) then
              print *, "ERROR: Failed to get BRDF-corrected SR from AERONET site: ", trim(asite), status
              return
            end if
            
            xnorm_412_f1 = ab412 / maero412(ii)
            xnorm_470_f1 = ab470 / maero470(ii)
            
            jj = sorted(i+1)
            asite = aero_sites(msiteindx(jj))
            status = get_aeronet_brdf_sr(trim(asite), month, ra, sa, vza, ndvi, stdv, ab412,  &
            &                            ab470, ab650, ac412, ac470, ac650, use_alternate_brdf, debug=dflag)
            if (status /= 0) then
              print *, "ERROR: Failed to get BRDF-corrected SR from AERONET site: ", trim(asite), status
              return
            end if
            
            xnorm_412_f2 = ab412 / maero412(jj)
            xnorm_470_f2 = ab470 / maero470(jj)
          
!           -- calculate AERONET site weights according to 650 values and adjust table refs.
            frac = (refsr650-maero650(ii)) / (maero650(jj)-maero650(ii))
            
            normfrac = frac*xnorm_412_f2 + (1.0-frac)*xnorm_412_f1
            sr412 = get_viirs_LER412(ilat, ilon)
            if (sr412 < -900.0) then
              status = -2
              return
            end if
            sr412 = sr412 * normfrac
            
            normfrac = frac*xnorm_470_f2 + (1.0-frac)*xnorm_470_f1
            sr470 = get_viirs_LER488(ilat, ilon)
            if (sr470 < -900.0) then
              status = -2
              return
            end if
            sr470 = sr470 * normfrac
            
            status = 0
            
            if (dflag) then
              print '(A,A,F11.6)', trim(func_name),", Pixel Ref 650 SR: ", refsr650
              print '(A,A,2(F11.6,1X))', trim(func_name),", Pixel Baseline SR, 412, 470: ", &
              & get_viirs_LER412(ilat, ilon), get_viirs_LER488(ilat, ilon)
              print '(4(A,1X))', trim(func_name),", interp sites: ", &
              & trim(aero_sites(msiteindx(ii))), trim(aero_sites(msiteindx(jj)))
              print '(A,A,3(F11.6,1X))', trim(func_name),", calcsr412, aerosr412: ", sr412, maero412(ii), maero412(jj)
              print '(A,A,3(F11.6,1X))', trim(func_name),", calcsr470, aerosr470: ", sr470, maero470(ii), maero470(jj)
              print '(A,A,3(F11.6,1X))', trim(func_name),", calcsr650, aerosr650: ", sr650, maero650(ii), maero650(jj)
              print '(A,A,2(F11.6,1X))', trim(func_name),", xnorm412: ", xnorm_412_f1, xnorm_412_f2
              print '(A,A,2(F11.6,1X))', trim(func_name),", xnorm470: ", xnorm_470_f1, xnorm_470_f2
            end if
            
            exit     !  jump out of loop, we're done!
          
          end if
        end do
                
!     -- no interpolation, use single AERONET site.
      else
        if (refsr650 <= minval(maero650)) then
          ii = sorted(1)  ! AERONET site w/ min. sr650 value
        else
          ii = sorted(m)  ! AERONET site w/ max. sr650 value
        end if
  
        asite = aero_sites(msiteindx(ii))
        status = get_aeronet_brdf_sr(trim(asite), month, ra, sa, vza, ndvi, stdv, ab412, ab470, &
        &                             ab650, ac412, ac470, ac650, use_alternate_brdf, debug=dflag)
        if (status /= 0) then
          print *, "ERROR: Failed ot get BRDF-corrected SR from AERONET site, single: ", trim(asite), status
          return
        end if

        ! correction to baseline adjustment factor for background aod, 10 January 2018 JLee
!        if ((gzone == 15) .and. bgaod > 0.0 .and. mbgaod(ii) > 0.0 &
!        &  .and. bgaod < mbgaod(ii)) then
!          if (ndvi < 0.3 .and. elev > 450) then
!!            select case (season)
!!              case (1)
!!                aod_corr_factor = min(0.1+(elev-300)*0.0015*amf,1.0)
!!              case (2)
!!                aod_corr_factor = min(0.1+(elev-300)*0.0010,0.5)
!!              case default
!!                aod_corr_factor = min(0.1+(elev-300)*0.0010,0.5)
!!            end select
!
!            select case (season)
!              case (1)
!                aod_corr_factor = 0.0
!              case (2)
!                aod_corr_factor = 0.0
!              case (3)
!                aod_corr_factor = 0.0
!              case (4)
!                aod_corr_factor = 0.0
!              case default
!                aod_corr_factor = 0.0
!            end select
!            xnorm_412_f1 = ab412 / (ac412+(maero412(ii)-ac412)*(aod_corr_factor*bgaod/mbgaod(ii)+(1.0-aod_corr_factor)))
!          else
!            xnorm_412_f1 = ab412 / (ac412+(maero412(ii)-ac412)*(0.1*bgaod/mbgaod(ii)+0.9))
!          end if
!        else
          xnorm_412_f1 = ab412 / maero412(ii)
!        end if

        sr412 = get_viirs_LER412(ilat, ilon)
        if (sr412 < -900.0) then
          status = -2
          return
        end if
        sr412 = sr412 * xnorm_412_f1
       
        ! correction to baseline adjustment factor for background aod, 10 January 2018 JLee
!        if ((gzone == 15) .and. bgaod > 0.0 .and. mbgaod(ii) > 0.0 &
!        &  .and. bgaod < mbgaod(ii)) then
!          if (ndvi < 0.3 .and. elev > 450) then
!!            select case (season)
!!              case (1)
!!                aod_corr_factor = min(0.1+(elev-300)*0.0015*amf,1.0)
!!              case(2)
!!                aod_corr_factor = min(0.1+(elev-300)*0.0010,0.5)
!!              case default
!!                aod_corr_factor = min(0.1+(elev-300)*0.0010,0.5)
!!            end select
!
!            select case (season)
!              case (1)
!                aod_corr_factor = 0.0
!              case (2)
!                aod_corr_factor = 0.0
!              case (3)
!                aod_corr_factor = 0.0
!              case (4)
!                aod_corr_factor = 0.0
!              case default
!                aod_corr_factor = 0.0
!            end select
!            xnorm_470_f1 = ab470 / (ac470+(maero470(ii)-ac470)*(aod_corr_factor*bgaod/mbgaod(ii)+(1.0-aod_corr_factor)))
!          else
!            xnorm_470_f1 = ab470 / (ac470+(maero470(ii)-ac470)*(0.1*bgaod/mbgaod(ii)+0.9))
!          end if
!        else
          xnorm_470_f1 = ab470 / maero470(ii)
!        end if

        sr470 = get_viirs_LER488(ilat, ilon)
        if (sr470 < -900.0) then
          status = -2
          return
        end if
        sr470 = sr470 * xnorm_470_f1
        
        status = 0
        
        if (dflag) then
          print '(A,A,F11.6)', trim(func_name), ", Pixel Ref 650 SR: ", refsr650
          print '(A,A,2(F11.6,1X))', trim(func_name), ", Pixel Baseline SR, 412, 470: ", &
          &  get_viirs_LER412(ilat, ilon),get_viirs_LER488(ilat, ilon)
          print '(3(A,1X))', trim(func_name), ", interp site: ", trim(aero_sites(msiteindx(ii)))
          print '(A,A,2(F11.6,1X))', trim(func_name), ", calcsr412, aerosr412: ", sr412, maero412(ii)
          print '(A,A,2(F11.6,1X))', trim(func_name), ", calcsr470, aerosr470: ", sr470, maero470(ii)
          print '(A,A,2(F11.6,1X))', trim(func_name), ", calcsr650, aerosr650 : ", sr650, maero650(ii)
          print '(A,A,F11.6)', trim(func_name), ", xnorm412: ", xnorm_412_f1
          print '(A,A,F11.6)', trim(func_name), ", xnorm470: ", xnorm_470_f1
          print '(A,A,2(F11.6,1X))', trim(func_name), ", aerobrdf412, 470: ", ab412, ab470
        end if
            
      end if
      
!   -- no AERONET sites (m==0), use table values.
    else
    
    sr412 = -999.0
    sr470 = -999.0
    sr650 = -999.0
     sr412 = get_LER412(ilat, ilon, ndvi, sa, ra, min_flag)    !get_viirs_modisbrdf_LER* ->get_LER*
     sr470 = get_LER470(ilat, ilon, ndvi, sa, ra, min_flag)    !7.5.2017 W.Kim
     sr650 = get_LER650(ilat, ilon, ndvi, sa, ra, min_flag)


         if (dflag) then
!          print '(A,A,14(F10.4))', trim(func_name), "lat, lon: ", lat, lon
!          print '(A,A,3(F10.4))', trim(func_name), "MBSR412, MBSR470, MBSR650: ", mb_sr412, mb_sr470, mb_sr650
!          print '(A,A,3(F10.4))', trim(func_name), "MSR412, MSR470, MSR650: ", get_modis_LER412(ilat,ilon), &
!          & get_modis_LER470(ilat,ilon), get_modis_LER650(ilat,ilon)
!          print '(A,A,3(F10.4))', trim(func_name), "VSR412, VSR488, VSR670: ", get_viirs_LER412(ilat,ilon), &
!          & get_viirs_LER488(ilat,ilon), get_viirs_LER670(ilat,ilon)
           print '(A,A,3(F10.4))', trim(func_name), "final SR412, SR470, SR650: ", sr412, sr470, sr650
         end if

         status = 1
    end if
    
!   -- final check
    if (sr412 < 0.0 .OR. sr470 < 0.0) then
      status = -1
    end if
      
    if (dflag) then
      print *, trim(func_name), ", table-based SR, 412, 470, 650: ", sr412, sr470, sr650
      print *, trim(func_name), ", final status: ", status
    end if
    
    return
          
  end function get_brdfcorr_sr

! -- returns surface reflectivity based on AERONET site (aero_site) BRDF and input
! --  geometry.  Helper function for get_bdrfcorr_sr(). Returns -1 on failure,
! --  otherwise 0.
  integer function get_aeronet_brdf_sr(aero_site, month, raa, sca, vza, ndvi, stdv, s412, s470,           &
  &                                     s650, c412, c470, c650, use_alternate_brdf, debug) result(status)

    implicit none
    
    character(len=50), parameter        ::  func_name = "get_aeronet_brdf_sr"
    
    integer, parameter                  ::  ndegs = 4
    
    character(len=*), intent(in)        ::  aero_site
    integer, intent(in)                 ::  month
    real, intent(in)                    ::  raa
    real, intent(in)                    ::  sca
    real, intent(in)                    ::  vza
    real, intent(in)                    ::  ndvi
    real, intent(in)                    ::  stdv
    real, intent(inout)                 ::  s412
    real, intent(inout)                 ::  s470
    real, intent(inout)                 ::  s650
    real, intent(inout)                 ::  c412
    real, intent(inout)                 ::  c470
    real, intent(inout)                 ::  c650
!    real, intent(inout)                 ::  ssa412
!    real, intent(inout)                 ::  ssa470
!    real, intent(inout)                 ::  ssa650
    logical, intent(in), optional       ::  use_alternate_brdf
    logical, intent(in), optional       ::  debug

    real, dimension(ndegs)              ::  co412
    real, dimension(ndegs)              ::  co470
    real, dimension(ndegs)              ::  co650
    
    integer                             ::  season
    
    logical                             :: fwd_scat
    logical                             :: dflag
    
    dflag = .false.
    if (present(debug)) then
      dflag = debug
    end if
    
    fwd_scat = .false.
    if (raa < 90.0) then
      fwd_scat = .true.
    end if
    
    s412 = -999.0
    s470 = -999.0
    s650 = -999.0

    c412 = -999.0
    c470 = -999.0
    c650 = -999.0
 
!   -- get season from month
    select case (month)
      case (12,1,2)
        season = 1
      case(3:5)
        season = 2
      case (6:8)
        season = 3
      case (9:11)
        season = 4
      case default
        print *, "ERROR: Invalid month specified: ", month
        status = -1
        return
    end select
    
    co412 = (/-999.0,-999.0,-999.0,-999.0/)
    co470 = (/-999.0,-999.0,-999.0,-999.0/)
    co650 = (/-999.0,-999.0,-999.0,-999.0/)
            
    select case (aero_site)
!     ------------------------------------
      case ("Banizoumbou")
        select case (season)
          case (1)
       if(ndvi >= 0.15) then
              co470 = (/1.02169058e01, 4.243027e-02, 1.54773501e-04, 0.0/)
              co412 = (/6.56567239, 1.46437509e-02, 0.0, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
       else
              co412 = (/5.05991244, 4.90682739e-02,0.0,0.0/)
              co470 = (/9.33658552, 7.10523577e-02, -2.60069034e-05,0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
          case (2)
            if (ndvi >= 0.12) then
              co412 = (/6.39744066, 4.00276042e-02, 0.0, 0.0/)
              co470 = (/1.04021434e01, 5.37401649e-02, 1.11642009e-04, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)

!             -- allow AERONET AOT<=1.0 rather than 0.5
!              co412 = (/6.35012418, 2.40388121e-02, 1.99004385e-04, 0.0/)
!              co470 = (/1.03240653e01, 4.68872309e-02, 1.41180318e-04, 0.0/)
!              co650 = (/0.0, 0.0, 0.0, 0.0/)
              
              
            else
              co412 = (/7.4842753, 0.0,0.0,0.0/)
              co470 = (/1.11423275e01, 8.99414273e-02, -7.65091516e-04, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
              
            ! -- all NDVI fits for NDVI < 0.12
              co412 = (/6.36294769, 5.40543846e-02, -1.71327907e-04, 0.0/)
              co470 = (/1.06517369e01, 6.01331952e-02, 0.0, 0.0/)
              
!             -- allow AERONET AOT<=1.0 rather than 0.5
!              co412 = (/5.72731182, 3.27134511e-02, 2.74284048e-04, 0.0/)
!              co470 = (/1.04276182e01, 8.12197464e-02, -4.61564311e-04, 0.0/)
!              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
            
          case (3)
!            if (ndvi >= 0.24) then
!              co412 = (/4.71101627, 6.64118803e-02, -5.38844041e-04, 0.0/)
!              co470 = (/7.65812386, 7.83366273e-02, -2.4482234e-04, 0.0/)
!              co650 = (/0.0, 0.0, 0.0, 0.0/)
!            else
!              co412 = (/4.83957621, 9.91643198e-02, -4.31830448e-04, 0.0/)
!              co470 = (/9.52056332, 1.18515428e-01, -8.79903107e-04, 0.0/)
!              co650 = (/0.0, 0.0, 0.0, 0.0/)

            if (ndvi >= 0.24) then
              co412 = (/4.68590202, 4.79492301e-02, 0.0, 0.0/)
              co470 = (/7.65714295, 7.81415694e-02, -2.40306130e-04, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/4.78835211, 8.16497137e-02, -3.62862493e-04, 0.0/)
              co470 = (/9.02297785, 1.21271052e-01, -1.09555283e-03, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)

              
!             -- test new BRDF for summer, Banizoumbou. Trying to limit low bias in june...
              co412 = (/4.33902797, 8.37586808e-02, -5.54927180e-04, 0.0/)
              co470 = (/9.25667223, 1.14446764e-01, -1.03215095e-03, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
              
            end if
            
          case (4)
            if (ndvi >= 0.21) then
              co412 = (/4.78430658, 3.2201042e-02, -1.62743787e-04, 6.68624319e-06/)
              co470 = (/7.8706363, 5.51922546e-02, 0.0, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/5.239379, 4.142642525e-02, 1.66142711e-05, 0.0/)
              co470 = (/9.09596194, 6.99428916e-02, -2.528929176e-04, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
          case default
            print *, "ERROR: Invalid season specified: ", season
            status = -1
            return
        end select
        
!     ------------------------------------
      case ("Kanpur")
        select case (season)
          case (1)
            co412 = (/4.5439484, 7.12450582e-2, 5.1565853e-4, 0.0/)
            co470 = (/5.1905870, 7.0479636e-2, -7.2117339e-4, 0.0/)
            co650 = (/0.0, 0.0, 0.0, 0.0/)
          case (2)
            co412 = (/4.0491380, 6.5488864e-2, -2.7752629e-4, 0.0/)
            co470 = (/5.5883717, 7.6432090e-2,  3.1404959e-5, 0.0/)
            co650 = (/0.0, 0.0, 0.0, 0.0/)
          case (3)
!           use fall BRDF, 26 January JLee, TEST
            co412 = (/4.7996805, 6.3598622e-2, -8.9915671e-6, 0.0/)
            co470 = (/5.6717516, 6.6228202e-2, 3.5839652e-4  , 0.0/)
            co650 = (/0.0, 0.0, 0.0, 0.0/)
!            original
!            co412 = (/4.76004180, -5.53951481e-02, 2.54996532e-03, 0.0/)
!            co470 = (/8.72009627, -3.85232353e-02, 2.63096171e-03, 0.0/)
!            co650 = (/0.0, 0.0, 0.0, 0.0/)
          case (4)
            co412 = (/4.7996805, 6.3598622e-2, -8.9915671e-6, 0.0/)
            co470 = (/5.6717516, 6.6228202e-2, 3.5839652e-4  , 0.0/)
            co650 = (/0.0, 0.0, 0.0, 0.0/)
          case default
            print *, "ERROR: Invalid season specified: ", season
            status = -1
            return
        end select
         
!     ------------------------------------
      case ("IER_Cinzana")
        select case (season)
          case (1)
            if (ndvi >= 0.2) then
         co412 = (/3.63552866, 3.80334634e-2, 0.0, 0.0/)
         co470 = (/7.1468792, 4.90084709e-2, -5.8487537e-4, 2.65330779e-5/)
         co650 = (/0.0, 0.0, 0.0, 0.0/)
       else
         co412 = (/3.50830056, 4.45759847e-2, 0.0, 0.0/)
         co470 = (/7.2185198, 7.14231475e-2, 0.0, 0.0/)
         co650 = (/0.0, 0.0, 0.0, 0.0/)
       end if
          case (2)
            if(ndvi >= 0.18) then
         co412 = (/8.03123187, 6.31014685e-2, -3.03999188e-4, 0.0/)
         co470 = (/8.00850678, 5.36508158e-2, 0.0, 0.0/)
         co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/4.78540468, 4.53379041e-2, 0.0, 0.0/)
              co470 = (/8.55690294, 6.26426162e-2, 0.0, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
          case (3)
            if (ndvi >= 0.3) then
              co412 = (/3.56969498, 1.48750348e-2, 0.0, 0.0/)
              co470 = (/5.218741, 3.68430117e-2, 0.0, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/4.1364437, 3.54993284e-2, 0.0, 0.0/)
              co470 = (/7.68177478, 5.48823787e-2, 0.0, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
            
          case (4)
           if (ndvi < 0.36) then
              co412 = (/3.3600304, 4.07091256e-2, -6.36564712e-4, 0.0/)
              co470 = (/6.00156726, 5.32555429e-2, -2.03417949e-4, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else if(ndvi < 0.48) then
              co412 = (/2.75902271, 5.86471497e-3, 8.9298557e-4, -2.84768722e-6/)
              co470 = (/4.49584493, 3.2333347e-2, 1.68247099e-3, -2.08146852e-5/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/2.54155592, 0.0, 0.0, 0.0/)
              co470 = (/4.08158909, 7.65889923e-3, 0.0, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
          case default
            print *, "ERROR: Invalid season specified: ", season
            status = -1
            return
        end select
        
!     ------------------------------------
      case ("Zinder_Airport")
        select case (season)
          case (1)
            if (ndvi >= 0.15) then
              co412 = (/5.74980283, 2.39050366e-02, 0.0, 0.0/)
              co470 = (/1.00530139e01, 5.17129972e-02, 0.0, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/5.86499038, 4.16328036e-02, 0.0, 0.0/)
              co470 = (/1.07113413e01, 6.94885771e-02, -2.12100636e-04, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
          case (2)
            co412 = (/5.48605593, 7.32715822e-02, -3.94126361e-04, 0.0/)
            co470 = (/1.02289147e01, 7.179364421e-02, 0.0, 0.0/)
            co650 = (/0.0, 0.0, 0.0, 0.0/)
          case (3)
            if (ndvi >= 0.24) then
              co412 = (/3.72175525, 1.33333812e-03, 0.0, 0.0/)
              co470 = (/6.28406012, 6.05189262e-02, 0.0, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
 
            else if (ndvi >= 0.15) then
              co412 = (/9.14958326, -6.22526837e-03, 0.0, 0.0/)
              co470 = (/1.22197489e01, 3.78885944e-02, 0.0, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
!              if (use_alternate_brdf) then
!                co412 = (/6.00999178, -2.42563433e-03, 0.0, 0.0 /)
!                co470 = (/9.32222748, 0.0, 0.0, 0.0/)
!              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/5.21423049, 8.21993929e-02, -3.75107523e-04, 0.0/)
              co470 = (/9.78066195, 8.69504404e-02, -3.28108722e-04, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
            
          case (4)
            if (ndvi >= 0.24) then
              co412 = (/2.07354907, 5.40878937e-02, 0.0, 0.0/)
              co470 = (/5.31006297, 8.37948079e-02, 0.0, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/5.05835071, 1.82023167e-02, -6.29515824e-04, 2.93686368e-05/)
              co470 = (/9.26108740, 6.81130675e-02, 0.0, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
          case default
            print *, "ERROR: Invalid season specified: ", season
            status = -1
            return
        end select
        
!     ------------------------------------
      case ("Beijing")
        select case (season)
          case (1)
            co412 = (/6.3014686,  4.1543979e-2, 1.1858985e-4, 0.0/)
            co470 = (/7.4044324,  8.1265848e-2, 4.9485414e-4, 0.0/)
            co650 = (/0.0, 0.0, 0.0, 0.0/)
          case (2)
            if (ndvi < 0.18) then
              co412 = (/5.89152474, 4.52651696e-2, 2.68647127e-4, 0.0/)
              co470 = (/7.30775670, 7.62461937e-2, -3.21229444e-4, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/5.30165348, 2.39439950e-2, 7.06604780e-4, 0.0/)
              co470 = (/6.28190845, 6.15060495e-2, 1.20706642e-4, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
          case (3)
            co412 = (/5.6957808, 5.3492403e-2, -2.1731312e-4, 0.0/)
            co470 = (/6.1378929, 7.4536939e-2, -3.7338370e-4, 0.0/)
            co650 = (/0.0, 0.0, 0.0, 0.0/)
          case (4)
            co412 = (/5.4468575, 4.9992086e-2, 1.0863964e-05, 0.0/)
            co470 = (/6.1302000, 6.7842888e-2, 6.9728565e-05, 0.0/)
            co650 = (/0.0, 0.0, 0.0, 0.0/)
          case default
            print *, "ERROR: Invalid season specified: ", season
            status = -1
            return
        end select

!     ------------------------------------
      case ("Hamim")
        select case (season)
          case (1)
            co412 = (/10.226094,     -0.14939743,    0.0064489043,  -7.3014348e-05/)
            co470 = (/14.875093,     0.046262226,   -0.0010471590,   1.1772507e-05/)
            co650 = (/34.573460,      0.39409698,    -0.017102505,   0.00021277064/)
 !           ssa412 = 0.94
 !           ssa470 = 0.96
 !           ssa650 = 0.995
          case (2)
            co412 = (/7.1646428,      0.15052042,   -0.0043741791,   4.4118414e-05/)
            co470 = (/12.347033,      0.22954843,   -0.0069846621,   7.1241796e-05/)
            co650 = (/34.616283,     0.042627358,   -0.0013109076,   2.2689966e-05/)
!            ssa412 = 0.94
!            ssa470 = 0.96
!            ssa650 = 0.995
          case (3)
            co412 = (/12.806989,     -0.46510965,     0.016412267,  -0.00016780296/)
            co470 = (/19.796796,     -0.62064603,     0.021296302,  -0.00021068886/)
            co650 = (/32.705724,      0.12309721,   -0.0029126864,   3.7611775e-05/)
!            ssa412 = 0.94
!            ssa470 = 0.96
!            ssa650 = 0.995
          case (4)
            co412 = (/8.4939626,     0.062900115,   -0.0016468367,   1.5106109e-05/)
            co470 = (/15.206425,     0.020553284,   -0.0013703752,   2.2283588e-05/)
            co650 = (/34.954706,      0.17635488,   -0.0058044940,   5.6107481e-05/)
!            ssa412 = 0.94
!            ssa470 = 0.96
!            ssa650 = 0.995
          case default
            print *, "ERROR: Invalid season specified: ", season
            status = -1
            return
        end select

!     ------------------------------------
      case ("Fresno_2")
        select case (season)
          case (1)
            co412 = (/5.80751073, 4.06728637e-02, 4.24590213e-05, 0.0/)
            co470 = (/7.20443297, 7.07862547e-02, 5.50556530e-04, 0.0/)
            co650 = (/0.0, 0.0, 0.0, 0.0/)
          case (2)
!            if (ndvi < 0.28) then
!              co412 = (/6.40431062, 3.72488428e-02, -1.66166477e-03, 3.28633059e-05/)
!              co470 = (/7.69437062, 6.82783655e-02, -6.49639866e-04, 7.40504460e-06/)
!              co650 = (/0.0, 0.0, 0.0, 0.0/)
!            else
              co412 = (/5.10554282, 3.12540398e-02, 3.70094580e-04, 0.0/)
              co470 = (/6.98775847, 6.37364825e-02, -5.27383741e-05, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
!            end if
          case (3)
!            if (ndvi < 0.28) then
!              co412 = (/5.58406544, 1.61289652e-02, 6.22273421e-04, 0.0/)
!              co470 = (/7.11813056, 4.77749714e-02, 4.23934048e-04, 0.0/)
!              co650 = (/0.0, 0.0, 0.0, 0.0/)
!            else
              co412 = (/4.76062810, 4.07850599e-02, 0.0, 0.0/)
              co470 = (/6.85966085, 5.99596507e-02, -8.44790522e-04, 2.47006074e-05/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
!            end if
          case (4)
!            if (ndvi < 0.27) then
!              co412 = (/4.63756258, 4.98512265e-02, 0.0, 0.0/)
!              co470 = (/6.85551799, 5.0385959e-02, -1.39476331e-04, 2.59945517e-05/)
!              co650 = (/0.0, 0.0, 0.0, 0.0/)
!            else
              co412 = (/4.53821328, 2.89784912e-02, 6.46309542e-04, 0.0/)
              co470 = (/6.80695142, 5.35553957e-02, 1.21293340e-04, 5.95530797e-06/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
!            end if
          case default
            print *, "ERROR: Invalid season specified: ", season
            status = -1
            return
        end select
        
      !     ------------------------------------
      case ("Fresno_GZ18")
        select case (season)
          case (1)
            co412 = (/5.80751073, 4.06728637e-02, 4.24590213e-05, 0.0/)
            co470 = (/7.20443297, 7.07862547e-02, 5.50556530e-04, 0.0/)
            co650 = (/0.0, 0.0, 0.0, 0.0/)
          case (2)
            if (ndvi < 0.28) then
              co412 = (/6.40431062, 3.72488428e-02, -1.66166477e-03, 3.28633059e-05/)
              co470 = (/7.69437062, 6.82783655e-02, -6.49639866e-04, 7.40504460e-06/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/5.10554282, 3.12540398e-02, 3.70094580e-04, 0.0/)
              co470 = (/6.98775847, 6.37364825e-02, -5.27383741e-05, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
          case (3)
            if (ndvi < 0.28) then
              co412 = (/5.58406544, 1.61289652e-02, 6.22273421e-04, 0.0/)
              co470 = (/7.11813056, 4.77749714e-02, 4.23934048e-04, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/4.76062810, 4.07850599e-02, 0.0, 0.0/)
              co470 = (/6.85966085, 5.99596507e-02, -8.44790522e-04, 2.47006074e-05/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
          case (4)
            if (ndvi < 0.27) then
              co412 = (/4.63756258, 4.98512265e-02, 0.0, 0.0/)
              co470 = (/6.85551799, 5.0385959e-02, -1.39476331e-04, 2.59945517e-05/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/4.53821328, 2.89784912e-02, 6.46309542e-04, 0.0/)
              co470 = (/6.80695142, 5.35553957e-02, 1.21293340e-04, 5.95530797e-06/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
          case default
            print *, "ERROR: Invalid season specified: ", season
            status = -1
            return
        end select

!     ------------------------------------
      case ("CCNY")
        select case (season)
          case (1)
             co412 = (/6.9337283, 4.4282749e-02, 1.4061984e-03, 0.0/)
             co470 = (/7.3485598, 6.9214518e-02, 1.8788394e-03, 0.0/)
             co650 = (/0.0, 0.0, 0.0, 0.0/)
          case (2)
             co412 = (/6.5373929, 2.6971498e-02, -7.8218577e-05, 0.0/)
             co470 = (/7.2815217, 4.3698874e-02, -2.7398479e-04, 0.0/)
             co650 = (/0.0, 0.0, 0.0, 0.0/)
          case (3)
            co412 = (/6.6297657, 3.1107265e-02, -1.2153919e-04, 0.0/)
            co470 = (/6.8517577, 5.2723244e-02, -4.9795072e-04, 0.0/)
            co650 = (/0.0, 0.0, 0.0, 0.0/)
          case (4)
            co412 = (/6.1744159, 1.5031182e-02, 0.0, 0.0/)
            co470 = (/6.7001040, 2.5827169e-02, 0.0, 0.0/)
            co650 = (/0.0, 0.0, 0.0, 0.0/)
          case default
            print *, "ERROR: Invalid season specified: ", season
            status = -1
            return
        end select
        
!     ------------------------------------
      case ("Agoufou")
        if (fwd_scat) then                      ! forward scattering
          select case (season)
            case (1)
                     co412 = (/4.79872574, -1.56512429e-02, -1.23860221e-03, 5.84017625e-05/)
                     co470 = (/8.05072532, 1.61294620e-02, -6.57200682e-04, 4.10211916e-05/)
                     co650 = (/0.0, 0.0, 0.0, 0.0/)
            case (2)
              co412 = (/5.34279489, 2.40838594e-03, 0.0, 0.0/)
              co470 = (/8.28211141, 2.88113903e-02, 0.0, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            case (3)
              if (ndvi >= 0.24) then
                co412 = (/3.26738286, 4.36827818e-02, -1.55895303e-04, -1.20224162e-05/)
                co470 = (/5.03384478, 7.29893383e-02, 1.31678743e-03, -4.31578699e-05/)
                co650 = (/0.0, 0.0, 0.0, 0.0/)
              else
                co412 = (/4.49283427, 4.04848382e-05, -3.09039795e-04, 1.65524662e-05/)
                co470 = (/7.34716100, 4.10633137e-02, 2.51289909e-04, -3.32283482e-06/)
                co650 = (/0.0, 0.0, 0.0, 0.0/)
              end if
            case (4)
              if (ndvi >= 0.21) then
                co412 = (/3.21406439, 2.84017859e-02, 2.25979209e-04, -1.35954941e-05/)
                co470 = (/5.32086992, 5.10973584e-02, 1.23563583e-03, -3.14731061e-05/)
                co650 = (/0.0, 0.0, 0.0, 0.0/)
              else if (ndvi >= 0.18) then
                co412 = (/4.03496081, -1.16119226e-02, -4.21409772e-04, 3.26077680e-05/)
                co470 = (/6.89540010, 3.44825524e-02, 4.25223284e-04, 0.0/)
                co650 = (/0.0, 0.0, 0.0, 0.0/)
              else
               co412 = (/4.68551647, 3.37375535e-03, 3.28094666e-04, 0.0/)
                co470 = (/7.97076041, 4.07903024e-02, 0.0, 0.0/)
                co650 = (/0.0, 0.0, 0.0, 0.0/)
              end if
            case default
              print *, "ERROR: Invalid season specified: ", season
              status = -1
              return
          end select
        
        else                                      ! backward scattering
               select case (season)
                  case (1)
                     co412 = (/4.79872574, -1.56512429e-02, -1.23860221e-03, 5.84017625e-05/)
                     co470 = (/8.05072532, 1.61294620e-02, -6.57200682e-04, 4.10211916e-05/)
                     co650 = (/0.0, 0.0, 0.0, 0.0/)
                  case (2)
                     co412 = (/5.31125838, 1.14257083e-02, 0.0, 0.0/)
                     co470 = (/8.30402936, 4.22170099e-02, 0.0, 0.0/)
                     co650 = (/0.0, 0.0, 0.0, 0.0/)
                  case (3)
                     if (ndvi >= 0.24) then
                        co412 = (/3.26738286, 4.36827818e-02, -1.55895303e-04, -1.20224162e-05/)
                        co470 = (/5.03384478, 7.29893383e-02, 1.31678743e-03, -4.31578699e-05/)
                        co650 = (/0.0, 0.0, 0.0, 0.0/)
                     else
                        co412 = (/4.49283427, 4.04848382e-05, -3.09039795e-04, 1.65524662e-05/)
                        co470 = (/7.34716100, 4.10633137e-02, 2.51289909e-04, -3.32283482e-06/)
                        co650 = (/0.0, 0.0, 0.0, 0.0/)
                     end if
                  case (4)
                     if (ndvi >= 0.21) then
                        co412 = (/3.21406439, 2.84017859e-02, 2.25979209e-04, -1.35954941e-05/)
                        co470 = (/5.32086992, 5.10973584e-02, 1.23563583e-03, -3.14731061e-05/)
                        co650 = (/0.0, 0.0, 0.0, 0.0/)
                     else if (ndvi >= 0.18) then
                        co412 = (/4.03496081, -1.16119226e-02, -4.21409772e-04, 3.26077680e-05/)
                        co470 = (/6.89540010, 3.44825524e-02, 4.25223284e-04, 0.0/)
                        co650 = (/0.0, 0.0, 0.0, 0.0/)
                     else
                        co412 = (/4.68551647, 3.37375535e-03, 3.28094666e-04, 0.0/)
                        co470 = (/7.97076041, 4.07903024e-02, 0.0, 0.0/)
                        co650 = (/0.0, 0.0, 0.0, 0.0/)
                     end if
                  case default
                     print *, "ERROR: Invalid season specified: ", season
                     status = -1
                     return
               end select
         end if

!     ------------------------------------
      case ("Tinga_Tingana")
        if (fwd_scat) then                      ! forward scattering
        select case (season)
          case (1)
            if (ndvi >= 0.13) then
              co412 = (/9.20413537, -5.17086422e-03, 1.02106845e-04, 0.0/)
              co470 = (/1.22306875e01, 4.82832263e-02, 3.68569587e-04, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/9.43803716, -3.72278076e-03 , 0.0, 0.0/)
              co470 = (/1.27377047e01, 8.35333189e-03, 0.0, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
          case (2)
            if (ndvi >= 0.18) then
              co412 = (/8.92737977, 1.76851898e-02, 0.0, 0.0/)
              co470 = (/12.634865, 0.028769372, -0.00041455473, -3.2556186e-05/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/9.37232462, -9.65154143e-03, -1.74009834e-04, 0.0/)
              co470 = (/12.634865, 0.028769372, -0.00041455473, -3.2556186e-05/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
          case (3)
            if (ndvi >= 0.20) then
              co412 = (/7.81846282, 3.15351932e-03, 0.0, 0.0/)
              co470 = (/1.08279281e01, 1.57698665e-02, 0.0, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else if (ndvi >= 0.14) then
              co412 = (/8.88977803, 1.51799616e-03, -7.23417787e-05, 0.0/)
              co470 = (/1.20676913e01, 2.91883388e-02, 0.0, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/9.08868817, -4.17799194e-02, -8.09302376e-04, 0.0/)
              co470 = (/1.25503323e01, -3.95329882e-02, -1.77305946e-03, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
          case (4)
            if (ndvi >= 0.14) then
              co412 = (/8.87184595, -1.57754324e-02, 0.0, 0.0/)
              co470 = (/1.19334546e01, 8.68690821e-03, 0.0, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/9.09342651, -4.17208574e-02, 0.0, 0.0/)
              co470 = (/1.23033719e01, -5.36893150e-02, -8.53259219e-04, 7.49332728e-05/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
          case default
            print *, "ERROR: Invalid season specified: ", season
            status = -1
            return
        end select

      else                                      ! backward scattering
        select case (season)
          case (1)
            if (ndvi >= 0.12) then
              co412 = (/9.41470858, -2.15844867e-03, -1.82266839e-04, 1.21977424e-05/)
              co470 = (/1.2501821e01, 4.02391785e-02, 2.08696099e-04, 1.11738774e-06/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/1.02358733e01, 7.77620010e-02, -1.98193828e-03, 0.0/)
              co470 = (/1.36183155e01, 4.05909269e-02, 0.0, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
          case (2)
            if (ndvi >= 0.18) then
              co412 = (/7.98223150, 1.40752752e-02, 2.69501591e-03, 0.0/)
              co470 = (/1.25564034e01, 5.50899365e-02, 0.0, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/9.42624421, 1.27763170e-02, 3.65495933e-04, -1.96877371e-05/)
              co470 = (/1.30602664e01, 7.57287787e-02, 4.73549830e-04, -6.40350281e-05/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
          case (3)
            if (ndvi >= 0.20) then
              co412 = (/7.89996147, -8.14332506e-03, -7.83342101e-04, 0.0/)
              co470 = (/1.13849147e01, 3.69927333e-02, -3.88390040e-04, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else if (ndvi >= 0.14) then
              co412 = (/8.64034181, 1.67726763e-03, 5.88861610e-04, 0.0/)
              co470 = (/1.24428324e01, 3.71238324e-02, -6.76748814e-04, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
              
            else
              co412 = (/9.41753918, -2.67338094e-02, -8.88155712e-04, 0.0/)
              co470 = (/1.34852983e01, 1.12436978e-02, -1.61943826e-03, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
          case (4)
            if (ndvi >= 0.14) then
              co412 = (/8.83070860, -1.51364763e-03, 8.54079578e-04, 0.0/)
              co470 = (/1.20956002e01, 3.68362395e-02, 7.40851587e-04, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/8.88610258, -2.01943966e-02, 1.61198594e-03, 4.14828495e-06/)
              co470 = (/1.20856901e01, 3.26791012e-02, 2.03660386e-03, -1.72052223e-05/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
          case default
            print *, "ERROR: Invalid season specified: ", season
            status = -1
            return
        end select
      end if
!     ------------------------------------
      case ("Moldova")
        select case (season)
          case (1)
            co412 = (/6.08251454, 8.6508708e-02, 1.4521957e-03, 0.0/)
            co470 = (/5.80924919, 6.4047214e-02, 1.4495468e-03, 0.0/)
            co650 = (/0.0, 0.0, 0.0, 0.0/)
          case (2)
            if (ndvi >= 0.3) then
              co412 = (/4.59092762, 3.22065937e-2, 1.72682251e-4, 0.0/)
              co470 = (/5.18390073, 5.80489830e-2, -2.1862814e-4, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/5.34520098, 3.28317599e-2, 3.17897070e-4, 0.0/)
              co470 = (/6.08505785, 5.89311646e-2, 2.59308378e-4, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
             end if
          case (3)
            co412 = (/4.41373796, 2.74747266e-2, -9.0896914e-5, 0.0/)
            co470 = (/4.84420545, 4.30808241e-2, -5.4624034e-5, 0.0/)
            co650 = (/0.0, 0.0, 0.0, 0.0/)
          case (4)
            if (ndvi >= 0.3) then
              co412 = (/4.45954687, 2.77637428e-2, 7.73159082e-4, 0.0/)
              co470 = (/5.06433779, 5.27772907e-2, 7.80966442e-4, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/5.22111926, 9.79821719e-2, 2.74956928e-3, 0.0/)
              co470 = (/5.63105884, 1.20658024e-1, 2.76971745e-3, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
          case default
            print *, "ERROR: Invalid season specified: ", season
            status = -1
            return
        end select
        
!     ------------------------------------
      case ("Modena")
        select case (season)
          case (1)
            co412 = (/4.8756241, 6.2534569e-2, 1.3572766e-3, 0.0/)
            co470 = (/5.4655822, 9.2120653e-2, 2.2203512e-3, 0.0/)
            co650 = (/0.0, 0.0, 0.0, 0.0/)
          case (2)
            co412 = (/5.0263816, 2.5373434e-2, 3.9146226e-4, 0.0/)
            co470 = (/5.4623850, 3.9852901e-2, 2.4863178e-4, 0.0/)
            co650 = (/0.0, 0.0, 0.0, 0.0/)
          case (3)
            co412 = (/5.1737469, 3.5232032e-2, -2.7832108e-4, 0.0/)
            co470 = (/5.4415591, 5.7405000e-2, -3.7084290e-4, 0.0/)
            co650 = (/0.0, 0.0, 0.0, 0.0/)
          case (4)
            co412 = (/3.57525283, 4.76106748e-2, 0.0, 0.0/)
            co470 = (/5.07540158, 6.07470043e-2, 0.0, 0.0/)
            co650 = (/0.0, 0.0, 0.0, 0.0/)
          case default
            print *, "ERROR: Invalid season specified: ", season
            status = -1
            return
        end select
        
!     ------------------------------------
      case ("Ispra")
        select case (season)
          case (1)
            co412 = (/2.94840397, 7.50469276e-02, -8.55631487e-03, 2.27740054e-04/)
            co470 = (/3.84084962, 2.46123068e-02, -3.93096482e-03, 1.39839185e-04/)
            co650 = (/5.32222028, 4.87768133e-02, 4.22486516e-03, -1.01546887e-04/)
 !           ssa412 = 0.96
 !           ssa470 = 0.96
 !           ssa650 = 0.995
          case (2)
            if (ndvi < 0.42) then
              co412 = (/2.13608185, -1.16229203e-02, 2.94401536e-03, -4.50006125e-05/)
              co470 = (/2.49471720, 2.33437868e-02, 3.68220345e-03, -7.55117909e-05/)
              co650 = (/3.45073310, 2.58147011e-01, -6.23027605e-03, 5.75888737e-05/)
            else
              co412 = (/2.23138624, -4.08180960e-02, 2.65231326e-03, -3.66271249e-05/)
              co470 = (/2.53458698, -5.68549547e-03, 1.36468352e-03, -1.8128715e-05/)
              co650 = (/3.3522585, 1.75290271e-02, 3.74108282e-04, -4.21350245e-06/)
            end if
!            ssa412 = 0.96
!            ssa470 = 0.96
!            ssa650 = 0.995
          case (3)
            co412 = (/5.85428528, -4.40189848e-01, 1.40438945e-02, -1.30337719e-04/)
            co470 = (/5.56715596, -3.51677373e-01, 1.13826478e-02, -1.01625742e-04/)
            co650 = (/5.32153704, -2.47300042e-01, 8.34454443e-03, -6.98915803e-05/)
!            ssa412 = 0.96
!            ssa470 = 0.96
!            ssa650 = 0.995
          case (4)
            co412 = (/9.78051816e-01, 3.83718088e-03, 2.03033445e-03, -3.88452812e-05/)
            co470 = (/1.34514727, 2.32911980e-02, 2.02551784e-03, -3.87492873e-05/)
            co650 = (/2.47427232, 4.49367275e-02, 1.42347421e-03, -3.36743928e-05/)
!            ssa412 = 0.96
!            ssa470 = 0.96
!            ssa650 = 0.995
          case default
            print *, "ERROR: Invalid season specified: ", season
            status = -1
            return
        end select

!     ------------------------------------
      case ("Saada")
        if (fwd_scat) then                      ! forward scattering
          select case (season)
            case (1)
              if (ndvi >= 0.3) then
                co412 = (/2.64678353, 0.0, 0.0, 0.0/)
                co470 = (/4.30869805, 2.23352643e-02, -3.08512532e-05, 0.0/)
                co650 = (/0.0, 0.0, 0.0, 0.0/)
              else
                co412 = (/2.89503206, -1.10265858e-02, 8.83157740e-05, 4.86510471e-05/)
                co470 = (/4.85292481, 1.80204315e-02, -7.41860898e-05, 4.92511685e-05/)
                co650 = (/0.0, 0.0, 0.0, 0.0/)
              end if
            case (2)
              if (ndvi < 0.32) then
                co412 = (/4.18631389, 2.5289862e-02, -9.62343375e-04, 3.39909450e-05/)
                co470 = (/5.32132814, 5.63754659e-02, -9.39866104e-04, 2.88986833e-05/)
                co650 = (/0.0, 0.0, 0.0, 0.0/)
              else if (ndvi >= 0.32 .AND. ndvi < 0.36) then
                co412 = (/4.06688221, -1.77611644e-02, 0.0, 0.0/)
                co470 = (/4.97204509, 8.48689442e-03, 7.25410501e-04, 0.0/)
                co650 = (/0.0, 0.0, 0.0, 0.0/)
              else
                co412 = (/3.42376327, -1.99413791e-02, 0.0, 0.0/)
                co470 = (/4.27929238, 2.86227327e-02, 0.0, 0.0/)
                co650 = (/0.0, 0.0, 0.0, 0.0/)
              end if
            case (3)
              if (ndvi < 0.27) then
                co412 = (/4.37627766, 8.01305335e-03, 8.12666253e-04, -1.00058657e-05/)
                co470 = (/5.52675332, 4.58968534e-02, 5.09676314e-04, -9.71283147e-06/)
                co650 = (/0.0, 0.0, 0.0, 0.0/)
              else
                co412 = (/4.31991803, -4.05160123e-02, 2.50853131e-03, -4.14415804e-05/)
                co470 = (/5.20304491, 3.31189616e-02, 2.11258050e-04, 0.0/)
                co650 = (/0.0, 0.0, 0.0, 0.0/)
              end if
            case (4)
              if (ndvi < 0.27) then
                co412 = (/3.90506839, 0.0, 0.0, 0.0/)
                co470 = (/4.55176618, 0.0, 0.0, 0.0/)
                co650 = (/0.0, 0.0, 0.0, 0.0/)
              else
                co412 = (/3.71875289, 1.75827423e-03, 0.0, 0.0/)
                co470 = (/5.24692412, 4.20630300e-02, 0.0, 0.0/)
                co650 = (/0.0, 0.0, 0.0, 0.0/)
              end if
            case default
              print *, "ERROR: Invalid season specified: ", season
              status = -1
              return
          end select
        else                                      ! backward scattering
          select case (season)
            case (1)
              if (ndvi >= 0.3) then
                co412 = (/2.40043582, -3.33151038e-02, 9.09508864e-04, 4.12504587e-05/)
                co470 = (/4.25935855, 9.65676859e-03, 2.92691096e-04, 2.94546746e-05/)
                co650 = (/0.0, 0.0, 0.0, 0.0/)
              else
                co412 = (/2.79446907, 0.0, 0.0, 0.0/)
                co470 = (/4.85292481, 1.80204315e-02, -7.41860898e-05, 4.92511685e-05/)
                co650 = (/0.0, 0.0, 0.0, 0.0/)
              end if
            case (2)
              if (ndvi < 0.32) then
                co412 = (/3.77066492, 1.58922076e-02, 3.03098949e-04, 0.0/)
                co470 = (/5.03028472, 5.31452578e-02, 1.26507640e-04, 0.0/)
                co650 = (/0.0, 0.0, 0.0, 0.0/)
              else if (ndvi >= 0.32 .AND. ndvi < 0.36) then
                co412 = (/4.15714313, -1.65945058e-02, -3.80237318e-04, 1.50147243e-05/)
                co470 = (/5.00549166, 1.15932385e-02, 6.83756479e-04, 0.0/)
                co650 = (/0.0, 0.0, 0.0, 0.0/)
              else
                co412 = (/3.09370523, -3.88747605e-02, 1.03242387e-03, 0.0/)
                co470 = (/4.35868847, 1.25746607e-02, -9.70954354e-04, 4.01565514e-05/)
                co650 = (/0.0, 0.0, 0.0, 0.0/)
              end if
            case (3)
              if (ndvi < 0.27) then
                co412 = (/4.42122008, -1.43446348e-03, -1.26690140e-04, 1.39839144e-05/)
                co470 = (/5.56874754, 3.92057334e-02, -2.26222735e-04, 1.30731833e-05/)
                co650 = (/0.0, 0.0, 0.0, 0.0/)
              else
                co412 = (/4.45603717, -3.60020814e-02, 8.76708868e-04, 0.0/)
                co470 = (/5.16502610, 2.56647090e-02, 6.43566423e-04, 0.0/)
                co650 = (/0.0, 0.0, 0.0, 0.0/)
              end if
            case (4)
              if (ndvi < 0.27) then
                co412 = (/3.16438602, 6.46839077e-02, 0.0, 0.0/)
                co470 = (/4.63876927, 9.09174541e-02, 0.0, 0.0/)
                co650 = (/0.0, 0.0, 0.0, 0.0/)
              else
                co412 = (/3.45291249, -2.23746051e-03, 6.18819075e-04, 0.0/)
                co470 = (/5.09438639, 3.96335945e-02, 3.59544361e-04, 0.0/)
                co650 = (/0.0, 0.0, 0.0, 0.0/)
              end if
            case default
              print *, "ERROR: Invalid season specified: ", season
              status = -1
              return
          end select
        endif

!     ------------------------------------
      case ("Palencia")
          select case (season)
            case (1)
              if (ndvi < 0.28) then
                co412 = (/2.68778125, 3.16549893e-02,0.0,0.0/)
                co470 = (/5.44559746,0.0,0.0,0.0/)
                co650 = (/0.0, 0.0, 0.0, 0.0/)
              else
                co412 = (/2.94527862, 0.0,0.0,0.0/)
                co470 = (/5.422391543, 7.49454787e-03,0.0,0.0/)
                co650 = (/0.0, 0.0, 0.0, 0.0/)
              end if
              
            case (2)
              if (ndvi >= 0.5) then
                co412 = (/2.11712331, 1.89785795e-02, -3.36702190e-04, 1.02791555e-05/)
                co470 = (/4.12074484, 3.94824004e-02, 0.0, 0.0/)
                co650 = (/0.0, 0.0, 0.0, 0.0/)
                
              else if (ndvi < 0.4) then
                co412 = (/3.38904321, 2.76987802e-02, 0.0, 0.0/)
                co470 = (/6.42545968, 5.19763334e-02, 0.0, 0.0/)
                co650 = (/0.0, 0.0, 0.0, 0.0/)
              else
                co412 = (/2.43331357, 3.25327465e-02, 0.0, 0.0/)
                co470 = (/4.76639599, 4.15582016e-02, 1.11252132e-04, 0.0/)
                co650 = (/0.0, 0.0, 0.0, 0.0/)
              end if

            case (3)
              if (ndvi >= 0.32) then
                co412 = (/2.84519284, 2.72442032e-2, 0.0, 0.0/)
                co470 = (/6.15825559, 3.47664076e-2, 1.9403065e-4, 0.0/)
                co650 = (/0.0, 0.0, 0.0, 0.0/)
                
                !co412 = (/2.04842972, 2.45804399e-02, 0.0, 0.0/)
                !co470 = (/5.86019014, 2.96439063e-02, 0.0, 0.0/)
                !co650 = (/0.0, 0.0, 0.0, 0.0/)
              else
                co412 = (/2.6482836, 4.0263797e-2, 0.0, 0.0/)
                co470 = (/7.06495618, 4.88471232e-2, -8.39825079e-5, 0.0/)
                co650 = (/0.0, 0.0, 0.0, 0.0/)
                
                !co412 = (/1.90045829, -2.51858297e-03, 8.11589166e-04, 0.0/)
                !co470 = (/7.05691697, 4.93936950e-02, 0.0, 0.0/)
                !co650 = (/0.0, 0.0, 0.0, 0.0/)
                
              end if
            
            case (4)
              co412 = (/2.97435558, 1.72237964e-02, 0.0, 0.0/)
              co470 = (/7.41527891, 5.88402071e-02, 7.96356631e-05, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
              
            case default
              print *, "ERROR: Invalid season specified: ", season
              status = -1
              return
          end select

!     ------------------------------------
      case ("Lecce_University")
          select case (season)
            case (1)
              co412 = (/3.66811317, 2.12308807e-02, 0.0, 0.0/)
              co470 = (/5.10650889, 7.12197835e-02, 6.15205093e-05, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
                            
            case (2)
              if (ndvi >= 0.45) then
                co412 = (/3.42708429, 2.01945894e-02, 5.70725971e-04, 0.0/)
                co470 = (/4.23230444, 3.76471806e-02, 6.98131497e-04, 0.0/)
                co650 = (/0.0, 0.0, 0.0, 0.0/)
              else if (ndvi >= 0.40 .AND. ndvi < 0.45) then
                co412 = (/4.82551344, 6.13242297e-03, 0.0, 0.0/)
                co470 = (/5.92910956, 5.63037114e-02, -9.05894219e-04, 0.0/)
                co650 = (/0.0, 0.0, 0.0, 0.0/)
              else
                co412 = (/4.22986020, 2.78424337e-02, -1.44379581e-04, 0.0/)
                co470 = (/5.51855341, 6.20406750e-02, -2.63676847e-04, 0.0/)
                co650 = (/0.0, 0.0, 0.0, 0.0/)
              end if
                            
            case (3)
              co412 = (/5.4920958, 4.2560058e-2, -3.3332266e-4, 0.0/)
              co470 = (/6.1583353, 6.8705510e-2, -2.7344898e-4, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
                                     
            case (4)
              co412 = (/5.3015222, 3.2603717e-2, 4.9840706e-4, 0.0/)
              co470 = (/5.8268639, 5.6895741e-2, 4.2632621e-4, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
              
            case default
              print *, "ERROR: Invalid season specified: ", season
              status = -1
              return
          end select

!     ------------------------------------
      case ("Carpentras")
          select case (season)
            case (1)
              if (ndvi >= 0.36) then
                co412 = (/3.63305224, 1.57411548e-2, 0.0, 0.0/)
                co470 = (/5.12339543, 3.27343565e-2, 0.0, 0.0/)
                co650 = (/0.0, 0.0, 0.0, 0.0/)
              else
                co412 = (/4.14705645, 4.21266106e-2, 0.0, 0.0/)
                co470 = (/6.03776651, 6.0694725e-2, 0.0, 0.0/)
                co650 = (/0.0, 0.0, 0.0, 0.0/)
              end if
  
!              -- 25km BRDF
!              co412 = (/3.65420405, 3.19272709e-02, 3.05220414e-04, 0.0/)
!              co470 = (/5.46812267, 8.72999543e-02, 1.51137682e-03, 0.0/)
!              co650 = (/0.0, 0.0, 0.0, 0.0/)
              
            case (2)
              if (ndvi < 0.3) then
                co412 = (/5.31754912, 5.6181342e-2, 0.0, 0.0/)
                co470 = (/7.21956575, 8.6234168e-2, 0.0, 0.0/)
                co650 = (/0.0, 0.0, 0.0, 0.0/)
              else if (ndvi < 0.42) then
                co412 = (/4.45713707, 2.95018243e-3, -9.6135173e-4, 4.91402021e-5/)
                co470 = (/5.88268514, 4.70497569e-2, 3.55853726e-5, 0.0/)
                co650 = (/0.0, 0.0, 0.0, 0.0/)
              else
                co412 = (/3.48800908, 3.96286071e-2, 0.0, 0.0/)
                co470 = (/5.29642181, 4.57011215e-2, -6.90453319e-4, 3.20776151e-5/)
                co650 = (/0.0, 0.0, 0.0, 0.0/)
              end if

!              -- 25km BRDF
!              if (ndvi >= 0.46) then
!                co412 = (/2.97987530, 2.27364122e-02, 9.88588305e-04, 0.0/)
!                co470 = (/4.54948280, 4.06782524e-02, 6.53794050e-04, 0.0/)
!                co650 = (/0.0, 0.0, 0.0, 0.0/)
!              else if (ndvi >= 0.4) then
!                co412 = (/3.42964905, -1.07598875e-02, 6.53505474e-04, 3.75382640e-05/)
!                co470 = (/4.74069593, 2.36714602e-02, 8.20597094e-04, 2.07381443e-05/)
!                co650 = (/0.0, 0.0, 0.0, 0.0/)
!              else if (ndvi >= 0.3) then
!                co412 = (/3.92563927, 3.23587748e-02, 1.11824371e-03, 0.0/)
!                co470 = (/5.51160066, 5.03020135e-02, 9.10344431e-04, 0.0/)
!                co650 = (/0.0, 0.0, 0.0, 0.0/)
!              else
!                co412 = (/5.19496351, 5.47740577e-02, 0.0, 0.0/)
!                co470 = (/6.90182111, 7.42751209e-02, 0.0, 0.0/)
!                co650 = (/0.0, 0.0, 0.0, 0.0/)
!              end if
                          
            case (3)
              if (ndvi < 0.42) then
                co412 = (/3.3539774, 4.88920328e-2, -1.7442563e-3, 4.09605495e-5/)
                co470 = (/5.11861399, 6.21094378e-2, 0.0, 0.0/)
                co650 = (/0.0, 0.0, 0.0, 0.0/)
              else
                co412 = (/2.94725072, 3.72651435e-2, 5.28386136e-4, 0.0/)
                co470 = (/4.75336791, 6.02706548e-2, -3.53602622e-4, 2.01191512e-5/)
                co650 = (/0.0, 0.0, 0.0, 0.0/)
              end if

!              -- 25km BRDF
!              if (ndvi >= 0.46) then
!                co412 = (/4.33893052, 4.70373358e-02, 2.55311732e-04, -4.89431042e-06/)
!                co470 = (/2.72831822, 1.99844060e-02, 4.82744505e-04, 0.0/)
!                co650 = (/0.0, 0.0, 0.0, 0.0/)
!              else
!                co412 = (/2.72831822, 1.99844060e-02, 4.82744505e-04, 0.0/)
!                co470 = (/4.75588609, 4.56984776e-02, 9.02843443e-05, 0.0/)
!                co650 = (/0.0, 0.0, 0.0, 0.0/)
!              end if
              
            case (4)
              co412 = (/3.13735391, 3.48937426e-2, 3.64366638e-4, 0.0/)
              co470 = (/4.56241718, 5.2711995e-2, 9.26548766e-4, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)

!              -- 25km BRDF
!              if (ndvi > 0.46) then
!                co412 = (/2.90669844, 2.65353619e-02, 4.88268055e-04, 0.0/)
!                co470 = (/4.17834040, 4.11956853e-02, 8.09480813e-04, 0.0/)
!                co650 = (/0.0, 0.0, 0.0, 0.0/)
!              else
!                co412 = (/3.11168018, 1.93021460e-02, 0.0, 0.0/)
!                co470 = (/4.68012390, 5.57148485e-02, 5.760506467e-04, -8.35696941e-06/)
!                co650 = (/0.0, 0.0, 0.0, 0.0/)
!              end if
              
            case default
              print *, "ERROR: Invalid season specified: ", season
              status = -1
              return
          end select
        
        case ("Trelew")                            ! backward scattering
          select case (season)
            case (1)
              if (ndvi < 0.2) then
                co412 = (/7.12231945, 4.85895624e-02, -1.28573607e-03, 2.25761306e-05/)
                co470 = (/9.29057722, 8.45987919e-02, -7.61137837e-04, 0.0/)
                co650 = (/0.0, 0.0, 0.0, 0.0/)
              else
                co412 = (/6.67352731, 2.31096853e-02, -1.07985947e-03, 5.12622592e-05/)
                co470 = (/8.92322447, 7.99348838e-02, 0.0, 0.0/)
                co650 = (/0.0, 0.0, 0.0, 0.0/)
              end if
            case (2)
              if (ndvi < 0.18) then
                co412 = (/6.79254216, 6.60335364e-02, 0.0, 0.0/)
                co470 = (/9.05448478, 1.01175023e-01, 1.31793002e-02, 0.0/)
                co650 = (/0.0, 0.0, 0.0, 0.0/)
              else
                co412 = (/6.63306271, 4.20332840e-02, 1.87758245e-04, 0.0/)
                co470 = (/8.71407292, 7.83112414e-02, 1.01125671e-03, 0.0/)
                co650 = (/0.0, 0.0, 0.0, 0.0/)
              end if
            case (3)
              if (ndvi < 0.15) then
                co412 = (/6.22897476, 0.0, 0.0, 0.0/)
                co470 = (/8.73022491, 0.0, 0.0, 0.0/)
                co650 = (/0.0, 0.0, 0.0, 0.0/)
              else
                co412 = (/5.83693596, 2.94875730e-02, 4.27587428e-04, 0.0/)
                co470 = (/8.12482274, 5.32465081e-02, 0.0, 0.0/)
                co650 = (/0.0, 0.0, 0.0, 0.0/)
              end if
            case (4)
              if (ndvi < 0.2) then
                co412 = (/9.20319079, 7.41965188e-02, -1.80650001e-04, 0.0/)
                co470 = (/6.17664662, 7.45518438e-02, -3.64408664e-04, -1.05957912e-05/)
                co650 = (/0.0, 0.0, 0.0, 0.0/)
              else
                co412 = (/8.20920358, 5.96572455e-02, 7.20029748e-04, 1.07035450e-05/)
                co470 = (/6.35028669, 3.72669462e-02, 1.39396118e-04, 1.09315712e-05/)
                co650 = (/0.0, 0.0, 0.0, 0.0/)
              end if
            case default
              print *, "ERROR: Invalid season specified: ", season
              status = -1
              return
          end select
      
      case ("Pune")
        select case (season)
          case (1)
            co412 = (/2.2107312, 1.9546884e-2, 1.0217331e-3, 0.0/)
            co470 = (/3.9707635, 4.5260501e-2, 8.2287074e-4, 0.0/)
            co650 = (/0.0, 0.0, 0.0, 0.0/)

          case (2)
            co412 = (/2.3410848, 1.1693201e-2, 9.9765075e-4, 0.0/)
            co470 = (/4.5157122, 5.7494184e-2, 2.3163434e-4, 0.0/)
            co650 = (/0.0, 0.0, 0.0, 0.0/)

          case (3)
          !Pune is unable to perform BRDF study at summer due to monsoon
            !use fall for summer, 26 January 2018, JLee TEST
            co412 = (/1.6509085, 3.2402477e-2, 1.5609563e-3, 0.0/)
            co470 = (/3.0505165, 4.3231390e-2, 1.4585233e-3, 0.0/)
            co650 = (/0.0, 0.0, 0.0, 0.0/)

          case (4)
            co412 = (/1.6509085, 3.2402477e-2, 1.5609563e-3, 0.0/)
            co470 = (/3.0505165, 4.3231390e-2, 1.4585233e-3, 0.0/)
            co650 = (/0.0, 0.0, 0.0, 0.0/)
            
          case default
            print *, "ERROR: Invalid season specified: ", season
            status = -1
            return
        end select

      case ("Evora")
        if (fwd_scat) then                      ! forward scattering
        select case (season)
          case (1)
            if (ndvi >= 0.54) then
              co412 = (/1.8666770, 0.026571400, 0.0, 0.0/)
              co470 = (/3.24662224, 2.49647745e-02, 0.0, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/2.1474521, 0.019420845, 0.0, 0.0/)
              co470 = (/3.8045176, 0.038118250, 0.00044871309, 1.4730596e-05/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
          case (2)
            if (ndvi >= 0.54) then
              co412 = (/2.67001676, -1.07772899e-02, 4.98010288e-04, -7.06071474e-06/)
              co470 = (/3.47056999, 1.75578420e-02, 1.07148973e-04, 7.60010333e-06/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else if (ndvi >= 0.42) then
              co412 = (/3.57946117, -2.10298778e-02, -4.11579021e-04, 1.55923969e-05/)
              co470 = (/4.34427890, 1.57998603e-02, -5.70677301e-04, 2.11628341e-05/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/3.49153138, -1.82641547e-02, 1.71132131e-03, -2.49950954e-05/)
              co470 = (/4.24504494, 4.49089153e-02, 1.81406770e-03, -3.88986042e-05/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
          case (3)
            if (ndvi >= 0.36) then
              co412 = (/3.51804408, -5.43557298e-03, -4.11664584e-05, 4.19463147e-06/)
              co470 = (/4.42020622, 3.56782224e-02, 2.32453523e-04, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else if (ndvi >= 0.32) then
              co412 = (/4.02658655, -3.41714839e-02, 1.49532934e-03, -1.20524324e-05/)
              co470 = (/5.29484758, 3.27714291e-02, 0.0, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else if (ndvi >= 0.27) then
              co412 = (/4.29981928, -1.11379081e-02, 1.05046171e-03, -1.30375107e-05/)
              co470 = (/5.64253900, 3.74786248e-02, 5.83243337e-04, -1.00991872e-05/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/4.86874172, 2.25252861e-02, 1.87141545e-04, 0.0/)
              co470 = (/6.16672976, 5.10490840e-02, 9.10921618e-05, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
          case (4)
            if (ndvi >= 0.47) then
              co412 = (/1.58754835, 1.20745361e-02, 0.0, 0.0/)
              co470 = (/2.98246568, 9.57936578e-03, 0.0, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else if (ndvi >= 0.36) then
              co412 = (/2.94407598, 1.56448599e-02, 0.0, 0.0/)
              co470 = (/4.64143550, 3.94157958e-02, 5.30955969e-05, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else if (ndvi >= 0.27) then
              co412 = (/3.55241655, 1.07045809e-02, 4.82439668e-04, 0.0/)
              co470 = (/5.21955006, 5.10611495e-02, 6.06966186e-04, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/4.93773812, -6.30823250e-04, 4.83307917e-04, 0.0/)
              co470 = (/6.54996856, 3.14618262e-02, 2.51836348e-04, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
          case default
            print *, "ERROR: Invalid season specified: ", season
            status = -1
            return
        end select
        
      else                                      ! backward scattering
        select case (season)
          case (1)
                  if (ndvi >= 0.54) then
              co412 = (/1.8666770, 0.026571400, 0.0, 0.0/)
              co470 = (/3.24662224, 2.49647745e-02, 0.0, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/2.1474521, 0.019420845, 0.0, 0.0/)
              co470 = (/3.8045176, 0.038118250, 0.00044871309, 1.4730596e-05/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
          case (2)
            if (ndvi >= 0.54) then
              co412 = (/2.67001676, -1.07772899e-02, 4.98010288e-04, -7.06071474e-06/)
              co470 = (/3.47056999, 1.75578420e-02, 1.07148973e-04, 7.60010333e-06/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else if (ndvi >= 0.42) then
              co412 = (/3.57946117, -2.10298778e-02, -4.11579021e-04, 1.55923969e-05/)
              co470 = (/4.34427890, 1.57998603e-02, -5.70677301e-04, 2.11628341e-05/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/3.49153138, -1.82641547e-02, 1.71132131e-03, -2.49950954e-05/)
              co470 = (/4.24504494, 4.49089153e-02, 1.81406770e-03, -3.88986042e-05/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
          case (3)
            if (ndvi >= 0.36) then
              co412 = (/3.51804408, -5.43557298e-03, -4.11664584e-05, 4.19463147e-06/)
              co470 = (/4.42020622, 3.56782224e-02, 2.32453523e-04, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else if (ndvi >= 0.32) then
              co412 = (/4.02658655, -3.41714839e-02, 1.49532934e-03, -1.20524324e-05/)
              co470 = (/5.29484758, 3.27714291e-02, 0.0, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else if (ndvi >= 0.27) then
              co412 = (/4.29981928, -1.11379081e-02, 1.05046171e-03, -1.30375107e-05/)
              co470 = (/5.64253900, 3.74786248e-02, 5.83243337e-04, -1.00991872e-05/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/4.48791268, 2.86111924e-02, 0.0, 0.0/)
              co470 = (/6.24638193, 6.48367744e-02, -1.46322122e-03, 2.77197095e-05/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
          case (4)
            if (ndvi >= 0.47) then
              co412 = (/1.58754835, 1.20745361e-02, 0.0, 0.0/)
              co470 = (/2.98246568, 9.57936578e-03, 0.0, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else if (ndvi >= 0.36) then
              co412 = (/2.94407598, 1.56448599e-02, 0.0, 0.0/)
              co470 = (/4.64143550, 3.94157958e-02, 5.30955969e-05, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else if (ndvi >= 0.27) then
              co412 = (/3.55241655, 1.07045809e-02, 4.82439668e-04, 0.0/)
              co470 = (/5.21955006, 5.10611495e-02, 6.06966186e-04, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/4.93773812, -6.30823250e-04, 4.83307917e-04, 0.0/)
              co470 = (/6.54996856, 3.14618262e-02, 2.51836348e-04, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
          case default
            print *, "ERROR: Invalid season specified: ", season
            status = -1
            return
        end select
      end if
      case ("Blida")
        if (fwd_scat) then                      ! forward scattering
        select case (season)
          case (1)
            if (ndvi >= 0.42) then
              co412 = (/1.32086397, 9.43028773e-03, 6.22359822e-04, 0.0/)
              co470 = (/2.85239261, 3.20907415e-02, 8.29035330e-04, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/1.36446568, 0.0, 0.0, 0.0/)
              co470 = (/2.92020512, -1.43458296e-02, -1.39739540e-03, 3.21024412e-05/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
!            co412 = (/1.91200337, 1.88742388e-02, 0.0, 0.0/)
!            co470 = (/3.51364747, 3.85886985e-02, 0.0, 0.0/)
!            co650 = (/0.0, 0.0, 0.0, 0.0/)
          case (2)
            co412 = (/2.29853905, 8.40506270e-03, 4.88300278e-04, 0.0/)
            co470 = (/3.44064029, 3.21979141e-02, -4.24620655e-04, 2.12964306e-05/)
            co650 = (/0.0, 0.0, 0.0, 0.0/)
!            if (ndvi >= 0.36) then
!              co412 = (/2.98730763, 1.52743195e-02, 3.42813724e-04, 0.0/)
!              co470 = (/4.10360363, 3.81184181e-02, -3.94054665e-04, 1.91764237e-05/)
!              co650 = (/0.0, 0.0, 0.0, 0.0/)
!            else
!              co412 = (/3.36923597, 8.91861898e-03, 5.41011164e-04, 0.0/)
!              co470 = (/4.31305027, 4.16934936e-02, 5.07485951e-04, 0.0/)
!              co650 = (/0.0, 0.0, 0.0, 0.0/)
!            end if
          case (3)
            if (ndvi >= 0.36) then
              co412 = (/3.51804408, -5.43557298e-03, -4.11664584e-05, 4.19463147e-06/)
              co470 = (/4.42020622, 3.56782224e-02, 2.32453523e-04, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else if (ndvi >= 0.30) then
              co412 = (/4.02658655, -3.41714839e-02, 1.49532934e-03, -1.20524324e-05/)
              co470 = (/5.29484758, 3.27714291e-02, 0.0, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/4.29981928, -1.11379081e-02, 1.05046171e-03, -1.30375107e-05/)
              co470 = (/5.64253900, 3.74786248e-02, 5.83243337e-04, -1.00991872e-05/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
!            if (ndvi >= 0.27) then
!              co412 = (/4.14926099, 1.42958399e-03, 4.64133459e-04, 0.0/)
!              co470 = (/5.02372861, 4.57270649e-02, -5.17735411e-05, 4.71527134e-06/)
!              co650 = (/0.0, 0.0, 0.0, 0.0/)
!            else
!              co412 = (/4.13229478, 2.66646676e-02, 0.0, 0.0/)
!              co470 = (/5.01590910, 6.48840906e-02, 0.0, 0.0/)
!              co650 = (/0.0, 0.0, 0.0, 0.0/)
!            end if
          case (4)
            if (ndvi >= 0.36) then
              co412 = (/1.31986294, -1.07513582e-02, 4.65319050e-04, 2.35774710e-05/)
              co470 = (/3.03892353, 2.84719270e-02, 4.80364309e-04, 2.19066757e-05/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/2.18453166, 3.51533148e-02, 0.0, 0.0/)
              co470 = (/3.68111291, 5.65877577e-02, 0.0, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
!            if (ndvi >= 0.27) then
!              co412 = (/2.21833952, 1.02488472e-02, 0.0, 0.0/)
!              co470 = (/4.01205578, 4.26871364e-02, -3.01909029e-05, 0.0/)
!              co650 = (/0.0, 0.0, 0.0, 0.0/)
!            else
!              co412 = (/3.84250107, 3.15218358e-02, 0.0, 0.0/)
!              co470 = (/5.07606928, 5.52702453e-02, 1.11978961e-04, 0.0/)
!              co650 = (/0.0, 0.0, 0.0, 0.0/)
!            end if
          case default
            print *, "ERROR: Invalid season specified: ", season
            status = -1
            return
        end select
        
      else                                      ! backward scattering
        select case (season)
          case (1)
            if (ndvi >= 0.42) then
              co412 = (/1.32086397, 9.43028773e-03, 6.22359822e-04, 0.0/)
              co470 = (/2.85239261, 3.20907415e-02, 8.29035330e-04, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/1.36446568, 0.0, 0.0, 0.0/)
              co470 = (/2.92020512, -1.43458296e-02, -1.39739540e-03, 3.21024412e-05/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
!            co412 = (/1.91200337, 1.88742388e-02, 0.0, 0.0/)
!            co470 = (/3.51364747, 3.85886985e-02, 0.0, 0.0/)
!            co650 = (/0.0, 0.0, 0.0, 0.0/)
          case (2)
            co412 = (/2.29853905, 8.40506270e-03, 4.88300278e-04, 0.0/)
            co470 = (/3.44064029, 3.21979141e-02, -4.24620655e-04, 2.12964306e-05/)
            co650 = (/0.0, 0.0, 0.0, 0.0/)
!            if (ndvi >= 0.36) then
!              co412 = (/2.98730763, 1.52743195e-02, 3.42813724e-04, 0.0/)
!              co470 = (/4.10360363, 3.81184181e-02, -3.94054665e-04, 1.91764237e-05/)
!              co650 = (/0.0, 0.0, 0.0, 0.0/)
!            else
!              co412 = (/3.36923597, 8.91861898e-03, 5.41011164e-04, 0.0/)
!              co470 = (/4.31305027, 4.16934936e-02, 5.07485951e-04, 0.0/)
!              co650 = (/0.0, 0.0, 0.0, 0.0/)
!            end if
          case (3)
            if (ndvi >= 0.36) then
              co412 = (/3.51804408, -5.43557298e-03, -4.11664584e-05, 4.19463147e-06/)
              co470 = (/4.42020622, 3.56782224e-02, 2.32453523e-04, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else if (ndvi >= 0.30) then
              co412 = (/4.02658655, -3.41714839e-02, 1.49532934e-03, -1.20524324e-05/)
              co470 = (/5.29484758, 3.27714291e-02, 0.0, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/4.29981928, -1.11379081e-02, 1.05046171e-03, -1.30375107e-05/)
              co470 = (/5.64253900, 3.74786248e-02, 5.83243337e-04, -1.00991872e-05/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
!            if (ndvi >= 0.27) then
!              co412 = (/4.14926099, 1.42958399e-03, 4.64133459e-04, 0.0/)
!              co470 = (/5.02372861, 4.57270649e-02, -5.17735411e-05, 4.71527134e-06/)
!              co650 = (/0.0, 0.0, 0.0, 0.0/)
!            else
!              co412 = (/4.13229478, 2.66646676e-02, 0.0, 0.0/)
!              co470 = (/5.01590910, 6.48840906e-02, 0.0, 0.0/)
!              co650 = (/0.0, 0.0, 0.0, 0.0/)
!            end if
          case (4)
            if (ndvi >= 0.36) then
              co412 = (/1.31986294, -1.07513582e-02, 4.65319050e-04, 2.35774710e-05/)
              co470 = (/3.03892353, 2.84719270e-02, 4.80364309e-04, 2.19066757e-05/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/2.18453166, 3.51533148e-02, 0.0, 0.0/)
              co470 = (/3.68111291, 5.65877577e-02, 0.0, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
!            if (ndvi >= 0.27) then
!              co412 = (/2.21833952, 1.02488472e-02, 0.0, 0.0/)
!              co470 = (/4.01205578, 4.26871364e-02, -3.01909029e-05, 0.0/)
!              co650 = (/0.0, 0.0, 0.0, 0.0/)
!            else
!              co412 = (/3.84250107, 3.15218358e-02, 0.0, 0.0/)
!              co470 = (/5.07606928, 5.52702453e-02, 1.11978961e-04, 0.0/)
!              co650 = (/0.0, 0.0, 0.0, 0.0/)
!            end if
          case default
            print *, "ERROR: Invalid season specified: ", season
            status = -1
            return
        end select
      end if
      case ("Blida_High")
        if (fwd_scat) then                      ! forward scattering
        select case (season)
          case (1)
            if (ndvi >= 0.42) then
              co412 = (/1.32086397, 9.43028773e-03, 6.22359822e-04, 0.0/)
              co470 = (/2.85239261, 3.20907415e-02, 8.29035330e-04, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/1.36446568, 0.0, 0.0, 0.0/)
              co470 = (/2.92020512, -1.43458296e-02, -1.39739540e-03, 3.21024412e-05/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
!            co412 = (/1.91200337, 1.88742388e-02, 0.0, 0.0/)
!            co470 = (/3.51364747, 3.85886985e-02, 0.0, 0.0/)
!            co650 = (/0.0, 0.0, 0.0, 0.0/)
          case (2)
            co412 = (/2.29853905, 8.40506270e-03, 4.88300278e-04, 0.0/)
            co470 = (/3.44064029, 3.21979141e-02, -4.24620655e-04, 2.12964306e-05/)
            co650 = (/0.0, 0.0, 0.0, 0.0/)
!            if (ndvi >= 0.36) then
!              co412 = (/2.98730763, 1.52743195e-02, 3.42813724e-04, 0.0/)
!              co470 = (/4.10360363, 3.81184181e-02, -3.94054665e-04, 1.91764237e-05/)
!              co650 = (/0.0, 0.0, 0.0, 0.0/)
!            else
!              co412 = (/3.36923597, 8.91861898e-03, 5.41011164e-04, 0.0/)
!              co470 = (/4.31305027, 4.16934936e-02, 5.07485951e-04, 0.0/)
!              co650 = (/0.0, 0.0, 0.0, 0.0/)
!            end if
          case (3)
            if (ndvi >= 0.36) then
              co412 = (/3.51804408, -5.43557298e-03, -4.11664584e-05, 4.19463147e-06/)
              co470 = (/4.42020622, 3.56782224e-02, 2.32453523e-04, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else if (ndvi >= 0.30) then
              co412 = (/4.02658655, -3.41714839e-02, 1.49532934e-03, -1.20524324e-05/)
              co470 = (/5.29484758, 3.27714291e-02, 0.0, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/4.29981928, -1.11379081e-02, 1.05046171e-03, -1.30375107e-05/)
              co470 = (/5.64253900, 3.74786248e-02, 5.83243337e-04, -1.00991872e-05/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
!            if (ndvi >= 0.27) then
!              co412 = (/4.14926099, 1.42958399e-03, 4.64133459e-04, 0.0/)
!              co470 = (/5.02372861, 4.57270649e-02, -5.17735411e-05, 4.71527134e-06/)
!              co650 = (/0.0, 0.0, 0.0, 0.0/)
!            else
!              co412 = (/4.13229478, 2.66646676e-02, 0.0, 0.0/)
!              co470 = (/5.01590910, 6.48840906e-02, 0.0, 0.0/)
!              co650 = (/0.0, 0.0, 0.0, 0.0/)
!            end if
          case (4)
            if (ndvi >= 0.36) then
              co412 = (/1.31986294, -1.07513582e-02, 4.65319050e-04, 2.35774710e-05/)
              co470 = (/3.03892353, 2.84719270e-02, 4.80364309e-04, 2.19066757e-05/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/2.18453166, 3.51533148e-02, 0.0, 0.0/)
              co470 = (/3.68111291, 5.65877577e-02, 0.0, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
!            if (ndvi >= 0.27) then
!              co412 = (/2.21833952, 1.02488472e-02, 0.0, 0.0/)
!              co470 = (/4.01205578, 4.26871364e-02, -3.01909029e-05, 0.0/)
!              co650 = (/0.0, 0.0, 0.0, 0.0/)
!            else
!              co412 = (/3.84250107, 3.15218358e-02, 0.0, 0.0/)
!              co470 = (/5.07606928, 5.52702453e-02, 1.11978961e-04, 0.0/)
!              co650 = (/0.0, 0.0, 0.0, 0.0/)
!            end if
          case default
            print *, "ERROR: Invalid season specified: ", season
            status = -1
            return
        end select
        
      else                                      ! backward scattering
        select case (season)
          case (1)
            if (ndvi >= 0.42) then
              co412 = (/1.32086397, 9.43028773e-03, 6.22359822e-04, 0.0/)
              co470 = (/2.85239261, 3.20907415e-02, 8.29035330e-04, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/1.36446568, 0.0, 0.0, 0.0/)
              co470 = (/2.92020512, -1.43458296e-02, -1.39739540e-03, 3.21024412e-05/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
!            co412 = (/1.91200337, 1.88742388e-02, 0.0, 0.0/)
!            co470 = (/3.51364747, 3.85886985e-02, 0.0, 0.0/)
!            co650 = (/0.0, 0.0, 0.0, 0.0/)
          case (2)
            co412 = (/2.29853905, 8.40506270e-03, 4.88300278e-04, 0.0/)
            co470 = (/3.44064029, 3.21979141e-02, -4.24620655e-04, 2.12964306e-05/)
            co650 = (/0.0, 0.0, 0.0, 0.0/)
!            if (ndvi >= 0.36) then
!              co412 = (/2.98730763, 1.52743195e-02, 3.42813724e-04, 0.0/)
!              co470 = (/4.10360363, 3.81184181e-02, -3.94054665e-04, 1.91764237e-05/)
!              co650 = (/0.0, 0.0, 0.0, 0.0/)
!            else
!              co412 = (/3.36923597, 8.91861898e-03, 5.41011164e-04, 0.0/)
!              co470 = (/4.31305027, 4.16934936e-02, 5.07485951e-04, 0.0/)
!              co650 = (/0.0, 0.0, 0.0, 0.0/)
!            end if
          case (3)
            if (ndvi >= 0.36) then
              co412 = (/3.51804408, -5.43557298e-03, -4.11664584e-05, 4.19463147e-06/)
              co470 = (/4.42020622, 3.56782224e-02, 2.32453523e-04, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else if (ndvi >= 0.30) then
              co412 = (/4.02658655, -3.41714839e-02, 1.49532934e-03, -1.20524324e-05/)
              co470 = (/5.29484758, 3.27714291e-02, 0.0, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/4.29981928, -1.11379081e-02, 1.05046171e-03, -1.30375107e-05/)
              co470 = (/5.64253900, 3.74786248e-02, 5.83243337e-04, -1.00991872e-05/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
!            if (ndvi >= 0.27) then
!              co412 = (/4.14926099, 1.42958399e-03, 4.64133459e-04, 0.0/)
!              co470 = (/5.02372861, 4.57270649e-02, -5.17735411e-05, 4.71527134e-06/)
!              co650 = (/0.0, 0.0, 0.0, 0.0/)
!            else
!              co412 = (/4.13229478, 2.66646676e-02, 0.0, 0.0/)
!              co470 = (/5.01590910, 6.48840906e-02, 0.0, 0.0/)
!              co650 = (/0.0, 0.0, 0.0, 0.0/)
!            end if
          case (4)
            if (ndvi >= 0.36) then
              co412 = (/1.31986294, -1.07513582e-02, 4.65319050e-04, 2.35774710e-05/)
              co470 = (/3.03892353, 2.84719270e-02, 4.80364309e-04, 2.19066757e-05/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/2.18453166, 3.51533148e-02, 0.0, 0.0/)
              co470 = (/3.68111291, 5.65877577e-02, 0.0, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
!            if (ndvi >= 0.27) then
!              co412 = (/2.21833952, 1.02488472e-02, 0.0, 0.0/)
!              co470 = (/4.01205578, 4.26871364e-02, -3.01909029e-05, 0.0/)
!              co650 = (/0.0, 0.0, 0.0, 0.0/)
!            else
!              co412 = (/3.84250107, 3.15218358e-02, 0.0, 0.0/)
!              co470 = (/5.07606928, 5.52702453e-02, 1.11978961e-04, 0.0/)
!              co650 = (/0.0, 0.0, 0.0, 0.0/)
!            end if
          case default
            print *, "ERROR: Invalid season specified: ", season
            status = -1
            return
        end select
      end if
      case ("Ilorin")
        select case (season)
          case (1)
!            if (ndvi >= 0.24) then
              co412 = (/1.95397318, 1.73386252E-02, -3.19877587E-04, 4.02574564E-06/)
              co470 = (/2.76241404, 5.21461588E-02, -4.94776417E-04, 5.45084985E-06/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
!            else
!              co412 = (/4.97022734, 1.39816544E-02, -6.20379323E-06, 4.10908585E-06/)
!              co470 = (/7.51990477, 5.62140108E-02, -2.65648653E-04, 0.0/)
!              co650 = (/0.0, 0.0, 0.0, 0.0/)
!            end if
          case (2)
!            if (ndvi >= 0.24) then
!              co412 = (/4.95397318, 1.73386252E-02, -3.19877587E-04, 4.02574564E-06/)
!              co470 = (/6.76241404, 5.21461588E-02, -4.94776417E-04, 5.45084985E-06/)
!              co650 = (/0.0, 0.0, 0.0, 0.0/)
              co412 = (/1.95397318, 1.73386252E-02, -3.19877587E-04, 4.02574564E-06/)
              co470 = (/2.76241404, 5.21461588E-02, -4.94776417E-04, 5.45084985E-06/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
!            else
!              co412 = (/4.97022734, 1.39816544E-02, -6.20379323E-06, 4.10908585E-06/)
!              co470 = (/7.51990477, 5.62140108E-02, -2.65648653E-04, 0.0/)
!              co650 = (/0.0, 0.0, 0.0, 0.0/)
!            end if
          case (3)
!            if (ndvi >= 0.24) then
!              co412 = (/4.95397318, 1.73386252E-02, -3.19877587E-04, 4.02574564E-06/)
!              co470 = (/6.76241404, 5.21461588E-02, -4.94776417E-04, 5.45084985E-06/)
!              co650 = (/0.0, 0.0, 0.0, 0.0/)
              co412 = (/1.95397318, 1.73386252E-02, -3.19877587E-04, 4.02574564E-06/)
              co470 = (/2.76241404, 5.21461588E-02, -4.94776417E-04, 5.45084985E-06/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
!            else
!              co412 = (/4.97022734, 1.39816544E-02, -6.20379323E-06, 4.10908585E-06/)
!              co470 = (/7.51990477, 5.62140108E-02, -2.65648653E-04, 0.0/)
!              co650 = (/0.0, 0.0, 0.0, 0.0/)
!            end if
          case (4)
!            if (ndvi >= 0.24) then
!              co412 = (/4.95397318, 1.73386252E-02, -3.19877587E-04, 4.02574564E-06/)
!              co470 = (/6.76241404, 5.21461588E-02, -4.94776417E-04, 5.45084985E-06/)
!              co650 = (/0.0, 0.0, 0.0, 0.0/)
              co412 = (/1.95397318, 1.73386252E-02, -3.19877587E-04, 4.02574564E-06/)
              co470 = (/2.76241404, 5.21461588E-02, -4.94776417E-04, 5.45084985E-06/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
!            else
!              co412 = (/4.97022734, 1.39816544E-02, -6.20379323E-06, 4.10908585E-06/)
!              co470 = (/7.51990477, 5.62140108E-02, -2.65648653E-04, 0.0/)
!              co650 = (/0.0, 0.0, 0.0, 0.0/)
!            end if
          case default
            print *, "ERROR: Invalid season specified: ", season
            status = -1
            return
        end select
      case ("Ilorin_Transition")
        select case (season)
          case (1)
!            if (ndvi >= 0.24) then
              co412 = (/1.95397318, 1.73386252E-02, -3.19877587E-04, 4.02574564E-06/)
              co470 = (/2.76241404, 5.21461588E-02, -4.94776417E-04, 5.45084985E-06/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
!            else
!              co412 = (/4.97022734, 1.39816544E-02, -6.20379323E-06, 4.10908585E-06/)
!              co470 = (/7.51990477, 5.62140108E-02, -2.65648653E-04, 0.0/)
!              co650 = (/0.0, 0.0, 0.0, 0.0/)
!            end if
          case (2)
!            if (ndvi >= 0.24) then
!              co412 = (/4.95397318, 1.73386252E-02, -3.19877587E-04, 4.02574564E-06/)
!              co470 = (/6.76241404, 5.21461588E-02, -4.94776417E-04, 5.45084985E-06/)
!              co650 = (/0.0, 0.0, 0.0, 0.0/)
              co412 = (/1.95397318, 1.73386252E-02, -3.19877587E-04, 4.02574564E-06/)
              co470 = (/2.76241404, 5.21461588E-02, -4.94776417E-04, 5.45084985E-06/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
!            else
!              co412 = (/4.97022734, 1.39816544E-02, -6.20379323E-06, 4.10908585E-06/)
!              co470 = (/7.51990477, 5.62140108E-02, -2.65648653E-04, 0.0/)
!              co650 = (/0.0, 0.0, 0.0, 0.0/)
!            end if
          case (3)
!            if (ndvi >= 0.24) then
!              co412 = (/4.95397318, 1.73386252E-02, -3.19877587E-04, 4.02574564E-06/)
!              co470 = (/6.76241404, 5.21461588E-02, -4.94776417E-04, 5.45084985E-06/)
!              co650 = (/0.0, 0.0, 0.0, 0.0/)
              co412 = (/1.95397318, 1.73386252E-02, -3.19877587E-04, 4.02574564E-06/)
              co470 = (/2.76241404, 5.21461588E-02, -4.94776417E-04, 5.45084985E-06/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
!            else
!              co412 = (/4.97022734, 1.39816544E-02, -6.20379323E-06, 4.10908585E-06/)
!              co470 = (/7.51990477, 5.62140108E-02, -2.65648653E-04, 0.0/)
!              co650 = (/0.0, 0.0, 0.0, 0.0/)
!            end if
          case (4)
!            if (ndvi >= 0.24) then
!              co412 = (/4.95397318, 1.73386252E-02, -3.19877587E-04, 4.02574564E-06/)
!              co470 = (/6.76241404, 5.21461588E-02, -4.94776417E-04, 5.45084985E-06/)
!              co650 = (/0.0, 0.0, 0.0, 0.0/)
              co412 = (/1.95397318, 1.73386252E-02, -3.19877587E-04, 4.02574564E-06/)
              co470 = (/2.76241404, 5.21461588E-02, -4.94776417E-04, 5.45084985E-06/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
!            else
!              co412 = (/4.97022734, 1.39816544E-02, -6.20379323E-06, 4.10908585E-06/)
!              co470 = (/7.51990477, 5.62140108E-02, -2.65648653E-04, 0.0/)
!              co650 = (/0.0, 0.0, 0.0, 0.0/)
!            end if
          case default
            print *, "ERROR: Invalid season specified: ", season
            status = -1
            return
        end select

      case ("SACOL")
        select case (season)
          case (1)
              co412 = (/1.42564317, 0.0, 0.0, 0.0/)
              co470 = (/5.32085036, 7.76166789e-02, 0.0, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
          case (2)
            if (ndvi < 0.22) then
              co412 = (/2.56745262, 3.14455912e-02, 0.0, 0.0/)
              co470 = (/5.72543787, 6.57833956e-02, 0.0, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/2.98578012, 9.14981794e-03, 2.78390281e-05, 0.0/)
              co470 = (/5.35045898, 5.95767492e-02, -2.40088902e-04, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
          case (3)
            if (ndvi < 0.32) then
              co412 = (/3.09389618, 3.43945319e-02, -1.73503979e-04, 0.0/)
              co470 = (/5.26831030, 7.26053989e-02, -4.36121233e-04, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else if (ndvi >= 0.32 .AND. ndvi < 0.38) then
              co412 = (/2.29526550, -2.85560957e-02, 1.16556608e-03, 0.0/)
              co470 = (/4.35768768, 2.04838826e-02, 7.12226936e-04, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/2.57989638, 3.85962492e-02, -1.23228526e-03, 0.0/)
              co470 = (/4.19697042, 6.43432125e-02, -1.15802337e-03, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
          case (4)
            if (ndvi < 0.20) then
              co412 = (/7.26239004e-01, 1.6566520e-02, 0.0, 0.0/)
              co470 = (/5.08897915, 2.07599476e-02, -1.86393057e-03, 2.40076652e-05/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else if (ndvi >= 0.2 .AND. ndvi < 0.3) then
              co412 = (/1.45776255, 7.95302789e-03, 0.0, 0.0/)
              co470 = (/4.23968258, 4.32695605e-02, 0.0, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/1.01734497, 0.0, 0.0, 0.0/)
              co470 = (/3.58458577, 4.92096024e-02, -1.69568606e-04, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
          case default
            print *, "ERROR: Invalid season specified: ", season
            status = -1
            return
        end select
      
      case ("Mexico_City")
        select case (season)
          case (1)
            co412 = (/1.5925121, 3.3497126e-3, 5.2599315e-4, 0.0/)
            co470 = (/5.1297925, 2.8796298e-2, 1.3826814e-3, 0.0/)
            co650 = (/0.0, 0.0, 0.0, 0.0/)

          case (2)
            co412 = (/2.0907425, -9.7213367e-3, 5.6728594e-4, 0.0/)
            co470 = (/5.4475574, 4.8689334e-2, 1.3585394e-4, 0.0/)
            co650 = (/0.0, 0.0, 0.0, 0.0/)

          case (3)
            if (ndvi < 0.2) then
              co412 = (/1.31131496, 1.90175933e-02, 0.0, 0.0/)
              co470 = (/4.28780511, 4.85799304e-02, 1.65298184e-03, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/3.19067520, 2.81951590e-02, 0.0, 0.0/)
              co470 = (/3.14833814, 1.05654939e-01, 0.0, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
            
          case (4)
            co412 = (/1.8093284, -9.0862879e-4, 2.2399699e-4, 0.0/)
            co470 = (/5.1765348, 3.8358608e-2, 9.4711253e-4, 0.0/)
            co650 = (/0.0, 0.0, 0.0, 0.0/)
            
          case default
            print *, "ERROR: Invalid season specified: ", season
            status = -1
            return
        end select
      
      case ("Jaipur")
        select case (season)
          case (1)
            if (ndvi < 0.2) then
              co412 = (/3.96991225, 5.96628950e-02, 2.13615454e-04, 0.0/)
              co470 = (/6.99850585, 8.36445599e-02, -1.10860486e-04, 6.74081297e-06/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else if (ndvi < 0.25) then
              co412 = (/3.833049805, 4.78326202e-02, -6.22253363e-04, 0.0/)
              co470 = (/6.31127139, 6.19777232e-02, -1.56997452e-04, 2.82474197e-05/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/4.66180293, 2.20314064e-02, 5.87633091e-04, -4.13864227e-06/)
              co470 = (/6.57270431, 5.7747869e-02, 7.34838329e-06, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
            
          case (2)
            if (ndvi < 0.18) then
              co412 = (/4.20196229, 5.29500002e-02, 3.20204802e-04, 0.0/)
              co470 = (/7.22026945, 9.03690963e-02, 0.0, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/4.65595362, 7.51832486e-02, -1.85215442e-03, 2.82980108e-05/)
              co470 = (/7.06658854, 8.85016889e-02, -8.29003477e-04, 1.30125894e-05/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
          case (3)
            if (ndvi < 0.2) then
              co412 = (/4.28651135, 9.9311289e-02, -5.34064043e-04, 0.0/)
              co470 = (/7.49201557, 1.37003165e-01, -1.41178674e-03, 8.95252657e-06/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/7.02860038, -4.65828642e-03, 0.0, 0.0/)
              co470 = (/7.89866675, 2.63625768e-02, 0.0, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
            
          case (4)
            if (ndvi < 0.22) then
              co412 = (/2.82799186, 3.51207193e-02, 0.0, 0.0/)
              co470 = (/6.03474644, 6.89646847e-02, 3.96325172e-04, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else if (ndvi < 0.25) then
              co412 = (/3.52929901, 5.18358835e-02, 0.0, 0.0/)
              co470 = (/5.94045545, 7.50816643e-02, 6.01703567e-04, -4.63811122e-06/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/4.70013291, 2.23865261e-02, -2.10276745e-05, 1.11426756e-05/)
              co470 = (/6.46922757, 3.76858031e-02, 7.48824729e-04, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
          case default
            print *, "ERROR: Invalid season specified: ", season
            status = -1
            return
        end select
      
      case ("NW_India_Desert")
        select case (season)
          case (1)
            if (ndvi < 0.2) then
              co412 = (/3.96991225, 5.96628950e-02, 2.13615454e-04, 0.0/)
              co470 = (/6.99850585, 8.36445599e-02, -1.10860486e-04, 6.74081297e-06/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else if (ndvi < 0.25) then
              co412 = (/3.833049805, 4.78326202e-02, -6.22253363e-04, 0.0/)
              co470 = (/6.31127139, 6.19777232e-02, -1.56997452e-04, 2.82474197e-05/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/4.66180293, 2.20314064e-02, 5.87633091e-04, -4.13864227e-06/)
              co470 = (/6.57270431, 5.7747869e-02, 7.34838329e-06, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
            
          case (2)
            if (ndvi < 0.18) then
              co412 = (/4.20196229, 5.29500002e-02, 3.20204802e-04, 0.0/)
              co470 = (/7.22026945, 9.03690963e-02, 0.0, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/4.65595362, 7.51832486e-02, -1.85215442e-03, 2.82980108e-05/)
              co470 = (/7.06658854, 8.85016889e-02, -8.29003477e-04, 1.30125894e-05/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
          case (3)
            if (ndvi < 0.2) then
              co412 = (/4.28651135, 9.9311289e-02, -5.34064043e-04, 0.0/)
              co470 = (/7.49201557, 1.37003165e-01, -1.41178674e-03, 8.95252657e-06/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/7.02860038, -4.65828642e-03, 0.0, 0.0/)
              co470 = (/7.89866675, 2.63625768e-02, 0.0, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
            
          case (4)
            if (ndvi < 0.22) then
              co412 = (/2.82799186, 3.51207193e-02, 0.0, 0.0/)
              co470 = (/6.03474644, 6.89646847e-02, 3.96325172e-04, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else if (ndvi < 0.25) then
              co412 = (/3.52929901, 5.18358835e-02, 0.0, 0.0/)
              co470 = (/5.94045545, 7.50816643e-02, 6.01703567e-04, -4.63811122e-06/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/4.70013291, 2.23865261e-02, -2.10276745e-05, 1.11426756e-05/)
              co470 = (/6.46922757, 3.76858031e-02, 7.48824729e-04, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
          case default
            print *, "ERROR: Invalid season specified: ", season
            status = -1
            return
        end select
        
      case ("Solar_Village")
        select case (season)
          case (1)
            co412 = (/9.57450477, 6.14694128e-02, 0.0, 0.0/)
            if (vza < 20.0) then
              co470 = (/1.60353001e01, 1.10436421e-02, 7.06173197e-04, 0.0/)
            else if (vza < 50.0) then
              co470 = (/1.50683392e01, 9.11896309e-02, 8.77940132e-04, -2.02723278e-05/)
            else
              co470 = (/1.58007075e01, 1.11295620e-01, 0.0, 0.0/)
            end if
            co650 = (/0.0, 0.0, 0.0, 0.0/)
            
          case (2)
            co412 = (/8.74736445, 6.51554815e-02, 0.0, 0.0/)
            co470 = (/1.53111360e01, 1.05595037e-01, -7.31500339e-04, 0.0/)
            co650 = (/0.0, 0.0, 0.0, 0.0/)
            
            if (raa < 50.0) then
              co412 = (/8.63300701, 1.29057864e-01, -2.95916295e-03, 1.83110100e-05/)
              co470 = (/1.46695788e01, 9.60421930e-02, -2.40220378e-03, 2.99138049e-05/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else if (raa >= 50.0 .AND. raa < 90.0) then
              co412 = (/1.01609630e01, 1.02434047e-01, -1.36452416e-03, 0.0/)
              co470 = (/1.51255831e01, 7.36233156e-02, 4.74327855e-04, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else if (raa >= 90.0 .AND. raa < 160.0) then
              co412 = (/1.09242704e01, 0.0, 0.0, 0.0/)
              co470 = (/1.95602923e01, -5.44866657e-02, 0.0, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/7.76882726, 7.90204491e-02, 0.0, 0.0/)
              co470 = (/1.65201198e01, 3.67393757e-02, 0.0, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
            
          case (3)
            co412 = (/9.22657792, 5.53831366e-02, 0.0, 0.0/)
            if (vza < 25.0) then
              co470 = (/1.39265362e01, 9.32481786e-02, 0.0, 0.0/)
            else if (vza < 45.0) then
              co470 = (/1.47453243e01, 9.60301995e-02, 0.0, 0.0/)
            else
              co470 = (/1.54562369e01, 9.40036429e-02, 0.0, 0.0/)
            endif
            
            co412 = (/8.9590389, 0.017604729, 0.00076547611, -2.1966118e-06/)
            co470 = (/14.904570, 0.053656442, 0.0010690852, -1.8405743e-05/)
            co650 = (/0.0, 0.0, 0.0, 0.0/)
            
            if (raa < 90.0) then
              co412 = (/9.1177595, 0.039785655, -0.00062878201, 2.1237891e-05/)
              co470 = (/14.843738, 0.024588279, 0.00021114164, 1.2466120e-05/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            else
              co412 = (/7.55891171, 8.05857229e-02, 0.0, 0.0/)
              co470 = (/1.66963898e01, 1.07268522e-02, 2.79875503e-04, 0.0/)
              co650 = (/0.0, 0.0, 0.0, 0.0/)
            end if
          case (4)
            co412 = (/9.81027873, 5.79203714e-02, 0.0, 0.0/)
            if (vza < 25.0) then
              co470 = (/1.50017355e01, 8.33530298e-02, 0.0, 0.0/)
            else if (vza < 50.0) then
              co470 = (/1.54093325e01, 9.33699483e-02, -8.75858189e-05, 1.76518228e-05/)
            else
              co470 = (/1.63141003e10, 1.14528713e-01, 0.0, 0.0/)
            end if

          case default
            print *, "ERROR: Invalid season specified: ", season
            status = -1
            return
        end select
 
        case ("Yuma")
        select case (season)
          case (1)
            co412 = (/8.88808, 4.0613965e-2, 1.1960322e-3, 0.0/)
            co470 = (/10.78728, 6.864285e-2, 1.15425508e-3, 0.0/)
            co650 = (/0.0, 0.0, 0.0, 0.0/)

          case (2)
            co412 = (/8.7442283, 2.53485488e-2, 4.54854215e-4,0.0/)
            co470 = (/10.872899, 5.6826757e-2, 1.35925198e-4, 0.0/)
            co650 = (/0.0, 0.0, 0.0, 0.0/)
              
          case (3)
            co412 = (/8.7692651, 3.3495399e-2, 3.7264769e-4,0.0/)
            co470 = (/10.92206559, 5.96282369e-2, 2.79797715e-4, 0.0/)
            co650 = (/0.0, 0.0, 0.0, 0.0/)
              
          case (4)
            co412 = (/8.9693997, 3.65029508e-2, 6.62722818e-4,0.0/)
            co470 = (/11.0841727, 5.91805818e-2, 7.60347456e-4, 0.0/)
            co650 = (/0.0, 0.0, 0.0, 0.0/)
               
          case default
            print *, "ERROR: Invalid season specified: ", season
            status = -1
            return
        end select

      case default
        print *, "ERROR: Invalid AERONET site specified.  No BRDF values."
        status = -1
        return
    end select

!   -- coefficients above were calculated using (sca-120) rather than sca.  Must do same here.
    s412 = co412(1) + co412(2)*(sca-120.0) + co412(3)*((sca-120.0)**2) + co412(4)*((sca-120.0)**3)
    s470 = co470(1) + co470(2)*(sca-120.0) + co470(3)*((sca-120.0)**2) + co470(4)*((sca-120.0)**3)
    s650 = co650(1) + co650(2)*(sca-120.0) + co650(3)*((sca-120.0)**2) + co650(4)*((sca-120.0)**3)

    c412 = co412(1)
    c470 = co470(1)
    c650 = co650(1)
 
    if (dflag) then
      print '(A,A,A,I4,I4,2(F11.6,1X))', trim(func_name), ", site, month, season, ndvi, scat: ", &
      &     aero_site, month, season, ndvi, sca
      print '(A,A,4(F11.6,1X))', trim(func_name), ", BRDF coeffs412: ", co412
      print '(A,A,4(F11.6,1X))', trim(func_name), ", BRDF coeffs470: ", co470
      print '(A,A,4(F11.6,1X))', trim(func_name), ", BRDF coeffs650: ", co650
      print '(A,A,3(F11.6,1x))', trim(func_name), ", BRDF SR, 412, s470, s650: ", s412, s470, s650
    end if
    
    status = 0
    return
    
  end function get_aeronet_brdf_sr


   ! -- Based on matching AERONET sites, return AOT value using appropriate AOT model.
   ! --  Returns 0 on success, otherwise -1.
   real function get_aot500(lat, lon, elev, sa, season, ndvi, gzone, lc_type, stdv02, &
      &       aot412_91, aot412_93, aot412_94, aot412_96, aot412_995, &
      &       aot470_91, aot470_92, aot470_93, aot470_94, aot470_95, aot470_96, aot470_995, &
      &       aot412_91_dust, aot412_93_dust, aot412_94_dust, aot412_96_dust, aot412_995_dust, &
      &       aot470_91_dust, aot470_92_dust, aot470_93_dust, aot470_94_dust, aot470_95_dust, &
      &       aot470_96_dust, aot470_995_dust, ae, status, debug,platform) result(aot500)

      implicit none
    
      character(len=20), parameter  ::  func_name = "get_aot500"
      character(len=*), intent(in)           ::  platform
      real, intent(in)          ::  lat
      real, intent(in)          ::  lon
      real, intent(in)          ::  elev            ! -- surface elevation
      real, intent(in)          ::  sa              ! -- scattering angle
      integer, intent(in)       ::  season
      real, intent(in)          ::  ndvi
      integer, intent(in)       ::  gzone
      integer, intent(in)       ::  lc_type
      real, intent(in)          ::  stdv02
      real, intent(in)          ::  aot412_91, aot412_93, aot412_94, aot412_96, aot412_995
      real, intent(in)          ::  aot470_91, aot470_92, aot470_93, aot470_94, aot470_95
      real, intent(in)          ::  aot470_96, aot470_995
      real, intent(in)          ::  aot412_91_dust, aot412_93_dust, aot412_94_dust, aot412_96_dust, aot412_995_dust
      real, intent(in)          ::  aot470_91_dust, aot470_92_dust, aot470_93_dust, aot470_94_dust, aot470_95_dust
      real, intent(in)          ::  aot470_96_dust, aot470_995_dust
      real, intent(in)          ::  ae              ! angstrom exponent
      integer, intent(inout)    ::  status
      logical, intent(in), optional ::  debug
    
      real                      ::  refsr650
      integer                   ::  ilat, ilon
      integer                   ::  m
    
      character(len=255)                  ::  asite
      real, dimension(:), allocatable     ::  maero412, maero470, maero650
      integer, dimension(:), allocatable  ::  msiteindx
      integer, dimension(:), allocatable  ::  sorted
      real                                ::  frac
      integer                             ::  i, ii, jj, cnt
    
      real                                ::  aot500_1, aot500_2
    
      logical                             ::  dflag
    
      dflag = .false.
      if (present(debug)) dflag = debug
    
      status = 0
      aot500 = -999.0
    
      !   -- convert geolocation into array indices.
      ilat = floor(lat*10.0) + 900 + 1
      ilon = floor(lon*10.0) + 1800 + 1
    
      if (ilat > 1800) ilat = 1800
      if (ilon > 3600) ilon = 3600
      if (dflag) print *, trim(func_name)//', lat, lon, ilat, ilon: ', lat, lon, ilat, ilon
    
      !   -- set up our reference surface reflectance
      refsr650 = brdf650(ilon,ilat)
        
      !   -- do we have an AERONET site in the same zone with the same land cover type?
      m = 0
      m = count(aero_zones == gzone .AND. aero_types == lc_type .AND. (elev < 500 .EQV. aero_elev < 500))

      !   -- create exception for china, europe, spain, morocco -- match by zone and elevation only
      if (gzone == 16 .OR. (gzone == 17 .OR. gzone == 2) .OR. gzone == 22) then
         m = count(aero_zones == gzone .AND. (elev < 500 .EQV. aero_elev < 500))
      end if
    
      !     -- create exception for Fresno Valley, Australia - match by region only.
      if (gzone == 18 .OR. gzone == 12 .OR. (gzone == 26 .OR. gzone == 27) .OR. gzone == 29) then
         m = count(aero_zones == gzone)
      end if

      !   -- create exception for high elevation Tibet/China zone - match by region
      !   only.
      if (gzone == 28) then
         m = count(aero_zones == gzone)
      end if

      !   -- create exception for Jaipur zone - match by region only.
      if (gzone == 20) then
         m = count(aero_zones == gzone)
      end if

      !   -- create exception for NW_India_Desert zone - match by region only.
      if (gzone == 30) then
         m = count(aero_zones == gzone)
      end if

      !   -- create exception for Pune - match by region only.
      if (gzone == 19) then
         m = count(aero_zones == gzone)
      end if

      !   -- create exception for Kanpur - match by region only, 9 January 2018 JLee
      if (gzone == 15) then
         m = count(aero_zones == gzone)
      end if

      !   -- create exception for Sahel, geozone =5, landcover=2 to ignore elevation.
      if ((gzone == 5 .AND. lc_type == 2) .OR. gzone == 1 .or. gzone == 13) then
         m = count(aero_zones == gzone .AND. aero_types == lc_type)
      end if

      !   -- create exception for Barren North America, geozone =31 to ignore
      !   elevation above 750m.
      if (gzone == 31 .AND. elev < 750) then
         m = count(aero_zones == gzone)
      end if
        
      if (m > 0) then
         !     -- allocate our arrays to store the matching AERONET data.
         !     -- no explicit deallocate() as these should automatically be
         !     -- deallocated when the function ends.
         if (allocated(maero412)) deallocate(maero412, stat=status)
         if (allocated(maero470)) deallocate(maero470, stat=status)
         if (allocated(maero650)) deallocate(maero650, stat=status)
         if (allocated(msiteindx)) deallocate(msiteindx, stat=status)
         if (allocated(sorted)) deallocate(sorted, stat=status)
         allocate(maero412(m), maero470(m), maero650(m), msiteindx(m), sorted(m), stat=status)
         if (status /= 0) then
            print *, "ERROR: Failed to allocate AERONET 650 SR match arrays: ", status
            return
         end if
      
         cnt = 0
    
         !     -- get and store base table SR values at each matching AERONET site at 412, 470, and 650.
         do i = 1, size(aero_sites)        ! i = AERONET site index
            select case (gzone)
               case (2, 16, 17, 22)   ! -- china, europe, spain, morocco, only match by zone, elevation
                  if (aero_zones(i) == gzone .AND. (elev < 500 .EQV. aero_elev(i) < 500)) then
                     cnt = cnt + 1
                     maero412(cnt) = aero_sr412(i,season)
                     maero470(cnt) = aero_sr470(i,season)
                     maero650(cnt) = aero_sr650(i,3)        ! < -- always use summer for 650nm to match refsr650
                     msiteindx(cnt) = i
              
                     if (dflag) then
                        print '(3(A,1X),I4,I4)', trim(func_name), ', matching site: ', trim(aero_sites(i)), aero_zones(i), aero_types(i)
                        print '(A,A,3(F11.6))', trim(func_name), ', AERONET Baseline SR, 412, 470, 650: ', &
                           &  maero412(cnt), maero470(cnt), maero650(cnt)
                     end if
                  end if
               case (18, 12, 20, 26, 27, 28, 29, 30, 31)                  ! Fresno Valley, Australia, only match by region
                  if (aero_zones(i) == gzone) then
                     cnt = cnt + 1
                     maero412(cnt) = aero_sr412(i,season)
                     maero470(cnt) = aero_sr470(i,season)
                     maero650(cnt) = aero_sr650(i,3)        ! < -- always use summer for 650nm to match refsr650
                     msiteindx(cnt) = i
              
                     if (dflag) then
                        print '(A,A,A,I4,I4)', trim(func_name), ', matching site: ', trim(aero_sites(i)), aero_zones(i), aero_types(i)
                        print '(A,A,3(F11.6,1X))', trim(func_name), ', AERONET Baseline SR: ', maero412(cnt), maero470(cnt), maero650(cnt)
                     end if
                  end if
               case (15, 19)          ! Pune, only match by region, added Kanpur and India high elevation 31 January 2018 JLee
                  if (aero_zones(i) == gzone) then
                     cnt = cnt + 1
                     maero412(cnt) = aero_sr412(i,season)
                     maero470(cnt) = aero_sr470(i,season)
                     maero650(cnt) = aero_sr650(i,3)        ! < -- always use summer for 650nm to match refsr650
                     msiteindx(cnt) = i
              
                     if (dflag) then
                        print '(A,A,A,I4,I4)', trim(func_name), ', matching site: ', trim(aero_sites(i)), aero_zones(i), aero_types(i)
                        print '(A,A,3(F11.6,1X))', trim(func_name), ', AERONET Baseline SR: ', maero412(cnt), maero470(cnt), maero650(cnt)
                     end if
                  end if
               case (1, 5, 10, 13)         ! N. Africa, Solar Villge (Saudi Arabia)
                  if (aero_zones(i) == gzone .AND. aero_types(i) == lc_type) then
                     cnt = cnt + 1
                     maero412(cnt) = aero_sr412(i,season)
                     maero470(cnt) = aero_sr470(i,season)
                     maero650(cnt) = aero_sr650(i,3)        ! < -- always use summer for 650nm to match refsr650
                     msiteindx(cnt) = i
              
                     if (dflag) then
                        print '(A,A,A,I4,I4)', trim(func_name), ', matching site: ', trim(aero_sites(i)), aero_zones(i), aero_types(i)
                        print '(A,A,3(F11.6,1X))', trim(func_name), ', AERONET Baseline SR: ', maero412(cnt), maero470(cnt), maero650(cnt)
                     end if
                  end if
               case default
                  if (aero_zones(i) == gzone .AND. aero_types(i) == lc_type .AND. (elev < 500 .EQV. aero_elev(i) < 500)) then
                     cnt = cnt + 1
                     maero412(cnt) = aero_sr412(i,season)
                     maero470(cnt) = aero_sr470(i,season)
                     maero650(cnt) = aero_sr650(i,3)        ! < -- always use summer for 650nm to match refsr650
                     msiteindx(cnt) = i
          
                     if (dflag) then
                        print '(3(A,1X),I4,I4)', trim(func_name), ', matching site: ', trim(aero_sites(i)), aero_zones(i), aero_types(i)
                        print '(A,A,3(F11.6))', trim(func_name), ', AERONET Baseline SR, 412, 470, 650: ', &
                           &  maero412(cnt), maero470(cnt), maero650(cnt)
                     end if
                  end if
            end select
         end do
      
         !     -- can we interpolate between zone's AERONET sites?
         call sortrx(m, maero650, sorted)
         if (refsr650 >= minval(maero650) .AND. refsr650 < maxval(maero650)) then

            !       -- find where refsr650 fits in maero650() and interpolate between the two sites
            do i = 1, m-1
               if (refsr650 >= maero650(sorted(i)) .AND. refsr650 < maero650(sorted(i+1))) then
                  ii = sorted(i)
                  asite = aero_sites(msiteindx(ii))
                  aot500_1 = get_aeronet_aot500(trim(asite), lat, lon, sa, season, ndvi, stdv02, &
                     &          aot412_91, aot412_93, aot412_94, aot412_96, aot412_995, &
                     &          aot470_91, aot470_92, aot470_93, aot470_94, aot470_95, aot470_96, aot470_995, &
                     &          aot412_91_dust, aot412_93_dust, aot412_94_dust, aot412_96_dust, aot412_995_dust, &
                     &          aot470_91_dust, aot470_92_dust, aot470_93_dust, aot470_94_dust, aot470_95_dust, &
                     &          aot470_96_dust, aot470_995_dust, ae, status,platform, debug=dflag)
                  if (status /= 0) then
                     print *, "ERROR: Failed to get AOT@500nm from AERONET site: ", trim(asite), status
                     return
                  end if
            
                  jj = sorted(i+1)
                  asite = aero_sites(msiteindx(jj))
                  aot500_2 = get_aeronet_aot500(trim(asite), lat, lon, sa, season, ndvi, stdv02, &
                     &          aot412_91, aot412_93, aot412_94, aot412_96, aot412_995, &
                     &          aot470_91, aot470_92, aot470_93, aot470_94, aot470_95, aot470_96, aot470_995, &
                     &          aot412_91_dust, aot412_93_dust, aot412_94_dust, aot412_96_dust, aot412_995_dust, &
                     &          aot470_91_dust, aot470_92_dust, aot470_93_dust, aot470_94_dust, aot470_95_dust, &
                     &          aot470_96_dust, aot470_995_dust, ae, status,platform, debug=dflag)
                  if (status /= 0) then
                     print *, "ERROR: Failed to get AOT@500nm from AERONET site: ", trim(asite), status
                     return
                  end if

                  !           -- calculate AERONET site weights according to 650 values and adjust AOT.
                  frac = (refsr650-maero650(ii)) / (maero650(jj)-maero650(ii))
            
                  aot500 = (1.0-frac)*aot500_1 + frac*aot500_2
            
                  if (dflag) then
                     print *, trim(func_name), ", Pixel Ref. SR, 650: ", refsr650
                     !              print *, trim(func_name), ", Pixel Baseline SR, 412, 470: ", xsfc412(ilon,ilat), xsfc470(ilon,ilat)
                     print *, trim(func_name), ', interp sites: ', trim(aero_sites(msiteindx(ii))), ' ', trim(aero_sites(msiteindx(jj)))
                     print *, trim(func_name), ', aot values, 412: ', aot412_91, aot412_93, aot412_94, aot412_96, aot412_995
                     print *, trim(func_name), ", aot values, 470: ", aot470_91, aot470_92, aot470_93, &
                        &                                                aot470_94, aot470_95, aot470_96, aot470_995
                     print *, trim(func_name), ', aot500: ', aot500
                  end if
            
                  exit     !  jump out of loop, we're done!
          
               end if
            end do
        
         !     -- no interpolation, use single AERONET site.
         else
            if (refsr650 <= minval(maero650)) then
               ii = sorted(1)  ! AERONET site w/ min. sr650 value
            else
               ii = sorted(m)  ! AERONET site w/ max. sr650 value
            end if
  
            asite = aero_sites(msiteindx(ii))
            aot500 = get_aeronet_aot500(trim(asite), lat, lon, sa, season, ndvi, stdv02, &
               &          aot412_91, aot412_93, aot412_94, aot412_96, aot412_995, &
               &          aot470_91, aot470_92, aot470_93, aot470_94, aot470_95, aot470_96, aot470_995, &
               &          aot412_91_dust, aot412_93_dust, aot412_94_dust, aot412_96_dust, aot412_995_dust, &
               &          aot470_91_dust, aot470_92_dust, aot470_93_dust, aot470_94_dust, aot470_95_dust, &
               &          aot470_96_dust, aot470_995_dust, ae, status, platform,debug=dflag)
            if (status /= 0) then
               print *, "ERROR: Failed to get AOT at 500nm over AERONET site, single: ", trim(asite), status
               return
            end if
        
            if (dflag) then
               print *, trim(func_name), ", Pixel Ref. SR, 650: ", refsr650
               !          print *, trim(func_name), ", Pixel Baseline SR, 412, 470: ", xsfc412(ilon,ilat), xsfc470(ilon,ilat)
               print *, trim(func_name), ', interp site: ', trim(aero_sites(msiteindx(ii)))
               print *, trim(func_name), ', aot values, 412: ', aot412_91, aot412_93, aot412_94, aot412_96, aot412_995
               print *, trim(func_name), ", aot values, 470: ", aot470_91, aot470_92, aot470_93, &
                  &                                                aot470_94, aot470_95, aot470_96, aot470_995
               print *, trim(func_name), ', aot500: ', aot500
            end if
            
         end if
      !   -- no matching site, no rules to select model.  Return error.
      else
         status = -1
         return
      end if

      return
            
   end function get_aot500
  
   ! -- returns AOT @ 500 nm over the AERONET site, aero_site for the given season.
   ! -- NOTE: thresholds below based on case studies over the each site and is specific to
   ! --        that site.
   real function get_aeronet_aot500(aero_site, lat, lon, sca, season, ndvi, stdv02, &
      &       aot412_91, aot412_93, aot412_94, aot412_96, aot412_995, &
      &       aot470_91, aot470_92, aot470_93, aot470_94, aot470_95, aot470_96, aot470_995, &
      &       aot412_91_dust, aot412_93_dust, aot412_94_dust, aot412_96_dust, aot412_995_dust, &
      &       aot470_91_dust, aot470_92_dust, aot470_93_dust, aot470_94_dust, aot470_95_dust, &
      &       aot470_96_dust, aot470_995_dust, ae, status,platform, debug) result(aot500)

      implicit none
    
      character(len=20), parameter            ::  func_name = "get_aeronet_aot500"
      character(len=*), intent(in)            ::  platform
      character(len=*), intent(in)            ::  aero_site
      real, intent(in)                        ::  sca           ! scattering angle
      real, intent(in)                        ::  lat, lon
      integer, intent(in)                     ::  season
      real, intent(in)                        ::  ndvi
      real, intent(in)                        ::  stdv02
      real, intent(in)                        ::  aot412_91, aot412_93, aot412_94, aot412_96, aot412_995
      real, intent(in)                        ::  aot470_91, aot470_92, aot470_93, aot470_94
      real, intent(in)                        ::  aot470_95, aot470_96, aot470_995
      real, intent(in)          ::  aot412_91_dust, aot412_93_dust, aot412_94_dust, aot412_96_dust, aot412_995_dust
      real, intent(in)          ::  aot470_91_dust, aot470_92_dust,aot470_93_dust, aot470_94_dust, aot470_95_dust
      real, intent(in)          ::  aot470_96_dust, aot470_995_dust
      real, intent(in)                        ::  ae            ! angstrom exponent
      integer, intent(inout)                  ::  status
      logical, intent(in), optional           ::  debug
    
      real          ::  aot412_92
      real          ::  dd
      real          ::  model_frac, model_frac2
      logical       ::  dflag
    
      aot500 = -999.0
      model_frac = 0.0
      model_frac2 = 0.0
      status = 0
    
      dflag = .false.
      if (present(debug)) dflag = debug

      select case (aero_site)
      
         case ("Agoufou")  !----------------------------------------------
            select case (season)
        
               !         -- winter
               case (1)
                  aot500 = aot412_94
                  if (aot412_94 >= 0.6) then
                     aot500 = (aot412_93 + aot412_91)/2.0
                  end if
            
               !         -- spring
               case (2)
                  aot500 = aot470_96
                  if (aot412_94 >= 0.6) then
                     aot500 =aot470_92
                  end if
                  if (aot500 >= 1.3) then
                     aot500 = (aot470_91 + aot470_92)/2.0
                  end if
            
               !         -- summer
               case (3)
                  aot500 = aot470_96
                  if (aot470_96 >= 0.5) then
                     aot500 = aot470_93
                  end if
                  if (aot470_96 >= 0.7) then
                     aot500 = aot470_91 * 1.1
                  end if
            
               !         -- fall
               case (4)
                  aot500 = aot412_94
                  if (aot412_94 >= 0.5) then
                     aot500 = (aot412_91 + aot412_93) / 2.0
                  end if
            
               !         -- default
               case default
                  print *, "ERROR: Invalid season specified: ", season
                  status = -1
                  return
        
            end select
        
         case ("IER_Cinzana")  !------------------------------------------
            select case (season)

               !         -- winter
               case (1)
                  aot500 = aot470_94
                  if (aot470_94 >= 0.8) then
                     aot500 = aot470_93
                  end if
                 !if (aot470_94 >= 0.7) then
                 !  aot500 = aot470_91
                 !end if

               !         -- spring
               case (2)
                  aot500 = aot470_96
                  if (aot500 > 0.6) then
                     aot500 = aot470_95
                  end if
                 !if (aot470_94 > 1.0) then
                 !  aot500 = aot470_93
                 !endif

               !         -- summer
               case (3)
                  aot500 = aot470_96
                  if (ndvi >= 0.3) then
                     aot500 = aot412_94
                  else
                     aot500 = aot470_995
                  end if

               !         -- fall
               case (4)
                  if (ndvi < 0.36) then
                     aot500 = aot412_94
                  else
                     aot500 = aot470_96
                  endif
          
               case default
                  print *, "ERROR: Invalid season specified: ", season
                  status = -1
                  return
          
            end select
        
         case ("Zinder_Airport")  !--------------------------------------
            select case (season)

               !         -- winter
               case (1)
                  aot500 = aot470_96
            
               !         -- spring
               case (2)
                  aot500 = aot412_94
            
               !         -- summer
               case (3)
                  aot500 = aot470_96
                  if (aot470_96 > 0.6 .AND. ndvi < 0.18) then
                     aot500 = aot470_93
                  end if
                  if (aot470_96 > 1.0 .AND. ndvi < 0.18) then
                     aot500 = aot470_92
                  end if
            
               !         -- fall
               case (4)
                  aot500 = aot412_94
            
               case default
                  print *, "ERROR: Invalid season specified: ", season
                  status = -1
                  return
            
            end select
        
         case ("Banizoumbou")  !------------------------------------------
            select case (season)

               !         -- winter
               case (1)
                  aot500 = aot412_94
                  if (aot500 < 0.4) then
                     aot500 = (aot412_96 + aot412_995) / 2.0
                  else
                     aot500 = aot412_93
                  end if
            
               !            if (aot470_94 > 0.6) then
               !              aot500 = (aot470_96 + aot470_94)/2.0
               !            end if
            
               !         -- spring
               case (2)
                  aot500 = aot412_93
                  !if (ndvi >= 0.12) then
                 !  aot500 = aot470_94
                 !end if
               !         -- summer
               case (3)
                  aot500 = (aot470_96 + aot470_995) / 2.0
                  if (aot470_96 > 0.7) then
                     aot500 = aot470_96
                  end if
            
               !         -- fall
               case (4)
                  aot500 = aot470_96
               !            if (ndvi > 0.24) then
               !              aot500 = aot412_94
               !            else
               !              aot500 = aot470_96
               !            end if
            
               !            if (aot470_94 > 0.4) then
               !              aot500 = aot412_93
               !            end if
               case default
                  print *, "ERROR: Invalid season specified: ", season
                  status = -1
                  return
            
            end select
          

         case ("Kanpur")  !------------------------------------------
            select case (season)
               !         -- winter
               case (1)
                  aot500 = aot412_96
                  if (platform .eq. 'AHI') aot500 = aot470_96
                  if (lon > 85) then !Dhaka University needs more absorbing aerosol model, longitudinal dependence for smooth transition
                     if (lon < 90) then
                        model_frac2 = 1.0-(lon-85.0)/5.0
                     else
                        model_frac2 = 0.0
                     end if

                     if (aot500 < 0.5) then
                        model_frac = 1.0
                     elseif (aot500 < 1.0) then
                        model_frac = 1.0-(aot500-0.5)/0.5*(1.0-model_frac2)
                     else
                        model_frac = model_frac2
                     end if
                     aot500 = aot412_96*model_frac+aot412_94*(1.0-model_frac)
                     if (platform .eq. 'AHI') aot500 = aot470_96*model_frac+aot470_94*(1.0-model_frac)
                  endif

               !         -- spring
               case (2)
                  aot500 = aot412_96
                  if (platform .eq. 'AHI') aot500 = aot470_96
                  if (aot500 < 0.6) then
                     model_frac = 1.0
                  elseif (aot500 < 1.2) then
                     model_frac = 1.0-(aot500-0.6)/0.6
                  else
                     model_frac = 0.0
                  end if
                  aot500 = aot412_96*model_frac+aot412_94*(1.0-model_frac)
                  if (platform .eq. 'AHI') aot500 = aot470_96*model_frac+aot470_94*(1.0-model_frac)

               !         -- summer
               case (3)
                  aot500 = aot412_96
                  if (platform .eq. 'AHI') aot500 = aot470_96
            
               !         -- fall
               case (4)
                  aot500 = aot412_96
                  if (platform .eq. 'AHI') aot500 = aot470_96
                  if (lon > 85) then
                     if (lon < 90) then
                        model_frac2 = 1.0-(lon-85.0)/5.0
                     else
                        model_frac2 = 0.0
                     end if

                     if (aot500 < 0.5) then
                        model_frac = 1.0
                     elseif (aot500 < 1.0) then
                        model_frac = 1.0-(aot500-0.5)/0.5*(1.0-model_frac2)
                     else
                        model_frac = model_frac2
                     end if
                     aot500 = aot412_96*model_frac+aot412_94*(1.0-model_frac)
                     if (platform .eq. 'AHI') aot500 = aot470_96*model_frac+aot470_94*(1.0-model_frac)
                  endif

               case default
                  print *, "ERROR: Invalid season specified: ", season
                  status = -1
                  return
            
            end select

         case ("Tinga_Tingana")  !----------------------------------------
            select case (season)

               !         -- all seasons
               case (1,2,3,4)
                  if (aot412_94 < 1.2) then
                     aot500 = aot412_995
                     if (platform .eq. 'AHI') aot500 = aot470_995
                  else
                     aot500 = aot412_94
                     if (platform .eq. 'AHI') aot500 = aot470_94
                  end if
               !            aot500 = aot412_995
               !            aot500 = aot412_94
               !            if (aot412_995 > 0.2 .AND. stdv02 < 0.002 .AND. ndvi < 0.1) then
               !              if (aot412_995 >= 0.3) then
               !                aot500 = aot412_94
               !              end if
               !              if (aot412_995 >= 0.2 .AND. aot412_995 < 0.3) then
               !                dd = (aot412_995 - 0.2) / 0.1
               !                aot500 = aot412_995 * (1.0-dd) + aot412_94*dd
               !              end if
               !            end if

               case default
                  print *, "ERROR: Invalid season specified: ", season
                  status = -1
                  return
            
            end select
      
         case ("GZ24_Only")  !----------------------------------------
            select case (season)
               !         -- winter, summer, fall
               case (1,3,4)
                  if (aot412_94 < 1.2) then
                     aot500 = aot412_995
                  else
                     aot500 = aot412_94
                  end if

               !         -- spring
               case (2)
                  if (aot412_94 < 0.5) then
                     aot500 = aot412_995
                  else
                     aot500 = aot412_94
                  endif
          
               case default
                  print *, "ERROR: Invalid season specified: ", season
                  status = -1
                  return
            end select
        
         case ("Fresno_2") !------------------------------------------------
            select case (season)

               !         -- winter
               case (1)
                  aot500 = aot470_96
                  if (aot500 < 0.5) then
                     model_frac = 1.0-(aot500-0.0)/0.5
                     aot500 = aot470_96*model_frac + aot470_94*(1.0-model_frac)
                  else
                     aot500 = aot470_94
                  endif

               !         -- spring
               case (2)
                  aot500 = aot470_995
                  if (aot500 < 0.5) then
                     model_frac = 1.0-(aot500-0.0)/0.5
                     aot500 = aot470_995*model_frac + aot470_94*(1.0-model_frac)
                  else
                     aot500 = aot470_94
                  endif

               !         -- summer
               case (3)
                  aot500 = aot470_96

               !         -- fall
               case (4)
                  aot500 = aot470_995
                  if (aot500 < 0.5) then
                     model_frac = 1.0-(aot500-0.0)/0.5
                     aot500 = aot470_995*model_frac + aot470_94*(1.0-model_frac)
                  else
                     aot500 = aot470_94
                  endif

               case default
                  print *, "ERROR: Invalid season specified: ", season
                  status = -1
                  return
            
            end select
        
         case ("Fresno_GZ18") !------------------------------------------------
            select case (season)

               !         -- winter
               case (1)
                  aot500 = aot470_96
            
               !         -- spring
               case (2)
                  aot500 = aot470_96
            
               !         -- summer
               case (3)
                  aot500 = aot470_995
                  if (aot470_995 > 0.3) then
                     aot500 = aot470_96 + aot470_995 / 2.0
                  end if
            
               !         -- fall
               case (4)
                  aot500 = aot470_96
            
               case default
                  print *, "ERROR: Invalid season specified: ", season
                  status = -1
                  return
            
            end select

         case ("CCNY") !------------------------------------------------
            select case (season)

               !         -- winter
               case (1)
                  aot500 = aot470_995
            
               !         -- spring
               case (2)
                  aot500 = aot412_995
            
               !         -- summer
               case (3)
                  aot500 = aot470_995
            
               !         -- fall
               case (4)
                  aot500 = aot470_96
            
               case default
                  print *, "ERROR: Invalid season specified: ", season
                  status = -1
                  return
            
            end select
        
         case ("Beijing")  !----------------------------------------------
            select case (season)
      
               !         -- winter
               case (1)
                  aot500 = aot470_96
                  if (aot470_96 > 0.5) then
                     aot500 = (aot470_94+aot470_92)/2.
                  end if

               !         -- spring
               case (2)
                  aot500 = aot470_96
                  if (ndvi < 0.18 .and. aot470_96 > 0.4) then
                     aot500 = (aot470_94+aot470_96)/2.
                  end if
 
                  if (aot470_96 > 0.6) then
                     aot500 = aot470_94
                  end if

               !           -- summer
               case (3)
                  aot500 = aot470_96

                  if (aot470_96 > 1.0) then
                     aot500 = (aot470_94+aot470_96)/2.
                  end if
            
               !           -- fall
               case (4)
                  aot500 = aot470_96
            
                  if (aot470_96 > 0.7) then
                     aot500 = aot470_94
                  end if

               case default
                  print *, "ERROR: Invalid season specified: ", season
                  status = -1
                  return
            
            end select
        
         case ("Hamim")  !----------------------------------------
            select case (season)

               !         -- all seasons
               case (1,2,3,4)
                  aot500 = aot412_94

               case default
                  print *, "ERROR: Invalid season specified: ", season
                  status = -1
                  return
            
            end select
      
         case ("Moldova")  !----------------------------------------
            select case (season)
               !         -- all seasons
               case (1)
                  aot500 = aot470_96
               case (2)
                  aot500 = aot412_995
                  if (platform .eq. 'AHI') aot500 = aot470_995
               case (3)
                  aot500 = aot470_96
               case (4)
                  aot500 = aot412_995
                  if (platform .eq. 'AHI') aot500 = aot470_995
               case default
                  print *, "ERROR: Invalid season specified: ", season
                  status = -1
                  return
            
            end select
        
         case ("Modena")   !----------------------------------------
            select case (season)
               case (1)
                  aot500 = aot470_96
               case (2)
                  aot500 = aot470_96
                  if (aot470_96 < 0.4) then
                     aot500 = aot470_995
                  end if
               case (3)
                  aot500 = aot470_96
                  if (aot470_96 < 0.4) then
                     aot500 = aot470_995
                  end if
               case (4)
                  aot500 = aot470_96
                  if (aot470_96 < 0.3) then
                     aot500 = aot470_995
                  end if
               case default
                  print *, "ERROR: Invalid season specified: ", season
                  status = -1
                  return
            end select
      
         case ("Ispra")    !----------------------------------------
            select case (season)
               case (1,2,3,4)
                  aot500 = aot470_96
          
               case default
                  print *, "ERROR: Invalid season specified: ", season
                  status = -1
                  return
            end select
      
         case ("Palencia") !----------------------------------------
            select case (season)
               case (1)
                  if (ndvi < 0.28) then
                     aot500 = aot412_995
                     if (platform .eq. 'AHI') aot500 = aot470_995
                  else
                     aot500 = aot470_995
                  endif
          
               case(2)
                  aot500 = aot470_995
          
               case(3)
                  aot500 = aot412_995
                  if (platform .eq. 'AHI') aot500 = aot470_995
                  if (aot412_995 > 0.2) then
                     aot500 = aot412_96
                     if (platform .eq. 'AHI') aot500 = aot470_96
                  end if
                  if (aot412_995 > 0.3) then
                     aot500 = aot412_94
                     if (platform .eq. 'AHI') aot500 = aot470_94
                  end if
            
               case (4)
                  aot500 = aot470_995
        
               case default
                  print *, "ERROR: Invalid season specified: ", season
                  status = -1
                  return
            
            end select
      
         case ("Saada")    !----------------------------------------
            select case (season)
               case (1,4)
                  aot500 = aot470_96
                  if (aot412_94 >= 0.4) then
                     aot500 = aot470_94
                  end if
                  if (aot412_94 >= 0.8) then
                     aot500 = aot470_92
                  end if
               case (2,3)
                  aot500 = aot470_96
                  if (aot412_94 >= 0.4) then
                     aot500 = aot470_94
                  end if
                  if (aot412_94 >= 0.8) then
                     aot500 = aot470_92
                  end if
                  if (stdv02 > 0.007) then
                     aot500 = aot470_995
                  end if
               case default
                  print *, "ERROR: Invalid season specified: ", season
                  status = -1
                  return
            end select
      
         case ("Solar_Village")    !----------------------------------------
            select case (season)
               case (1)
                  aot500 = aot412_93
                  if (aot412_93 > 0.5) then
                     aot500 = (aot412_93 + aot412_91) / 2.0
                  end if
            
               case (2)
                  aot500 = aot412_93
                 ! if (aot412_94 > 0.4) then
               !               aot500 = aot412_93
               !             end if
               !             if (aot412_94 > 0.8) then
               !               aot500 = (aot412_91 + aot412_93) / 2.0
               !             end if
               case (3)
                  aot500 = aot412_94
            
               case (4)
                  aot500 = aot412_96
            
               case default
                  print *, "ERROR: Invalid season specified: ", season
                  status = -1
                  return
            end select
      
         case ("Lecce_University")    !----------------------------------------
            select case (season)
               case (1)
                  aot500 = aot470_96
               case (2)
                  if (ndvi > 0.4) then
                     aot500 = aot412_94
                     if (platform .eq. 'AHI') aot500 = aot470_94
                  else
                     aot500 = (aot470_96 + aot470_995) / 2.0
                  end if
          
               case (3)
                  aot500 = aot412_995
                  if (platform .eq. 'AHI') aot500 = aot470_995
                  if (aot412_995 < 0.2) then
                     aot500 = aot412_96
                     if (platform .eq. 'AHI') aot500 = aot470_96
                  end if
            
               case (4)
                  aot500 = aot470_96
          
               case default
                  print *, "ERROR: Invalid season specified: ", season
                  status = -1
                  return
            end select
      
         case ("Carpentras")    !----------------------------------------
            select case (season)
               case (1)
                  !aot500 = aot470_96
                  aot500 = aot412_93
                  if (platform .eq. 'AHI') aot500 = aot470_93
               case (3,4)
                  !aot500 = aot470_96
                  aot500 = aot412_995
                  if (platform .eq. 'AHI') aot500 = aot470_995
               case (2)
                  if(ndvi < 0.3) then
                     !aot500 = aot412_94
                     aot500 = aot412_96
                     if (platform .eq. 'AHI') aot500 = aot470_96
                  else
                     aot500 = aot470_96
                  end if
               case default
                  print *, "ERROR: Invalid season specified: ", season
                  status = -1
                  return
            end select
        
         case ("Trelew")  !----------------------------------------
            select case (season)

               !         -- all seasons
               case(1)
                  aot500 = aot470_96
               case (2,3)
                  aot500 = aot470_96
               case (4)
                  if (ndvi < 0.2) then
                     aot500 = aot412_94
                  else
                     aot500 = aot470_96
                  end if
            
               case default
                  print *, "ERROR: Invalid season specified: ", season
                  status = -1
                  return
            
            end select
        
         case ("Pune")  !----------------------------------------
            select case (season)

               !         -- winter
               case (1)
                  aot500 = aot412_96
                  if (platform .eq. 'AHI') aot500 = aot470_96
                  if (aot500 < 0.5) then
                     model_frac = 1.0
                  elseif (aot500 < 1.0) then
                     model_frac = 1.0-(aot500-0.5)/0.5
                  else
                     model_frac = 0.0
                  end if
                  aot500 = aot412_96*model_frac+aot412_94*(1.0-model_frac)
                  if (platform .eq. 'AHI') aot500 = aot470_96*model_frac+aot470_94*(1.0-model_frac)
               !         -- spring
               case (2)
                  aot500 = aot412_96
                  if (platform .eq. 'AHI') aot500 = aot470_96
                  if (aot500 < 0.6) then
                     model_frac = 1.0
                  elseif (aot500 < 1.2) then
                     model_frac = 1.0-(aot500-0.6)/0.6
                  else
                     model_frac = 0.0
                  end if
                  aot500 = aot412_96*model_frac+aot412_94*(1.0-model_frac)
                  if (platform .eq. 'AHI') aot500 = aot470_96*model_frac+aot470_94*(1.0-model_frac)
               !         -- summer
               case (3)
                  aot500 = aot412_96
                  if (platform .eq. 'AHI') aot500 = aot470_96
               !         -- fall
               case (4)
                  aot500 = aot412_96
                  if (platform .eq. 'AHI') aot500 = aot470_96
                  if (lon > 85) then
                     if (lon < 90) then
                        model_frac2 = 1.0-(lon-85.0)/5.0
                     else
                        model_frac2 = 0.0
                     end if

                     if (aot500 < 0.5) then
                        model_frac = 1.0
                     elseif (aot500 < 1.0) then
                        model_frac = 1.0-(aot500-0.5)/0.5*(1.0-model_frac2)
                     else
                        model_frac = model_frac2
                     end if
                     aot500 = aot412_96*model_frac+aot412_94*(1.0-model_frac)
                     if (platform .eq. 'AHI') aot500 = aot470_96*model_frac+aot470_94*(1.0-model_frac)
                  endif
           
               case default
                  print *, "ERROR: Invalid season specified: ", season
                  status = -1
                  return
            
            end select

         case ("Evora")  !----------------------------------------
            select case (season)
               !         -- winter
               case (1)
                  aot500 = aot470_96

               !         -- spring
               case (2)
                  aot500 = aot412_94
                  if (platform .eq. 'AHI') aot500 = aot470_94
               !         -- summer
               case (3)
                  aot500 = aot470_995
            
               !         -- fall
               case (4)
                  aot500 = aot412_995
                  if (platform .eq. 'AHI') aot500 = aot470_995

               case default
                  print *, "ERROR: Invalid season specified: ", season
                  status = -1
                  return
            
            end select
         case ("Blida")  !----------------------------------------
            select case (season)
               !         -- winter
               case (1)
                  aot500 = aot470_96

               !         -- spring
               case (2)
                  aot500 = aot470_96
            
               !         -- summer
               case (3)
                  aot500 = aot470_96
            
               !         -- fall
               case (4)
                  aot500 = aot470_96

               case default
                  print *, "ERROR: Invalid season specified: ", season
                  status = -1
                  return
            
            end select
         case ("Blida_High")  !----------------------------------------
            select case (season)
               !         -- winter
               case (1)
                  aot500 = aot470_96

               !         -- spring
               case (2)
                  aot500 = aot470_96
            
               !         -- summer
               case (3)
                  aot500 = aot470_96
            
               !         -- fall
               case (4)
                  aot500 = aot470_96

               case default
                  print *, "ERROR: Invalid season specified: ", season
                  status = -1
                  return
            
            end select
         case ("Ilorin")  !------------------------------------------
            select case (season)

               !         -- winter
               case (1)
                  aot412_92 = (aot412_91+aot412_93)/2.0
                  aot500 = aot412_92

               !            aot500 = aot412_94
               !            if (aot500 < 0.5) then
               !              model_frac = 1.0
               !            elseif (aot500 < 1.2) then
               !              model_frac = 1.0-(aot500-0.5)/0.7
               !            else
               !              model_frac = 0.0
               !            end if
               !            aot500 = aot412_94*model_frac+aot412_91*(1.0-model_frac)

               !            if (aot412_94 > 0.6) then
               !              aot500 = aot412_91
               !            else if (aot412_94 > 0.4) then
               !              aot500 = aot412_94
               !            else
               !              aot500 = aot412_96
               !            endif
            
               !         -- spring
               case (2)
                  aot500 = aot412_94
                  if (aot500 < 0.5) then
                     model_frac = 1.0
                  elseif (aot500 < 1.0) then
                     model_frac = 1.0-(aot500-0.5)/0.5
                  else
                     model_frac = 0.0
                  end if
                  aot500 = aot412_94*model_frac+aot412_93*(1.0-model_frac)

               !            aot500 = aot412_94
               !            if (aot412_94 > 1.0) then
               !              aot500 = aot412_93
               !            end if
               !            if (aot412_94 > 1.5) then
               !              aot500 = aot412_91
               !            end if
            
               !         -- summer
               case (3)
                  aot500 = aot412_96

               !         -- fall
               case (4)
                  aot500 = aot412_94
               !            aot500 = aot412_94
               !            if (aot412_94 > 0.5) then
               !              aot500 = aot412_93
               !            end if
               !            if (aot412_94 > 0.8) then
               !              aot500 = aot412_91
               !            end if
               case default
                  print *, "ERROR: Invalid season specified: ", season
                  status = -1
                  return
            
            end select
         case ("Ilorin_Transition")  !------------------------------------------
            select case (season)

               !         -- winter
               case (1)
                  aot500 = aot412_96
                  if (aot500 < 0.5) then
                     model_frac = 1.0
                  elseif (aot500 < 1.0) then
                     model_frac = 1.0-(aot500-0.5)/0.5
                  else
                     model_frac = 0.0
                  end if
                  aot500 = aot412_96*model_frac+aot412_94*(1.0-model_frac)

               !            if (aot412_94 > 0.5) then
               !              aot500 = aot412_94
               !            else
               !              aot500 = aot412_96
               !            end if

               !         -- spring
               case (2)
                  aot500 = aot412_995
                  if (aot412_94 < 0.5) then
                     model_frac = 1.0
                  elseif (aot412_94 < 1.0) then
                     model_frac = 1.0-(aot500-0.5)/0.5
                  else
                     model_frac = 0.0
                  end if
                  aot500 = aot412_995*model_frac+aot412_96*(1.0-model_frac)

               !            aot500 = aot412_995
               !            if (aot412_94 > 1.0) then
               !              aot500 = aot412_96
               !            end if
               !            if (aot412_94 > 1.5) then
               !              aot500 = aot412_96
               !            end if
            
               !         -- summer
               case (3)
                  aot500 = aot412_995

               !         -- fall
               case (4)
                  aot500 = aot412_995
                  if (aot412_94 < 0.5) then
                     model_frac = 1.0
                  elseif (aot412_94 < 1.0) then
                     model_frac = 1.0-(aot500-0.5)/0.5
                  else
                     model_frac = 0.0
                  end if
                  aot500 = aot412_995*model_frac+aot412_96*(1.0-model_frac)

               !            aot500 = aot412_995
               !            if (aot412_94 > 0.5) then
               !              aot500 = aot412_96
               !            end if
               !            if (aot412_94 > 0.8) then
               !              aot500 = aot412_94
               !            end if
               case default
                  print *, "ERROR: Invalid season specified: ", season
                  status = -1
                  return
            
            end select
        
         case ("SACOL") !----------------------------------------
            select case (season)
               case (1,2,3,4)
                  aot500 = aot470_995

               case default
                  print *, "ERROR: Invalid season specified: ", season
                  status = -1
                  return
            
            end select
      
         case ("Mexico_City") !----------------------------------------
            select case (season)
               case (1)
                  aot500 = aot412_96
          
               case (2)
                  aot500 = aot412_995
                        
               case (3)
                  aot500 = aot470_96

               case (4)
                  aot500 = aot412_96

               case default
                  print *, "ERROR: Invalid season specified: ", season
                  status = -1
                  return
            
            end select
        
         case ("Jaipur") !----------------------------------------
            select case (season)
               case (1)
                  aot500 = aot412_96
                  if (platform .eq. 'AHI') aot500 = aot470_96
          
               case (2)
                  aot500 = aot412_96
                  if (platform .eq. 'AHI') aot500 = aot470_96
            
               case(3)
                  aot500 = aot412_96
                  if (platform .eq. 'AHI') aot500 = aot470_96

               case (4)
                  aot500 = aot412_96
                  if (platform .eq. 'AHI') aot500 = aot470_96
            
               case default
                  print *, "ERROR: Invalid season specified: ", season
                  status = -1
                  return
            
            end select
        
         case ("NW_India_Desert") !----------------------------------------
            select case (season)
               case (1)
                  aot500 = aot412_94
                  if (platform .eq. 'AHI') aot500 = aot470_94
          
               case (2)
                  aot500 = aot412_94
                  if (platform .eq. 'AHI') aot500 = aot470_94
                        
               case(3)
                  aot500 = aot412_94
                  if (platform .eq. 'AHI') aot500 = aot470_94

               case (4)
                  aot500 = aot412_96
                  if (platform .eq. 'AHI') aot500 = aot470_96
            
               case default
                  print *, "ERROR: Invalid season specified: ", season
                  status = -1
                  return
            
            end select

         case ("Yuma") !----------------------------------------
            select case (season)
               case (1)
                  aot500 = aot412_995
          
               case (2)
                  aot500 = aot470_995
                        
               case (3)
                  aot500 = aot412_995

               case (4)
                  aot500 = aot412_995
            
               case default
                  print *, "ERROR: Invalid season specified: ", season
                  status = -1
                  return
            end select

         case default  !--------------------------------------------------
            print *, "ERROR: No data for AERONET site specified: ", aero_site
            status = -1
            aot500 = -999.0
            return

      end select      ! -- AERONET site
    
      if (dflag) then
         print *, trim(func_name), ', AERONET site: ', aero_site
         print *, trim(func_name), ', season, ndvi, stdv02: ', season, ndvi, stdv02
         print '(A,A,12(F10.4,1x))', trim(func_name), ', aot values: ', aot412_91, aot412_93, aot412_94,  &
            &                                   aot412_96, aot412_995, aot470_91, aot470_92, aot470_93, &
            &                                   aot470_94, aot470_95, aot470_96, aot470_995
         print *, trim(func_name), ', ae: ', ae
         print *, trim(func_name), ', aot500: ', aot500
      end if
      return
    
   end function get_aeronet_aot500

  
   ! --  Find granule limits and set LER offsets, and allocate space for the tables
   integer function set_limits(locedge, lat, long) RESULT(status)

      integer, intent(in)  :: locedge(2)
      real, intent(in)     :: lat(locedge(1),locedge(2)), long(locedge(1),locedge(2))

      integer              :: checkvariable
      integer              :: i, j
      character (len=256)  :: msg
      real                 :: eastedge, westedge
      status = 0
    
      !   -- if processing extracts, if region of interest is within bounds of the granule,
      !   -- but not actually contained in the granule quadralateral, lat and lon values may be
      !   -- set to fill values (-999). Check this condition and fail.
      !   -- @TODO - utilize the status variable and check where called.
      !    if (minval(lat) < -900.0 .OR. minval(long) < -900.0) then
      !      msg = "Min lat or min lon is fill value. Failing."
      !      call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'set_limits')
      !    end if
    
      eastedge = -999.0
      westedge =  999.0
      dateline = 0
      if (minval(long, long > -900.0) < -175.0 .and. maxval(long, long > -900.0) > 175.0) then
         eastedge = 180.0
         westedge = -180.0
         do i=1, locedge(1)
            do j=1, locedge(2)
               if (long(i,j) < -900.0) cycle    ! skip undefined
               if (long(i,j) > 0.0 .and. long(i,j) < eastedge) eastedge = long(i,j)
               if (long(i,j) < 0.0 .and. long(i,j) > westedge) westedge = long(i,j)
            enddo
         enddo
         LERstart(1) = 10*(180+floor(eastedge)-1)
         if (LERstart(1) .le. 0) LERstart(1) = 1
         dateline = 3600 - LERstart(1)
         LERedge(1) = 10*(180+(floor(westedge)+2)) + dateline
      else
         LERstart(1) = 10*(180+(floor(minval(long, long > -900.0))-1))
         if (LERstart(1) .le. 0) LERstart(1) = 1
         LERedge(1) = 10*(180+(floor(maxval(long, long > -900.0))+2)) - LERstart(1)
         if (LERedge(1)+LERstart(1) > 3600) LERedge(1) = 3600 - LERstart(1)
      endif
    
      LERstart(2) = 10*(90+(floor(minval(lat, lat > -900.0))-1))
      if (LERstart(2) .le. 0) LERstart(2) = 1
      LERedge(2) = 10*(90+(floor(maxval(lat, lat > -900.0))+2)) - LERstart(2)
      if (LERedge(2)+LERstart(2) > 1800) LERedge(2) = 1800 - LERstart(2)
        
      if (allocated(gref412_all)) deallocate(gref412_all)
      allocate (gref412_all(LERedge(1),LERedge(2)), stat = checkvariable)
      if ( checkvariable /= 0 ) goto 90

      if (allocated(gref470_all)) deallocate(gref470_all)
      allocate (gref470_all(LERedge(1),LERedge(2)), stat = checkvariable)
      if ( checkvariable /= 0 ) goto 90

      if (allocated(gref650_all)) deallocate(gref650_all)
      allocate (gref650_all(LERedge(1),LERedge(2)), stat = checkvariable)
      if ( checkvariable /= 0 ) goto 90

      if (allocated(gref412_fwd)) deallocate(gref412_fwd)
      allocate (gref412_fwd(LERedge(1),LERedge(2)), stat = checkvariable)
      if ( checkvariable /= 0 ) goto 90

      if (allocated(gref470_fwd)) deallocate(gref470_fwd)
      allocate (gref470_fwd(LERedge(1),LERedge(2)), stat = checkvariable)
      if ( checkvariable /= 0 ) goto 90

      if (allocated(gref650_fwd)) deallocate(gref650_fwd)
      allocate (gref650_fwd(LERedge(1),LERedge(2)), stat = checkvariable)
      if ( checkvariable /= 0 ) goto 90

      !if (allocated(gref412_bkd)) deallocate(gref412_bkd)
      !allocate (gref412_bkd(LERedge(1),LERedge(2)), stat = checkvariable)
      !if ( checkvariable /= 0 ) goto 90

      !if (allocated(gref470_bkd)) deallocate(gref470_bkd)
      !allocate (gref470_bkd(LERedge(1),LERedge(2)), stat = checkvariable)
      !if ( checkvariable /= 0 ) goto 90

      !if (allocated(gref650_bkd)) deallocate(gref650_bkd)
      !allocate (gref650_bkd(LERedge(1),LERedge(2)), stat = checkvariable)
      !if ( checkvariable /= 0 ) goto 90

      if (allocated(gref865_all)) deallocate(gref865_all)
      allocate (gref865_all(LERedge(1),LERedge(2)), stat = checkvariable)
      if ( checkvariable /= 0 ) goto 90

      !if (allocated(coefs412_all_tp)) deallocate(coefs412_all_tp)
      !allocate (coefs412_all_tp(LERedge(1),LERedge(2),4,3,2), stat = checkvariable)
      !if ( checkvariable /= 0 ) goto 90

      !if (allocated(coefs412_fwd_tp)) deallocate(coefs412_fwd_tp)
      !allocate (coefs412_fwd_tp(LERedge(1),LERedge(2),4,3,2), stat = checkvariable)
      !if ( checkvariable /= 0 ) goto 90

      if (allocated(coefs412_all)) deallocate(coefs412_all)
      allocate (coefs412_all(LERedge(1),LERedge(2),4,3), stat = checkvariable)
      if ( checkvariable /= 0 ) goto 90

      if (allocated(coefs412_fwd)) deallocate(coefs412_fwd)
      allocate (coefs412_fwd(LERedge(1),LERedge(2),4,3), stat = checkvariable)
      if ( checkvariable /= 0 ) goto 90

      !if (allocated(coefs470_all_tp)) deallocate(coefs470_all_tp)
      !allocate (coefs470_all_tp(LERedge(1),LERedge(2),4,3,2), stat = checkvariable)
      !if ( checkvariable /= 0 ) goto 90

      !if (allocated(coefs470_fwd_tp)) deallocate(coefs470_fwd_tp)
      !allocate (coefs470_fwd_tp(LERedge(1),LERedge(2),4,3,2), stat = checkvariable)
      !if ( checkvariable /= 0 ) goto 90

      if (allocated(coefs470_all)) deallocate(coefs470_all)
      allocate (coefs470_all(LERedge(1),LERedge(2),4,3), stat = checkvariable)
      if ( checkvariable /= 0 ) goto 90

      if (allocated(coefs470_fwd)) deallocate(coefs470_fwd)
      allocate (coefs470_fwd(LERedge(1),LERedge(2),4,3), stat = checkvariable)
      if ( checkvariable /= 0 ) goto 90

      !if (allocated(coefs650_all_tp)) deallocate(coefs650_all_tp)
      !allocate (coefs650_all_tp(LERedge(1),LERedge(2),4,3,2), stat = checkvariable)
      !if ( checkvariable /= 0 ) goto 90

      !if (allocated(coefs650_fwd_tp)) deallocate(coefs650_fwd_tp)
      !allocate (coefs650_fwd_tp(LERedge(1),LERedge(2),4,3,2), stat = checkvariable)
      !if ( checkvariable /= 0 ) goto 90

      if (allocated(coefs650_all)) deallocate(coefs650_all)
      allocate (coefs650_all(LERedge(1),LERedge(2),4,3), stat = checkvariable)
      if ( checkvariable /= 0 ) goto 90

      if (allocated(coefs650_fwd)) deallocate(coefs650_fwd)
      allocate (coefs650_fwd(LERedge(1),LERedge(2),4,3), stat = checkvariable)
      if ( checkvariable /= 0 ) goto 90
    
      !   -- allocate arrays for VIIRS, all-angle surface database
      if (allocated(vgref412_all)) deallocate(vgref412_all)
      allocate (vgref412_all(LERedge(1),LERedge(2)), stat = checkvariable)
      if ( checkvariable /= 0 ) goto 90

      if (allocated(vgref488_all)) deallocate(vgref488_all)
      allocate (vgref488_all(LERedge(1),LERedge(2)), stat = checkvariable)
      if ( checkvariable /= 0 ) goto 90

      if (allocated(vgref670_all)) deallocate(vgref670_all)
      allocate (vgref670_all(LERedge(1),LERedge(2)), stat = checkvariable)
      if ( checkvariable /= 0 ) goto 90

      goto 100
    
90 continue
   print *, "ERROR: Unable to allocate coefficient array: ", status
   return
    
100 continue
    return

 end function set_limits

 ! --  Find granule limits and set offsets for 2.2 um surface database, and allocate space for the
 ! tables, 0.06 degree resolution
 integer function set_limits6(locedge, lat, long) RESULT(status)

    integer, intent(in)  :: locedge(2)
    real, intent(in)     :: lat(locedge(1),locedge(2)), long(locedge(1),locedge(2))

    integer              :: checkvariable
    integer              :: i, j
    character (len=256)  :: msg
    real                 :: eastedge, westedge
    status = 0

    !   -- if processing extracts, if region of interest is within bounds of the
    !   granule,
    !   -- but not actually contained in the granule quadralateral, lat and lon
    !   values may be
    !   -- set to fill values (-999). Check this condition and fail.
    !   -- @TODO - utilize the status variable and check where called.
    !    if (minval(lat) < -900.0 .OR. minval(long) < -900.0) then
    !      msg = "Min lat or min lon is fill value. Failing."
    !      call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'set_limits')
    !    end if

    eastedge = -999.0
    westedge =  999.0
    dateline6 = 0
    if (minval(long, long > -900.0) < -175.0 .and. maxval(long, long > -900.0) > 175.0) then
       eastedge = 180.0
       westedge = -180.0
       do i=1, locedge(1)
          do j=1, locedge(2)
             if (long(i,j) < -900.0) cycle    ! skip undefined
             if (long(i,j) > 0.0 .and. long(i,j) < eastedge) eastedge = long(i,j)
             if (long(i,j) < 0.0 .and. long(i,j) > westedge) westedge = long(i,j)
          enddo
       enddo
       LERstart6(1) = (180+(floor(eastedge)-1))/0.06
       if (LERstart6(1) .le. 0) LERstart6(1) = 1
       dateline6 = 6000 - LERstart6(1)
       LERedge6(1) = (180+(floor(westedge)+2))/0.06 + dateline6
    else
       LERstart6(1) = (180+(floor(minval(long, long > -900.0))-1))/0.06
       if (LERstart6(1) .le. 0) LERstart6(1) = 1
       LERedge6(1) = (180+(floor(maxval(long, long > -900.0))+2))/0.06 - LERstart6(1)
       if (LERedge6(1)+LERstart6(1) > 6000) LERedge6(1) = 6000 - LERstart6(1)
    endif

    LERstart6(2) = (90+(floor(minval(lat, lat > -900.0))-1))/0.06
    if (LERstart6(2) .le. 0) LERstart6(2) = 1
    LERedge6(2) = (90+(floor(maxval(lat, lat > -900.0))+2))/0.06 - LERstart6(2)
    if (LERedge6(2)+LERstart6(2) > 3000) LERedge6(2) = 3000 - LERstart6(2)

    if (allocated(swir_coeffs412)) deallocate(swir_coeffs412)
    allocate (swir_coeffs412(LERedge6(1),LERedge6(2),3), stat = checkvariable)
    if ( checkvariable /= 0 ) goto 90

    if (allocated(swir_coeffs470)) deallocate(swir_coeffs470)
    allocate (swir_coeffs470(LERedge6(1),LERedge6(2),3), stat = checkvariable)
    if ( checkvariable /= 0 ) goto 90

    if (allocated(swir_stderr412)) deallocate(swir_stderr412)
    allocate (swir_stderr412(LERedge6(1),LERedge6(2)), stat = checkvariable)
    if ( checkvariable /= 0 ) goto 90

    if (allocated(swir_stderr470)) deallocate(swir_stderr470)
    allocate (swir_stderr470(LERedge6(1),LERedge6(2)), stat = checkvariable)
    if ( checkvariable /= 0 ) goto 90

    if (allocated(swir_min)) deallocate(swir_min)
    allocate (swir_min(LERedge6(1),LERedge6(2)), stat = checkvariable)
    if ( checkvariable /= 0 ) goto 90

    if (allocated(swir_max)) deallocate(swir_max)
    allocate (swir_max(LERedge6(1),LERedge6(2)), stat = checkvariable)
    if ( checkvariable /= 0 ) goto 90

    goto 100

90 continue
   print *, "ERROR: Unable to allocate 2.2 um surface database array: ", status
   return

100 continue
    return

 end function set_limits6

 ! -- Load surface LER coefficient tables.
 integer function load_hdfLER(lut_file, season) RESULT(status)

    !   include 'hdf.f90'
    !   include 'dffunc.f90'
    use netcdf
    USE OCIUAAER_Config_Module

    character(len=255), intent(in)    ::  lut_file
    integer, intent(in)               ::  season

    integer             :: start2(3), stride2(3), edges2(3)
    integer             :: start4(5), stride4(5), edge4(5)
      
    character(len=255)    ::  sds_name
    character(len=255)    ::  dset_name
    character(len=255)    ::  test_name
    character(len=255)    ::  group_name

    integer               ::  nc_id
    integer               ::  dim_id
    integer               ::  dset_id
    integer               ::  grp_id
    integer             ::  sd_id, sds_index, sds_id

    start2 = (/LERstart(1),LERstart(2),season/)
    edges2 = (/LERedge(1),LERedge(2),1/)
    stride2 = (/1,1,1/)
      
    start4 = (/LERstart(1),LERstart(2),1,1,season/)
    edge4 = (/LERedge(1),LERedge(2),4,3,1/)
    stride4 = (/1,1,1,1,1/)

    test_name = trim(lut_file)
    status = nf90_open(test_name, nf90_nowrite, nc_id)
    if (status /= NF90_NOERR) then
       print *, "ERROR: Failed to open deepblue lut_nc4 file: ", status
       return
    end if

    group_name = 'surface_coeffs'
    status = nf90_inq_ncid(nc_id, group_name, grp_id)
    if (status /= NF90_NOERR) then
       print *, "ERROR: Failed to get ID of group "//trim(group_name)//": ", status
       return
    end if

    sds_name = "412_fwd"
    status = readLER5(start4, edge4, stride4, sds_name, grp_id, coefs412_fwd)
    if (status < 0) goto 90

    sds_name = "412_all"
    status = readLER5(start4, edge4, stride4, sds_name, grp_id, coefs412_all)
    if (status < 0) goto 90

    sds_name = "470_fwd"
    status = readLER5(start4, edge4, stride4, sds_name, grp_id, coefs470_fwd)
    if (status < 0) goto 90

    sds_name = "470_all"
    status = readLER5(start4, edge4, stride4, sds_name, grp_id, coefs470_all)
    if (status < 0) goto 90

    sds_name = "650_fwd"
    status = readLER5(start4, edge4, stride4, sds_name, grp_id, coefs650_fwd)
    if (status < 0) goto 90

    sds_name = "650_all"
    status = readLER5(start4, edge4, stride4, sds_name, grp_id, coefs650_all)
    if (status < 0) goto 90

    status = nf90_close(nc_id)
    if (status /= NF90_NOERR) then
       print *, "ERROR: Failed to close lut_nc4 file: ", status
       return
    end if

    !   -- MODIS, all-angle surface database
    status = nf90_open(trim(lut_file), nf90_nowrite, nc_id)
    if (status /= NF90_NOERR) then
       print *, "ERROR: Failed to open deepblue lut_nc4 file: ", status
       return
    end if

    group_name = 'modis_surface'
    status = nf90_inq_ncid(nc_id, group_name, grp_id)
    if (status /= NF90_NOERR) then
       print *, "ERROR: Failed to get ID of group "//trim(group_name)//": ", status
       return
    end if

    sds_name = "412_all"
    status = readLER2(start2, edges2, stride2, sds_name, grp_id, gref412_all)
    if (status < 0) goto 90

    sds_name = "412_fwd"
    status = readLER2(start2, edges2, stride2, sds_name, grp_id, gref412_fwd)
    if (status < 0) goto 90

    sds_name = "470_all"
    status = readLER2(start2, edges2, stride2, sds_name, grp_id, gref470_all)
    if (status < 0) goto 90

    sds_name = "470_fwd"
    status = readLER2(start2, edges2, stride2, sds_name, grp_id, gref470_fwd)
    if (status < 0) goto 90

    sds_name = "650_all"
    status = readLER2(start2, edges2, stride2, sds_name, grp_id, gref650_all)
    if (status < 0) goto 90

    sds_name = "650_fwd"
    status = readLER2(start2, edges2, stride2, sds_name, grp_id, gref650_fwd)
    if (status < 0) goto 90

    sds_name = "865_all_all"
    status = readLER2(start2, edges2, stride2, sds_name, grp_id, gref865_all)
    if (status < 0) goto 90

    status = nf90_close(nc_id)
    if (status /= NF90_NOERR) then
       print *, "ERROR: Failed to close lut_nc4 file: ", status
       return
    end if

    !   -- VIIRS, all-angle surface database
    status = nf90_open(trim(lut_file), nf90_nowrite, nc_id)
    if (status /= NF90_NOERR) then
       print *, "ERROR: Failed to open deepblue lut_nc4 file: ", status
       return
    end if

    group_name = 'viirs_surface'
    status = nf90_inq_ncid(nc_id, group_name, grp_id)
    if (status /= NF90_NOERR) then
       print *, "ERROR: Failed to get ID of group "//trim(group_name)//": ", status
       return
    end if

    sds_name = "412_all"
    status = readLER2(start2, edges2, stride2, sds_name, grp_id, vgref412_all)
    if (status < 0) goto 90

    sds_name = "488_all"
    status = readLER2(start2, edges2, stride2, sds_name, grp_id, vgref488_all)
    if (status < 0) goto 90
    
    sds_name = "670_all"
    status = readLER2(start2, edges2, stride2, sds_name, grp_id, vgref670_all)
    if (status < 0) goto 90
    
    status = nf90_close(nc_id)
    if (status /= NF90_NOERR) then
       print *, "ERROR: Failed to close lut_nc4 file: ", status
       return
    end if
    
    ! Close up shop and go home
    goto 100

90 continue
   print *, "Error reading "//trim(sds_name)//" from file "//trim(lut_file)
   return
100 continue

 end function load_hdfLER

 ! -- Load 2.2 um surface database
 integer function load_swir_coeffs(lut_file, season) result(status) !jlee added 05/16/2017
    use netcdf
    USE OCIUAAER_Config_Module
    implicit none

    character(len=255), intent(in)  ::  lut_file
    integer, intent(in)   ::  season

    ! HDF vars
    character(len=255)    ::  sds_name
    character(len=255)    ::  dset_name
    character(len=255)    ::  attr_name
    character(len=255)    ::  group_name

    integer               ::  nc_id
    integer               ::  dim_id
    integer               ::  dset_id
    integer               ::  grp_id
    integer               ::  sd_id, sds_index, sds_id
    integer, dimension(3) ::  start2, stride2, edges2
    integer, dimension(4) ::  start3, stride3, edges3

    status = -1

    start2 = (/LERstart6(1),LERstart6(2),season/)
    edges2 =  (/LERedge6(1),LERedge6(2),1/)
    stride2 =(/1,1,1/)

    start3 = (/LERstart6(1),LERstart6(2),1,season/)
    edges3 =  (/LERedge6(1),LERedge6(2),3,1/)
    stride3 =(/1,1,1,1/)

    status = nf90_open(cfg%db_nc4, nf90_nowrite, nc_id)
    if (status /= NF90_NOERR) then
       print *, "ERROR: Failed to open deepblue lut_nc4 file: ", status
       return
    end if

    group_name = 'swir_vis_surface_coeffs'
    status = nf90_inq_ncid(nc_id, group_name, grp_id)
    if (status /= NF90_NOERR) then
       print *, "ERROR: Failed to get ID of group "//trim(group_name)//": ", status
       return
    end if

    sds_name = 'coeffs_2250_to_412'
    status = readSWIR3(start3, edges3, stride3, sds_name, grp_id, swir_coeffs412)
    if (status < 0) goto 90
    !============================================================================
    sds_name = 'stderr_412'
    status = readSWIR2(start2, edges2, stride2, sds_name, grp_id, swir_stderr412)
    if (status < 0) goto 90
    !============================================================================
    sds_name = 'coeffs_2250_to_488'
    status = readSWIR3(start3, edges3, stride3, sds_name, grp_id, swir_coeffs470)
    if (status < 0) goto 90
    !============================================================================
    sds_name = 'stderr_488'
    status = readSWIR2(start2, edges2, stride2, sds_name, grp_id, swir_stderr470)
    if (status < 0) goto 90
    !============================================================================
    sds_name = 'min_2250_for_488'
    status = readSWIR2(start2, edges2, stride2, sds_name, grp_id, swir_min)
    if (status < 0) goto 90
    !============================================================================
    sds_name = 'max_2250_for_488'
    status = readSWIR2(start2, edges2, stride2, sds_name, grp_id, swir_max)
    if (status < 0) goto 90
    !============================================================================

    status = nf90_close(nc_id)
    if (status /= NF90_NOERR) then
       print *, "ERROR: Failed to close lut_nc4 file: ", status
       return
    end if

    goto 100

90 continue
   print *, "Error reading "//trim(sds_name)//" from file "//trim(lut_file)
   return

100 continue
  
 end function load_swir_coeffs
   

 ! Retrieve swir vs. vis coefficients
 function get_swir_coeffs412(latidx, lonidx)
    integer, intent(in) :: latidx, lonidx
    integer :: i,j
    real    :: get_swir_coeffs412(3)

    if (dateline6 .eq. 0 .OR. lonidx .gt. LERstart6(1)) then
       i = lonidx - LERstart6(1)
    else
       i = lonidx + dateline6
    end if
    j = latidx - LERstart6(2)

    get_swir_coeffs412 = swir_coeffs412(i,j,:)

 end function

 ! Retrieve swir vs. vis coefficients
 function get_swir_coeffs470(latidx, lonidx)
    integer, intent(in) :: latidx, lonidx
    integer :: i,j
    real    :: get_swir_coeffs470(3)

    if (dateline6 .eq. 0 .OR. lonidx .gt. LERstart6(1)) then
       i = lonidx - LERstart6(1)
    else
       i = lonidx + dateline6
    end if
    j = latidx - LERstart6(2)

    get_swir_coeffs470 = swir_coeffs470(i,j,:)

 end function

 ! Retrieve swir vs. vis stderr
 real function get_swir_stderr412(latidx, lonidx) result(stderr)
    integer, intent(in) :: latidx, lonidx
    integer :: i,j

    if (dateline6 .eq. 0 .OR. lonidx .gt. LERstart6(1)) then
       i = lonidx - LERstart6(1)
    else
       i = lonidx + dateline6
    end if
    j = latidx - LERstart6(2)

    stderr = swir_stderr412(i,j)

 end function get_swir_stderr412

 ! Retrieve swir vs. vis stderr
 real function get_swir_stderr470(latidx, lonidx) result(stderr)
    integer, intent(in) :: latidx, lonidx
    integer :: i,j

    if (dateline6 .eq. 0 .OR. lonidx .gt. LERstart6(1)) then
       i = lonidx - LERstart6(1)
    else
       i = lonidx + dateline6
    end if
    j = latidx - LERstart6(2)

    stderr = swir_stderr470(i,j)

 end function get_swir_stderr470

 ! Retrieve swir vs. vis min-max range
 function get_swir_range(latidx, lonidx)
    integer, intent(in) :: latidx, lonidx
    integer :: i,j
    real    :: get_swir_range(2)

    if (dateline6 .eq. 0 .OR. lonidx .gt. LERstart6(1)) then
       i = lonidx - LERstart6(1)
    else
       i = lonidx + dateline6
    end if
    j = latidx - LERstart6(2)

    get_swir_range = (/swir_min(i,j),swir_max(i,j)/)

 end function

 ! Retrieve LER865 value
 real function get_LER865(latidx, lonidx) RESULT(LER)
    integer, intent(in) :: latidx, lonidx

    integer :: i,j

    if (dateline .eq. 0 .OR. lonidx .gt. LERstart(1)) then
       i = lonidx - LERstart(1)
    else
       i = lonidx + dateline
    end if
    j = latidx - LERstart(2)

    LER = gref865_all(i,j)

 end function get_LER865

 ! -- Retrieve LER412 value from surface reflectivity coefficient table for pixel at
 ! --  latidx, lonidx in table.
 real function get_LER412(latidx, lonidx, NDVI, scatAngle, relaz,min_flag) RESULT(LER)
    integer, intent(in) :: latidx, lonidx
    real, intent(in)    :: NDVI, scatAngle, relaz
    integer, intent(inout) :: min_flag

    integer :: i,j,nidx
    real    :: coefs(4), acoefs(4), ncoefs(4), mcoefs(4), tLER

    LER = -999.0
    if (dateline .eq. 0 .OR. lonidx .gt. LERstart(1)) then
       i = lonidx - LERstart(1)
    else
       i = lonidx + dateline
    end if
    j = latidx - LERstart(2)

    !    print *,"in get_LER412"
    !    print *,"latidx, lonidx = ",latidx, lonidx
    !    print *,"i,j = ",i,j
    !    print *,"LERstart = ",LERstart
    !    print *,"LERedge = ",LERedge
    
    if (NDVI < NDVI1_CUTOFF) then
       nidx = 1
    elseif (NDVI < NDVI2_CUTOFF) then
       nidx = 2
    else
       nidx = 3
    endif

    if (relaz < 90.0) then
       !ncoefs(:) = coefs412_fwd_tp(i,j,:,nidx,1)
       !acoefs(:) = coefs412_fwd_tp(i,j,:,1,1)
       !mcoefs(:) = coefs412_fwd_tp(i,j,:,1,2)
       ncoefs(:) = coefs412_fwd(i,j,:,nidx)
       !acoefs(:) = coefs412_fwd(i,j,:,1,1)
       !mcoefs(:) = coefs412_fwd(i,j,:,1,2)
       tLER = vgref412_all(i,j)!if there is no BRDF value, use VIIRS min ref (Jul 2017 W.KIM)
    else
       !ncoefs(:) = coefs412_all_tp(i,j,:,nidx,1)
       !acoefs(:) = coefs412_all_tp(i,j,:,1,1)
       !mcoefs(:) = coefs412_all_tp(i,j,:,1,2)
       ncoefs(:) = coefs412_all(i,j,:,nidx)
       !acoefs(:) = coefs412_all(i,j,:,1,1)
       !mcoefs(:) = coefs412_all(i,j,:,1,2)
       tLER = vgref412_all(i,j)!if there is no BRDF value, use VIIRS min ref (Jul 2017 W.KIM)
    endif

    !print *,"relaz = ",relaz,NDVI
       !coefs(:) = ncoefs(:)
      !LER = coefs(1) + scatAngle*(coefs(2) + scatAngle*(coefs(3) + scatAngle*coefs(4)))
    !print *,"ncoefs = ",ncoefs, LER
       !coefs(:) = acoefs(:)
      !LER = coefs(1) + scatAngle*(coefs(2) + scatAngle*(coefs(3) + scatAngle*coefs(4)))
    !print *,"acoefs = ",acoefs, LER
       !coefs(:) = mcoefs(:)
      !LER = coefs(1) + scatAngle*(coefs(2) + scatAngle*(coefs(3) + scatAngle*coefs(4)))
    !print *,"mcoefs = ",mcoefs, LER
    !print *,"tLER = ",tLER

    if (maxval(ncoefs) > 0.0) then
       coefs(:) = ncoefs(:)
    !    elseif (maxval(acoefs) > 0.0) then
    !       coefs(:) = acoefs(:)
    !    elseif (maxval(mcoefs) > 0.0) then
    !       coefs(:) = mcoefs(:)
    else
       coefs(:) = -999.0
    endif

    if (maxval(coefs) > 0.0) then
       LER = coefs(1) + scatAngle*(coefs(2) + scatAngle*(coefs(3) + scatAngle*coefs(4)))
    endif

    if (LER < 0.0) then
       LER = tLER
       min_flag = 1
    endif

 end function get_LER412

 ! -- Retrieve LER470 value from surface reflectivity coefficient table for pixel at
 ! --  latidx, lonidx in table.
 real function get_LER470(latidx, lonidx, NDVI, scatAngle, relaz, min_flag) RESULT(LER)
    integer, intent(in) :: latidx, lonidx
    real, intent(in)    :: NDVI, scatAngle, relaz
    integer, intent(inout) :: min_flag

    integer :: i,j,nidx
    real    :: coefs(4), acoefs(4), ncoefs(4), mcoefs(4), tLER

    LER = -999.0
    if (dateline .eq. 0 .OR. lonidx .gt. LERstart(1)) then
       i = lonidx - LERstart(1)
    else
       i = lonidx + dateline
    end if
    j = latidx - LERstart(2)

    !print *,"in get_LER470"
    !print *,"latidx, lonidx = ",latidx, lonidx
    !print *,"i,j = ",i,j
    !print *,"LERstart = ",LERstart
    !print *,"LERedge = ",LERedge

    if (NDVI < NDVI1_CUTOFF) then
       nidx = 1
    elseif (NDVI < NDVI2_CUTOFF) then
       nidx = 2
    else
       nidx = 3
    endif

    acoefs(:) = -999.0
    mcoefs(:) = -999.0
    ncoefs(:) = -999.0
    if (relaz < 90.0) then
       !ncoefs(:) = coefs470_fwd_tp(i,j,:,nidx,1)
       !acoefs(:) = coefs470_fwd_tp(i,j,:,1,1)
       !coefs(:) = coefs470_fwd_tp(i,j,:,1,2)
       ncoefs(:) = coefs470_fwd(i,j,:,nidx)
       !acoefs(:) = coefs470_fwd(i,j,:,1,1)
       !mcoefs(:) = coefs470_fwd(i,j,:,1,2)
       tLER = vgref488_all(i,j)!if there is no BRDF value, use VIIRS min ref (Jul 2017 W.KIM)
    else
       !ncoefs(:) = coefs470_all_tp(i,j,:,nidx,1)
       !acoefs(:) = coefs470_all_tp(i,j,:,1,1)
       !mcoefs(:) = coefs470_all_tp(i,j,:,1,2)
       ncoefs(:) = coefs470_all(i,j,:,nidx)
       !acoefs(:) = coefs470_all(i,j,:,1,1)
       !mcoefs(:) = coefs470_all(i,j,:,1,2)
       tLER = vgref488_all(i,j)!if there is no BRDF value, use VIIRS min ref (Jul 2017 W.KIM)
    endif

    !print *,"relaz = ",relaz,NDVI
       !coefs(:) = ncoefs(:)
      !LER = coefs(1) + scatAngle*(coefs(2) + scatAngle*(coefs(3) + scatAngle*coefs(4)))
    !print *,"ncoefs = ",ncoefs, LER
       !coefs(:) = acoefs(:)
      !LER = coefs(1) + scatAngle*(coefs(2) + scatAngle*(coefs(3) + scatAngle*coefs(4)))
    !print *,"acoefs = ",acoefs, LER
       !coefs(:) = mcoefs(:)
      !LER = coefs(1) + scatAngle*(coefs(2) + scatAngle*(coefs(3) + scatAngle*coefs(4)))
    !print *,"mcoefs = ",mcoefs, LER
    !print *,"tLER = ",tLER
    
    if (maxval(ncoefs) > 0.0) then
       coefs(:) = ncoefs(:)
    elseif (maxval(acoefs) > 0.0) then
       coefs(:) = acoefs(:)
    elseif (maxval(mcoefs) > 0.0) then
       coefs(:) = mcoefs(:)
    else
       coefs(:) = -999.0
    endif

    if (maxval(coefs) > 0.0) then
       LER = coefs(1) + scatAngle*(coefs(2) + scatAngle*(coefs(3) + scatAngle*coefs(4)))
    endif

    if (LER < 0.0) then
       LER = tLER
       min_flag = 1
    endif

 end function get_LER470

 ! -- Retrieve LER650 value from surface reflectivity coefficient table for pixel at
 ! --  latidx, lonidx in table.
 real function get_LER650(latidx, lonidx, NDVI, scatAngle, relaz,min_flag) RESULT(LER)
    integer, intent(in) :: latidx, lonidx
    real, intent(in)    :: NDVI, scatAngle, relaz
    integer, intent(inout) :: min_flag

    integer :: i,j,nidx
    real    :: coefs(4), acoefs(4), ncoefs(4), mcoefs(4), tLER

    LER = -999.0
    if (dateline .eq. 0 .OR. lonidx .gt. LERstart(1)) then
       i = lonidx - LERstart(1)
    else
       i = lonidx + dateline
    end if
    j = latidx - LERstart(2)

    !    print *,"in get_LER650"
    !    print *,"latidx, lonidx = ",latidx, lonidx
    !    print *,"i,j = ",i,j
    !    print *,"LERstart = ",LERstart
    !    print *,"LERedge = ",LERedge
    !    print *,"dateline = ", dateline
    
    if (NDVI < NDVI1_CUTOFF) then
       nidx = 1
    elseif (NDVI < NDVI2_CUTOFF) then
       nidx = 2
    else
       nidx = 3
    endif

    if (relaz < 90.0) then
       !ncoefs(:) = coefs650_fwd_tp(i,j,:,nidx,1)
       !acoefs(:) = coefs650_fwd_tp(i,j,:,1,1)
       !mcoefs(:) = coefs650_fwd_tp(i,j,:,1,2)
       ncoefs(:) = coefs650_fwd(i,j,:,nidx)
       !acoefs(:) = coefs650_fwd(i,j,:,1)
       !mcoefs(:) = coefs650_fwd(i,j,:,1)
       tLER = vgref670_all(i,j)!if there is no BRDF value, use VIIRS min ref (Jul 2017 W.KIM)
    else
       !ncoefs(:) = coefs650_all_tp(i,j,:,nidx,1)
       !acoefs(:) = coefs650_all_tp(i,j,:,1,1)
       !mcoefs(:) = coefs650_all_tp(i,j,:,1,2)
       ncoefs(:) = coefs650_all(i,j,:,nidx)
        !acoefs(:) = coefs650_all(i,j,:,1,1)
        !mcoefs(:) = coefs650_all(i,j,:,1,2)
       tLER = vgref670_all(i,j)!if there is no BRDF value, use VIIRS min ref (Jul 2017 W.KIM)
    endif

    !print *,"relaz = ",relaz,NDVI
       !coefs(:) = ncoefs(:)
      !LER = coefs(1) + scatAngle*(coefs(2) + scatAngle*(coefs(3) + scatAngle*coefs(4)))
    !print *,"ncoefs = ",ncoefs, LER
       !coefs(:) = acoefs(:)
      !LER = coefs(1) + scatAngle*(coefs(2) + scatAngle*(coefs(3) + scatAngle*coefs(4)))
    !print *,"acoefs = ",acoefs, LER
       !coefs(:) = mcoefs(:)
      !LER = coefs(1) + scatAngle*(coefs(2) + scatAngle*(coefs(3) + scatAngle*coefs(4)))
    !print *,"mcoefs = ",mcoefs, LER
    !print *,"tLER = ",tLER

    if (maxval(ncoefs) > 0.0) then
       coefs(:) = ncoefs(:)
    !    elseif (maxval(acoefs) > 0.0) then
    !       coefs(:) = acoefs(:)
    !    elseif (maxval(mcoefs) > 0.0) then
    !       coefs(:) = mcoefs(:)
    else
       coefs(:) = -999.0
    endif

    if (maxval(coefs) > 0.0) then
       LER = coefs(1) + scatAngle*(coefs(2) + scatAngle*(coefs(3) + scatAngle*coefs(4)))
    endif

    if (LER < 0.0) then
       LER = tLER
       min_flag = 1
    endif

 end function get_LER650

 ! Read in LER tables
 integer function readLER2(start, edge, stride, sds_name, grp_id, outref) RESULT(status)

    !   include 'hdf.f90'
    !   include 'dffunc.f90'
    use netcdf

    implicit none
    
    integer, dimension(3), intent(in)   ::  start, edge, stride
    integer, intent(in)                 ::  grp_id
    real, intent(out)                   ::  outref(edge(1),edge(2))

    ! HDF vars
    character(len=255)    ::  sds_name
    character(len=255)    ::  dset_name
    character(len=255)    ::  attr_name
    character(len=255)    ::  group_name

    integer               ::  nc_id
    integer               ::  dim_id
    integer               ::  dset_id
    integer                           ::  sds_index, sds_id 
    integer, dimension(3)             ::  tmpedge, tmpstart
    real, dimension(:,:), allocatable ::  tmpout

    dset_name = sds_name
    status = nf90_inq_varid(grp_id, dset_name, dset_id)
    if (status /= NF90_NOERR) then
       print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
       return
    end if

    if (dateline .eq. 0) then
       status = nf90_get_var(grp_id, dset_id, outref, start=start, &
          stride=stride, count=edge)
       if (status /= NF90_NOERR) then
          print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
          return
       end if
    else
       ! The granule straddles the dateline, so we need to make an accommodation
       tmpedge(:) = edge(:)
       tmpedge(1) = dateline
       allocate(tmpout(tmpedge(1), tmpedge(2)), stat=status)
       if (status /= 0) then
          print *, "ERROR: Unable to allocate tmpedge: ", status
          return
       end if
       status = nf90_get_var(grp_id, dset_id, tmpout, start=start, &
          stride=stride, count=tmpedge)
       if (status /= NF90_NOERR) then
          print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
          return
       end if
       outref(1:dateline, :) = tmpout

       deallocate(tmpout, stat=status)
       if (status /= 0) then
          print *, "Failed to deallocate tmpout: ", status
          return
       end if

       tmpstart(:) = start(:)
       tmpstart(1) = 1
       tmpedge(1) = edge(1) - dateline
       allocate(tmpout(tmpedge(1), tmpedge(2)), stat=status)
       if (status /= 0) then
          print *, "ERROR: Unable to allocate tmpedge: ", status
          return
       end if
       status = nf90_get_var(grp_id, dset_id, tmpout, start=tmpstart, &
          stride=stride, count=tmpedge)
       if (status /= NF90_NOERR) then
          print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
          return
       end if
       outref(dateline+1:edge(1),:) = tmpout
      
       deallocate(tmpout, stat=status)
       if (status /= 0) then
          print *, "Failed to deallocate tmpout: ", status
          return
       end if
    end if
    
    return
 end function readLER2

 ! Read in LER tables
 integer function readLER5(start, edge, stride, sds_name, grp_id, outref) RESULT(status)

    !   include 'hdf.f90'
    !   include 'dffunc.f90'
    use netcdf

    implicit none
    
    integer, dimension(:), intent(in)       ::  start, edge, stride
    character(len=255), intent(in)          ::  sds_name
    integer, intent(in)                     ::  grp_id
    real, dimension(:,:,:,:), intent(inout) ::  outref
    
    ! HDF vars
    character(len=255)    ::  dset_name

    integer               ::  nc_id
    integer               ::  dim_id
    integer               ::  dset_id
    integer               ::  sds_index, sds_id 
    integer, dimension(5) ::  tmpedge, tmpstart
    real, dimension(:,:,:,:), allocatable ::  tmpout
    character(len=255)    ::  tmp_name
    integer               ::  rank, ntype, nattrs
    integer, dimension(5) ::  dims
    
    status = -1
    
    dset_name = sds_name
    status = nf90_inq_varid(grp_id, dset_name, dset_id)
    if (status /= NF90_NOERR) then
       print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
       return
    end if
    
    if (dateline .eq. 0) then

       status = nf90_get_var(grp_id, dset_id, outref, start=start, &
          stride=stride, count=edge)
       if (status /= NF90_NOERR) then
          print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
          return
       end if
    else
       ! The granule straddles the dateline, so we need to make an accommodation
       tmpedge(:) = edge(:)
       tmpedge(1) = dateline
       allocate(tmpout(tmpedge(1), tmpedge(2), tmpedge(3), tmpedge(4)), stat=status)
       if (status /= 0) then
          print *, "ERROR: Unable to allocate tmpedge: ", status
          return
       end if
       status = nf90_get_var(grp_id, dset_id, tmpout, start=start, &
          stride=stride, count=tmpedge)
       if (status /= NF90_NOERR) then
          print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
          return
       end if
       outref(1:dateline, :, :, :) = tmpout(:,:,:,:)

       deallocate(tmpout, stat=status)
       if (status /= 0) then
          print *, "Failed to deallocate tmpout: ", status
          return
       end if

       tmpstart(:) = start(:)
       tmpstart(1) = 1
       tmpedge(1) = edge(1) - dateline
       allocate(tmpout(tmpedge(1), tmpedge(2), tmpedge(3), tmpedge(4)), stat=status)
       if (status /= 0) then
          print *, "ERROR: Unable to allocate tmpedge: ", status
          return
       end if
       status = nf90_get_var(grp_id, dset_id, tmpout, start=tmpstart, &
          stride=stride, count=tmpedge)
       if (status /= NF90_NOERR) then
          print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
          return
       end if
       outref(dateline+1:edge(1),:,:,:) = tmpout(:,:,:,:)

       deallocate(tmpout, stat=status)
       if (status /= 0) then
          print *, "Failed to deallocate tmpout: ", status
          return
       end if
      
    end if
    status = 0
    return

 end function readLER5

 ! Read in 2.2 um surface database
 integer function readSWIR2(start, edge, stride, sds_name, grp_id, outref) RESULT(status)

    !   include 'hdf.f90'
    !   include 'dffunc.f90'
    use netcdf

    implicit none

    integer, dimension(3), intent(in)   ::  start, edge, stride
    integer, intent(in)                 ::  grp_id
    real, intent(out)                   ::  outref(edge(1),edge(2))

    ! HDF vars
    character(len=255)    ::  sds_name
    character(len=255)    ::  dset_name
    character(len=255)    ::  attr_name
    character(len=255)    ::  group_name

    integer               ::  nc_id
    integer               ::  dim_id
    integer               ::  dset_id
    integer                           ::  sds_index, sds_id
    integer, dimension(3)             ::  tmpedge, tmpstart
    real, dimension(:,:), allocatable ::  tmpout

    dset_name = sds_name
    status = nf90_inq_varid(grp_id, dset_name, dset_id)
    if (status /= NF90_NOERR) then
       print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
       return
    end if

    if (dateline6 .eq. 0) then
       status = nf90_get_var(grp_id, dset_id, outref, start=start, &
          stride=stride, count=edge)
       if (status /= NF90_NOERR) then
          print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
          return
       end if
    else
       ! The granule straddles the dateline6, so we need to make an accommodation
       tmpedge(:) = edge(:)
       tmpedge(1) = dateline6
       allocate(tmpout(tmpedge(1), tmpedge(2)), stat=status)
       if (status /= 0) then
          print *, "ERROR: Unable to allocate tmpedge: ", status
          return
       end if
       status = nf90_get_var(grp_id, dset_id, tmpout, start=start, &
          stride=stride, count=tmpedge)
       if (status /= NF90_NOERR) then
          print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
          return
       end if
       outref(1:dateline6, :) = tmpout

       deallocate(tmpout, stat=status)
       if (status /= 0) then
          print *, "Failed to deallocate tmpout: ", status
          return
       end if

       tmpstart(:) = start(:)
       tmpstart(1) = 1
       tmpedge(1) = edge(1) - dateline6
       allocate(tmpout(tmpedge(1), tmpedge(2)), stat=status)
       if (status /= 0) then
          print *, "ERROR: Unable to allocate tmpedge: ", status
          return
       end if
       status = nf90_get_var(grp_id, dset_id, tmpout, start=tmpstart, &
          stride=stride, count=tmpedge)
       if (status /= NF90_NOERR) then
          print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
          return
       end if
       outref(dateline6+1:edge(1),:) = tmpout

       deallocate(tmpout, stat=status)
       if (status /= 0) then
          print *, "Failed to deallocate tmpout: ", status
          return
       end if
    end if
    return

 end function readSWIR2

 ! Read in 2.2 um surface database
 integer function readSWIR3(start3, edges3, stride3, sds_name, grp_id, outref) RESULT(status)

    !   include 'hdf.f90'
    !   include 'dffunc.f90'
    use netcdf

    implicit none

    integer, dimension(4), intent(in)     ::  start3, edges3, stride3
    integer, intent(in)                   ::  grp_id
    real, dimension(:,:,:), intent(inout) ::  outref
    ! HDF vars
    character(len=255)    ::  sds_name
    character(len=255)    ::  dset_name
    character(len=255)    ::  attr_name
    character(len=255)    ::  group_name

    integer               ::  nc_id
    integer               ::  dim_id
    integer               ::  dset_id
    integer               ::  sds_index, sds_id
    integer, dimension(4) ::  tmpedge, tmpstart
    real, dimension(:,:,:), allocatable ::  tmpout
    character(len=255)    ::  tmp_name
    integer               ::  rank, ntype, nattrs
    integer, dimension(4) ::  dims

    status = -1

    dset_name = sds_name
    status = nf90_inq_varid(grp_id, dset_name, dset_id)
    if (status /= NF90_NOERR) then
       print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
       return
    end if

    if (dateline6 .eq. 0) then
       status = nf90_get_var(grp_id, dset_id, outref, start=start3, &
          stride=stride3, count=edges3)
       if (status /= NF90_NOERR) then
          print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
          return
       end if
    else
       ! The granule straddles the dateline6, so we need to make an accommodation
       tmpedge(:) = edges3(:)
       tmpedge(1) = dateline6
       allocate(tmpout(tmpedge(1), tmpedge(2), tmpedge(3)), stat=status)
       if (status /= 0) then
          print *, "ERROR: Unable to allocate tmpedge: ", status
          return
       end if
       status = nf90_get_var(grp_id, dset_id, tmpout, start=start3, &
          stride=stride3, count=tmpedge)
       if (status /= NF90_NOERR) then
          print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
          return
       end if
       outref(1:dateline6, :, :) = tmpout

       deallocate(tmpout, stat=status)
       if (status /= 0) then
          print *, "Failed to deallocate tmpout: ", status
          return
       end if

       tmpstart(:) = start3(:)
       tmpstart(1) = 1
       tmpedge(1) = edges3(1) - dateline6
       allocate(tmpout(tmpedge(1), tmpedge(2), tmpedge(3)), stat=status)
       if (status /= 0) then
          print *, "ERROR: Unable to allocate tmpedge: ", status
          return
       end if
       status = nf90_get_var(grp_id, dset_id, tmpout, start=tmpstart, &
          stride=stride3, count=tmpedge)
       if (status /= NF90_NOERR) then
          print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
          return
       end if
       outref(dateline6+1:edges3(1),:,:) = tmpout

       deallocate(tmpout, stat=status)
       if (status /= 0) then
          print *, "Failed to deallocate tmpout: ", status
          return
       end if
    end if
    return

    dset_name = sds_name
    status = nf90_inq_varid(grp_id, dset_name, dset_id)
    if (status /= NF90_NOERR) then
       print *, "ERROR: Failed to get ID of dataset "//trim(dset_name)//": ", status
       return
    end if

    if (dateline6 .eq. 0) then
       status = nf90_get_var(grp_id, dset_id, tmpout, start=start3, &
          stride=stride3, count=edges3)
       if (status /= NF90_NOERR) then
          print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
          return
       end if
    else
       ! The granule straddles the dateline6, so we need to make an accommodation
       tmpedge(:) = edges3(:)
       tmpedge(1) = dateline6
       allocate(tmpout(tmpedge(1), tmpedge(2), tmpedge(3)), stat=status)
       if (status /= 0) then
          print *, "ERROR: Unable to allocate tmpedge: ", status
          return
       end if
       status = nf90_get_var(grp_id, dset_id, tmpout, start=start3, &
          stride=stride3, count=tmpedge)
       if (status /= NF90_NOERR) then
          print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
          return
       end if
       outref(1:dateline6, :, :) = tmpout

       deallocate(tmpout, stat=status)
       if (status /= 0) then
          print *, "Failed to deallocate tmpout: ", status
          return
       end if

       tmpstart(:) = start3(:)
       tmpstart(1) = 0
       tmpedge(1) = edges3(1) - dateline6
       allocate(tmpout(tmpedge(1), tmpedge(2), tmpedge(3)), stat=status)
       if (status /= 0) then
          print *, "ERROR: Unable to allocate tmpedge: ", status
          return
       end if
       status = nf90_get_var(grp_id, dset_id, tmpout, start=tmpstart, &
          stride=stride3, count=tmpedge)
       if (status /= NF90_NOERR) then
          print *, "ERROR: Failed to read dataset "//trim(dset_name)//": ", status
          return
       end if
       outref(dateline6+1:edges3(1),:,:) = tmpout

       deallocate(tmpout, stat=status)
       if (status /= 0) then
          print *, "Failed to deallocate tmpout: ", status
          return
       end if
    end if
    return
 end function readSWIR3
  
 integer function latlon_to_index_LER(lat, lon, ilat, ilon) result(status)
    implicit none
    
    real, intent(in)        ::  lat
    real, intent(in)        ::  lon
    integer, intent(inout)  ::  ilat
    integer, intent(inout)  ::  ilon
    
    status = 0
    if (lat > 90.0 .OR. lat < -90.0) then
       print *, "ERROR: Invalid latitude specified: ", lat
       status = -1
       return
    end if
    if (lon > 180.0 .OR. lon < -180.0) then
       print *, "ERROR: Invalid longitude specified: ", lon
       status = -1
       return
    end if
    
    ilat = (lat + 90.0) * 10.0 + 1
    if (ilat > 1800)  ilat = 1800
    if (ilat < 1)     ilat = 1

    ilon = (lon + 180.0) * 10.0 + 1
    if (ilon > 3600)  ilon = 3600
    if (ilon < 1)     ilon = 1
    
    return
  
 end function latlon_to_index_LER
  
 real function get_viirs_LER412(latidx, lonidx) result(ref)
    implicit none
    
    integer, intent(in)   ::  latidx
    integer, intent(out)  ::  lonidx

    integer               ::  i, j
    
    if (dateline .eq. 0 .OR. lonidx .gt. LERstart(1)) then
       i = lonidx - LERstart(1)
    else
       i = lonidx + dateline
    end if
    j = latidx - LERstart(2)
    
    ref = vgref412_all(i,j)
    return
    
 end function get_viirs_LER412
  
 real function get_viirs_LER488(latidx, lonidx) result(ref)
    implicit none
    
    integer, intent(in)   ::  latidx
    integer, intent(out)  ::  lonidx

    integer               ::  i, j
    
    if (dateline .eq. 0 .OR. lonidx .gt. LERstart(1)) then
       i = lonidx - LERstart(1)
    else
       i = lonidx + dateline
    end if
    j = latidx - LERstart(2)
    
    ref = vgref488_all(i,j)
    return
    
 end function get_viirs_LER488
  
 real function get_viirs_LER670(latidx, lonidx) result(ref)
    implicit none
    
    integer, intent(in)   ::  latidx
    integer, intent(out)  ::  lonidx

    integer               ::  i, j
    
    if (dateline .eq. 0 .OR. lonidx .gt. LERstart(1)) then
       i = lonidx - LERstart(1)
    else
       i = lonidx + dateline
    end if
    j = latidx - LERstart(2)
    
    ref = vgref670_all(i,j)
    return
    
 end function get_viirs_LER670
  
 real function get_modis_LER412(latidx, lonidx) result(ref)
    implicit none
    
    integer, intent(in)   ::  latidx
    integer, intent(out)  ::  lonidx
    
    integer               ::  i, j

    if (dateline .eq. 0 .OR. lonidx .gt. LERstart(1)) then
       i = lonidx - LERstart(1)
    else
       i = lonidx + dateline
    end if
    j = latidx - LERstart(2)
    
    ref = gref412_all(i,j)
    return
    
 end function get_modis_LER412
  
 real function get_modis_LER470(latidx, lonidx) result(ref)
    implicit none
    
    integer, intent(in)   ::  latidx
    integer, intent(out)  ::  lonidx
    
    integer               ::  i, j
    
    if (dateline .eq. 0 .OR. lonidx .gt. LERstart(1)) then
       i = lonidx - LERstart(1)
    else
       i = lonidx + dateline
    end if
    j = latidx - LERstart(2)
    
    ref = gref470_all(i,j)
    return
    
 end function get_modis_LER470
  
 real function get_modis_LER650(latidx, lonidx) result(ref)
    implicit none
    
    integer, intent(in)   ::  latidx
    integer, intent(out)  ::  lonidx

    integer               ::  i, j

    if (dateline .eq. 0 .OR. lonidx .gt. LERstart(1)) then
       i = lonidx - LERstart(1)
    else
       i = lonidx + dateline
    end if
    j = latidx - LERstart(2)
    
    ref = gref650_all(i,j)
    return
    
 end function get_modis_LER650
  
 real function get_modis_LER865(latidx, lonidx) result(ref)
    implicit none
    
    integer, intent(in)   ::  latidx
    integer, intent(out)  ::  lonidx
    
    integer               ::  i, j

    if (dateline .eq. 0 .OR. lonidx .gt. LERstart(1)) then
       i = lonidx - LERstart(1)
    else
       i = lonidx + dateline
    end if
    j = latidx - LERstart(2)
    
    ref = gref865_all(i,j)
    return
    
 end function get_modis_LER865
  
 real function get_viirs_modisbrdf_LER412(ilat, ilon, ndvi, sa, ra) result(ref)
    implicit none
    
    integer, intent(in)   ::  ilat
    integer, intent(out)  ::  ilon
    real, intent(in)      ::  ndvi
    real, intent(in)      ::  sa
    real, intent(in)      ::  ra
      
    real                  ::  mb_sr412
    real                  ::  m_sr412
    real                  ::  v_sr412
    integer               ::  min_flag
      
    ref = -999.0
      
    mb_sr412 =  get_LER412(ilat, ilon, ndvi, sa, ra,min_flag)
    if (mb_sr412 < -900.0) then
       !        print *, "ERROR: Undefined MODIS BRDF-corrected surface reflectance found: ", mb_sr412
       return
    end if
            
    m_sr412 = get_modis_LER412(ilat,ilon)
    if (m_sr412 < -900.0) then
       !        print *, "ERROR: Undefined MODIS all-angle surface reflectance found: ", m_sr412
       return
    end if
      
    v_sr412 = get_viirs_LER412(ilat,ilon)
    if (v_sr412 < -900.0) then
       !        print *, "ERROR: Undefined VIIRS all-angle surface reflectance found: ", v_sr412
       return
    end if
      
    ref = v_sr412 * (mb_sr412 / m_sr412)
    return
            
 end function get_viirs_modisbrdf_LER412
  
 real function get_viirs_modisbrdf_LER488(ilat, ilon, ndvi, sa, ra) result(ref)
    implicit none
    
    integer, intent(in)   ::  ilat
    integer, intent(out)  ::  ilon
    real, intent(in)      ::  ndvi
    real, intent(in)      ::  sa
    real, intent(in)      ::  ra
      
    real                  ::  mb_sr470
    real                  ::  m_sr470
    real                  ::  v_sr488
    integer               ::  min_flag
      
    ref = -999.0
      
    mb_sr470 =  get_LER470(ilat, ilon, ndvi, sa, ra, min_flag)
    if (mb_sr470 < -900.0) then
       !        print *, "ERROR: Undefined MODIS BRDF-corrected surface reflectance found: ", mb_sr470
       return
    end if
            
    m_sr470 = get_modis_LER470(ilat,ilon)
    if (m_sr470 < -900.0) then
       !        print *, "ERROR: Undefined MODIS all-angle surface reflectance found: ", m_sr470
       return
    end if
      
    v_sr488 = get_viirs_LER488(ilat,ilon)
    if (v_sr488 < -900.0) then
       !        print *, "ERROR: Undefined VIIRS all-angle surface reflectance found: ", v_sr488
       return
    end if
      
    ref = v_sr488 * (mb_sr470 / m_sr470)
      
    return
            
 end function get_viirs_modisbrdf_LER488
  
 real function get_viirs_modisbrdf_LER670(ilat, ilon, ndvi, sa, ra) result(ref)
    implicit none
    
    integer, intent(in)   ::  ilat
    integer, intent(out)  ::  ilon
    real, intent(in)      ::  ndvi
    real, intent(in)      ::  sa
    real, intent(in)      ::  ra
      
    real                  ::  mb_sr650
    real                  ::  m_sr650
    real                  ::  v_sr670
    integer               ::  min_flag
      
    ref = -999.0
      
    mb_sr650 =  get_LER650(ilat, ilon, ndvi, sa, ra, min_flag)
    if (mb_sr650 < -900.0) then
       !        print *, "ERROR: Undefined MODIS BRDF-corrected surface reflectance found: ", mb_sr650
       return
    end if
            
    m_sr650 = get_modis_LER650(ilat,ilon)
    if (m_sr650 < -900.0) then
       !        print *, "ERROR: Undefined MODIS all-angle surface reflectance found: ", m_sr650
       return
    end if
      
    v_sr670 = get_viirs_LER670(ilat,ilon)
    if (v_sr670 < -900.0) then
       !        print *, "ERROR: Undefined VIIRS all-angle surface reflectance found: ", v_sr670
       return
    end if
      
    ref = v_sr670 * (mb_sr650 / m_sr650)
      
    return
            
 end function get_viirs_modisbrdf_LER670
  
  
 integer function get_geographic_zone(lat, lon, status) result(zone)
    implicit none
    
    real, intent(in)    ::  lat
    real, intent(in)    ::  lon
      
    integer   ::  ilat, ilon
    integer   ::  status
       
    status = -1
    
    if (lat > 90.0 .OR. lat < -90.0) then
       print *, "ERROR: Invalid latitude specified: ", lat
       status = -1
       return
    end if
    if (lon > 180.0 .OR. lon < -180.0) then
       print *, "ERROR: Invalid longitude specified: ", lon
       status = -1
       return
    end if
      
    ilat = (lat + 90.0) * 10 + 1
    if (ilat > 1800)  ilat = 1800
    if (ilat < 1)     ilat = 1

    ilon = (lon + 180.0) * 10 + 1
    if (ilon > 3600)  ilon = 3600
    if (ilon < 1)     ilon = 1
    
    zone = terrain_flag_new(ilon,ilat)
    status = 0
    
    return
    
 end function get_geographic_zone


 integer function get_sfc_elev_std(lat, lon, status) result(elev_std)
    implicit none

    real, intent(in)    ::  lat
    real, intent(in)    ::  lon

    integer   ::  ilat, ilon
    integer   ::  status

    status = -1

    if (lat > 90.0 .OR. lat < -90.0) then
       print *, "ERROR: Invalid latitude specified: ", lat
       status = -1
       return
    end if
    if (lon > 180.0 .OR. lon < -180.0) then
       print *, "ERROR: Invalid longitude specified: ", lon
       status = -1
       return
    end if

    ilat = (lat + 90.0) * 10 + 1
    if (ilat > 1800)  ilat = 1800
    if (ilat < 1)     ilat = 1

    ilon = (lon + 180.0) * 10 + 1
    if (ilon > 3600)  ilon = 3600
    if (ilon < 1)     ilon = 1

    elev_std = sfc_elev_std(ilon,ilat)
    status = 0

    return

 end function get_sfc_elev_std

 real function get_background_aod(lat, lon, season, status) result(aod)
    implicit none

    real, intent(in)    ::  lat
    real, intent(in)    ::  lon
    integer, intent(in) ::  season

    integer   ::  ilat, ilon
    integer   ::  status

    status = -1

    if (lat > 90.0 .OR. lat < -90.0) then
       print *, "ERROR: Invalid latitude specified: ", lat
       status = -1
       return
    end if
    if (lon > 180.0 .OR. lon < -180.0) then
       print *, "ERROR: Invalid longitude specified: ", lon
       status = -1
       return
    end if

    ilat = floor(lat + 90.0) + 1
    if (ilat > 180)  ilat = 180
    if (ilat < 1)    ilat = 1

    ilon = floor(lon + 180.0) + 1
    if (ilon > 360)  ilon = 360
    if (ilon < 1)    ilon = 1

    aod = bg_aod(ilon,ilat)
    status = 0

    return

 end function get_background_aod

 subroutine check_status(status, str)
    integer, intent (in) :: status
    character(len=*), intent (in)    :: str

    if(status /= 0) then
       print *, str, status
    end if
 end subroutine check_status

 end module modis_surface
