      subroutine modis(cbdoy, xscan, scan, cxscan, cscan, time, lat, lon, sza, vza, raa,                    &
      &   sca, m01_newgc, m03_newgc, m04_newgc, m05_newgc, m07_newgc, m08_newgc, m10_newgc, m11_newgc,      &
      &   m01_nogc, m03_nogc, m05_nogc, m11_oldgc, ls_mask, ndvi, elev, ps,                          &
      &   lbp_mask, obp_mask, smoke_mask, smoke_ae_mask,pyrocb_mask, high_alt_smoke_mask, n_total_pixels,   &
      &   output_file, windsp, winddir, oz, wv, ler412, ler488, ler670, qdf412, qdf488, qdf670, bathy, chl, &
      &   vcc, bt11, btd11, hires, cell_resolution,nrt,platform,dtlat, dtlon,DBDTaot,dtaod, dt_cldmsk,dt_qa,&
      &   uvaod,dtspec,dtfmf,dbdtfmf,lcld_mask, uv1, uv2, dbdt_refl)
!-----------------------------------------------------------------------
! !F90
!
! !DESCRIPTION:
!
!    The modis subroutine drives the processing of MODIS data through
!    the Deep Blue algorithm.  modis sequentially steps through the
!    intermediate data file, passing each datapoint first through
!    a reflectivity preprocessor, and then through the Deep Blue 
!    algorithm itself (find_v).  As the processed data is passed out 
!    from find_v, it is sorted and binned.  Once all the data has been
!    processed, it is averaged and then written out to the MOD04 data 
!    product.
!
! !INPUT PARAMETERS:
!
!    Type       Name             Description
!    ====       ====             ===========
!    INTEGER*4  fh               file handle for the intermediate 
!                                 data file
!    BYTE       QAFlags          3-D array of bit flags from the 
!                                 MOD04 product
!    INTEGER*4  dims             array of SDS dimensions from the 
!                                 MOD04 product
!
! !OUTPUT PARAMETERS:  none
!
! !REVISION HISTORY:
!
!    Initial Version by Jeremy Warner   12/01/2006
!    Updated by Clare Salustro          05/02/2008
!
! !TEAM-UNIQUE HEADER:
!
!    This software is developed by the Deep Blue Science Team
!    for the National Aeronautics and Space Administration,
!    Goddard Space Flight Center, under contract NAS5-02041.
!
! !REFERENCES AND CREDITS
!
! !DESIGN NOTES:
!
!
!   Externals:
!
!     MODIS_W_GENERIC            (MODIS_39500.f)
!
!   Functions:
!
!     extract_data
!     read_data
!     total
!     find_v
!     write_db
!
! !END
!-----------------------------------------------------------------------
      use modis_surface, only:    get_geographic_zone
                                                             
      use viirs_db_utils, only:   create_viirs_l2,      &
                                  viirs_calib_correction, &
                                  create_dbdt_aot550,&
                                  create_dbdt_pace,&
                                  extrapolate_uv
      
      use viirs_ler_luts, only:   realbuf
                         
      use viirs_ocean_aero, only: run_ocean_retrieval,  &
                                  viirs_ocean_cell_output, &
                                  check_status
     
      use seawifs_surface_pressure, only: get_elevation ! jlee added 05/16/2017
      use read_l1b, only:         latitude, longitude
      implicit none
      
      integer, parameter  ::  i64 = selected_int_kind(18)
      
      character(len=*)    ::  output_file, platform
      integer             ::  scan, xscan
      integer             ::  cxscan, cscan
      
      integer             ::  cbdoy
      real(kind=8), dimension(scan)   ::  time          ! scan times
      real, dimension(xscan,scan)  :: lat, lon
      real, dimension(xscan,scan)  ::  sza, vza, raa, sca
      real, dimension(xscan,scan)  ::  m01_newgc, m03_newgc, m04_newgc, m05_newgc, m07_newgc, m08_newgc
      real, dimension(xscan,scan)  ::  m10_newgc, m11_newgc, uv1, uv2
      real, dimension(xscan,scan)  ::  windsp, winddir, oz, wv,dum
      real, dimension(8, xscan,scan)  ::  oreflc
      real, dimension(10, xscan,scan) ::  lreflc

      real, dimension(xscan,scan)  ::  m11_oldgc
      real, dimension(xscan,scan)  ::  m01_nogc, m03_nogc, m05_nogc

      integer, dimension(xscan,scan) ::  ls_mask, lbp_mask, obp_mask
      integer, dimension(xscan,scan) ::  smoke_mask, lcld_mask
      integer, dimension(xscan,scan) ::  smoke_ae_mask, pyrocb_mask, high_alt_smoke_mask      
            
      real, dimension(xscan,scan)  ::  dstar, ndvi, elev, ps
      integer, dimension(cxscan,cscan)      :: n_total_pixels
      
      real, dimension(xscan,scan)  ::  turbid_res      ! temp debug
      real, dimension(xscan,scan)  ::  ler412, ler488, ler670
      real, dimension(xscan,scan)  ::  qdf412, qdf488, qdf670
      integer, dimension(xscan,scan)  ::  bathy ! bathymetry
      real, dimension(xscan,scan)  ::  chl ! Climatological Chl
      real, dimension(xscan,scan)  ::  bt11, btd11


      integer(kind=4) :: i, j, k, m
      integer         :: i1, i2, j1, j2, ii, jj
      integer         :: ls_flag,cell_resolution,nrt
      integer(kind=4) :: idim, jdim
      integer(kind=4) :: iidx, jidx,iidx2, jidx2
      integer(kind=4) :: r412_count(cxscan,cscan)
      integer(kind=4) :: t650_count(cxscan,cscan)
      integer         :: naot550_avg(cxscan,cscan), noreflc(cxscan,cscan), nlreflc(cxscan,cscan)
      integer(kind=2) :: nae_avg(cxscan,cscan)
      integer(kind=2) :: naot412_avg(cxscan,cscan)
      integer(kind=2) :: naot488_avg(cxscan,cscan)
      integer(kind=2) :: naot670_avg(cxscan,cscan)
      integer(kind=2) :: nssa_avg(cxscan,cscan,3)
      integer(kind=2) :: nsr_avg(cxscan,cscan,3), nref_avg(cxscan,cscan,3)
      real(kind=4)    :: aot550_best(cxscan,cscan), oaot550_best(cxscan,cscan),DBDTaot(cxscan,cscan,3)
      real(kind=4)    :: oae_best(cxscan,cscan),ofmf_best(cxscan,cscan)
      real(kind=4)    :: oreflc_mean(8,cxscan,cscan), lreflc_mean(10,cxscan,cscan),nref_sum(10,cxscan,cscan)
      real(kind=4)    :: noref_sum(8,cxscan,cscan)      
      real            :: uvaod(cxscan,cscan,2),dtspec(cxscan,cscan,9)
      real            :: dtlat(cxscan,cscan), dtlon(cxscan,cscan),dtfmf(cxscan,cscan),dbdtfmf(cxscan,cscan)
      real            :: dtaod(cxscan,cscan,4),dbdt_refl(cxscan,cscan,9)
   
      
      integer(kind=4) :: flags(4)
      integer(kind=4) :: aerosol_type(cxscan,cscan), combined_type(cxscan,cscan)

      real            ::  ws_avg(cxscan,cscan), wd_avg(cxscan,cscan)
      real            ::  wv_avg(cxscan,cscan), oz_avg(cxscan,cscan), elev_avg(cxscan,cscan), chl_avg(cxscan,cscan)
      real            ::  elev_avg_land(cxscan,cscan), elev_avg_ocean(cxscan,cscan)
      
      integer         ::  smoke_count(cxscan,cscan)
      integer         ::  smoke_ae_count(cxscan,cscan), pyrocb_count(cxscan,cscan), high_alt_smoke_count(cxscan,cscan) 
      
      integer(kind=2) :: alg_cnt(cxscan,cscan,3) 
      integer         ::  alg_flag(cxscan,cscan),alg_flag_old(cxscan,cscan),alg_flag_old2(cxscan,cscan)
      integer         ::  db_cnt(cxscan,cscan),db_cnt2(cxscan,cscan),db_cnt3(cxscan,cscan),db_cnt_old(cxscan,cscan)      
      integer         ::  veg_cnt(cxscan,cscan),veg_cnt2(cxscan,cscan),veg_cnt_old(cxscan,cscan)    
      integer(kind=4) :: conf_flag(cxscan,cscan), usefulness(cxscan,cscan), aero_type(cxscan,cscan), rcond(cxscan,cscan)
      integer(kind=4) :: oconf_flag(cxscan,cscan),DBDTqa(cxscan,cscan),DBDTflag(cxscan,cscan)
      integer(kind=2) :: edgecnt(cxscan,cscan)
      integer(kind=2) :: veg_bflag_cnt(cxscan,cscan)
      integer         :: dt_cldmsk(xscan,scan)
      real            :: dt_qa(cxscan,cscan)
      
      real(kind=4) :: outbuf(22)
      real(kind=4) :: tmpvg(7)
      real(kind=4) :: var
      real(kind=4) :: aot550_avg(cxscan,cscan), ae_avg(cxscan,cscan),  ae_best(cxscan,cscan),oaot550_avg(cxscan,cscan)
      real(kind=4) :: caot550_avg(cxscan,cscan), caot550_best(cxscan,cscan)
      real(kind=4) :: cae_avg(cxscan,cscan), cae_best(cxscan,cscan)
      real(kind=4) :: aot550_list(cxscan,cscan,100)
      real(kind=4) :: sd550(cxscan,cscan)
      real(kind=4) :: aot412_avg(cxscan,cscan)
      real(kind=4) :: aot488_avg(cxscan,cscan)
      real(kind=4) :: aot670_avg(cxscan,cscan)
      real(kind=4) :: ssa_avg(cxscan,cscan,3)
      real(kind=4) :: sr_avg(cxscan,cscan,3), ref_avg(cxscan,cscan,3)
      real(kind=4) :: sfc_typ(cxscan,cscan)
      real(kind=4) :: alpha,  alpha_550, xxrat,  dd
      real(kind=4) :: alpha2, xxrat2, alpha2_550
      real(kind=4) :: refbuf(3)
      real(kind=4) :: lat_sav(cxscan,cscan), lon_sav(cxscan,cscan)
      real(kind=4) :: sza_sav(cxscan,cscan), vza_sav(cxscan,cscan)
      real(kind=4) :: raa_sav(cxscan,cscan), sca_sav(cxscan,cscan)
!       real(kind=4) :: dstar_avg(cxscan,cscan)
      real(kind=4) :: btd11_avg(cxscan,cscan)
      real(kind=4) :: ratio_avg(cxscan,cscan)
      real(kind=4) :: ndvi_avg(cxscan,cscan) 
      real(kind=4) :: rcndvi_avg(cxscan,cscan) 
      integer      :: nrcndvi_avg(cxscan,cscan)
      integer      :: dateline(cxscan,cscan)
      integer(kind=2) ::  nll_sav(cxscan,cscan),nll_sav_land(cxscan,cscan),nll_sav_ocean(cxscan,cscan)
      integer      ::  cls_mask(cxscan,cscan)
    
      real(kind=8)    ::  time_avg(cscan)
      integer         ::  ntime(cscan)
        
      real, dimension(:,:), allocatable     ::  ls_sav

      type(viirs_ocean_cell_output)           ::  o_out
      type(viirs_ocean_cell_output), dimension(:), allocatable  ::  model_results

      real(kind=4) :: cl_flag  ! MJ test20110727
      integer, dimension(cxscan,cscan)      :: gzone
      integer         ::  status, doy, cnt, cnt2, cnt3
            
      real            ::  px_elev ! jlee added 05/16/2017 
      type(viirs_calib_correction) :: vcc   
      logical         ::  hires, px_debug          
      integer         ::  nagg
      
!     ---common parameter      
      common  /xday/ doy

!      print *, 'in modis.f90...'

!     -- aggregation threshold
      nagg = 9
      if (hires) nagg = 0

!     -- initialize stuff
      
!      print *, 'start initialization...'
      naot412_avg(:,:) = 0
      naot488_avg(:,:) = 0
      naot670_avg(:,:) = 0
      aot412_avg(:,:) = 0.0
      aot488_avg(:,:) = 0.0
      aot670_avg(:,:) = 0.0
      ref_avg(:,:,:)  = 0.0
      nssa_avg(:,:,:) = 0
      ssa_avg(:,:,:)  = 0.0
      sr_avg(:,:,:)   = 0.0
      nsr_avg(:,:,:)  = 0
      nref_avg(:,:,:) = 0
      
      nae_avg(:,:) = 0
      ae_avg(:,:) = 0.0
      ae_best(:,:) = -999.0
      naot550_avg(:,:) = 0
      aot550_avg(:,:) = 0.0
      caot550_avg(:,:) = -999.0
      caot550_best(:,:) = -999.0
      aot550_best(:,:) = -999.0
      oaot550_best(:,:) = -999.0
      oae_best(:,:) = -999.0
      ofmf_best(:,:) = -999.0
      r412_count(:,:) = 0
      t650_count(:,:) = 0
      sd550(:,:) = -999.
      sfc_typ(:,:) = -999.
      aot550_list(:,:,:) = -999.0
      DBDTaot(:,:,:) = -999.
      DBDTqa(:,:)=0
      DBDTflag(:,:)=0
      
      sza_sav(:,:) = 0.0
      vza_sav(:,:) = 0.0
      raa_sav(:,:) = 0.0
      sca_sav(:,:) = 0.0
!       dstar_avg(:,:) = 0.0  
      btd11_avg(:,:) = 0.0
      ratio_avg(:,:) = 0.0
      ndvi_avg(:,:) = 0.0
      rcndvi_avg(:,:) = 0.0
      nrcndvi_avg(:,:) = 0
      lat_sav(:,:) = 0.0
      lon_sav(:,:) = 0.0
      nll_sav(:,:) = 0
      nll_sav_land(:,:) = 0
      nll_sav_ocean(:,:) = 0
      gzone(:,:) = -999.0
      time_avg(:) = 0
      ntime(:)    = 0
      dateline(:,:) = 0     
      noreflc(:,:) = 0
      nlreflc(:,:) = 0 
      oreflc_mean(:,:,:) = 0
      lreflc_mean(:,:,:) = 0
 
      ws_avg(:,:)   = 0.0
      wd_avg(:,:)   = 0.0
      wv_avg(:,:)   = 0.0
      oz_avg(:,:)   = 0.0
      elev_avg(:,:)   = 0.0
      elev_avg_land(:,:)  = 0.0
      elev_avg_ocean(:,:) = 0.0
      
      chl_avg(:,:)   = 0.0

!     -- initialize algorithm flag arrays
      alg_cnt(:,:,:)  = 0
      alg_flag(:,:)   = -999
      alg_flag_old(:,:)= -999
      alg_flag_old2(:,:)= -999
      db_cnt(:,:)     = 0
      db_cnt_old(:,:) = 0
      veg_cnt(:,:)    = 0
      veg_cnt_old(:,:)= 0
      db_cnt2(:,:)    = 0
      db_cnt3(:,:)    = 0
      veg_cnt2(:,:)   = 0      
      conf_flag(:,:)  = 0
      oconf_flag(:,:) = 0
      usefulness(:,:) = 0
      aero_type(:,:)  = 0
      rcond(:,:)      = 0
      
      veg_bflag_cnt(:,:) = 0
      smoke_count(:,:)   = 0
      
!      print *, 'end initialization.'

      allocate(ls_sav(xscan,scan), stat=status)
      if (status /= 0) then
        print *, "ERROR: Failed to allocate array for land sea flag: ", status
        stop
      end if
      ls_sav(:,:) = -999
      
!-------------------------------------------------------------------------------------------
! ocean
!-------------------------------------------------------------------------------------------
!     -- perform ocean retrievals if we have ocean pixels.
!      print *, 'skip ocean retrieval....'

!     -- perform calibration correction for ocean retrievals
!       if (platform .eq. 'VIIRS') then
!         where (m01_newgc > -900.0) m01_newgc = m01_newgc*vcc%m01
!         where (m08_newgc > -900.0) m08_newgc = m08_newgc*vcc%m08   
!         where (m03_newgc > -900.0) m03_newgc = m03_newgc*vcc%m03
!         where (m04_newgc > -900.0) m04_newgc = m04_newgc*vcc%m04
!         where (m05_newgc > -900.0) m05_newgc = m05_newgc*vcc%m05
!         where (m07_newgc > -900.0) m07_newgc = m07_newgc*vcc%m07
!         where (m10_newgc > -900.0) m10_newgc = m10_newgc*vcc%m10
!         where (m11_newgc > -900.0) m11_newgc = m11_newgc*vcc%m11
!       end if
      
      allocate(model_results(4), stat=status)
      if (status /= 0) then
        print *, 'ERROR: Failed to allocate array for all model results: ', status
        return
      end if
      turbid_res(:,:) = -999.0
		
      o_out = run_ocean_retrieval(lat, lon, sza, vza, raa, m03_newgc,     &
      &         m04_newgc, m05_newgc, m07_newgc, m08_newgc, m10_newgc, m11_newgc, windsp, ls_mask,     &
      &         obp_mask, status, turbid_res, bathy, chl, .false., model_results,elev,cell_resolution,&
      &         nrt,platform,cxscan,cscan)

      if (status /= 0) then
        print *, "ERROR: Ocean retrieval failed: ", status
        stop
      end if
!      print *, 'end ocean retrieval.'
 
    !create spectral radiance variable
!       if (platform .eq. 'VIIRS') then
        dum = m01_newgc
        where (obp_mask == 1 .or. ls_mask < 2) dum = -999.
        oreflc(1,:,:) = dum
        
        dum = m08_newgc
        where (obp_mask == 1 .or. ls_mask < 2) dum = -999.
        oreflc(6,:,:) = dum
!       end if

      dum = m03_newgc
      where (obp_mask == 1 .or. ls_mask < 2) dum = -999.
      oreflc(2,:,:) = dum
      
!       if (platform .ne. 'GOES') then
        dum = m04_newgc
        where (obp_mask == 1 .or. ls_mask < 2) dum = -999.
        oreflc(3,:,:) = dum
!       end if

      dum = m05_newgc
      where (obp_mask == 1 .or. ls_mask < 2) dum = -999.
      oreflc(4,:,:) = dum
 
      dum = m07_newgc
      where (obp_mask == 1 .or. ls_mask < 2) dum = -999.
      oreflc(5,:,:) = dum
      
      dum = m10_newgc
      where (obp_mask == 1 .or. ls_mask < 2) dum = -999.
      oreflc(7,:,:) = dum

      dum = m11_newgc
      where (obp_mask == 1 .or. ls_mask < 2) dum = -999.
      oreflc(8,:,:) = dum      
      
!-------------------------------------------------------------------------------------------
! land
!-------------------------------------------------------------------------------------------
!      print *, 'start DB land retrieval....'

!     -- perform land retrievals

!     -- revert calibration correction for land retrievals
!       if (platform .eq. 'VIIRS') then
!         where (m01_newgc > -900.0) m01_newgc = m01_newgc/vcc%m01
!         where (m08_newgc > -900.0) m08_newgc = m08_newgc/vcc%m08
!         where (m03_newgc > -900.0) m03_newgc = m03_newgc/vcc%m03
!         where (m04_newgc > -900.0) m04_newgc = m04_newgc/vcc%m04
!         where (m05_newgc > -900.0) m05_newgc = m05_newgc/vcc%m05
!         where (m07_newgc > -900.0) m07_newgc = m07_newgc/vcc%m07
!         where (m10_newgc > -900.0) m10_newgc = m10_newgc/vcc%m10
!         where (m11_newgc > -900.0) m11_newgc = m11_newgc/vcc%m11        
!       end if 

    !create spectral radiance variable
!       if (platform .eq. 'VIIRS') then
        dum = m01_newgc
        where (lbp_mask == 1 .or. ls_mask >= 2) dum = -999.
        lreflc(1,:,:) = dum
        dum = m08_newgc
        where (lbp_mask == 1 .or. ls_mask >= 2) dum = -999.
        lreflc(6,:,:) = dum
!       end if 
   
      dum = m03_newgc
      where (lbp_mask == 1 .or. ls_mask >= 2) dum = -999.
      lreflc(2,:,:) = dum

!       if (platform .ne. 'GOES') then      
        dum = m04_newgc
        where (lbp_mask == 1 .or. ls_mask >= 2) dum = -999.
        lreflc(3,:,:) = dum
!       end if

      dum = m05_newgc
      where (lbp_mask == 1 .or. ls_mask >= 2) dum = -999.
      lreflc(4,:,:) = dum
 
      dum = m07_newgc
      where (lbp_mask == 1 .or. ls_mask >= 2) dum = -999.
      lreflc(5,:,:) = dum
      
      dum = m10_newgc
      where (lbp_mask == 1 .or. ls_mask >= 2) dum = -999.
      lreflc(7,:,:) = dum

      dum = m11_newgc
      where (lbp_mask == 1 .or. ls_mask >= 2) dum = -999.
      lreflc(8,:,:) = dum      

      dum = uv1
      where (lbp_mask == 1 .or. ls_mask >= 2) dum = -999.
      lreflc(9,:,:) = dum    

      dum = uv2
      where (lbp_mask == 1 .or. ls_mask >= 2) dum = -999.
      lreflc(10,:,:) = dum   
                         
!       print *, 'start modis.f loop :', xscan, scan
      do jdim = 1, scan 
        do idim = 1, xscan   

        !reflectance aggregation
        iidx = (idim-1)/cell_resolution + 1
        jidx = (jdim-1)/cell_resolution + 1 
           
        if (iidx >cxscan .or. jidx>cscan) then 
          if (idim ==xscan  .and. jdim==scan) go to 120
          cycle
        end if
        !for aggregated pixel level debugging  
!         px_debug = .false.
!         if (iidx == 397 .and. jidx == 47) px_debug = .true.
!         if (px_debug) print*,oz(idim, jdim) 
          
!         if (oreflc(2,idim,jdim) .ge. -900) then
! !           if (platform .eq. 'VIIRS') then
!             oreflc_mean(1,iidx,jidx)  = oreflc_mean(1,iidx,jidx) + oreflc(1,idim,jdim)
!             oreflc_mean(6,iidx,jidx)  = oreflc_mean(6,iidx,jidx) + oreflc(6,idim,jdim)
! !           end if
! !           if (platform .ne. 'GOES') then
!             oreflc_mean(3,iidx,jidx)  = oreflc_mean(3,iidx,jidx) + oreflc(3,idim,jdim)  
! !           end if 
!           oreflc_mean(2,iidx,jidx)  = oreflc_mean(2,iidx,jidx) + oreflc(2,idim,jdim)        
!           oreflc_mean(4,iidx,jidx)  = oreflc_mean(4,iidx,jidx) + oreflc(4,idim,jdim)
!           oreflc_mean(5,iidx,jidx)  = oreflc_mean(5,iidx,jidx) + oreflc(5,idim,jdim)
!           oreflc_mean(7,iidx,jidx)  = oreflc_mean(7,iidx,jidx) + oreflc(7,idim,jdim)
!           oreflc_mean(8,iidx,jidx)  = oreflc_mean(8,iidx,jidx) + oreflc(8,idim,jdim)
!           noreflc(iidx,jidx)        = noreflc(iidx,jidx) + 1
!         end if
!         if (lreflc(2,idim,jdim) .ge. -900) then
! !           if (platform .eq. 'VIIRS') then        
!             lreflc_mean(1,iidx,jidx)  = lreflc_mean(1,iidx,jidx) + lreflc(1,idim,jdim)
!             lreflc_mean(6,iidx,jidx)  = lreflc_mean(6,iidx,jidx) + lreflc(6,idim,jdim)
! !           end if
! !           if (platform .ne. 'GOES') then
!             lreflc_mean(3,iidx,jidx)  = lreflc_mean(3,iidx,jidx) + lreflc(3,idim,jdim)         
! !           end if
!           lreflc_mean(2,iidx,jidx)  = lreflc_mean(2,iidx,jidx) + lreflc(2,idim,jdim)
!           lreflc_mean(4,iidx,jidx)  = lreflc_mean(4,iidx,jidx) + lreflc(4,idim,jdim)
!           lreflc_mean(5,iidx,jidx)  = lreflc_mean(5,iidx,jidx) + lreflc(5,idim,jdim)
!           lreflc_mean(7,iidx,jidx)  = lreflc_mean(7,iidx,jidx) + lreflc(7,idim,jdim)
!           lreflc_mean(8,iidx,jidx)  = lreflc_mean(8,iidx,jidx) + lreflc(8,idim,jdim)
!           nlreflc(iidx,jidx)        = nlreflc(iidx,jidx) + 1
!         end if
                
!          if (abs(31.5-lat(idim,jdim)) > 0.025 .OR. abs(71.75-lon(idim,jdim)) >  0.025) cycle 
        if (lbp_mask(idim,jdim) == 1) cycle    ! skip bad pixels, cloudy, etc
        if (ls_mask(idim,jdim) >= 2) cycle      ! skip ocean pixels, taken care of above.
        if (vza(idim,jdim) > 90. .or. sza(idim,jdim) > 90. .or. vza(idim,jdim) < 0. .or. sza(idim,jdim) < 0.) cycle
        doy = cbdoy
        ls_flag = ls_mask(idim,jdim)
        if (platform .eq. 'VIIRS') then
          tmpvg = (/m05_newgc(idim,jdim), m03_newgc(idim,jdim), m01_newgc(idim,jdim), m11_oldgc(idim,jdim), &
            &         m07_newgc(idim,jdim), m08_newgc(idim,jdim), m11_oldgc(idim,jdim)/)
        end if 
        if (platform .eq. 'AHI' .or. platform .eq. 'GOES' ) then
          tmpvg = (/m05_newgc(idim,jdim), m03_newgc(idim,jdim),-999., m11_oldgc(idim,jdim), &
            &         m07_newgc(idim,jdim), -999.,  m11_oldgc(idim,jdim)/)    
        end if       
        cl_flag = -999.0  ! keep this -999.0 otherwise crazy cloud checks will be turned on below...
!     -- move read_data operations here.
        if (platform .eq. 'VIIRS') then  
          if (m05_newgc(idim,jdim) < 0.0 .OR. m03_newgc(idim,jdim) < 0.0 .OR. m01_newgc(idim,jdim) < 0.0 &
          & .or. m11_oldgc(idim,jdim) <=0.0) cycle
        end if 
        if (platform .eq. 'AHI' .or. platform .eq. 'GOES' ) then 
          if (m05_newgc(idim,jdim) < 0.0 .OR. m03_newgc(idim,jdim) < 0.0 ) cycle
        end if
     
!     Parse the input data into a form deepblue is comfortable with

        realbuf(1)  = lat(idim, jdim)
        realbuf(2)  = lon(idim, jdim)
        realbuf(3)  = sza(idim,jdim)
        realbuf(4)  = vza(idim,jdim)
        realbuf(5)  = raa(idim,jdim)
        realbuf(6)  = ler670(idim,jdim)
        realbuf(7)  = m05_nogc(idim,jdim)
        realbuf(8)  = m03_newgc(idim,jdim)
        realbuf(9)  = m03_nogc(idim,jdim)

        realbuf(12) = m11_oldgc(idim,jdim)
        realbuf(13) = -999.0  
        realbuf(14) = ndvi(idim,jdim)
!         realbuf(15) = bt11(idim,jdim)
!         realbuf(16) = dstar(idim,jdim)
        realbuf(17) = cl_flag

        realbuf(19) = qdf488(idim,jdim)
        realbuf(20) = qdf670(idim,jdim)
        realbuf(21) = ps(idim,jdim)/1013.25 
        realbuf(22) = ler412(idim,jdim)
        realbuf(23) = ler488(idim,jdim)
        realbuf(24) = -999.0
        realbuf(25) = -999.0
        realbuf(26) = -999.0
      
!         if (platform .eq. 'VIIRS') then 
           realbuf(10) = m01_newgc(idim,jdim)
           realbuf(11) = m01_nogc(idim,jdim)
           realbuf(18) = qdf412(idim,jdim)  
!         end if    
      
        outbuf(:)       = -999.0
        
        if (realbuf(1) .lt. -900.0 .or. realbuf(2) .lt. -900.0) cycle
 
        ls_sav(idim,jdim) = ls_flag

        iidx = (idim-1)/cell_resolution + 1
        jidx = (jdim-1)/cell_resolution + 1
      
  !     Extract reflectance from realbuf (are ref*100)
  !     refbuf(1) = realbuf(22)  ! original; LER at 412nm (MJ)
  !     refbuf(2) = realbuf(23)  ! original; LER at 470nm (MJ)
  !     refbuf(3) = realbuf(6)   ! original
      
        refbuf(1) = 3.14159*(realbuf(11))*100.  ! MJ test TOA refl(20100929), B8
        refbuf(2) = 3.14159*(realbuf(9))*100.  ! MJ test TOA refl(20101105), B3
        refbuf(3) = 3.14159*(realbuf(7))*100.  ! MJ test TOA refl(20101105), B1
!         print *, lon(idim, jdim),lat(idim, jdim),idim, jdim, elev(idim,jdim)
        call find_v_viirs(realbuf, tmpvg, outbuf, flags, elev(idim,jdim), ls_flag, &
              & windsp(idim,jdim), wv(idim,jdim), platform)
        !if (iidx==85.and.jidx==394) print *, outbuf(1:3)     
 
  !     -- record algorithm path for this pixel
        if (outbuf(7) .gt. -900.0) then 
          select case (int(outbuf(18)))
  !           case (0)      ! Deep Blue
  !             db_cnt(iidx,jidx) = db_cnt(iidx,jidx) + 1
  !           case (10)      ! Vegetated
  !             veg_cnt(iidx,jidx) = veg_cnt(iidx,jidx) + 1
            case (1)      ! Deep Blue AERONET brdf
              db_cnt(iidx,jidx) = db_cnt(iidx,jidx) + 1
            case (2)      ! Deep Blue BRDF from table
              db_cnt2(iidx,jidx) = db_cnt2(iidx,jidx) + 1  
            case (3)      ! Deep Blue min SR from db
              db_cnt3(iidx,jidx) = db_cnt3(iidx,jidx) + 1    
            case (11)      ! Vegetated
              veg_cnt(iidx,jidx) = veg_cnt(iidx,jidx) + 1
            case (12)      ! 2.2 micron
              veg_cnt2(iidx,jidx) = veg_cnt2(iidx,jidx) + 1  
            case (-999)   ! no retrieval
            case default
!               print *, "ERROR: Invalid algorithm flag: ", outbuf(18)
              cycle
          end select
        endif
        !old_algorithm flag
        if (outbuf(7) .gt. -900.0) then
          select case (int(outbuf(22)))
            case (0)      ! Deep Blue
              db_cnt_old(iidx,jidx) = db_cnt_old(iidx,jidx) + 1
            case (1)      ! Vegetated
              veg_cnt_old(iidx,jidx) = veg_cnt_old(iidx,jidx) + 1
            case (-999)   ! no retrieval
            case default
!               print *, "ERROR: Invalid algorithm flag: ", outbuf(22)
              cycle
          end select
        endif      
      
  !     Extract and average output
        if (outbuf(1) .ge. 0.0 .and. outbuf(1) .le. 5.0) then
          aot412_avg(iidx,jidx)  = aot412_avg(iidx,jidx) + outbuf(1)
          naot412_avg(iidx,jidx) = naot412_avg(iidx,jidx) + 1
        end if
        if (outbuf(2) .ge. 0.0 .and. outbuf(2) .le. 5.0) then
          aot488_avg(iidx,jidx)  = aot488_avg(iidx,jidx) + outbuf(2)
          naot488_avg(iidx,jidx) = naot488_avg(iidx,jidx) + 1
        end if
        if (outbuf(3) .ge. 0.0 .and. outbuf(3) .le. 5.0) then
          aot670_avg(iidx,jidx)  = aot670_avg(iidx,jidx) + outbuf(3)
          naot670_avg(iidx,jidx) = naot670_avg(iidx,jidx) + 1
        end if
      
        do i=1,3
           if (refbuf(i) .gt. 0.0) then 
             ref_avg(iidx,jidx,i)  = ref_avg(iidx,jidx,i) + refbuf(i)
             nref_avg(iidx,jidx,i) = nref_avg(iidx,jidx,i) + 1
           endif
           if (outbuf(i+3) .ge. 0.0 .and. outbuf(i+3) .le. 5.0) then
              ssa_avg(iidx,jidx,i)  = ssa_avg(iidx,jidx,i) + outbuf(i+3)
              nssa_avg(iidx,jidx,i) = nssa_avg(iidx,jidx,i) + 1
           endif
       enddo
      
        if (outbuf(7) .gt. -900.0) then
           aot550_avg(iidx,jidx) = aot550_avg(iidx,jidx) + outbuf(7)
           naot550_avg(iidx,jidx) = naot550_avg(iidx,jidx) + 1
           aot550_list(iidx,jidx,naot550_avg(iidx,jidx)) = outbuf(7)         
!            dstar_avg(iidx,jidx) = dstar_avg(iidx,jidx) + dstar(idim,jdim)
!            btd11_avg(iidx,jidx) = btd11_avg(iidx,jidx) + btd11(idim,jdim)
           if (outbuf(21) .gt. -900.0)   &
           ratio_avg(iidx,jidx) = ratio_avg(iidx,jidx) + outbuf(21)   !  LER412/488nm ratio
           ndvi_avg(iidx,jidx) = ndvi_avg(iidx,jidx) + ndvi(idim,jdim)
           
           lreflc_mean(1,iidx,jidx)  = lreflc_mean(1,iidx,jidx) + lreflc(1,idim,jdim)                 
           lreflc_mean(2,iidx,jidx)  = lreflc_mean(2,iidx,jidx) + lreflc(2,idim,jdim)
           lreflc_mean(3,iidx,jidx)  = lreflc_mean(3,iidx,jidx) + lreflc(3,idim,jdim)  
           lreflc_mean(4,iidx,jidx)  = lreflc_mean(4,iidx,jidx) + lreflc(4,idim,jdim)
           lreflc_mean(5,iidx,jidx)  = lreflc_mean(5,iidx,jidx) + lreflc(5,idim,jdim)
           lreflc_mean(6,iidx,jidx)  = lreflc_mean(6,iidx,jidx) + lreflc(6,idim,jdim)          
           lreflc_mean(7,iidx,jidx)  = lreflc_mean(7,iidx,jidx) + lreflc(7,idim,jdim)
           lreflc_mean(8,iidx,jidx)  = lreflc_mean(8,iidx,jidx) + lreflc(8,idim,jdim)
           lreflc_mean(9,iidx,jidx)  = lreflc_mean(9,iidx,jidx) + lreflc(9,idim,jdim)
           lreflc_mean(10,iidx,jidx) = lreflc_mean(10,iidx,jidx) + lreflc(10,idim,jdim)
           nlreflc(iidx,jidx)        = nlreflc(iidx,jidx) + 1          
        endif

        if (outbuf(8) .gt. -900.0) then
           ae_avg(iidx,jidx) = ae_avg(iidx,jidx) + outbuf(8)
           nae_avg(iidx,jidx) = nae_avg(iidx,jidx) + 1
        endif

        if (outbuf(9) .gt. 0.12) r412_count(iidx,jidx) = r412_count(iidx,jidx)+1
        if (outbuf(10) .gt. 0.) t650_count(iidx,jidx) = t650_count(iidx,jidx)+1  !  out-of-bound flag

        if (outbuf(9) .gt. 0.0) then
  !         if (outbuf(18) == 0) then    ! SR412 is only defined for DB retrievals.
          if (outbuf(22) == 0) then    ! SR412 is only defined for DB retrievals.
           sr_avg(iidx,jidx,1) = sr_avg(iidx,jidx,1) + outbuf(9)/100
           nsr_avg(iidx,jidx,1) = nsr_avg(iidx,jidx,1) + 1
          end if
        endif
        if (outbuf(11) .gt. 0.0) then
           sr_avg(iidx,jidx,2) = sr_avg(iidx,jidx,2) + outbuf(11)/100
           nsr_avg(iidx,jidx,2) = nsr_avg(iidx,jidx,2) + 1
        endif
        if (outbuf(12) .gt. 0.0) then
           sr_avg(iidx,jidx,3) = sr_avg(iidx,jidx,3) + outbuf(12)/100
           nsr_avg(iidx,jidx,3) = nsr_avg(iidx,jidx,3) + 1
        endif

  !     -- count number of smoke_mask pixels in this cell and save.
        if (outbuf(14) > 0.0 .or. smoke_mask(idim,jdim) == 1) then
!         if (outbuf(14) > 0.0) then
          smoke_count(iidx,jidx) = smoke_count(iidx,jidx) + 1
        end if

  !     -- count number of pixels triggering the vegetated bright surface flag. 
        if (outbuf(19) > 0.5) then 
          veg_bflag_cnt(iidx,jidx) = veg_bflag_cnt(iidx,jidx) + 1
        end if

        if (outbuf(20) > -900.0) then
          nrcndvi_avg(iidx,jidx) = nrcndvi_avg(iidx,jidx) + 1
          rcndvi_avg(iidx,jidx) = rcndvi_avg(iidx,jidx) + outbuf(20)
        end if
      
  !     Want values for central pixel (or nearby) for three variables
        if (hires .NEQV. .true.) then
         if (mod(idim,8) .lt. 1. .or. mod(jdim,8) .lt. 1. .or. &
                 mod(idim,8) .gt. 6.5 .or. mod(jdim,8) .gt. 6.5) cycle

         if (mod(idim,8) .le. 4 .and. mod(jdim,8) .le. 4) then
         
            if (outbuf(15) .gt. -900.) then
               sfc_typ(iidx,jidx) = outbuf(15)
            endif
         endif
      
         if (mod(idim,8) .gt. 4 .or. mod(jdim,8) .gt. 4) then
            if (outbuf(15) .gt. -900. .and. sfc_typ(iidx,jidx) .lt. -900.) then
               sfc_typ(iidx,jidx) = outbuf(15)
            endif
         endif
        endif    

! End of the main loop - that's pretty short

        end do  !idim
120   end do    !jdim
    
!       print *, "End of pixel processing."

! -- find and flag cells that actually cross the dateline.

  if (maxval(longitude) > 179.0 .AND. minval(longitude,mask=longitude .gt. -900) < -179.0) then
    do j = 1, cscan
      do i = 1, cxscan
        iidx = (i-1)*cell_resolution + 1
        iidx2= iidx+cell_resolution-1
        if (iidx2 > xscan) iidx2 = xscan
        jidx = (j-1)*cell_resolution + 1
        jidx2= jidx+cell_resolution-1
        if (jidx2 > scan) jidx2 = scan
        if ((maxval(longitude(iidx:iidx2,jidx:jidx2),mask=longitude(iidx:iidx2,jidx:jidx2) .gt. -900)) - &
        & (minval(longitude(iidx:iidx2,jidx:jidx2),mask=longitude(iidx:iidx2,jidx:jidx2) .gt. -900)) &
        & > 355.0) then
          dateline(i,j) = 1
        end if
      end do
    end do
  endif

! -- average scan times on their own since it's only 1D
!   do i = 1, size(time)
!     iidx = (i-1)/8 + 1
!     if (platform .eq. 'VIIRS') time_avg(iidx) = time_avg(iidx) + time(i)
!     if (platform .eq. 'AHI' .or. platform .eq. 'GOES' ) time_avg(iidx) = time_avg(iidx) + 3.
! 
!     ntime(iidx) = ntime(iidx) + 1
!   end do
!   where (ntime > 0)
!     time_avg = time_avg / ntime
!   end where

!   -- calculate average lat/lon for each 10x10 pixel cell. Record number of smoke mask pixels as well.    
!    smoke_count(:,:) = 0
    smoke_ae_count(:,:) = 0
		pyrocb_count(:,:)     = 0				
		high_alt_smoke_count(:,:) = 0		
    do i = 1, size(latitude,1)
      do j = 1, size(latitude,2)

!     -- skip undefined lat/lon values.
      if (latitude(i,j) < -900.0 .OR. longitude(i,j) < -900.0) cycle
      
      iidx = (i-1)/cell_resolution + 1
      jidx = (j-1)/cell_resolution + 1

      if (iidx >cxscan .or. jidx>cscan) then 
        if (idim ==xscan  .and. jdim==scan) go to 130
        cycle
      end if

      lat_sav(iidx,jidx) = lat_sav(iidx,jidx) + latitude(i,j)

      !for aggregated pixel level debugging  
!       px_debug = .false.
!       if (iidx == 397 .and. jidx == 47) px_debug = .true.
!       if (px_debug) print*, oz(i,j), i, j
      
!     -- convert longitudes to all positive values if cell stradles the dateline, otherwise
!     -- values will cancel out.
      if (dateline(iidx,jidx) == 1 .AND. longitude(i,j) < 0.0) then
        lon_sav(iidx,jidx) = lon_sav(iidx,jidx) + longitude(i,j) + 360.0
      else
        lon_sav(iidx,jidx) = lon_sav(iidx,jidx) + longitude(i,j)
      end if
          
      sza_sav(iidx,jidx) = sza_sav(iidx,jidx) + sza(i,j)
      vza_sav(iidx,jidx) = vza_sav(iidx,jidx) + vza(i,j)
      raa_sav(iidx,jidx) = raa_sav(iidx,jidx) + raa(i,j)
      sca_sav(iidx,jidx) = sca_sav(iidx,jidx) + sca(i,j)
      
      ws_avg(iidx,jidx) = ws_avg(iidx,jidx) + windsp(i,j)
      wd_avg(iidx,jidx) = wd_avg(iidx,jidx) + winddir(i,j)
      oz_avg(iidx,jidx) = oz_avg(iidx,jidx) + oz(i,j)
      wv_avg(iidx,jidx) = wv_avg(iidx,jidx) + wv(i,j)
      elev_avg(iidx,jidx) = elev_avg(iidx,jidx) + elev(i,j)
!       chl_avg(iidx,jidx) = chl_avg(iidx,jidx) + chl(i,j)

      nll_sav(iidx,jidx) = nll_sav(iidx,jidx) + 1

!     count land and ocean elevation separately
      if (ls_mask(i,j) < 2) then   !land pixel
        elev_avg_land(iidx,jidx) = elev_avg_land(iidx,jidx) + elev(i,j)
        nll_sav_land(iidx,jidx) = nll_sav_land(iidx,jidx) + 1
      endif
     
      if (ls_mask(i,j) >= 2 .and. ls_mask(i,j) /= 5 .and. elev(i,j)< 10 ) then   !ocean pixel
        elev_avg_ocean(iidx,jidx) = elev_avg_ocean(iidx,jidx) + elev(i,j)
        nll_sav_ocean(iidx,jidx) = nll_sav_ocean(iidx,jidx) + 1
      endif  
                
!     -- count number of smoke_mask pixels in this cell and save.
      if (smoke_mask(i,j) == 1) then 
        smoke_count(iidx,jidx) = smoke_count(iidx,jidx) + 1
      end if
      
      if (smoke_ae_mask(i,j) == 1) then 
        smoke_ae_count(iidx,jidx) = smoke_ae_count(iidx,jidx) + 1
      end if
      
      if (pyrocb_mask(i,j) == 1) then 
        pyrocb_count(iidx,jidx) = pyrocb_count(iidx,jidx) + 1
      end if
      
      if (high_alt_smoke_mask(i,j) == 1) then 
        high_alt_smoke_count(iidx,jidx) = high_alt_smoke_count(iidx,jidx) + 1
      end if      
      
      end do
130 end do
    print *, 'out'    

    where (nll_sav > 0)
      lat_sav = lat_sav / nll_sav
      lon_sav = lon_sav / nll_sav
      sza_sav = sza_sav / nll_sav
      vza_sav = vza_sav / nll_sav
      raa_sav = raa_sav / nll_sav
      sca_sav = sca_sav / nll_sav     
      ws_avg = ws_avg / nll_sav
      wd_avg = wd_avg / nll_sav
      oz_avg = oz_avg / nll_sav
      wv_avg = wv_avg / nll_sav
      elev_avg = elev_avg / nll_sav
!       chl_avg = 10.**(chl_avg / nll_sav) ! Convert averaged value from log10 to linear for output diag
    elsewhere
      lat_sav = -999.0
      lon_sav = -999.0
      sza_sav = -999.0
      vza_sav = -999.0
      raa_sav = -999.0
      sca_sav = -999.0
      ws_avg = -999.0
      wd_avg = -999.0
      oz_avg = -999.0
      wv_avg = -999.0
      elev_avg = -999.0
    end where

!land elevation average
    where (nll_sav_land > 0)
      elev_avg_land = elev_avg_land / nll_sav_land
    elsewhere
      elev_avg_land = -999.0
    end where
!ocean elevation average
    where (nll_sav_ocean > 0)
      elev_avg_ocean = elev_avg_ocean / nll_sav_ocean
    elsewhere
      elev_avg_ocean = -999.0
    end where    
    
!   -- back out our longitude conversion to positive values where needed.
    if (count(dateline > 0) > 0) then
      where (lon_sav > 180.0)
        lon_sav = lon_sav - 360.0
      end where
    end if
      
! Compute AOT and SSA averages
    where (naot550_avg > nagg)
      aot550_avg = aot550_avg / naot550_avg
    elsewhere
      aot550_avg = -999.0
      naot550_avg = 0
    end where
    
    where (naot550_avg > nagg .AND. naot412_avg > 0)  
      aot412_avg = aot412_avg / naot412_avg
    elsewhere
      aot412_avg = -999.0
    end where
    
    where (naot550_avg > nagg .AND. naot488_avg > 0)  
      aot488_avg = aot488_avg / naot488_avg
    elsewhere
      aot488_avg = -999.0
    end where

    where (naot550_avg > nagg .AND. naot670_avg > 0)  
      aot670_avg = aot670_avg / naot670_avg
    elsewhere
      aot670_avg = -999.0
    end where
    
    where (nref_avg > nagg .AND. ref_avg > 0.0)
      ref_avg = ref_avg/nref_avg
    elsewhere
      ref_avg = -999.0
    end where
    
    where (nssa_avg > nagg)
      ssa_avg = ssa_avg / nssa_avg
    elsewhere
      ssa_avg = -999.0
    end where
    
    where (nsr_avg > nagg)
      sr_avg = sr_avg / nsr_avg
    elsewhere
      sr_avg = -999.0
    end where
    
    where (nae_avg > nagg) 
      ae_avg = ae_avg / nae_avg
    elsewhere
      ae_avg = -999.0
    end where
       
    where (naot550_avg > nagg) 
!       dstar_avg = dstar_avg / naot550_avg
      btd11_avg = btd11_avg / naot550_avg
      ratio_avg = ratio_avg / naot550_avg
      ndvi_avg = ndvi_avg / naot550_avg 
    elsewhere
!       dstar_avg = -999.0
      btd11_avg = -999.0
      ratio_avg = -999.0
      ndvi_avg = -999.0 
    end where

    where (nrcndvi_avg > 0)
      rcndvi_avg = rcndvi_avg / nrcndvi_avg
    elsewhere
      rcndvi_avg = -999.0
    end where
  
! aggregated reflectance average
    noref_sum(1,:,:)= noreflc
    noref_sum(2,:,:)= noreflc
    noref_sum(3,:,:)= noreflc
    noref_sum(4,:,:)= noreflc
    noref_sum(5,:,:)= noreflc
    noref_sum(6,:,:)= noreflc
    noref_sum(7,:,:)= noreflc
    noref_sum(8,:,:)= noreflc

    where (noref_sum > 0)
      oreflc_mean  = oreflc_mean/noref_sum
    elsewhere
      oreflc_mean = -999.0
    end where
    
    nref_sum(1,:,:)= nlreflc
    nref_sum(2,:,:)= nlreflc
    nref_sum(3,:,:)= nlreflc
    nref_sum(4,:,:)= nlreflc
    nref_sum(5,:,:)= nlreflc
    nref_sum(6,:,:)= nlreflc
    nref_sum(7,:,:)= nlreflc
    nref_sum(8,:,:)= nlreflc
    nref_sum(9,:,:)= nlreflc
    nref_sum(10,:,:)= nlreflc  
    
    where (nref_sum > nagg)
      lreflc_mean  = lreflc_mean/nref_sum
    elsewhere
      lreflc_mean = -999.0
    end where    
      
!   -- set cell algorithm flag over land and ocean.
    do j = 1, cscan
      do i = 1, cxscan
        !alg_flag_old (0=DB 1=Veg 2=mix)
        if (db_cnt_old(i,j) > 0 .OR. veg_cnt_old(i,j) > 0) then 
          if (db_cnt_old(i,j) > 0 .AND. veg_cnt_old(i,j) <= 0) then
            alg_flag_old(i,j) = 0
          else if (db_cnt_old(i,j) <= 0 .AND. veg_cnt_old(i,j) > 0) then
            alg_flag_old(i,j) = 1
          else
            alg_flag_old(i,j) = 2
          end if
        end if
        !alg_flag (1= AERONET BRDF, 2=BRDF coeff, 3=min table, 11= veg, 12, 2.2 micron, 20=mix)       
         if (db_cnt(i,j) > 0 .OR. veg_cnt(i,j) > 0 .OR. db_cnt2(i,j) > 0 .OR. &
            veg_cnt2(i,j) > 0 .OR. db_cnt3(i,j) > 0) then 
            alg_flag(i,j) = 20
            if (db_cnt(i,j) > 0 .AND. veg_cnt(i,j) <= 0 .AND. db_cnt2(i,j) <= 0 &
                    .AND. veg_cnt2(i,j) <= 0 .AND. db_cnt3(i,j) <= 0) alg_flag(i,j) = 1
            if (db_cnt(i,j) <= 0 .AND. veg_cnt(i,j) > 0 .AND. db_cnt2(i,j) <= 0 &
                    .AND. veg_cnt2(i,j) <= 0 .AND. db_cnt3(i,j) <= 0) alg_flag(i,j) = 11
            if (db_cnt(i,j) <= 0 .AND. veg_cnt(i,j) <= 0 .AND. db_cnt2(i,j) > 0 &
                    .AND. veg_cnt2(i,j) <= 0 .AND. db_cnt3(i,j) <= 0) alg_flag(i,j) = 2
            if (db_cnt(i,j) <= 0 .AND. veg_cnt(i,j) <= 0 .AND. db_cnt2(i,j) <= 0 &
                    .AND. veg_cnt2(i,j) > 0 .AND. db_cnt3(i,j) <= 0) alg_flag(i,j) = 12
            if (db_cnt(i,j) <= 0 .AND. veg_cnt(i,j) <= 0 .AND. db_cnt2(i,j) <= 0 &
                    .AND. veg_cnt2(i,j) <= 0 .AND. db_cnt3(i,j) > 0) alg_flag(i,j) = 3                   
         end if
        !alg_flag_old2 (0=DB, 1=Veg, 2=2.2 micron, 3=mix)    
        if (alg_flag(i,j) > -1) then 
          if (alg_flag(i,j) == 1) alg_flag_old2(i,j) = 0
          if (alg_flag(i,j) == 2) alg_flag_old2(i,j) = 0
          if (alg_flag(i,j) == 3) alg_flag_old2(i,j) = 0
          if (alg_flag(i,j) == 11) alg_flag_old2(i,j) = 1
          if (alg_flag(i,j) == 12) alg_flag_old2(i,j) = 2
          if (alg_flag(i,j) == 20) then 
            alg_flag_old2(i,j) = 3
            if (veg_cnt(i,j) <= 0 .AND. veg_cnt2(i,j) <= 0) alg_flag_old2(i,j) = 0
          end if
        end if    
         
      enddo          
    enddo
 
! Calc Standard Deviation of tau
! TODO: Replace with Welford, single-pass algorithm
      do j = 1, cscan
        do i = 1, cxscan
          
          if (naot550_avg(i,j) .gt. nagg .and. aot550_avg(i,j) .gt. -900.0) then
            var = 0.0       
            do m=1,naot550_avg(i,j)
              var = var + (aot550_avg(i,j) - aot550_list(i,j,m))**2
            enddo
            var = var/(naot550_avg(i,j)-1)
            sd550(i,j) = sqrt(var)
          endif
          
         enddo
      enddo   
    
! Calculate tau550 and alpha
      do j = 1, cscan
        do i = 1, cxscan   

	        if (lat_sav(i,j) > -900.0 .AND. lon_sav(i,j) > -900.0) then
            gzone(i,j) = get_geographic_zone(lat_sav(i,j), lon_sav(i,j), status)          
            if (status /= 0) then
              print *, "ERROR: Failed to get geographic zone: ", lat_sav(i,j), lon_sav(i,j), status
              cycle
            end if
          end if

!           -- calculate AE for 550nm
            alpha_550  = -999.
            alpha2_550 = -999.
            xxrat  = (aot412_avg(i,j)) / (aot488_avg(i,j))
            xxrat2 = (aot670_avg(i,j)) / (aot488_avg(i,j))
            if (xxrat < 0.0 .AND. xxrat2 < 0.0) then
              aot550_avg(i,j) = -999.0
              aot412_avg(i,j) = -999.0
              aot488_avg(i,j) = -999.0
              aot670_avg(i,j) = -999.0
              ae_avg(i,j)     = -999.0
              cycle ! abort if ratios < 0.0, fill values in aot.
            end if
            if (aot412_avg(i,j) < -900.0 .AND. aot488_avg(i,j) < -900.0 .AND. aot670_avg(i,j) < -900.0) then
              aot550_avg(i,j) = -999.0
              ae_avg(i,j)     = -999.0
              cycle ! abort if ratios < 0.0, fill values in aot.
            end if


            if (aot412_avg(i,j) .gt. 0 .and. xxrat .gt. 0.0) then
               dd     = alog(412./490.)
               alpha_550  = alog(xxrat)
               alpha_550  = -1.*alpha_550/dd
               if (aot488_avg(i,j) .lt. 0.2) alpha_550  = 1.0
               if (alpha_550 .gt. 1.8) alpha_550  = 1.8
            endif
            if (aot670_avg(i,j) .gt. 0 .and. xxrat2 .gt. 0.0) then
               dd     = alog(670./490.)
               alpha2_550 = alog(xxrat2)
               alpha2_550 = -1.*alpha2_550/dd
               if (aot488_avg(i,j) .lt. 0.2) alpha2_550  = 1.0
               if (alpha2_550 .gt. 1.8) alpha2_550  = 1.8
            endif
            if (alpha_550 .lt. -900.0 .and. alpha2_550 .gt. -900.0) alpha_550 = alpha2_550
            
            if (aot412_avg(i,j) .lt. 0 .and. aot488_avg(i,j) .lt. 0) then
               if (ae_avg(i,j) .gt. -900.0) alpha_550 = ae_avg(i,j)
            endif
            if (aot550_avg(i,j) .gt. 1.8 .and. ae_avg(i,j) .gt. -900.0) then
               alpha_550 = ae_avg(i,j)
            endif
            if (alpha_550 .lt. -0.4 .and. alpha_550 .gt. -900.0) alpha_550 = -0.4
            if (alpha_550 .gt. 1.8) alpha_550  = 1.8
            
!           -- case of thin dust over very bright surface. Reset AE to 1.0.
            if (aot488_avg(i,j) .lt. 0.7 .and. sr_avg(i,j,3) .ge. 0.251 .and. &
                alpha_550 .gt. 1.0) alpha_550  = 1.0
            
!           -- now calculate AE for the rest of the bands
            alpha  = -999.
            alpha2 = -999.
            xxrat  = (aot412_avg(i,j)) / (aot488_avg(i,j))
            xxrat2 = (aot670_avg(i,j)) / (aot488_avg(i,j))
            
            if (aot412_avg(i,j) .gt. 0 .and. xxrat .gt. 0.0) then
               dd     = alog(412./490.)
               alpha  = alog(xxrat)
               alpha  = -1.*alpha/dd
               if (aot488_avg(i,j) .lt. 0.2) alpha  = 1.5
               if (alpha .gt. 1.8) alpha  = 1.8
            endif
            if (aot670_avg(i,j) .gt. 0 .and. xxrat2 .gt. 0.0) then
               dd     = alog(670./490.)
               alpha2 = alog(xxrat2)
               alpha2 = -1.*alpha2/dd
               if (aot488_avg(i,j) .lt. 0.2) alpha2  = 1.5
               if (alpha2 .gt. 1.8) alpha2  = 1.8
            endif
            if (alpha .lt. -900.0 .and. alpha2 .gt. -900.0) alpha = alpha2

            if (aot412_avg(i,j) .lt. 0 .and. aot488_avg(i,j) .lt. 0) then
               if (ae_avg(i,j) .gt. -900.0) alpha = ae_avg(i,j)
            endif
            if (aot550_avg(i,j) .gt. 1.8 .and. ae_avg(i,j) .gt. -900.0) then
               alpha = ae_avg(i,j)
            endif
            if (alpha .lt. 0.0 .and. alpha .gt. -900.0) alpha = 0.0
            if (alpha .gt. 1.8) alpha  = 1.8

!           -- case of thin dust over very bright surface. Reset AE to 1.0.
            if (aot488_avg(i,j) .lt. 0.7 .and. sr_avg(i,j,3) .ge. 0.251 .and. &
                alpha .gt. 1.0) then
                  alpha  = 1.0
            end if
            
!           -- D* says dust. Set AE to dust. Only over N.Africa and Arabia.
!             if ((alpha > -900.0 .AND. dstar_avg(i,j) > 1.02) .AND. &
            if ((alpha > -900.0 ) .AND. &
            &   ((gzone(i,j) >= 1 .AND. gzone(i,j) <= 5) .OR. (gzone(i,j) == 26 .OR. gzone(i,j) == 27))) then
              alpha = 0.0
            end if
            
!           -- if AOT is maxed out for all channels, AE will be 0.0. But what if
!           -- it's smoke? Use D* to detect smoke (D* < 1.06) and set AE accordingly.
!           -- also reset depending on smoke_count.
!             if (aot412_avg(i,j) .gt. 3.47 .AND. dstar_avg(i,j) .lt. 1.06) then
!               alpha = 1.8
!             endif
            if (aot550_avg(i,j) .ge. 0.7) then
              if (smoke_count(i,j) .gt. 1 .OR. high_alt_smoke_count(i,j) .gt. 1) alpha = 1.8			
            else 
              if (smoke_count(i,j) .gt. 5 .OR. high_alt_smoke_count(i,j) .gt. 5) alpha = 1.8			
            end if
            
!           -- reset AE if smoke detected and AOT is high. AOT threshold needed
!           -- otherwise clear vegetated regions will be reset.
            if (smoke_ae_count(i,j) > 0 .AND. aot550_avg(i,j) > 0.2) then
              alpha = 1.8
            end if

!           -- reset AE if high LER 412/488nm ratio is detected
            if (ratio_avg(i,j) >= 0.95 .AND. alpha < 0.7 .AND. aot550_avg(i,j) > 0.2) then
              alpha = ratio_avg(i,j) + 0.6
            endif
            
!           -- apply AE to AOT and get everything consistent.
            if (aot550_avg(i,j) .gt. 0.0 .and. alpha .gt. -900.) then 
               dd     = 550./500.
               dd     = dd**(-1.*alpha_550)
               aot550_avg(i,j)  = aot550_avg(i,j) * dd
               dd     = 412./550.
               dd     = dd**(-1.*alpha)
               aot412_avg(i,j)  = aot550_avg(i,j) * dd
               dd     = 490./550.
               dd     = dd**(-1.*alpha)
               aot488_avg(i,j)  = aot550_avg(i,j) * dd
               dd     = 670./550.
               dd     = dd**(-1.*alpha)
               aot670_avg(i,j)  = aot550_avg(i,j) * dd
            else
              aot412_avg(i,j) = -999.0
              aot488_avg(i,j) = -999.0
              aot670_avg(i,j) = -999.0
              aot550_avg(i,j)  = -999.0
              naot550_avg(i,j) = 0
            endif
            
            if (alpha .lt. 0.0) then
               ae_avg(i,j) = -9999
            else
               ae_avg(i,j) = alpha
            endif

         enddo
         
      enddo

!     -- reset high AOT's to 5.0
      where (aot412_avg > 5.0)  aot412_avg    = 5.0
      where (aot488_avg > 5.0)  aot488_avg    = 5.0
      where (aot670_avg > 5.0)  aot670_avg    = 5.0
      where (aot550_avg > 5.0)  aot550_avg  = 5.0
    
!     -- cloud edge screening
!     -- skip if smoke has been detected and we're in an applicable zone.
!     -- skip is D* is high enough (indicating aerosols) in that cell.
      do k = 1, 3
        do j = 2, cscan - 1
          do i = 2, cxscan - 1
            
            if (platform .eq. 'VIIRS') then 
              if (gzone(i,j) == 13 .OR. gzone(i,j) == 18 .OR. gzone(i,j) == 31) cycle     ! ConUS, Fresno, Yuma
              if (gzone(i,j) == 14) cycle                           ! S. America
              if (gzone(i,j) == 29) cycle                           ! Mexico City
              if (gzone(i,j) == 15 .and. smoke_count(i,j) > 1 .and. naot550_avg(i,j) > 20) cycle    !  India 
              if (gzone(i,j) == 12 .and. smoke_count(i,j) > 1 .and. naot550_avg(i,j) > 20) cycle    !  Australia

            end if

!             if (gzone(i,j) >= 1 .AND. gzone(i,j) <= 11) then
!               if (dstar_avg(i,j) .gt. 1.06) cycle
!             else
!               if (dstar_avg(i,j) .gt. 1.1) cycle
!             end if
            
            if (k == 1) then
              if (aot412_avg(i,j) .gt. 1.2 .and. &
                  abs(aot412_avg(i-1,j)-aot412_avg(i,j)) .gt. 0.8 .and. &
                  abs(aot412_avg(i,j)-aot412_avg(i+1,j)) .gt. 0.8) then 
                aot412_avg(i,j) = -999.0
                aot488_avg(i,j) = -999.0
                aot670_avg(i,j) = -999.0
                aot550_avg(i,j) = -999.0
                naot550_avg(i,j) = 0          
                ae_avg(i,j) = -999.0
              endif
              if (aot412_avg(i,j) .gt. 1.2 .and. &
                  abs(aot412_avg(i,j+1)-aot412_avg(i,j)) .gt. 0.8 .and. &
                  abs(aot412_avg(i,j)-aot412_avg(i,j-1)) .gt. 0.8) then 
                aot412_avg(i,j) = -999.0
                aot488_avg(i,j) = -999.0
                aot670_avg(i,j) = -999.0
                aot550_avg(i,j) = -999.0
                naot550_avg(i,j) = 0
                ae_avg(i,j) = -999.0

              end if
             
            end if

            if (k == 2) then
              if (aot488_avg(i,j) .gt. 1.2 .and. &
                  abs(aot488_avg(i-1,j)-aot488_avg(i,j)) .gt. 0.8 .and. &
                  abs(aot488_avg(i,j)-aot488_avg(i+1,j)) .gt. 0.8) then 
                aot412_avg(i,j) = -999.0
                aot488_avg(i,j) = -999.0
                aot670_avg(i,j) = -999.0
                aot550_avg(i,j) = -999.0
                naot550_avg(i,j) = 0          
                ae_avg(i,j) = -999.0
              end if
              if (aot488_avg(i,j) .gt. 1.2 .and. &
                  abs(aot488_avg(i,j+1)-aot488_avg(i,j)) .gt. 0.8 .and. &
                  abs(aot488_avg(i,j)-aot488_avg(i,j-1)) .gt. 0.8) then 
                aot412_avg(i,j) = -999.0
                aot488_avg(i,j) = -999.0
                aot670_avg(i,j) = -999.0
                aot550_avg(i,j) = -999.0
                naot550_avg(i,j) = 0
                ae_avg(i,j) = -999.0
              end if
             
            end if
            
            if (k == 3) then
              if (aot670_avg(i,j) .gt. 1.2 .and. &
                  abs(aot670_avg(i-1,j)-aot670_avg(i,j)) .gt. 0.8 .and. &
                  abs(aot670_avg(i,j)-aot670_avg(i+1,j)) .gt. 0.8) then 
                aot412_avg(i,j) = -999.0
                aot488_avg(i,j) = -999.0
                aot670_avg(i,j) = -999.0
                aot550_avg(i,j) = -999.0
                naot550_avg(i,j) = 0          
                ae_avg(i,j) = -999.0
              end if
              if (aot670_avg(i,j) .gt. 1.2 .and. &
                  abs(aot670_avg(i,j+1)-aot670_avg(i,j)) .gt. 0.8 .and. &
                  abs(aot670_avg(i,j)-aot670_avg(i,j-1)) .gt. 0.8) then 
                aot412_avg(i,j) = -999.0
                aot488_avg(i,j) = -999.0
                aot670_avg(i,j) = -999.0
                aot550_avg(i,j) = -999.0
                naot550_avg(i,j) = 0
                ae_avg(i,j) = -999.0
              end if
             
            end if
            
          enddo
        enddo
      enddo
     
      do j = 2, cscan - 1
        do i = 2, cxscan - 1
          
          if (platform .eq. 'VIIRS') then
            if (gzone(i,j) == 13 .OR. gzone(i,j) == 18 .OR. gzone(i,j) == 31) cycle     ! ConUS, Fresno, Yuma
            if (gzone(i,j) == 14) cycle                           ! S. America
            if (gzone(i,j) == 29) cycle                           ! Mexico City
            if (gzone(i,j) == 15 .and. smoke_count(i,j) > 1 .and. naot550_avg(i,j) > 20) cycle    !  India
            if (gzone(i,j) == 12 .and. smoke_count(i,j) > 1 .and. naot550_avg(i,j) > 20) cycle    !  Australia
          end if
          
!           if (gzone(i,j) >= 1 .AND. gzone(i,j) <= 11) then
!             if (dstar_avg(i,j) .gt. 1.06) cycle
!           else
!             if (dstar_avg(i,j) .gt. 1.1) cycle
!           end if  
            
          if (aot550_avg(i,j) .gt. 1.2 .and. &
              abs(aot550_avg(i-1,j)-aot550_avg(i,j)) .gt. 0.8 .and. &
              abs(aot550_avg(i,j)-aot550_avg(i+1,j)) .gt. 0.8) then
            aot550_avg(i,j) = -999.
            aot412_avg(i,j) = -999.0
            aot488_avg(i,j) = -999.0
            aot670_avg(i,j) = -999.0
            naot550_avg(i,j) = 0
            ae_avg(i,j) = -999.0

          endif
          if (aot550_avg(i,j) .gt. 1.2 .and. &
              abs(aot550_avg(i,j+1)-aot550_avg(i,j)) .gt. 0.8 .and. &
              abs(aot550_avg(i,j)-aot550_avg(i,j-1)) .gt. 0.8) then
            aot550_avg(i,j) = -999.
            aot412_avg(i,j) = -999.0
            aot488_avg(i,j) = -999.0
            aot670_avg(i,j) = -999.0
            naot550_avg(i,j) = 0
            ae_avg(i,j) = -999.0
          endif
          
          if (platform .ne. 'VIIRS') then
           if (aot550_avg(i,j) .gt. 0.8 .and. &
               abs(aot550_avg(i-1,j)-aot550_avg(i,j)) .gt. 10. .and. &
               abs(aot550_avg(i,j)-aot550_avg(i+1,j)) .gt. 10.) then
             aot550_avg(i,j) = -999.
             aot412_avg(i,j) = -999.0
             aot488_avg(i,j) = -999.0
             aot670_avg(i,j) = -999.0
             naot550_avg(i,j) = 0
             ae_avg(i,j) = -999.0

           endif
           if (aot550_avg(i,j) .gt. 0.8 .and. &
               abs(aot550_avg(i,j+1)-aot550_avg(i,j)) .gt. 10. .and. &
               abs(aot550_avg(i,j)-aot550_avg(i,j-1)) .gt. 10.) then
             aot550_avg(i,j) = -999.
             aot412_avg(i,j) = -999.0
             aot488_avg(i,j) = -999.0
             aot670_avg(i,j) = -999.0
             naot550_avg(i,j) = 0
             ae_avg(i,j) = -999.0
           endif   
          endif
                 
        enddo
      enddo

! Replace all variables with fill value if aot550_avg undefined  
      do j = 1, cscan
        do i = 1, cxscan
        
					 if (aot550_avg(i,j) .lt. -900.0) then
							aot412_avg(i,j) = -999.0
              aot488_avg(i,j) = -999.0
              aot670_avg(i,j) = -999.0
              naot412_avg(i,j) = 0   
              naot488_avg(i,j) = 0   
              naot670_avg(i,j) = 0   

							do k = 1, 3
								 ssa_avg(i,j,k) = -999.0
								 sr_avg(i,j,k)  = -999.0
								 ref_avg(i,j,k) = -999.0
							enddo ! k
							ae_avg(i,j) = -999.0
							sd550(i,j)  = -999.0
							
! 							dstar_avg(i,j) = -999.0
              btd11_avg(i,j) = -999.0
              ndvi_avg(i,j) = -999.0 !jlee added
							aot550_avg(i,j) = -999.0
							naot550_avg(i,j) = 0
							alg_flag(i,j) = -999
							alg_flag_old(i,j) = -999			
							alg_flag_old2(i,j) = -999				
					 endif
					 
         enddo ! i
      enddo ! j

    edgecnt(:,:) = 0
    do j = 2, cscan-1
      do i = 2, cxscan-1
        edgecnt(i,j) = 0
        
        if (aot550_avg(i,j) .gt. -1 .AND. aot550_avg(i,j) .lt. 1.0) then

          if (abs(aot550_avg(i,j)-aot550_avg(i+1,j)) .gt. 0.4) then
            edgecnt(i,j)=edgecnt(i,j)+1
          endif
          if (abs(aot550_avg(i,j)-aot550_avg(i-1,j)) .gt. 0.4) then
            edgecnt(i,j)=edgecnt(i,j)+1
          endif
          if (abs(aot550_avg(i,j)-aot550_avg(i,j+1)) .gt. 0.4) then
            edgecnt(i,j)=edgecnt(i,j)+1
          endif
          if (abs(aot550_avg(i,j)-aot550_avg(i,j-1)) .gt. 0.4) then
            edgecnt(i,j)=edgecnt(i,j)+1
          endif
          
          if (abs(aot550_avg(i,j)-aot550_avg(i+1,j+1)) .gt. 0.4) then
            edgecnt(i,j)=edgecnt(i,j)+1
          endif
          if (abs(aot550_avg(i,j)-aot550_avg(i+1,j-1)) .gt. 0.4) then
            edgecnt(i,j)=edgecnt(i,j)+1
          endif
          if (abs(aot550_avg(i,j)-aot550_avg(i-1,j+1)) .gt. 0.4) then
            edgecnt(i,j)=edgecnt(i,j)+1
          endif
          if (abs(aot550_avg(i,j)-aot550_avg(i-1,j-1)) .gt. 0.4) then
            edgecnt(i,j)=edgecnt(i,j)+1
          endif
        else if (aot550_avg(i,j) .gt. -1.0 .AND. aot550_avg(i,j) .ge. 1.0) then
          if (abs(aot550_avg(i,j)-aot550_avg(i+1,j))/aot550_avg(i,j) .gt. 0.3) then
            edgecnt(i,j)=edgecnt(i,j)+1
          endif
          if (abs(aot550_avg(i,j)-aot550_avg(i-1,j))/aot550_avg(i,j) .gt. 0.3) then
            edgecnt(i,j)=edgecnt(i,j)+1
          endif
          if (abs(aot550_avg(i,j)-aot550_avg(i,j+1))/aot550_avg(i,j) .gt. 0.3) then
            edgecnt(i,j)=edgecnt(i,j)+1
          endif
          if (abs(aot550_avg(i,j)-aot550_avg(i,j-1))/aot550_avg(i,j) .gt. 0.3) then
            edgecnt(i,j)=edgecnt(i,j)+1
          endif
          
          if (abs(aot550_avg(i,j)-aot550_avg(i+1,j+1))/aot550_avg(i,j) .gt. 0.3) then
            edgecnt(i,j)=edgecnt(i,j)+1
          endif
          if (abs(aot550_avg(i,j)-aot550_avg(i+1,j-1))/aot550_avg(i,j) .gt. 0.3) then
            edgecnt(i,j)=edgecnt(i,j)+1
          endif
          if (abs(aot550_avg(i,j)-aot550_avg(i-1,j+1))/aot550_avg(i,j) .gt. 0.3) then
            edgecnt(i,j)=edgecnt(i,j)+1
          endif
          if (abs(aot550_avg(i,j)-aot550_avg(i-1,j-1))/aot550_avg(i,j) .gt. 0.3) then
            edgecnt(i,j)=edgecnt(i,j)+1
          endif
        end if
      end do
    end do

! Set QA Flags
  
    do j = 1, cscan
      do i = 1, cxscan
       
          flags(1) = 0
          flags(2) = 0
          flags(3) = 0
          flags(4) = 0
            
          if (naot550_avg(i,j) > 0) then
            flags(1) = 1
            flags(2) = 2
            call set_qa_land(cxscan,naot550_avg(i,j),sd550(i,j), &
                      aot550_avg(i,j),edgecnt(i,j),0, &
                      smoke_count(i,j), high_alt_smoke_count(i,j), gzone(i,j),flags(2), platform)
            if (ae_avg(i,j).le.0.5) flags(3)= 1
            if (ae_avg(i,j).gt.1.2) flags(3)= 2
            if (r412_count(i,j) .gt. 20) flags(4) = 1
            if (naot488_avg(i,j) .lt. 50) flags(4) = 2
            if (100.*t650_count(i,j)/naot550_avg(i,j) .gt. 50.) flags(4) = 3     !Set out of bound cell   

            !apply outbound flag to quality flag (Jan 2018 WKim)
            if (flags(2) == 2 .or. flags(2) == 3) then
              if (flags(4) == 3) then
                flags(2) = 1
                !very high AOD case over N_Africa and Australia
                if ((gzone(i,j) < 6 .and. 0 < gzone(i,j)) .or. gzone(i,j) == 12) then   
                  flags(2) = 2
                endif    
              endif
            endif
            
            conf_flag(i,j)  = flags(2)
            usefulness(i,j) = flags(1)
            aero_type(i,j) = flags(3)
            rcond(i,j)    = flags(4)          

          endif
          if (o_out%aot550(i,j) > 0) then
            call set_qa_ocean(lat_sav(i,j), 0.0, o_out%aot550(i,j), o_out%npixels(i,j),o_out%aot550_stdv(i,j), &
            & n_total_pixels(i,j), o_out%ss(i,j), o_out%alg_flag(i,j), oconf_flag(i,j), platform)
          end if                            
    
       enddo
    enddo
 
!   -- reset QA if majority of pixels in cell have were retrieved by the veg code but
!   -- were over bright surfaces (bflag > 0).
    where (aot550_avg > 0.0 .AND. veg_bflag_cnt/real(naot550_avg) > 0.5)
      conf_flag = 1
    end where

!     -- reset QA to 1 if cell is adjacent to an empty (cloudy) pixel, i.e QA=0.
!     -- set usefulness flag to 0 if QA==0,1.  then save QA to MYD04 file.
    edgecnt(:,:) = 0
    do j = 1, cscan
      do i = 1, cxscan

!        -- skip cells over N.Africa with Dstar > 1.2 -- indicate dust. Don't reset.
!        -- Otherwise, use 1.06 threshold elsewhere.
!         if ((gzone(i,j) >= 1 .AND. gzone(i,j) <= 5) .OR. (gzone(i,j) == 26 .OR. gzone(i,j) == 27)) then
!           if (dstar_avg(i,j) > 1.2) cycle
!         else
!           if (dstar_avg(i,j) >= 1.06) cycle
!         endif

!       -- do not reset QA to 1 if smoke pixels are present in zones 13, 18 (ConUS)
        if (gzone(i,j) == 13 .OR. gzone(i,j) == 18 .OR. gzone(i,j) == 29) then 
          if (smoke_count(i,j) > 20) cycle ! 21 February 2018 JLee  threshod changed from 0 to 20, preliminary
        end if
        
        if (gzone(i,j) == 15 .and. smoke_count(i,j) > 1 .and. naot550_avg(i,j) > 20) cycle    ! India
        if (gzone(i,j) == 12 .and. smoke_count(i,j) > 1 .and. naot550_avg(i,j) > 20) cycle    ! Australia


        ! 16 May 2017 JLee
        ! to turn off the cloud edge test for high elevation region in geozone 30;
        ! cloud edge test removed too many data points because AOD from swir-vis
        ! surface database is sparse. We might be able to turn the filter on if
        ! surface elevation is correctly handled in the retrieval algorithm and
        ! the surface database (i.e. Rayleigh correction). We found that the
        ! swir-vis surface relationship can be derived much more accurately if
        ! considering the surface elevation, such that more AOD data can be
        ! retrieved.
        px_elev = -999.0
        px_elev = get_elevation(lat_sav(i,j), lon_sav(i,j), status)
        if (status /= 0) then
          print *, "ERROR: Failed to get elevation for pixel, ", i, j, lat_sav(i,j), lon_sav(i,j), status
        end if

        if (gzone(i,j) == 30 .AND. px_elev > 750.0) cycle 
        ! end jlee added 
       
        i1 = max(i-1,1)
        i2 = min(i+1,cxscan)
        j1 = max(j-1,1)
        j2 = min(j+1,cscan)
        if (conf_flag(i,j) /= 0) then
          cnt = count(conf_flag(i1:i2,j1:j2) == 0)
          if (cnt > 0) then
            conf_flag(i,j) = 1
          end if
        end if
      
      end do
    end do

    do j = 1, cscan
      do i = 1, cxscan
      if (gzone(i,j) == 12 .or. gzone(i,j) == 15) then
        i1 = max(i-2,1)
        i2 = min(i+2,cxscan)
        j1 = max(j-2,1)
        j2 = min(j+2,cscan)
        cnt  = count(smoke_count(i1:i2,j1:j2) > 1)
        cnt2 = count(naot550_avg(i1:i2,j1:j2) > 20)
          if (cnt > 1 .and. cnt2 > 16 .and. conf_flag(i,j) == 1) then
            conf_flag(i,j) = 2
          end if
        cnt3 = count(conf_flag(i1:i2,j1:j2) > 1)
          if (aot550_avg(i,j) > 0.3 .and. conf_flag(i,j) > 1 .and. cnt > 1 .and. cnt3 < 11) then
            conf_flag(i,j) = 1
          end if
      endif

      end do
    end do


!   -- reset confidence flag to 1 if pyrocb pixels are present.
    where (pyrocb_count > 0 .and. aot550_avg > -1) conf_flag = 1

!   -- update usefulness flags
    where (conf_flag == 1) usefulness = 0
    
!   -- copy AOT550 value to best estimate array if QA=3 
    where (conf_flag == 3 .OR. conf_flag == 2) 
      aot550_best = aot550_avg
      ae_best = ae_avg
    end where
    !where(oconf_flag == 3)
    !updated for NRT retrieval. Standard Ocean retrieval only have oconf_flag 1 or 3
    !sunglint area in NRT retrieval will be added as QA=2 and these pixels will be included as best estimate
    where(oconf_flag >= 2)
      oaot550_best = o_out%aot550
      oae_best = o_out%ae
      ofmf_best = o_out%fmf
    end where

!   -- create combined land/ocean datasets
    cls_mask(:,:) = -999
    do j = 1, cscan
      do i = 1, cxscan
        i1 = (i-1)*cell_resolution  + 1
        i2 = i1 + cell_resolution -1
        if (i2 > xscan) i2 = xscan
        j1 = (j-1)*cell_resolution  + 1
        j2 = j1 + cell_resolution -1
        if (j2 > scan) j2 = scan  
        
        if (real(n_total_pixels(i,j))==0.0) cycle
        if (count(ls_mask(i1:i2,j1:j2) < 2 .AND. ls_mask(i1:i2,j1:j2) > -900) / real(n_total_pixels(i,j)) > 0.5) then          
          cls_mask(i,j) = 1     ! land
        else
          cls_mask(i,j) = 0     ! water
        end if
      end do
    end do
   
    ! Combined land and ocean data sets.
    ! Land
    where (cls_mask == 1)
      caot550_avg = aot550_avg
      caot550_best = aot550_best
      cae_avg = ae_avg
      cae_best = ae_best
    end where
    ! Ocean
    where (cls_mask == 0)
      caot550_avg = o_out%aot550
      caot550_best = oaot550_best
      cae_avg = o_out%ae
      cae_best = oae_best
    end where

!   -- sync up coverage of Rayleigh-corrected NDVI with AOT550.
    where (aot550_avg < -900.0)
      rcndvi_avg = -999.0
    end where

!     -- create aerosol type SDS
      aerosol_type(:,:) = -999
      where (aot550_avg > -900)                                 aerosol_type = 5 ! mixed, default
      where (ae_avg > 1.2 .and. aot550_avg > 0.4)               aerosol_type = 4 ! non-smoke fine mode
      where (smoke_count > 0 .and. aot550_avg > -900)           aerosol_type = 1 ! smoke
      where (high_alt_smoke_count > 0 .and. aot550_avg > -900)  aerosol_type = 2 ! high altitude smoke
      where (pyrocb_count > 0 .and. aot550_avg > -900)          aerosol_type = 3 ! pyrocumulonimbus clouds
      where (aot550_avg > -900 .and. aot550_avg < 0.2)          aerosol_type = 6 ! background
!       where (dstar_avg > 1.1 .and. aot550_avg > -900)           aerosol_type = 0 ! dust
      where (ae_avg < 0.1 .and. ratio_avg < 0.78 .and. aot550_avg > 0.4) aerosol_type = 0    ! dust
      where (ae_avg < 0.5 .and. gzone > 0 .and. gzone < 12 .and. aot550_avg > 0.2)     aerosol_type = 0    ! dust
      where (gzone == 15 .and. ae_avg <= 1.2 .and. aerosol_type == 1)      aerosol_type = 5 ! reset to mix type

!     -- create combined aerosol type SDS      
      combined_type(:,:) = -999
      combined_type      = aerosol_type
      where (o_out%model_flag == 1)         combined_type  = 0 ! ocean dust -> Dust
      where (o_out%model_flag == 3)         combined_type  = 6 ! ocean maritime -> background
      where (o_out%model_flag == 4)         combined_type  = 5 ! ocean mixed -> mixed
      where (o_out%model_flag == 2)         combined_type  = 7 ! fine dominated
      where (caot550_avg < -900)            combined_type = -999  

!   match reflectance coverage and all QA AOD coverage
      oaot550_avg   = o_out%aot(:,:,2)
      nref_sum(1,:,:)= oaot550_avg
      nref_sum(2,:,:)= oaot550_avg
      nref_sum(3,:,:)= oaot550_avg
      nref_sum(4,:,:)= oaot550_avg
      nref_sum(5,:,:)= oaot550_avg
      nref_sum(6,:,:)= oaot550_avg
      nref_sum(7,:,:)= oaot550_avg
      nref_sum(8,:,:)= oaot550_avg
      
      where (nref_sum > -900.)
        oreflc_mean  = oreflc_mean
      elsewhere
!         oreflc_mean = -999.0
      end where
    
      nref_sum(1,:,:)= aot550_avg
      nref_sum(2,:,:)= aot550_avg
      nref_sum(3,:,:)= aot550_avg
      nref_sum(4,:,:)= aot550_avg
      nref_sum(5,:,:)= aot550_avg
      nref_sum(6,:,:)= aot550_avg
      nref_sum(7,:,:)= aot550_avg
      nref_sum(8,:,:)= aot550_avg
      where (nref_sum > -900.)
        lreflc_mean  = lreflc_mean
      elsewhere
!         lreflc_mean = -999.0
      end where    
!     lat_sav = dtlat
!     lon_sav = dtlon


    status = create_dbdt_pace(lat_sav, lon_sav, aot488_avg, aot550_avg, aot670_avg, conf_flag,aerosol_type, &
    &                    dtaod,dt_qa, DBDTaot,DBDTqa, DBDTflag,dtspec,lreflc_mean,dtfmf,dbdtfmf,dt_cldmsk,&
    &                    lcld_mask,lat, lon,dbdt_refl, sza_sav )    
    uvaod(:,:,:)  = -999.
    status = extrapolate_uv(DBDTaot,uvaod)    

!    print *, 'output file: ', output_file

!    status = create_viirs_l2(output_file, time_avg, lat_sav, lon_sav, sza_sav, vza_sav, raa_sav,  &
!    &          sca_sav, aot550_avg, aot550_best, naot550_avg, aot412_avg, aot488_avg, aot670_avg, sd550, &
!    &          ae_avg, ae_best, ssa_avg, conf_flag, alg_flag, o_out%alg_flag, o_out%aot,o_out%aot550, oaot550_best,   &
!    &          o_out%ae, oae_best, o_out%fmf, ofmf_best, o_out%aot550_stdv, o_out%npixels,    &
!    &          o_out%ss, oconf_flag, o_out%model_flag, DBDTaot(:,:,2), DBDTaot(:,:,2), cae_avg, cae_best, n_total_pixels, &
!    &          ws_avg, oz_avg, wv_avg, wd_avg, ndvi_avg, rcndvi_avg, ndvi, sr_avg, dstar_avg, btd11_avg, &
!    &          turbid_res, elev_avg,  aerosol_type, combined_type,elev_avg_land, &
!    &          elev_avg_ocean, smoke_count,oreflc_mean, lreflc_mean(1:8,:,:),hires, alg_flag_old,alg_flag_old2,platform, &
!    &          o_out%fmf,DBDTaot,DBDTqa, DBDTflag,dbdt_refl,dbdtfmf)
!     &          o_out%opt_err)


!    if (status /= 0) then
!      print *, "ERROR: Failed to create VIIRS DB L2 output file: ", status
!      return
!    end if

    return
  end
